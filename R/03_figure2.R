# =============================================================================
# 03_figure2.R  --  cross-species validation and forest plot
#
# Computes, does not hard-code. No GEOquery, no XML.
#
# Discovery column is read from results/concentration_response.csv (verified).
# Human column   : GSE244740, DESeq2, clozapine highest concentration vs vehicle.
# Rat column     : GSE59905 / GPL5426, limma, clozapine heart vs dose-0 control.
#
# The two validation sections are independent and each is wrapped so that a
# failure in one still lets the other finish and write its CSV. Read the
# diagnostics printed before each fit; they are there to be checked, not skipped.
#
#   Rscript 03_figure2.R 2>&1 | tee logs/fig2.log
# =============================================================================

WORKDIR <- getwd()
DATA <- file.path(WORKDIR, "data"); RES <- file.path(WORKDIR, "results")
FIGS <- file.path(WORKDIR, "figures")
for (d in c(RES, FIGS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
options(timeout = 1800)

lib <- Sys.getenv("R_LIBS_USER", unset = file.path(path.expand("~"), "Rlibs", "R-4.6.0"))
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))

GENES <- c("TNFRSF12A", "EGLN3", "VEGFA")
# Ensembl ids for the three genes. Printed back with their row sums below so a
# wrong id shows up as an absent or empty row rather than as a silent zero.
ENS <- c(TNFRSF12A = "ENSG00000006327",
         EGLN3     = "ENSG00000129521",
         VEGFA     = "ENSG00000112715")

hr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

# ------------------------------------------------ series matrix helpers -----
sm_lines <- function(path) readLines(gzfile(path), warn = FALSE)

sm_meta <- function(ln) {
  b <- grep("^!series_matrix_table_begin", ln)
  if (length(b)) ln <- ln[seq_len(b[1] - 1)]
  ln <- ln[grepl("^!Sample_", ln)]
  p <- strsplit(ln, "\t")
  list(keys = vapply(p, function(x) sub("^!", "", x[1]), character(1)),
       vals = lapply(p, function(x) gsub('^"|"$', "", x[-1])))
}

# Pull a characteristic by its "prefix:" label. GEO characteristic columns are
# ragged - the same variable appears in different positions for different
# samples - so position indexing is not safe here.
sm_char <- function(m, prefix) {
  ch <- m$vals[m$keys == "Sample_characteristics_ch1"]
  if (!length(ch)) return(character(0))
  out <- rep(NA_character_, length(ch[[1]]))
  pat <- paste0("^", prefix, "\\s*:\\s*")
  for (k in seq_along(ch)) {
    hit <- grepl(pat, ch[[k]], ignore.case = TRUE) & is.na(out)
    out[hit] <- trimws(sub(pat, "", ch[[k]][hit], ignore.case = TRUE))
  }
  out
}

sm_field <- function(m, key) {
  i <- which(m$keys == key)
  if (!length(i)) return(character(0))
  m$vals[[i[1]]]
}

sm_expr <- function(ln) {
  b <- grep("^!series_matrix_table_begin", ln)[1]
  e <- grep("^!series_matrix_table_end", ln)[1]
  if (is.na(b)) stop("no expression table in this series matrix")
  if (is.na(e)) e <- length(ln) + 1
  df <- read.delim(text = paste(ln[(b + 1):(e - 1)], collapse = "\n"),
                   row.names = 1, check.names = FALSE)
  as.matrix(df)
}

human <- NULL
rat   <- NULL

# ===================================================== HUMAN: GSE244740 ======
hr("GSE244740  --  human iPSC-CM validation")

human <- tryCatch({
  suppressPackageStartupMessages(library(DESeq2))

  sm <- file.path(DATA, "GSE244740", "GSE244740_series_matrix.txt.gz")
  ct <- file.path(DATA, "GSE244740", "GSE244740_processed_data_counts.txt.gz")
  stopifnot(file.exists(sm), file.exists(ct))

  m    <- sm_meta(sm_lines(sm))
  desc <- sm_field(m, "Sample_description")
  trt  <- sm_char(m, "treatment")
  ttl  <- sm_field(m, "Sample_title")

  cat("Samples:", length(desc), " with a treatment label:", sum(!is.na(trt)), "\n")

  # what does clozapine look like in this series?
  cz <- grep("clozapine", trt, ignore.case = TRUE)
  cat("\nClozapine entries:", length(cz), "\n")
  if (!length(cz)) stop("no clozapine samples in GSE244740; inspect unique(trt)")
  cat("Distinct clozapine treatment labels:\n")
  for (s in sort(unique(trt[cz]))) cat("   ", s, "  n=", sum(trt == s, na.rm = TRUE), "\n", sep = "")

  # plate membership matters: compare within the plate clozapine sits on
  plate_of <- sub("^.*,\\s*(.*plate[^,]*),.*$", "\\1", ttl, ignore.case = TRUE)
  cat("\nPlates carrying clozapine:", paste(unique(plate_of[cz]), collapse = " | "), "\n")

  # vehicle / control labels present on those plates
  on_plate <- plate_of %in% unique(plate_of[cz])
  ctrl_lab <- unique(trt[on_plate & grepl("DMSO|vehicle|control|untreated",
                                          trt, ignore.case = TRUE)])
  cat("Control labels on those plates:\n")
  if (!length(ctrl_lab)) {
    cat("   NONE FOUND. All labels on these plates:\n")
    for (s in head(sort(unique(trt[on_plate])), 40)) cat("     ", s, "\n", sep = "")
    stop("could not identify vehicle wells; pick the label from the list above")
  }
  for (s in ctrl_lab) cat("   ", s, "  n=", sum(trt == s & on_plate, na.rm = TRUE), "\n", sep = "")

  # highest clozapine concentration: Con8 ... Con1, Con8 assumed highest
  cz_lab  <- sort(unique(trt[cz]))
  top_lab <- cz_lab[which.max(as.integer(sub("^.*Con", "", cz_lab)))]
  cat("\nUsing clozapine label:", top_lab, "\n")
  cat("NOTE: this assumes ConN numbering runs low to high. If Con1 is the\n",
      "highest dose in this series, change top_lab and rerun.\n", sep = "")

  keep_desc <- desc[(trt == top_lab | (trt %in% ctrl_lab & on_plate)) & !is.na(trt)]
  grp <- ifelse(trt[desc %in% keep_desc] == top_lab, "Clozapine", "VEH")

  cat("\nReading counts (large file, this takes a minute) ...\n")
  cts <- read.delim(gzfile(ct), row.names = 1, check.names = FALSE)
  cat("Counts matrix:", nrow(cts), "genes x", ncol(cts), "columns\n")

  present <- intersect(keep_desc, colnames(cts))
  cat("Wells matched to count columns:", length(present), "of", length(keep_desc), "\n")
  if (length(present) < 4) stop("too few matched wells; check description-to-column mapping")

  md <- data.frame(row.names = present,
                   grp = factor(ifelse(trt[match(present, desc)] == top_lab,
                                       "Clozapine", "VEH"),
                                levels = c("VEH", "Clozapine")))
  cat("Design:", sum(md$grp == "Clozapine"), "clozapine vs",
      sum(md$grp == "VEH"), "vehicle\n")

  cat("\nTarget genes present in the counts matrix:\n")
  for (g in names(ENS)) {
    hit <- grep(paste0("^", ENS[g]), rownames(cts))
    cat(sprintf("  %-10s %-16s rows=%d  total counts=%s\n", g, ENS[g], length(hit),
                if (length(hit)) format(sum(cts[hit, present, drop = FALSE])) else "-"))
  }

  dds <- DESeqDataSetFromMatrix(round(cts[, present, drop = FALSE]), md, ~ grp)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds, quiet = TRUE)
  r <- as.data.frame(results(dds, contrast = c("grp", "Clozapine", "VEH")))
  r$ens <- sub("\\..*$", "", rownames(r))

  out <- do.call(rbind, lapply(names(ENS), function(g) {
    i <- which(r$ens == ENS[g])
    if (!length(i)) return(NULL)
    data.frame(symbol = g, logFC = r$log2FoldChange[i[1]],
               adj.P.Val = r$padj[i[1]], dataset = "Independent hiPSC-CM")
  }))
  cat("\nHuman validation result:\n"); print(out, row.names = FALSE)
  write.csv(out, file.path(RES, "validation_human.csv"), row.names = FALSE)
  out
}, error = function(e) { cat("\nHUMAN SECTION FAILED: ", conditionMessage(e), "\n", sep = ""); NULL })

# ======================================================= RAT: GSE59905 =======
hr("GSE59905 / GPL5426  --  rat heart in vivo")

rat <- tryCatch({
  suppressPackageStartupMessages(library(limma))

  sm <- file.path(DATA, "GSE59905", "GSE59905-GPL5426_series_matrix.txt.gz")
  stopifnot(file.exists(sm))

  ln <- sm_lines(sm)
  m  <- sm_meta(ln)
  cmp <- sm_char(m, "compound"); dose <- sm_char(m, "dose")
  tim <- sm_char(m, "time");     veh  <- sm_char(m, "vehicle")
  cat("Samples:", length(cmp), "\n")

  is_cz <- !is.na(cmp) & grepl("clozapine", cmp, ignore.case = TRUE)
  cat("\nClozapine arrays:", sum(is_cz), "\n")
  if (!sum(is_cz)) stop("no clozapine arrays found")
  cat("  doses :", paste(unique(dose[is_cz]), collapse = ", "), "\n")
  cat("  times :", paste(unique(tim[is_cz]),  collapse = ", "), "\n")
  cat("  vehicle:", paste(unique(veh[is_cz]), collapse = ", "), "\n")

  cz_veh <- unique(veh[is_cz]); cz_tim <- unique(tim[is_cz])
  is_ctl <- !is.na(dose) & grepl("^0\\s*mg", dose) &
            !is.na(veh) & veh %in% cz_veh &
            !is.na(tim) & tim %in% cz_tim
  cat("Vehicle controls matched on vehicle and timepoint:", sum(is_ctl), "\n")
  if (sum(is_ctl) < 3) stop("too few matched controls")

  cat("\nReading expression table ...\n")
  ex <- sm_expr(ln)
  cat("Probes:", nrow(ex), " arrays:", ncol(ex), "\n")

  sel <- is_cz | is_ctl
  grp <- factor(ifelse(is_cz[sel], "Clozapine", "Control"),
                levels = c("Control", "Clozapine"))
  cat("Design:", sum(grp == "Clozapine"), "clozapine vs", sum(grp == "Control"), "control\n")

  fit <- eBayes(lmFit(ex[, sel, drop = FALSE], model.matrix(~ grp)))
  tt  <- topTable(fit, coef = 2, number = Inf)

  # probe -> symbol from the GPL5426 annotation
  ann <- file.path(DATA, "GPL5426.annot.gz")
  if (!file.exists(ann)) {
    cat("\nFetching GPL5426 annotation ...\n")
    download.file("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL5nnn/GPL5426/annot/GPL5426.annot.gz",
                  ann, mode = "wb", quiet = TRUE)
  }
  al <- readLines(gzfile(ann), warn = FALSE)
  b <- grep("^!platform_table_begin", al)[1]; e <- grep("^!platform_table_end", al)[1]
  if (is.na(e)) e <- length(al) + 1
  tab <- read.delim(text = paste(al[(b + 1):(e - 1)], collapse = "\n"),
                    check.names = FALSE, quote = "")
  sym_col <- grep("^Gene symbol$|^Symbol$|GENE_SYMBOL", names(tab),
                  ignore.case = TRUE, value = TRUE)[1]
  cat("Annotation rows:", nrow(tab), " symbol column:", sym_col, "\n")
  map <- setNames(toupper(tab[[sym_col]]), as.character(tab[[1]]))

  tt$symbol <- map[rownames(tt)]
  out <- tt[!is.na(tt$symbol) & tt$symbol %in% GENES, c("symbol", "logFC", "adj.P.Val")]
  out <- out[order(out$symbol, out$adj.P.Val), ]
  out <- out[!duplicated(out$symbol), ]
  out$dataset <- "Rat heart in vivo"
  cat("\nRat result:\n"); print(out, row.names = FALSE)
  write.csv(out, file.path(RES, "validation_rat.csv"), row.names = FALSE)
  out
}, error = function(e) { cat("\nRAT SECTION FAILED: ", conditionMessage(e), "\n", sep = ""); NULL })

# ============================================================ FOREST =========
hr("Figure 2")

cr <- read.csv(file.path(RES, "concentration_response.csv"))
disc <- cr[abs(cr$conc - 10) < 1e-6 & cr$gene %in% GENES, c("gene", "lfc", "padj")]
names(disc) <- c("symbol", "logFC", "adj.P.Val"); disc$dataset <- "Discovery hiPSC-CM"

d <- do.call(rbind, Filter(Negate(is.null), list(disc, human, rat)))
cat("Datasets in the plot:", paste(unique(d$dataset), collapse = " | "), "\n\n")
print(d, row.names = FALSE)
write.csv(d, file.path(RES, "figure2_data.csv"), row.names = FALSE)

cat("\nManuscript values for comparison:\n",
    "  TNFRSF12A  -2.14 / -1.66 / -1.48\n",
    "  EGLN3      -2.12 / -1.70 / -0.64\n",
    "  VEGFA      -1.75 / -0.46 n.s. / +0.08 n.s.\n", sep = "")

suppressPackageStartupMessages(library(ggplot2))
d$symbol  <- factor(d$symbol, levels = rev(GENES))
d$dataset <- factor(d$dataset, levels = c("Discovery hiPSC-CM",
                                          "Independent hiPSC-CM", "Rat heart in vivo"))
d$sig <- !is.na(d$adj.P.Val) & d$adj.P.Val < 0.05

fig2 <- ggplot(d, aes(logFC, symbol, colour = dataset, shape = sig)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_point(size = 3, position = position_dodge(width = 0.6), fill = "white", stroke = 0.7) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c("Not significant", "Adjusted p < 0.05"), name = NULL) +
  scale_colour_brewer(palette = "Dark2", name = NULL, drop = FALSE) +
  labs(x = expression(log[2]~"fold-change"), y = NULL,
       caption = "Filled points meet adjusted p < 0.05. Rat values are limma moderated estimates.") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "bold.italic"),
        legend.position = "bottom", legend.box = "vertical",
        plot.caption = element_text(hjust = 0.5, size = 8.5, colour = "grey30"))

ggsave(file.path(FIGS, "Fig2_forest.png"), fig2, width = 7.5, height = 4.2, dpi = 600)
ggsave(file.path(FIGS, "Fig2_forest.pdf"), fig2, width = 7.5, height = 4.2)
writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo_figure2.txt"))
cat("\nWrote figures/Fig2_forest.png and results/figure2_data.csv\n")
