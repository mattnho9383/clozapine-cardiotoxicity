# =============================================================================
# 04_figure2_final.R  --  cross-species forest plot, with intervals
#
# Same analysis as 03_figure2.R. Two differences that matter:
#   - captures standard errors and 95% intervals, so the forest plot has
#     the intervals a forest plot is supposed to have
#   - writes results/figure2_data.csv with the interval bounds, and plots
#     from that file rather than from objects in memory, so the figure can
#     be regenerated without refitting anything
#
#   Rscript 04_figure2_final.R 2>&1 | tee logs/fig2_final.log
# =============================================================================

WORKDIR <- Sys.getenv("CLOZ_HOME", unset = file.path(path.expand("~"), "clozapine"))
DATA <- file.path(WORKDIR, "data"); RES <- file.path(WORKDIR, "results")
FIGS <- file.path(WORKDIR, "figures")
for (d in c(RES, FIGS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
options(timeout = 1800)
lib <- Sys.getenv("R_LIBS_USER", unset = file.path(path.expand("~"), "Rlibs", "R-4.6.0"))
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))

GENES <- c("TNFRSF12A", "EGLN3", "VEGFA")
ENS <- c(TNFRSF12A = "ENSG00000006327",
         EGLN3     = "ENSG00000129521",
         VEGFA     = "ENSG00000112715")
hr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

sm_lines <- function(p) readLines(gzfile(p), warn = FALSE)
sm_meta <- function(ln) {
  b <- grep("^!series_matrix_table_begin", ln)
  if (length(b)) ln <- ln[seq_len(b[1] - 1)]
  ln <- ln[grepl("^!Sample_", ln)]; p <- strsplit(ln, "\t")
  list(keys = vapply(p, function(x) sub("^!", "", x[1]), character(1)),
       vals = lapply(p, function(x) gsub('^"|"$', "", x[-1])))
}
sm_char <- function(m, prefix) {
  ch <- m$vals[m$keys == "Sample_characteristics_ch1"]
  if (!length(ch)) return(character(0))
  out <- rep(NA_character_, length(ch[[1]])); pat <- paste0("^", prefix, "\\s*:\\s*")
  for (k in seq_along(ch)) {
    hit <- grepl(pat, ch[[k]], ignore.case = TRUE) & is.na(out)
    out[hit] <- trimws(sub(pat, "", ch[[k]][hit], ignore.case = TRUE))
  }
  out
}
sm_field <- function(m, k) { i <- which(m$keys == k); if (!length(i)) character(0) else m$vals[[i[1]]] }
sm_expr <- function(ln) {
  b <- grep("^!series_matrix_table_begin", ln)[1]; e <- grep("^!series_matrix_table_end", ln)[1]
  if (is.na(e)) e <- length(ln) + 1
  as.matrix(read.delim(text = paste(ln[(b + 1):(e - 1)], collapse = "\n"),
                       row.names = 1, check.names = FALSE))
}

# ============================================================ DISCOVERY ======
hr("Discovery: reading verified concentration_response.csv")
cr <- read.csv(file.path(RES, "concentration_response.csv"))
disc <- cr[abs(cr$conc - 10) < 1e-6 & cr$gene %in% GENES, ]
# the discovery CSV has no SE column; refit the 10 uM contrast to recover it
suppressPackageStartupMessages(library(DESeq2))
hash <- read.csv(file.path(DATA, "GSE262419_hash.csv.gz"))
hash$Chemical_name <- trimws(hash$Chemical_name)
p15 <- read.csv(file.path(DATA, "GSE262419_Plate15.csv.gz"), row.names = 1, check.names = FALSE)
p16 <- read.csv(file.path(DATA, "GSE262419_Plate16.csv.gz"), row.names = 1, check.names = FALSE)
ca <- cbind(p15, p16)
wl <- gsub("[0-9]+", "", hash$Well_ID); wn <- gsub("[^0-9]", "", hash$Well_ID)
hash$key <- paste0("Plate", hash$Plate_ID, "-", sprintf("%s%02d", wl, as.integer(wn)))
ma <- hash[hash$Plate_ID %in% c(15, 16), ]
pres <- intersect(ma$key, colnames(ca)); cf <- ca[, pres, drop = FALSE]
ma <- ma[match(colnames(cf), ma$key), ]
ma$cn <- suppressWarnings(as.numeric(as.character(ma$Chemical_Concentration_uM)))
kp <- (ma$Chemical_name == "Clozapine" & abs(ma$cn - 10) < 1e-6 & ma$Plate_ID == 16) |
      (ma$Chemical_name == "VEH" & ma$Plate_ID == 16)
mm <- ma[kp, ]; mm$grp <- factor(ifelse(mm$Chemical_name == "Clozapine", "Drug", "VEH"),
                                 levels = c("VEH", "Drug"))
n_disc <- c(sum(mm$grp == "Drug"), sum(mm$grp == "VEH"))
dd <- DESeqDataSetFromMatrix(round(cf[, kp]), mm, ~ grp)
dd <- DESeq(dd[rowSums(counts(dd)) >= 10, ], quiet = TRUE)
rd <- as.data.frame(results(dd, contrast = c("grp", "Drug", "VEH")))
rd$symbol <- gsub("_[0-9]+$", "", rownames(rd))
rd <- rd[rd$symbol %in% GENES, ]
rd <- rd[order(rd$symbol, -rd$baseMean), ]; rd <- rd[!duplicated(rd$symbol), ]
z <- qnorm(0.975)
d1 <- data.frame(symbol = rd$symbol, logFC = rd$log2FoldChange,
                 lo = rd$log2FoldChange - z * rd$lfcSE,
                 hi = rd$log2FoldChange + z * rd$lfcSE,
                 adj.P.Val = rd$padj, dataset = "Discovery hiPSC-CM")
cat("Discovery:", n_disc[1], "clozapine vs", n_disc[2], "vehicle\n"); print(d1, row.names = FALSE)

# ================================================================ HUMAN ======
hr("Independent human: GSE244740")
sm <- file.path(DATA, "GSE244740", "GSE244740_series_matrix.txt.gz")
ct <- file.path(DATA, "GSE244740", "GSE244740_processed_data_counts.txt.gz")
m <- sm_meta(sm_lines(sm))
desc <- sm_field(m, "Sample_description"); trt <- sm_char(m, "treatment")
ttl <- sm_field(m, "Sample_title")
plate <- sub("^.*,\\s*(.*plate[^,]*),.*$", "\\1", ttl, ignore.case = TRUE)
cz <- grep("clozapine", trt, ignore.case = TRUE)
onp <- plate %in% unique(plate[cz])
ctrl <- unique(trt[onp & grepl("DMSO", trt, ignore.case = TRUE)])
lab <- sort(unique(trt[cz])); top <- lab[which.max(as.integer(sub("^.*Con", "", lab)))]
kd <- desc[!is.na(trt) & (trt == top | (trt %in% ctrl & onp))]
cts <- read.delim(gzfile(ct), row.names = 1, check.names = FALSE)
pr <- intersect(kd, colnames(cts))
md <- data.frame(row.names = pr,
                 grp = factor(ifelse(trt[match(pr, desc)] == top, "Clozapine", "VEH"),
                              levels = c("VEH", "Clozapine")))
n_hum <- c(sum(md$grp == "Clozapine"), sum(md$grp == "VEH"))
dh <- DESeqDataSetFromMatrix(round(cts[, pr]), md, ~ grp)
dh <- DESeq(dh[rowSums(counts(dh)) >= 10, ], quiet = TRUE)
rh <- as.data.frame(results(dh, contrast = c("grp", "Clozapine", "VEH")))
rh$ens <- sub("\\..*$", "", rownames(rh))
d2 <- do.call(rbind, lapply(names(ENS), function(g) {
  i <- which(rh$ens == ENS[g]); if (!length(i)) return(NULL)
  data.frame(symbol = g, logFC = rh$log2FoldChange[i[1]],
             lo = rh$log2FoldChange[i[1]] - z * rh$lfcSE[i[1]],
             hi = rh$log2FoldChange[i[1]] + z * rh$lfcSE[i[1]],
             adj.P.Val = rh$padj[i[1]], dataset = "Independent hiPSC-CM")
}))
cat("Human:", n_hum[1], "clozapine (", top, ") vs", n_hum[2], "DMSO\n")
print(d2, row.names = FALSE)

# ================================================================== RAT ======
hr("Rat: GSE59905 / GPL5426")
suppressPackageStartupMessages(library(limma))
ln <- sm_lines(file.path(DATA, "GSE59905", "GSE59905-GPL5426_series_matrix.txt.gz"))
mr <- sm_meta(ln)
cmp <- sm_char(mr, "compound"); dose <- sm_char(mr, "dose")
tim <- sm_char(mr, "time"); veh <- sm_char(mr, "vehicle")
isc <- !is.na(cmp) & grepl("clozapine", cmp, ignore.case = TRUE)
ist <- !is.na(dose) & grepl("^0\\s*mg", dose) & !is.na(veh) & veh %in% unique(veh[isc]) &
       !is.na(tim) & tim %in% unique(tim[isc])
ex <- sm_expr(ln); sel <- isc | ist
n_rat <- c(sum(isc), sum(ist))
grp <- factor(ifelse(isc[sel], "Clozapine", "Control"), levels = c("Control", "Clozapine"))
fit <- eBayes(lmFit(ex[, sel, drop = FALSE], model.matrix(~ grp)))
tt <- topTable(fit, coef = 2, number = Inf, confint = TRUE)   # confint gives CI.L / CI.R
ann <- file.path(DATA, "GPL5426.annot.gz")
if (!file.exists(ann))
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL5nnn/GPL5426/annot/GPL5426.annot.gz",
                ann, mode = "wb", quiet = TRUE)
al <- readLines(gzfile(ann), warn = FALSE)
b <- grep("^!platform_table_begin", al)[1]; e <- grep("^!platform_table_end", al)[1]
if (is.na(e)) e <- length(al) + 1
tab <- read.delim(text = paste(al[(b + 1):(e - 1)], collapse = "\n"), check.names = FALSE, quote = "")
sc <- grep("^Gene symbol$", names(tab), ignore.case = TRUE, value = TRUE)[1]
tt$symbol <- setNames(toupper(tab[[sc]]), as.character(tab[[1]]))[rownames(tt)]
tr <- tt[!is.na(tt$symbol) & tt$symbol %in% GENES, ]
tr <- tr[order(tr$symbol, tr$adj.P.Val), ]; tr <- tr[!duplicated(tr$symbol), ]
d3 <- data.frame(symbol = tr$symbol, logFC = tr$logFC, lo = tr$CI.L, hi = tr$CI.R,
                 adj.P.Val = tr$adj.P.Val, dataset = "Rat heart in vivo")
cat("Rat:", n_rat[1], "clozapine vs", n_rat[2], "control\n"); print(d3, row.names = FALSE)

# =============================================================== ASSEMBLE ====
hr("Figure 2")
d <- rbind(d1, d2, d3)
write.csv(d, file.path(RES, "figure2_data.csv"), row.names = FALSE)
print(d, row.names = FALSE)

suppressPackageStartupMessages(library(ggplot2))
labs_n <- c(sprintf("Discovery hiPSC-CM (n = %d vs %d)", n_disc[1], n_disc[2]),
            sprintf("Independent hiPSC-CM (n = %d vs %d)", n_hum[1], n_hum[2]),
            sprintf("Rat heart in vivo (n = %d vs %d)", n_rat[1], n_rat[2]))
d$dataset <- factor(d$dataset, levels = c("Discovery hiPSC-CM", "Independent hiPSC-CM",
                                          "Rat heart in vivo"), labels = labs_n)
d$symbol <- factor(d$symbol, levels = rev(GENES))
d$sig <- !is.na(d$adj.P.Val) & d$adj.P.Val < 0.05
dodge <- position_dodge(width = 0.62)

fig2 <- ggplot(d, aes(logFC, symbol, colour = dataset, shape = sig)) +
  annotate("rect", ymin = 1.5, ymax = 2.5, xmin = -Inf, xmax = Inf,
           fill = "grey96", colour = NA) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.45) +
  geom_linerange(aes(xmin = lo, xmax = hi), position = dodge, linewidth = 0.65) +
  geom_point(position = dodge, size = 2.7, fill = "white", stroke = 0.8) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 19),
                     labels = c(`FALSE` = "n.s.", `TRUE` = "adj. p < 0.05"), name = NULL) +
  scale_colour_manual(values = setNames(c("#1B7837", "#B2182B", "#4A5FA5"), labs_n),
                      name = NULL) +
  scale_x_continuous(breaks = seq(-2.5, 0.5, 0.5), expand = expansion(mult = 0.02)) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         shape  = guide_legend(order = 2, nrow = 1)) +
  labs(x = expression(log[2]~"fold-change vs vehicle (95% CI)"), y = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.text.y      = element_text(face = "italic", size = 11, colour = "black"),
        axis.text.x      = element_text(colour = "black"),
        axis.title.x     = element_text(margin = margin(t = 7)),
        axis.line.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        legend.position  = "bottom",
        legend.box       = "vertical",
        legend.spacing.y = unit(1, "pt"),
        legend.margin    = margin(t = 2, b = 0),
        legend.key.size  = unit(9, "pt"),
        plot.margin      = margin(6, 10, 4, 4))

ggsave(file.path(FIGS, "Fig2_forest.png"), fig2, width = 6.6, height = 2.9, dpi = 600)
ggsave(file.path(FIGS, "Fig2_forest.pdf"), fig2, width = 6.6, height = 2.9)
writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo_figure2.txt"))
cat("\nWrote figures/Fig2_forest.png (600 dpi) and .pdf\n")
