# =============================================================================
# clozapine_figure1.R  --  everything in one file
#
# Installs what it needs, downloads the data from GEO, runs the clozapine
# concentration-response, writes the values to CSV and saves Figure 1.
#
# Mac / RStudio:            open this file, click Source.
# Great Lakes OnDemand:     start an RStudio session, open this file, Source.
# Terminal either place:    Rscript clozapine_figure1.R
#
# Nothing to edit. First run takes a while because DESeq2 compiles and the
# plate files download. Later runs skip both.
# =============================================================================

WORKDIR <- getwd()
DATA    <- file.path(WORKDIR, "data")
RES     <- file.path(WORKDIR, "results")
FIGS    <- file.path(WORKDIR, "figures")
for (d in c(WORKDIR, DATA, RES, FIGS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

cat("Working in:", WORKDIR, "\n\n")

# ------------------------------------------------------------- 1. packages --
need <- function(pkg, bioc = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  cat("Installing", pkg, "...\n")
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
need("ggplot2")
need("DESeq2", bioc = TRUE)

suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })

# ------------------------------------------------------------- 2. get data --
BASE  <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262419/suppl"
FILES <- c("GSE262419_hash.csv.gz",
           "GSE262419_Plate15.csv.gz",
           "GSE262419_Plate16.csv.gz")

for (f in FILES) {
  dest <- file.path(DATA, f)
  if (file.exists(dest) && file.size(dest) > 1000) {
    cat("have", f, "\n"); next
  }
  cat("downloading", f, "...\n")
  ok <- tryCatch({
    download.file(file.path(BASE, f), dest, mode = "wb", quiet = TRUE); TRUE
  }, error = function(e) FALSE)

  if (!ok || !file.exists(dest) || file.size(dest) < 1000) {
    # filenames on GEO occasionally differ from what is expected; show what
    # is actually there rather than failing silently
    cat("\nCould not fetch", f, "\nFiles actually present in that GEO directory:\n")
    listing <- tryCatch(readLines(paste0(BASE, "/")), error = function(e) character(0))
    hits <- unique(unlist(regmatches(listing,
              gregexpr("GSE262419[A-Za-z0-9_.\\-]*", listing))))
    if (length(hits)) print(hits) else
      cat("  (could not read the directory; open ",
          "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262419 ",
          "and download the hash file plus plates 15 and 16 by hand into\n  ",
          DATA, "\n", sep = "")
    stop("stopping: data not available")
  }
}
cat("\n")

# --------------------------------------------------------------- 3. load ----
hash <- read.csv(file.path(DATA, "GSE262419_hash.csv.gz"))
hash$Chemical_name <- trimws(hash$Chemical_name)

p15 <- read.csv(file.path(DATA, "GSE262419_Plate15.csv.gz"), row.names = 1, check.names = FALSE)
p16 <- read.csv(file.path(DATA, "GSE262419_Plate16.csv.gz"), row.names = 1, check.names = FALSE)
counts_all <- cbind(p15, p16)

wl <- gsub("[0-9]+", "", hash$Well_ID)
wn <- gsub("[^0-9]", "", hash$Well_ID)
hash$match_key <- paste0("Plate", hash$Plate_ID, "-",
                         sprintf("%s%02d", wl, as.integer(wn)))

meta_all <- hash[hash$Plate_ID %in% c(15, 16), ]
present  <- intersect(meta_all$match_key, colnames(counts_all))
counts_full <- counts_all[, present, drop = FALSE]
meta_all <- meta_all[match(colnames(counts_full), meta_all$match_key), ]

# concentration is sometimes stored as text; coerce once and match on tolerance
meta_all$conc_num <- suppressWarnings(
  as.numeric(as.character(meta_all$Chemical_Concentration_uM)))

cat("Loaded", nrow(counts_full), "probes x", ncol(counts_full), "wells\n\n")

GENES <- c("TNFRSF12A","EGLN3","VEGFA","SH3D21","PRSS45P",
           "HIST1H1B","HIST1H4I","HIST2H4B","NUSAP1")
UP    <- c("HIST1H1B","HIST1H4I","HIST2H4B","NUSAP1")
CONCS <- c(0.2, 1, 10)
MIN_EXPR <- 100

# ---------------------------------------------------------- 4. diagnostic ---
# Read this before trusting anything downstream. Exact float comparison on the
# concentration column silently returns zero wells, which is the most likely
# reason the earlier dose loop came back empty at the two lower doses.
cat("Clozapine concentrations present in the metadata:\n")
print(sort(unique(meta_all$conc_num[meta_all$Chemical_name == "Clozapine"])))
cat("\nWells matched per requested concentration:\n")
for (cc in CONCS) {
  n <- sum(meta_all$Chemical_name == "Clozapine" & abs(meta_all$conc_num - cc) < 1e-6)
  cat(sprintf("  %6.2f uM : %2d clozapine wells%s\n", cc, n,
              if (n == 0) "   <-- nothing matched" else ""))
}
cat("\n")

# ------------------------------------------------------------ 5. DE helper --
deseq_drug <- function(drug, conc = 10) {
  plate <- unique(meta_all$Plate_ID[meta_all$Chemical_name == drug])[1]
  keep <- (meta_all$Chemical_name == drug &
             abs(meta_all$conc_num - conc) < 1e-6 &
             meta_all$Plate_ID == plate) |
          (meta_all$Chemical_name == "VEH" & meta_all$Plate_ID == plate)
  mm <- meta_all[keep, ]
  n_drug <- sum(mm$Chemical_name == drug)
  if (n_drug < 2) stop(sprintf("only %d wells for %s at %.2f uM", n_drug, drug, conc))
  ct <- round(counts_full[, keep, drop = FALSE])
  mm$grp <- factor(ifelse(mm$Chemical_name == drug, "Drug", "VEH"),
                   levels = c("VEH", "Drug"))
  cat(sprintf("  %s %5.2f uM: %2d treated vs %2d vehicle (plate %s)\n",
              drug, conc, n_drug, sum(mm$grp == "VEH"), plate))
  d <- DESeqDataSetFromMatrix(ct, mm, ~ grp)
  d <- d[rowSums(counts(d)) >= 10, ]
  d <- DESeq(d, quiet = TRUE)
  r <- as.data.frame(results(d, contrast = c("grp", "Drug", "VEH")))
  r$symbol <- gsub("_[0-9]+$", "", rownames(r))
  r
}

# --------------------------------------------------------- 6. dose loop -----
cat("Running dose loop:\n")
dose_tab <- do.call(rbind, lapply(CONCS, function(cc) {
  r <- deseq_drug("Clozapine", conc = cc)
  # no expression filter here: that rule decides DEG status, not extraction
  sub <- r[r$symbol %in% GENES, c("symbol","log2FoldChange","padj","baseMean")]
  # TempO-Seq probes carry numeric suffixes; keep the best-expressed per symbol
  sub <- sub[order(sub$symbol, -sub$baseMean), ]
  sub <- sub[!duplicated(sub$symbol), ]
  data.frame(conc = cc, gene = sub$symbol, lfc = sub$log2FoldChange,
             padj = sub$padj, baseMean = sub$baseMean, stringsAsFactors = FALSE)
}))

grid <- expand.grid(gene = GENES, conc = CONCS, stringsAsFactors = FALSE)
d <- merge(grid, dose_tab, by = c("gene","conc"), all.x = TRUE)
d$tested <- !is.na(d$padj)                                   # NA != null
d$sig    <- d$tested & d$padj < 0.05 & d$baseMean >= MIN_EXPR
d <- d[order(d$gene, d$conc), ]

cat("\n"); print(d, row.names = FALSE)
cat("\nCells with no result at all :", sum(is.na(d$lfc)),
    "\nCells DESeq2 filtered (NA)  :", sum(!d$tested), "\n\n")
write.csv(d, file.path(RES, "concentration_response.csv"), row.names = FALSE)

# ------------------------------------- 7. check against the manuscript ------
manuscript_10uM <- c(TNFRSF12A = -2.14, EGLN3 = -2.12, VEGFA = -1.75,
                     SH3D21 = -1.44, PRSS45P = -1.39, HIST1H1B = 1.66,
                     HIST1H4I = 1.25, HIST2H4B = 1.24, NUSAP1 = 1.13)
chk <- d[abs(d$conc - 10) < 1e-6, c("gene","lfc")]
chk$manuscript <- manuscript_10uM[as.character(chk$gene)]
chk$diff <- round(chk$lfc - chk$manuscript, 3)
cat("10 uM values vs those quoted in Section 3.5:\n")
print(chk, row.names = FALSE)
cat("\nTNFRSF12A, EGLN3 and VEGFA were the confirmed three and should match.\n",
    "Large differences in the other six mean the prose describing the lower\n",
    "concentrations needs updating too, not just the figure.\n\n", sep = "")

# ---------------------------------------------------------------- 8. plot ---
d$dir  <- ifelse(d$gene %in% UP, "Upregulated at 10 \u00b5M", "Downregulated at 10 \u00b5M")
d$gene <- factor(d$gene, levels = GENES)
rng  <- range(d$lfc, na.rm = TRUE)
ylim <- c(min(-2.6, floor(rng[1])), max(2.3, ceiling(rng[2])))

fig1 <- ggplot(d[!is.na(d$lfc), ], aes(conc, lfc, group = gene, colour = dir)) +
  annotate("rect", xmin = 1.05, xmax = 1.80, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.55) +
  annotate("text", x = 1.37, y = ylim[2] - 0.25, label = "therapeutic\nrange",
           size = 2.7, colour = "grey35", lineheight = 0.9) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  geom_point(aes(shape = sig), size = 2.4, fill = "white", stroke = 0.7) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c("Not significant", "Adjusted p < 0.05"), name = NULL) +
  scale_colour_manual(values = c("Downregulated at 10 \u00b5M" = "#B2182B",
                                 "Upregulated at 10 \u00b5M"   = "#2166AC"), name = NULL) +
  scale_x_log10(breaks = CONCS, labels = c("0.2","1","10")) +
  scale_y_continuous(limits = ylim, breaks = seq(-4, 4, 1)) +
  facet_wrap(~ gene, ncol = 3) +
  labs(x = expression("Clozapine concentration ("*mu*"M, log scale)"),
       y = expression(log[2]~"fold-change vs vehicle"),
       caption = paste(
         "Shaded band shows the reported therapeutic trough range (1.05-1.80 \u00b5M).",
         "Filled points meet adjusted p < 0.05 and mean normalized expression >= 100.",
         sep = "\n")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold.italic", size = 9.5),
        strip.background = element_rect(fill = "grey93"),
        axis.title       = element_text(face = "bold"),
        legend.position  = "bottom", legend.box = "vertical",
        legend.margin    = margin(t = -4),
        plot.caption     = element_text(hjust = 0.5, size = 8.5, colour = "grey30"))

ggsave(file.path(FIGS, "Fig1_concentration.png"), fig1, width = 8.5, height = 7, dpi = 600)
ggsave(file.path(FIGS, "Fig1_concentration.pdf"), fig1, width = 8.5, height = 7)

cat("Wrote:\n  ", file.path(FIGS, "Fig1_concentration.png"),
    "\n  ", file.path(RES,  "concentration_response.csv"), "\n", sep = "")
