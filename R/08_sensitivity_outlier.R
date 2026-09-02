# =============================================================================
# 08_sensitivity_outlier.R
#
# Plate16-B02 is a vehicle well sitting at PC2 = 30 while every other well on
# that plate falls between -5 and +8. It is in the vehicle group for both the
# clozapine and risperidone contrasts, so an aberrant control could inflate
# within-group variance and cost power. A reviewer can reasonably ask whether
# the clozapine null is partly an artifact of a bad control.
#
# This refits both contrasts with and without that well and reports whether
# the DEG list changes. Either answer is useful: unchanged gives you a
# sensitivity analysis to cite, changed tells you now rather than at review.
#
#   Rscript 08_sensitivity_outlier.R 2>&1 | tee logs/sensitivity.log
# =============================================================================

WORKDIR <- getwd()
DATA <- file.path(WORKDIR, "data"); RES <- file.path(WORKDIR, "results")
lib <- Sys.getenv("R_LIBS_USER", unset = file.path(path.expand("~"), "Rlibs", "R-4.6.0"))
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(DESeq2))

OUTLIER  <- "Plate16-B02"
MIN_EXPR <- 100
hr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

hash <- read.csv(file.path(DATA, "GSE262419_hash.csv.gz"))
hash$Chemical_name <- trimws(hash$Chemical_name)
p15 <- read.csv(file.path(DATA, "GSE262419_Plate15.csv.gz"), row.names = 1, check.names = FALSE)
p16 <- read.csv(file.path(DATA, "GSE262419_Plate16.csv.gz"), row.names = 1, check.names = FALSE)
ca <- cbind(p15, p16)
wl <- gsub("[0-9]+", "", hash$Well_ID); wn <- gsub("[^0-9]", "", hash$Well_ID)
hash$key <- paste0("Plate", hash$Plate_ID, "-", sprintf("%s%02d", wl, as.integer(wn)))
meta <- hash[hash$Plate_ID %in% c(15, 16), ]
pres <- intersect(meta$key, colnames(ca)); cf <- ca[, pres, drop = FALSE]
meta <- meta[match(colnames(cf), meta$key), ]
meta$cn <- suppressWarnings(as.numeric(as.character(meta$Chemical_Concentration_uM)))

# =================================================================== QC ======
hr(paste("Why is", OUTLIER, "aberrant?"))
v16 <- meta$key[meta$Plate_ID == 16 & meta$Chemical_name == "VEH"]
qc <- data.frame(
  well     = v16,
  lib_size = colSums(cf[, v16, drop = FALSE]),
  detected = colSums(cf[, v16, drop = FALSE] > 0),
  max_probe_share = round(apply(cf[, v16, drop = FALSE], 2,
                                function(x) max(x) / sum(x)), 4))
qc$flag <- ifelse(qc$well == OUTLIER, "  <-- outlier", "")
print(qc[order(-qc$lib_size), ], row.names = FALSE)
cat("\nMedian library size of the other vehicles: ",
    format(median(qc$lib_size[qc$well != OUTLIER]), big.mark = ","), "\n",
    OUTLIER, " library size: ", format(qc$lib_size[qc$well == OUTLIER], big.mark = ","), "\n",
    "\nA much smaller library or fewer detected probes means a technical failure,\n",
    "which is the cleanest justification for excluding it. Similar depth with a\n",
    "different expression profile is a weaker case and should stay a sensitivity\n",
    "analysis rather than an exclusion.\n", sep = "")

# ============================================================== REFIT ========
fit_drug <- function(drug, drop = NULL) {
  keep <- (meta$Chemical_name == drug & abs(meta$cn - 10) < 1e-6 & meta$Plate_ID == 16) |
          (meta$Chemical_name == "VEH" & meta$Plate_ID == 16)
  if (!is.null(drop)) keep <- keep & !(meta$key %in% drop)
  mm <- meta[keep, ]
  mm$grp <- factor(ifelse(mm$Chemical_name == drug, "Drug", "VEH"), levels = c("VEH", "Drug"))
  d <- DESeqDataSetFromMatrix(round(cf[, keep, drop = FALSE]), mm, ~ grp)
  d <- DESeq(d[rowSums(counts(d)) >= 10, ], quiet = TRUE)
  r <- as.data.frame(results(d, contrast = c("grp", "Drug", "VEH")))
  r$symbol <- gsub("_[0-9]+$", "", rownames(r))
  # TempO-Seq probes carry numeric suffixes; keep the best-expressed
  # probe per symbol, matching 01_figure1_concentration.R
  r <- r[order(r$symbol, -r$baseMean), ]
  r <- r[!duplicated(r$symbol), ]
  attr(r, "n") <- c(sum(mm$grp == "Drug"), sum(mm$grp == "VEH"))
  r
}
degs <- function(r) sort(unique(r$symbol[!is.na(r$padj) & r$padj < 0.05 &
                                         !is.na(r$baseMean) & r$baseMean >= MIN_EXPR]))

for (drug in c("Clozapine", "Risperidone")) {
  hr(paste(drug, "with and without", OUTLIER))
  a <- fit_drug(drug)
  b <- fit_drug(drug, drop = OUTLIER)
  ga <- degs(a); gb <- degs(b)
  cat(sprintf("all vehicles   : %2d vs %2d wells | %d DEGs\n", attr(a,"n")[1], attr(a,"n")[2], length(ga)))
  cat(sprintf("outlier dropped: %2d vs %2d wells | %d DEGs\n", attr(b,"n")[1], attr(b,"n")[2], length(gb)))

  if (identical(ga, gb)) {
    cat("\nIdentical gene list. The result does not depend on that well.\n")
  } else {
    cat("\nLost when dropped : ", paste(setdiff(ga, gb), collapse = ", "), "\n",
        "Gained when dropped: ", paste(setdiff(gb, ga), collapse = ", "), "\n", sep = "")
  }

  if (length(ga) || length(gb)) {
    u <- union(ga, gb)
    cmp <- data.frame(
      gene   = u,
      lfc_all  = round(a$log2FoldChange[match(u, a$symbol)], 3),
      padj_all = signif(a$padj[match(u, a$symbol)], 3),
      lfc_drop  = round(b$log2FoldChange[match(u, b$symbol)], 3),
      padj_drop = signif(b$padj[match(u, b$symbol)], 3))
    cat("\n"); print(cmp, row.names = FALSE)
    write.csv(cmp, file.path(RES, paste0("sensitivity_", tolower(drug), ".csv")),
              row.names = FALSE)
  }
}

hr("Interpretation")
cat("If clozapine still gives the same nine genes, report this in the Methods or\n",
    "Limitations as a sensitivity analysis: the result is robust to exclusion of\n",
    "an atypical vehicle well. That converts a reviewer question into a strength.\n\n",
    "If the list changes, decide on the QC evidence above whether the well is a\n",
    "technical failure. Exclude it only with a stated reason, and apply the same\n",
    "rule to every contrast, not just this one.\n", sep = "")
