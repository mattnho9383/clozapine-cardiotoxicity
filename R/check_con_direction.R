#!/usr/bin/env Rscript
# Which end of the Con1-Con8 series is the high dose?
# A dose-response gives far more DEGs at the top concentration.
suppressPackageStartupMessages({library(DESeq2)})
sm <- "data/GSE244740/GSE244740_series_matrix.txt.gz"
ct <- "data/GSE244740/GSE244740_processed_data_counts.txt.gz"

ln  <- readLines(gzfile(sm), warn = FALSE)
tit <- strsplit(sub("^.Sample_title\t", "", grep("^.Sample_title", ln, value = TRUE)), "\t")[[1]]
tit <- gsub('"', '', tit)
ids <- sub(".*\\[(.*)\\].*", "\\1", tit)
cnt <- read.delim(gzfile(ct), row.names = 1, check.names = FALSE)

veh <- ids[grepl("DMSO, Con0, Screening plate 2", tit, fixed = TRUE)]
veh <- veh[veh %in% colnames(cnt)]

run <- function(lab) {
  drug <- ids[grepl(paste0("Clozapine, ", lab, ","), tit, fixed = TRUE)]
  drug <- drug[drug %in% colnames(cnt)]
  use  <- c(drug, veh)
  m <- data.frame(row.names = use,
                  grp = factor(c(rep("Drug", length(drug)), rep("VEH", length(veh))),
                               levels = c("VEH", "Drug")))
  d <- DESeqDataSetFromMatrix(round(cnt[, use, drop = FALSE]), m, ~ grp)
  d <- DESeq(d[rowSums(counts(d)) >= 10, ], quiet = TRUE)
  r <- results(d, contrast = c("grp", "Drug", "VEH"))
  n <- sum(!is.na(r$padj) & r$padj < 0.05 & r$baseMean >= 100)
  tn <- r[grep("^TNFRSF12A", rownames(r)), ]
  cat(sprintf("%-5s  n=%d vs %d veh | DEGs %5d | TNFRSF12A lfc %7.3f  padj %.3g\n",
              lab, length(drug), length(veh), n,
              if (nrow(tn)) tn$log2FoldChange[1] else NA,
              if (nrow(tn)) tn$padj[1] else NA))
}
for (l in paste0("Con", 1:8)) run(l)
