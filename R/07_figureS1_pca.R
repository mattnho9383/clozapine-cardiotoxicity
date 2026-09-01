# =============================================================================
# 07_figureS1_pca.R  --  Figure S1, PCA computed within each plate
#
# Replaces the version in 06, which ran a single PCA across both plates and
# then facetted it. That made PC1 largely a plate effect, so the percentages
# in the axis labels described between-plate variance rather than drug
# response, and haloperidol appeared not to separate despite 95 DEGs.
#
# Here each plate gets its own PCA, matching the plate-matched design of the
# DESeq2 contrasts. Percentages appear in the panel headings because they
# now differ between plates.
#
#   Rscript 07_figureS1_pca.R 2>&1 | tee logs/figS1.log
# =============================================================================

WORKDIR <- Sys.getenv("CLOZ_HOME", unset = file.path(path.expand("~"), "clozapine"))
DATA <- file.path(WORKDIR, "data"); RES <- file.path(WORKDIR, "results")
FIGS <- file.path(WORKDIR, "figures")
lib <- Sys.getenv("R_LIBS_USER", unset = file.path(path.expand("~"), "Rlibs", "R-4.6.0"))
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })

NTOP <- 500   # top variable probes, matching plotPCA's default; stated in the caption

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

DRUGS <- c("Haloperidol", "Clozapine", "Risperidone")

pca_plate <- function(pl) {
  sel <- meta$Plate_ID == pl & meta$Chemical_name %in% c("VEH", DRUGS) &
         (meta$Chemical_name == "VEH" | abs(meta$cn - 10) < 1e-6)
  if (sum(sel) < 4) return(NULL)
  mm <- meta[sel, ]
  dds <- DESeqDataSetFromMatrix(round(cf[, sel, drop = FALSE]),
                                data.frame(row.names = mm$key, drug = mm$Chemical_name), ~ 1)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  mat <- assay(vst(dds, blind = TRUE))
  rv  <- apply(mat, 1, var)
  keep <- order(rv, decreasing = TRUE)[seq_len(min(NTOP, length(rv)))]
  p <- prcomp(t(mat[keep, , drop = FALSE]))
  pv <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)

  df <- data.frame(PC1 = p$x[, 1], PC2 = p$x[, 2], drug = mm$Chemical_name,
                   well = mm$key, plate = pl, stringsAsFactors = FALSE)
  df$panel <- sprintf("Plate %d  (PC1 %.0f%%, PC2 %.0f%%)", pl, pv[1], pv[2])

  # flag wells far from the plate centroid, which usually means a failed well
  cen <- c(median(df$PC1), median(df$PC2))
  dist <- sqrt((df$PC1 - cen[1])^2 + (df$PC2 - cen[2])^2)
  out <- df[dist > 3 * mad(dist) + median(dist), ]
  if (nrow(out)) {
    cat("\nPlate ", pl, " - wells far from centroid:\n", sep = "")
    print(out[, c("well", "drug", "PC1", "PC2")], row.names = FALSE)
  }
  cat("\nPlate ", pl, ": ", sum(sel), " wells, PC1 ", pv[1], "%, PC2 ", pv[2], "%\n", sep = "")
  df
}

d <- do.call(rbind, lapply(c(15, 16), pca_plate))
d$drug <- factor(d$drug, levels = c("VEH", DRUGS),
                 labels = c("Vehicle", "Haloperidol", "Clozapine", "Risperidone"))
write.csv(d, file.path(RES, "figureS1_pca.csv"), row.names = FALSE)

figS1 <- ggplot(d, aes(PC1, PC2, colour = drug, shape = drug)) +
  geom_point(size = 2.7, stroke = 0.9) +
  scale_colour_manual(values = c(Vehicle = "grey50", Haloperidol = "#D55E00",
                                 Clozapine = "#0072B2", Risperidone = "#009E73"), name = NULL) +
  scale_shape_manual(values = c(Vehicle = 4, Haloperidol = 17,
                                Clozapine = 16, Risperidone = 15), name = NULL) +
  facet_wrap(~ panel, scales = "free") +
  labs(x = "PC1", y = "PC2",
       caption = paste0("Principal components computed separately within each plate, on the ",
                        NTOP, " most variable probes\nof variance-stabilised counts. ",
                        "Drug wells at 10 \u00b5M with plate-matched vehicle.")) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        strip.background = element_rect(fill = "grey93"),
        strip.text       = element_text(face = "bold", size = 9),
        plot.caption     = element_text(hjust = 0, size = 7.4, colour = "grey30",
                                        lineheight = 1.2))

ggsave(file.path(FIGS, "FigS1_pca.png"), figS1, width = 6.8, height = 3.7, dpi = 600)
ggsave(file.path(FIGS, "FigS1_pca.pdf"), figS1, width = 6.8, height = 3.7)
cat("\nWrote figures/FigS1_pca.png at 600 dpi\n")
cat("Now check: does haloperidol separate on plate 15, and do clozapine and\n",
    "risperidone stay mixed with vehicle on plate 16? If clozapine separates,\n",
    "cut the sentence in Section 3.4 rather than publishing a figure against it.\n", sep = "")
