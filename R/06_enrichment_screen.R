# =============================================================================
# 06_enrichment_screen.R
#
# Covers four manuscript items in one pass over the discovery data:
#   Table 3   hypergeometric enrichment against three MSigDB GO:BP sets
#   Sec 3.2   DEG counts (expected: haloperidol 95, clozapine 9, risperidone 0)
#   Sec 3.4   21-compound cardiac enrichment screen
#   Fig S1    PCA by plate, for the "clozapine did not separate" claim
#
# Writes:
#   results/deg_counts.csv
#   results/table3_enrichment.csv
#   results/drug_screen.csv
#   results/plate_universes.csv
#   figures/FigS1_pca.png / .pdf
#   results/sessionInfo_enrichment.txt
#
#   Rscript 06_enrichment_screen.R 2>&1 | tee logs/enrichment.log
#
# Every computed value is printed next to the value the manuscript claims, so
# agreement or disagreement is visible without cross-referencing by hand.
# =============================================================================

WORKDIR <- Sys.getenv("CLOZ_HOME", unset = file.path(path.expand("~"), "clozapine"))
DATA <- file.path(WORKDIR, "data"); RES <- file.path(WORKDIR, "results")
FIGS <- file.path(WORKDIR, "figures")
for (d in c(RES, FIGS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
lib <- Sys.getenv("R_LIBS_USER", unset = file.path(path.expand("~"), "Rlibs", "R-4.6.0"))
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))

need <- setdiff(c("DESeq2", "ggplot2", "msigdbr"), rownames(installed.packages()))
if (length(need))
  stop("Missing: ", paste(need, collapse = ", "),
       "\nInstall on a login node:\n  R -e 'install.packages(\"msigdbr\", lib=\"", lib, "\")'")
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2); library(msigdbr) })

MIN_EXPR <- 100
hr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

# ------------------------------------------------------------------ load ----
hash <- read.csv(file.path(DATA, "GSE262419_hash.csv.gz"))
hash$Chemical_name <- trimws(hash$Chemical_name)
p15 <- read.csv(file.path(DATA, "GSE262419_Plate15.csv.gz"), row.names = 1, check.names = FALSE)
p16 <- read.csv(file.path(DATA, "GSE262419_Plate16.csv.gz"), row.names = 1, check.names = FALSE)
counts_all <- cbind(p15, p16)

wl <- gsub("[0-9]+", "", hash$Well_ID); wn <- gsub("[^0-9]", "", hash$Well_ID)
hash$key <- paste0("Plate", hash$Plate_ID, "-", sprintf("%s%02d", wl, as.integer(wn)))
meta <- hash[hash$Plate_ID %in% c(15, 16), ]
pres <- intersect(meta$key, colnames(counts_all))
cf   <- counts_all[, pres, drop = FALSE]
meta <- meta[match(colnames(cf), meta$key), ]
meta$cn <- suppressWarnings(as.numeric(as.character(meta$Chemical_Concentration_uM)))
cat("Loaded", nrow(cf), "probes x", ncol(cf), "wells\n")

# --------------------------------------------------------- DE helper --------
deseq_drug <- function(drug, conc = 10, quiet = TRUE) {
  pl <- unique(meta$Plate_ID[meta$Chemical_name == drug])[1]
  if (is.na(pl)) return(NULL)
  keep <- (meta$Chemical_name == drug & abs(meta$cn - conc) < 1e-6 & meta$Plate_ID == pl) |
          (meta$Chemical_name == "VEH" & meta$Plate_ID == pl)
  mm <- meta[keep, ]
  if (sum(mm$Chemical_name == drug) < 2) return(NULL)
  mm$grp <- factor(ifelse(mm$Chemical_name == drug, "Drug", "VEH"), levels = c("VEH", "Drug"))
  d <- DESeqDataSetFromMatrix(round(cf[, keep, drop = FALSE]), mm, ~ grp)
  d <- DESeq(d[rowSums(counts(d)) >= 10, ], quiet = TRUE)
  r <- as.data.frame(results(d, contrast = c("grp", "Drug", "VEH")))
  r$symbol <- gsub("_[0-9]+$", "", rownames(r))
  attr(r, "n") <- c(sum(mm$grp == "Drug"), sum(mm$grp == "VEH"))
  attr(r, "plate") <- pl
  r
}

# universe = probes meeting the expression threshold in that comparison
universe_of <- function(r) unique(r$symbol[!is.na(r$baseMean) & r$baseMean >= MIN_EXPR])
degs_of <- function(r) unique(r$symbol[!is.na(r$padj) & r$padj < 0.05 &
                                       !is.na(r$baseMean) & r$baseMean >= MIN_EXPR])

# ============================================================ DEG COUNTS =====
hr("DEG counts  (adjusted p < 0.05 and baseMean >= 100)")
MAIN <- c("Haloperidol", "Clozapine", "Risperidone")
EXPECT_DEG <- c(Haloperidol = 95, Clozapine = 9, Risperidone = 0)

fits <- list()
deg_rows <- list()
for (dr in MAIN) {
  r <- deseq_drug(dr)
  if (is.null(r)) { cat(dr, ": no wells\n"); next }
  fits[[dr]] <- r
  u <- universe_of(r); g <- degs_of(r)
  cat(sprintf("%-13s plate %-3s %2d vs %2d wells | universe %5d | DEGs %3d  (manuscript %d)%s\n",
              dr, attr(r, "plate"), attr(r, "n")[1], attr(r, "n")[2],
              length(u), length(g), EXPECT_DEG[dr],
              if (length(g) == EXPECT_DEG[dr]) "  MATCH" else "  <-- DIFFERS"))
  deg_rows[[dr]] <- data.frame(drug = dr, plate = attr(r, "plate"),
                               n_treated = attr(r, "n")[1], n_vehicle = attr(r, "n")[2],
                               universe = length(u), n_deg = length(g),
                               manuscript = EXPECT_DEG[dr])
}
deg_tab <- do.call(rbind, deg_rows)
write.csv(deg_tab, file.path(RES, "deg_counts.csv"), row.names = FALSE)

cat("\nClozapine DEGs:\n"); print(sort(degs_of(fits[["Clozapine"]])))

# plate universes: CHANGES.md back-calculated ~1598 (plate 15) and ~960 (plate 16)
uni <- data.frame(drug = names(fits),
                  plate = sapply(fits, function(r) attr(r, "plate")),
                  universe = sapply(fits, function(r) length(universe_of(r))))
write.csv(uni, file.path(RES, "plate_universes.csv"), row.names = FALSE)
cat("\nPlate universes (Anand will ask why these differ):\n"); print(uni, row.names = FALSE)

# ========================================================== ENRICHMENT =======
hr("Table 3  --  hypergeometric enrichment against MSigDB GO:BP")

SETNAMES <- c("GOBP_CARDIAC_MUSCLE_CONTRACTION",
              "GOBP_SARCOMERE_ORGANIZATION",
              "GOBP_CARDIAC_MUSCLE_TISSUE_DEVELOPMENT")
m <- tryCatch(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"),
              error = function(e) msigdbr(species = "Homo sapiens",
                                          category = "C5", subcategory = "GO:BP"))
SETS <- lapply(SETNAMES, function(s) unique(m$gene_symbol[m$gs_name == s]))
names(SETS) <- SETNAMES
cat("Gene set sizes in MSigDB:", paste(sprintf("%s=%d", SETNAMES, lengths(SETS)),
                                       collapse = ", "), "\n\n")

enrich <- function(r, set) {
  u <- universe_of(r); g <- degs_of(r)
  su <- intersect(set, u)                       # eligible members of the set
  hit <- intersect(g, su)
  expct <- length(su) * length(g) / length(u)
  p <- if (!length(su) || !length(g)) 1 else
       phyper(length(hit) - 1, length(su), length(u) - length(su), length(g), lower.tail = FALSE)
  list(eligible = length(su), observed = length(hit), expected = expct,
       p = p, genes = sort(hit))
}

t3 <- list()
for (sn in SETNAMES) {
  cat("--", sn, "\n")
  for (dr in names(fits)) {
    e <- enrich(fits[[dr]], SETS[[sn]])
    cat(sprintf("   %-13s eligible %3d | observed %2d | expected %5.2f | p = %.4g%s\n",
                dr, e$eligible, e$observed, e$expected, e$p,
                if (length(e$genes)) paste0("  [", paste(e$genes, collapse = ", "), "]") else ""))
    t3[[length(t3) + 1]] <- data.frame(gene_set = sn, drug = dr, eligible = e$eligible,
                                       observed = e$observed, expected = round(e$expected, 2),
                                       p = e$p,
                                       genes = paste(e$genes, collapse = "; "))
  }
  cat("\n")
}
t3 <- do.call(rbind, t3)
write.csv(t3, file.path(RES, "table3_enrichment.csv"), row.names = FALSE)

cat("Manuscript claims for comparison:\n",
    "  CONTRACTION  haloperidol 8 obs vs 2.20 exp, p = 0.0011, eligible 37\n",
    "               clozapine   0 of 32, p = 1\n",
    "  SARCOMERE    haloperidol 4 obs vs 1.19 exp, p = 0.027, eligible 20\n",
    "  DEVELOPMENT  haloperidol p = 0.096; clozapine VEGFA only, p = 0.30\n", sep = "")

# ========================================================= DRUG SCREEN =======
hr("Section 3.4  --  21-compound cardiac enrichment screen")
CANDIDATES <- c("Haloperidol","Droperidol","Domperidone","Risperidone","Dofetilide",
                "D,1 Sotalol","Quinidine","Bepridil","Cisapride","Terfenadine",
                "Astemizole","Ibutilide","Azimilide","Vernakalant","Ranolazine",
                "Verapamil","Nifedipine","Disopyramide","Mexiletine","Nitrendipine","Atenolol")
gs <- SETS$GOBP_CARDIAC_MUSCLE_CONTRACTION

screen <- do.call(rbind, lapply(CANDIDATES, function(dr) {
  r <- if (!is.null(fits[[dr]])) fits[[dr]] else tryCatch(deseq_drug(dr), error = function(e) NULL)
  if (is.null(r)) { cat(sprintf("  %-14s not present on plates 15/16\n", dr)); return(NULL) }
  e <- enrich(r, gs)
  data.frame(drug = dr, plate = attr(r, "plate"), n_deg = length(degs_of(r)),
             eligible = e$eligible, observed = e$observed,
             expected = round(e$expected, 2), p = e$p)
}))
screen <- screen[order(screen$p), ]
print(screen, row.names = FALSE)
write.csv(screen, file.path(RES, "drug_screen.csv"), row.names = FALSE)
cat("\nManuscript reports significant: ibutilide 0.0001, dofetilide 0.0013,\n",
    "astemizole 0.0071, haloperidol 0.012, bepridil 0.021, verapamil 0.037;\n",
    "droperidol, domperidone and risperidone all 0 DEGs.\n", sep = "")

# ============================================================= FIG S1 ========
hr("Figure S1  --  PCA by plate")
sel <- meta$Chemical_name %in% c("VEH", MAIN) &
       (meta$Chemical_name == "VEH" | abs(meta$cn - 10) < 1e-6)
mm <- meta[sel, ]
dp <- DESeqDataSetFromMatrix(round(cf[, sel, drop = FALSE]),
                             data.frame(row.names = mm$key,
                                        drug = mm$Chemical_name,
                                        plate = factor(mm$Plate_ID)), ~ 1)
dp <- dp[rowSums(counts(dp)) >= 10, ]
vs <- vst(dp, blind = TRUE)
pc <- plotPCA(vs, intgroup = c("drug", "plate"), returnData = TRUE)
pv <- round(100 * attr(pc, "percentVar"))
pc$plate <- paste("Plate", pc$plate)
pc$drug  <- factor(pc$drug, levels = c("VEH", MAIN),
                   labels = c("Vehicle", "Haloperidol", "Clozapine", "Risperidone"))

figS1 <- ggplot(pc, aes(PC1, PC2, colour = drug, shape = drug)) +
  geom_point(size = 2.6, stroke = 0.8) +
  scale_colour_manual(values = c(Vehicle = "grey45", Haloperidol = "#D55E00",
                                 Clozapine = "#0072B2", Risperidone = "#009E73"), name = NULL) +
  scale_shape_manual(values = c(Vehicle = 4, Haloperidol = 17, Clozapine = 16,
                                Risperidone = 15), name = NULL) +
  facet_wrap(~ plate, scales = "free") +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"),
       caption = "Variance-stabilised counts, all probes. Drug wells at 10 \u00b5M against plate-matched vehicle.") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.background = element_rect(fill = "grey93"),
        strip.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey30"))

ggsave(file.path(FIGS, "FigS1_pca.png"), figS1, width = 6.8, height = 3.6, dpi = 600)
ggsave(file.path(FIGS, "FigS1_pca.pdf"), figS1, width = 6.8, height = 3.6)
cat("Wrote figures/FigS1_pca.png\n")
cat("\nCheck by eye: haloperidol should separate from vehicle on its plate,\n",
    "clozapine and risperidone should not. If clozapine DOES separate, the\n",
    "sentence in Section 3.4 is wrong and must be cut, not illustrated.\n", sep = "")

writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo_enrichment.txt"))
hr("Done")
