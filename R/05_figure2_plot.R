# =============================================================================
# 05_figure2_plot.R  --  Figure 2, plotting only
#
# Reads results/figure2_data.csv and plots. Fits nothing, downloads nothing,
# so it runs in seconds and you can iterate on the design without refitting.
# This is also the arrangement the repository needs: figures generated from
# result files, not from values typed into the plotting script.
#
#   Rscript 05_figure2_plot.R
#
# Changes from 04:
#   - Okabe-Ito colourblind-safe palette; the previous green/red pair is the
#     hardest combination for red-green deficiency, which affects ~8% of men
#   - shape encodes dataset redundantly with colour, so the figure survives
#     greyscale printing and photocopying
#   - axis limits padded so no interval is clipped at the panel edge
#   - sample sizes and the adjusted/unadjusted distinction moved to the caption,
#     which also stops the legend overflowing the panel width
# =============================================================================

WORKDIR <- Sys.getenv("CLOZ_HOME", unset = file.path(path.expand("~"), "clozapine"))
RES  <- file.path(WORKDIR, "results"); FIGS <- file.path(WORKDIR, "figures")
lib <- Sys.getenv("R_LIBS_USER", unset = file.path(path.expand("~"), "Rlibs", "R-4.6.0"))
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(ggplot2))

GENES <- c("TNFRSF12A", "EGLN3", "VEGFA")
LEV   <- c("Discovery hiPSC-CM", "Independent hiPSC-CM", "Rat heart in vivo")

# Okabe-Ito: distinguishable under deuteranopia and protanopia
PAL <- c("Discovery hiPSC-CM"   = "#009E73",   # bluish green
         "Independent hiPSC-CM" = "#D55E00",   # vermillion
         "Rat heart in vivo"    = "#0072B2")   # blue
SHP <- c("Discovery hiPSC-CM"   = 21,          # circle
         "Independent hiPSC-CM" = 22,          # square
         "Rat heart in vivo"    = 24)          # triangle

d <- read.csv(file.path(RES, "figure2_data.csv"), stringsAsFactors = FALSE)
stopifnot(all(c("symbol", "logFC", "lo", "hi", "adj.P.Val", "dataset") %in% names(d)))
cat("Rows:", nrow(d), " datasets:", paste(unique(d$dataset), collapse = " | "), "\n")

d$sig   <- !is.na(d$adj.P.Val) & d$adj.P.Val < 0.05
d$fillc <- ifelse(d$sig, PAL[d$dataset], "white")   # open symbol when not significant
d$dataset <- factor(d$dataset, levels = LEV)
d$symbol  <- factor(d$symbol,  levels = rev(GENES))

# flag any point whose interval and adjusted p disagree, so it is a deliberate
# choice rather than something a reviewer notices first
mismatch <- d[!d$sig & (d$lo > 0 | d$hi < 0), ]
if (nrow(mismatch)) {
  cat("\nNOTE: interval excludes zero but FDR-adjusted p >= 0.05 for:\n")
  print(mismatch[, c("symbol", "dataset", "logFC", "lo", "hi", "adj.P.Val")], row.names = FALSE)
  cat("The caption states that intervals are unadjusted. Keep that sentence.\n\n")
}

pad  <- 0.18
xlim <- c(min(d$lo, na.rm = TRUE) - pad, max(d$hi, na.rm = TRUE) + pad)
dodge <- position_dodge(width = 0.60)

cap <- paste(
  "Discovery n = 3 vs 12 wells; independent human n = 3 vs 22 wells; rat n = 9 vs 32 arrays.",
  "Bars are unadjusted 95% intervals; open symbols indicate FDR-adjusted p \u2265 0.05.",
  sep = "\n")

fig2 <- ggplot(d, aes(logFC, symbol,
                      colour = dataset, shape = dataset, fill = fillc)) +
  annotate("rect", ymin = 1.5, ymax = 2.5, xmin = -Inf, xmax = Inf,
           fill = "grey96", colour = NA) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.45) +
  geom_linerange(aes(xmin = lo, xmax = hi), position = dodge, linewidth = 0.6) +
  geom_point(position = dodge, size = 2.6, stroke = 0.85) +
  scale_fill_identity(guide = "none") +
  scale_colour_manual(values = PAL, name = NULL, drop = FALSE) +
  scale_shape_manual(values = SHP, name = NULL, drop = FALSE) +
  scale_x_continuous(limits = xlim, breaks = seq(-3, 1, 0.5),
                     expand = expansion(mult = 0)) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(fill = PAL, linetype = 0)),
         shape  = guide_legend(nrow = 1, override.aes = list(fill = PAL, linetype = 0))) +
  labs(x = expression(log[2]~"fold-change vs vehicle"), y = NULL, caption = cap) +
  theme_classic(base_size = 10) +
  theme(axis.text.y     = element_text(face = "italic", size = 10.5, colour = "black"),
        axis.text.x     = element_text(colour = "black"),
        axis.title.x    = element_text(margin = margin(t = 6)),
        axis.line.y     = element_blank(),
        axis.ticks.y    = element_blank(),
        legend.position = "bottom",
        legend.margin   = margin(t = 0, b = 0),
        legend.key.size = unit(10, "pt"),
        legend.text     = element_text(size = 8.5),
        plot.caption    = element_text(hjust = 0, size = 7.4, colour = "grey30",
                                       margin = margin(t = 6), lineheight = 1.15),
        plot.margin     = margin(6, 8, 4, 4))

ggsave(file.path(FIGS, "Fig2_forest.png"), fig2, width = 6.8, height = 3.1, dpi = 600)
ggsave(file.path(FIGS, "Fig2_forest.pdf"), fig2, width = 6.8, height = 3.1)
cat("Wrote figures/Fig2_forest.png and .pdf at 600 dpi\n")
