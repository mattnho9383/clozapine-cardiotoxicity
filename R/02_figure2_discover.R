# =============================================================================
# 02_figure2_discover.R
#
# Inspects the two validation datasets so the Figure 2 analysis can be written
# against what is actually in the files. Computes nothing, fits nothing,
# writes no figure. Base R only: no GEOquery, no XML.
#
#   Rscript 02_figure2_discover.R 2>&1 | tee logs/fig2_discover.log
# =============================================================================

WORKDIR <- getwd()
DATA    <- file.path(WORKDIR, "data")
dir.create(DATA, recursive = TRUE, showWarnings = FALSE)
options(timeout = 1800)

hr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

# Parse the !Sample_* metadata block of a GEO series matrix without GEOquery.
sm_meta <- function(path) {
  ln <- readLines(gzfile(path), warn = FALSE)
  stop_at <- grep("^!series_matrix_table_begin", ln)
  if (length(stop_at)) ln <- ln[seq_len(stop_at[1] - 1)]
  ln <- ln[grepl("^!Sample_", ln)]
  parts <- strsplit(ln, "\t")
  keys  <- vapply(parts, function(p) sub("^!", "", p[1]), character(1))
  vals  <- lapply(parts, function(p) gsub('^"|"$', "", p[-1]))
  list(keys = keys, vals = vals)
}

# Print each metadata field with its distinct values, truncated.
sm_report <- function(m, max_show = 8) {
  for (i in seq_along(m$keys)) {
    v <- m$vals[[i]]
    u <- unique(v)
    if (length(u) <= 1 && !grepl("characteristic|title", m$keys[i])) next
    cat(sprintf("\n  %-34s  n=%d  distinct=%d\n", m$keys[i], length(v), length(u)))
    show <- head(u, max_show)
    for (s in show) cat("      ", substr(s, 1, 90), "\n", sep = "")
    if (length(u) > max_show) cat("       ... and ", length(u) - max_show, " more\n", sep = "")
  }
}

# ===================================================== GSE244740, human ======
hr("GSE244740  -- human iPSC-CM validation (files already local)")

sm244 <- file.path(DATA, "GSE244740", "GSE244740_series_matrix.txt.gz")
ct244 <- file.path(DATA, "GSE244740", "GSE244740_processed_data_counts.txt.gz")

if (!file.exists(sm244)) {
  cat("MISSING:", sm244, "\n")
} else {
  m <- sm_meta(sm244)
  cat("Metadata fields with more than one distinct value:\n")
  sm_report(m)

  ga <- m$vals[[which(m$keys == "Sample_geo_accession")[1]]]
  cat("\n  Samples in series matrix:", length(ga), "\n")
  cat("  First few GSM ids:", paste(head(ga, 4), collapse = ", "), "\n")
}

if (!file.exists(ct244)) {
  cat("MISSING:", ct244, "\n")
} else {
  cat("\nCounts file header (first 12 column names):\n")
  h <- readLines(gzfile(ct244), n = 1)
  cn <- strsplit(h, "\t")[[1]]
  cat("  total columns:", length(cn), "\n")
  for (x in head(cn, 12)) cat("      ", substr(gsub('"', "", x), 1, 70), "\n", sep = "")

  cat("\nFirst 3 data rows, first 5 columns:\n")
  d <- read.delim(gzfile(ct244), nrows = 3, check.names = FALSE)
  print(d[, seq_len(min(5, ncol(d)))])

  cat("\nDo the counts columns look like GSM ids, or like sample titles?\n")
  cat("  columns matching ^GSM :", sum(grepl("^GSM", gsub('"', "", cn))), "\n")
}

# ====================================================== GSE59905, rat ========
hr("GSE59905  -- rat DrugMatrix validation (needs downloading)")

MBASE <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59905/matrix/"
cat("Listing", MBASE, "\n")
listing <- tryCatch(readLines(MBASE, warn = FALSE), error = function(e) {
  cat("  could not read directory:", conditionMessage(e), "\n"); character(0)
})
mfiles <- unique(unlist(regmatches(listing,
            gregexpr("GSE59905[A-Za-z0-9_.\\-]*series_matrix\\.txt\\.gz", listing))))

if (!length(mfiles)) {
  cat("  No series matrix files found. Open the accession page by hand:\n",
      "  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59905\n")
} else {
  cat("  Series matrix files available:\n")
  for (f in mfiles) cat("      ", f, "\n", sep = "")

  dir.create(file.path(DATA, "GSE59905"), showWarnings = FALSE)
  for (f in mfiles) {
    dest <- file.path(DATA, "GSE59905", f)
    if (file.exists(dest) && file.size(dest) > 1000) { cat("  have ", f, "\n", sep = ""); next }
    cat("  downloading ", f, " ...\n", sep = "")
    tryCatch(download.file(paste0(MBASE, f), dest, mode = "wb", quiet = TRUE),
             error = function(e) cat("    FAILED: ", conditionMessage(e), "\n", sep = ""))
  }

  # Which matrix actually contains clozapine heart samples?
  for (f in mfiles) {
    dest <- file.path(DATA, "GSE59905", f)
    if (!file.exists(dest) || file.size(dest) < 1000) next
    cat("\n---- ", f, " ----\n", sep = "")
    m <- tryCatch(sm_meta(dest), error = function(e) NULL)
    if (is.null(m)) { cat("  could not parse\n"); next }

    ga <- m$vals[[which(m$keys == "Sample_geo_accession")[1]]]
    cat("  samples:", length(ga), "\n")

    ch <- m$vals[m$keys == "Sample_characteristics_ch1"]
    flat <- unlist(ch)
    n_cloz <- sum(grepl("clozapine", flat, ignore.case = TRUE))
    n_heart <- sum(grepl("heart", flat, ignore.case = TRUE))
    cat("  characteristic entries mentioning clozapine:", n_cloz, "\n")
    cat("  characteristic entries mentioning heart    :", n_heart, "\n")

    if (n_cloz > 0) {
      cat("  >>> THIS MATRIX HAS CLOZAKINE SAMPLES <<<\n")
      cat("  distinct characteristic fields:\n")
      for (k in seq_along(ch)) {
        u <- unique(ch[[k]])
        cat(sprintf("    field %d: %d distinct\n", k, length(u)))
        for (s in head(u, 6)) cat("        ", substr(s, 1, 80), "\n", sep = "")
        if (length(u) > 6) cat("        ... and ", length(u) - 6, " more\n", sep = "")
      }
      pl <- m$vals[[which(m$keys == "Sample_platform_id")[1]]]
      cat("  platform:", paste(unique(pl), collapse = ", "), "\n")
    }
  }
}

hr("Done. Paste this whole log back.")
