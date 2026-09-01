#!/usr/bin/env Rscript
# 09_faers_tables.R
# Builds Table 1 (disproportionality), Table 2 (outcomes) and Table S2
# (case characteristics) from the ten deduplicated OpenVigil extractions.
#
# Inputs   data/faers/*_DEDUPLICATED.xlsx   ten workbooks, one per MedDRA PT
#          data/faers/faers_2x2_cells.csv   contingency cells (see note below)
# Outputs  results/faers_table1.csv, table2_outcomes.csv, tableS2_demographics.csv
#
# Every parsing rule was validated against the published tables before being
# written here; the assertions at the foot re-check them on each run.
#
# NOTE ON THE 2x2 CELLS. Cells A and B are computed from the workbooks. Cells
# C and D (event counts across all other drugs) come from OpenVigil queries run
# through its web interface and cannot be derived from these files, so they are
# supplied as a documented input.

suppressPackageStartupMessages({ library(readxl) })

DATA_DIR    <- "data/faers"
RESULTS_DIR <- "results"
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

CLOZ_CASES  <- 103331
TOTAL_CASES <- 14104742

EVENTS <- c(
  "myocarditis_8_30_26_DEDUPLICATED.xlsx"                    = "Myocarditis",
  "left_atrial_enlargement_8_30_26_DEDUPLICATED.xlsx"        = "Left atrial enlargement",
  "left_ventricular_hypertrophy_8_30_26_DEDUPLICATED.xlsx"   = "Left ventricular hypertrophy",
  "echocardiogram_abnormal_8_30_DEDUPLICATED.xlsx"           = "Echocardiogram abnormal",
  "bundle_branch_block_right_8_30_26_DEDUPLICATED.xlsx"      = "Bundle branch block right",
  "Congestive_cardiomyopathy_8_30_26_DEDUPLICATED.xlsx"      = "Congestive cardiomyopathy",
  "Tachycardia_clz_DEDUPLICATED.xlsx"                        = "Tachycardia",
  "electrocardiogram_QT_prolonged_8_30_26_DEDUPLICATED.xlsx" = "Electrocardiogram QT prolonged",
  "cardiomyopathy_8_30_26_DEDUPLICATED.xlsx"                 = "Cardiomyopathy",
  "ejection_fraction_decreased_8_30_DEDUPLICATED.xlsx"       = "Ejection fraction decreased"
)

OUTCOME_LEVELS <- c(
  "Death", "Life-Threatening", "Hospitalization - Initial or Prolonged",
  "Disability", "Required Intervention to Prevent Permanent Impairment/Damage",
  "Congenital Anomaly")
split_outcome <- function(x) trimws(unlist(strsplit(as.character(x), "/ ", fixed = TRUE)))

COUNTRY_MAP <- c(
  AU = "Australia",      AUSTRALIA        = "Australia",
  US = "United States",  `UNITED STATES`  = "United States",
  CA = "Canada",         CANADA           = "Canada",
  GB = "United Kingdom", `UNITED KINGDOM` = "United Kingdom",
  DE = "Germany",        GERMANY          = "Germany",
  `NULL` = "Not specified", NONE          = "Not specified")
COUNTRY_LEVELS <- c("Australia","United States","Canada","United Kingdom",
                    "Germany","Other","Not specified")
map_country <- function(x) {
  out <- unname(COUNTRY_MAP[toupper(trimws(as.character(x)))])
  out[is.na(out)] <- "Other"; out
}

AGE_LEVELS <- c("<18","18-44","45-64","65+","Unknown")
band_age <- function(x) {
  a <- suppressWarnings(as.numeric(as.character(x)))
  ifelse(is.na(a),"Unknown", ifelse(a<18,"<18", ifelse(a<45,"18-44",
    ifelse(a<65,"45-64","65+"))))
}

ATTRIB_LEVELS <- c("Primary suspect","Secondary suspect","Interacting / concomitant")
map_attrib <- function(x) {
  r <- tolower(as.character(x))
  ifelse(grepl("primary",r),"Primary suspect",
    ifelse(grepl("secondary",r),"Secondary suspect","Interacting / concomitant"))
}

pct <- function(n, d) sprintf("%d (%.1f%%)", n, 100 * n / d)

read_event <- function(fname) {
  path <- file.path(DATA_DIR, fname)
  if (!file.exists(path)) stop("missing input: ", path)
  sheets <- excel_sheets(path)
  dedup  <- sheets[grepl("^Dedup", sheets)]
  if (length(dedup) != 1L) stop("expected one 'Deduplicated' sheet in ", fname)
  df <- read_excel(path, sheet = dedup, col_types = "text", .name_repair = "minimal")
  df <- df[!is.na(df$Case_id) & df$Case_id != "", , drop = FALSE]
  if (anyDuplicated(df$Case_id)) stop("duplicate Case_id in ", fname)
  df
}

cat("Reading", length(EVENTS), "workbooks from", DATA_DIR, "\n")
dat <- lapply(names(EVENTS), read_event)
names(dat) <- unname(EVENTS)
for (e in names(dat)) cat(sprintf("  %-32s %6d cases\n", e, nrow(dat[[e]])))

cells <- read.csv(file.path(DATA_DIR, "faers_2x2_cells.csv"), stringsAsFactors = FALSE)

table1 <- do.call(rbind, lapply(names(dat), function(e) {
  A <- as.numeric(nrow(dat[[e]]))
  row <- cells[cells$event == e, ]
  if (nrow(row) != 1L) stop("no 2x2 row for ", e)
  if (row$A_cloz_event != A)
    stop(sprintf("2x2 cell A (%d) disagrees with workbook (%d) for %s",
                 row$A_cloz_event, A, e))
  B <- CLOZ_CASES - A
  C <- as.numeric(row$C_nocloz_event)
  D <- (TOTAL_CASES - CLOZ_CASES) - C
  N <- A + B + C + D
  ror    <- (A * D) / (B * C)
  se     <- sqrt(1/A + 1/B + 1/C + 1/D)
  prr    <- (A / (A + B)) / (C / (C + D))
  chi2   <- N * (abs(A * D - B * C) - N / 2)^2 /
              ((A + B) * (C + D) * (A + C) * (B + D))
  expect <- (A + B) * (A + C) / N
  ic     <- log2((A + 0.5) / (expect + 0.5))
  icsd   <- sqrt((1/log(2)^2) * ((N - A + 0.5) / ((A + 0.5) * (1 + N + 1)) +
                                 (N - expect + 0.5) / ((expect + 0.5) * (1 + N + 1))))
  data.frame(event = e, n_cases = A, rate_pct = round(100 * A / CLOZ_CASES, 3),
    ROR = round(ror, 2),
    ROR_lo = round(exp(log(ror) - 1.96 * se), 2),
    ROR_hi = round(exp(log(ror) + 1.96 * se), 2),
    PRR = round(prr, 2), chi2_yates = round(chi2, 1),
    IC = round(ic, 2), IC025 = row$IC025_openvigil,
    IC975 = row$IC975_openvigil, stringsAsFactors = FALSE)
}))
table1 <- table1[order(-table1$ROR), ]
write.csv(table1, file.path(RESULTS_DIR, "faers_table1.csv"), row.names = FALSE)

table2 <- do.call(rbind, lapply(names(dat), function(e) {
  d <- dat[[e]]; n <- nrow(d)
  labs <- lapply(d$Outcome, split_outcome)
  counts <- vapply(OUTCOME_LEVELS, function(L)
    sum(vapply(labs, function(v) L %in% v, TRUE)), 1L)
  out <- data.frame(Event = e, `Total cases` = n, check.names = FALSE,
                    stringsAsFactors = FALSE)
  for (i in seq_along(OUTCOME_LEVELS)) out[[OUTCOME_LEVELS[i]]] <- pct(counts[i], n)
  out
}))
names(table2)[names(table2) == "Hospitalization - Initial or Prolonged"] <- "Hospitalisation"
names(table2)[names(table2) ==
  "Required Intervention to Prevent Permanent Impairment/Damage"] <- "Required intervention"
write.csv(table2, file.path(RESULTS_DIR, "table2_outcomes.csv"), row.names = FALSE)

s2 <- list()
s2[["Total unique cases"]] <- list(Total = vapply(dat, nrow, 1L))
s2[["Suspect attribution"]] <- lapply(setNames(ATTRIB_LEVELS, ATTRIB_LEVELS), function(L)
  vapply(dat, function(d) sum(map_attrib(d$`Role code`) == L), 1L))
s2[["Sex"]] <- lapply(setNames(c("Male","Female","Unknown / not specified"),
                               c("Male","Female","Unknown / not specified")), function(L) {
  vapply(dat, function(d) {
    g <- toupper(trimws(as.character(d$Gender)))
    if (L == "Male") sum(g %in% c("M","MALE"))
    else if (L == "Female") sum(g %in% c("F","FEMALE"))
    else sum(!g %in% c("M","MALE","F","FEMALE"))
  }, 1L)
})
s2[["Age (years)"]] <- lapply(setNames(AGE_LEVELS, AGE_LEVELS), function(L)
  vapply(dat, function(d) sum(band_age(d$`Age in report`) == L), 1L))
s2[["Reporter country"]] <- lapply(setNames(COUNTRY_LEVELS, COUNTRY_LEVELS), function(L)
  vapply(dat, function(d) sum(map_country(d$`Reporter country`) == L), 1L))

rows <- list()
for (blk in names(s2)) {
  first <- TRUE
  for (sub in names(s2[[blk]])) {
    v <- s2[[blk]][[sub]]
    cells_txt <- if (blk == "Total unique cases") as.character(v)
                 else mapply(pct, v, vapply(dat, nrow, 1L))
    rows[[length(rows) + 1L]] <- c(if (first) blk else "", if (blk == "Total unique cases") "" else sub, cells_txt)
    first <- FALSE
  }
}
s2df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(s2df) <- c("Category", "Subgroup", names(dat))
write.csv(s2df, file.path(RESULTS_DIR, "tableS2_demographics.csv"), row.names = FALSE)

stopifnot(
  nrow(dat[["Myocarditis"]]) == 1558L,
  nrow(dat[["Tachycardia"]]) == 2072L,
  sum(map_attrib(dat[["Myocarditis"]]$`Role code`) == "Primary suspect") == 1097L,
  sum(toupper(trimws(dat[["Myocarditis"]]$Gender)) %in% c("M","MALE")) == 1059L,
  sum(band_age(dat[["Myocarditis"]]$`Age in report`) == "18-44") == 840L,
  sum(vapply(dat[["Myocarditis"]]$Outcome,
             function(x) "Death" %in% split_outcome(x), TRUE)) == 122L,
  sum(vapply(dat[["Myocarditis"]]$Outcome,
             function(x) "Hospitalization - Initial or Prolonged" %in% split_outcome(x),
             TRUE)) == 852L,
  abs(table1$ROR[table1$event == "Myocarditis"] - 35.45) < 0.02,
  abs(table1$IC[table1$event == "Myocarditis"]  -  4.79) < 0.02,
  all(table1$ROR_lo > 1), all(table1$PRR > 2),
  all(table1$chi2_yates > 4), all(table1$IC025 > 0)
)

cat("\nAll assertions passed.\n")
cat("Wrote faers_table1.csv, table2_outcomes.csv, tableS2_demographics.csv to",
    RESULTS_DIR, "\n")
