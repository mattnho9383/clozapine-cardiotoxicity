#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 00_download_geo.R  --  fetch the public GEO inputs. Run once before 01-08.
#
# These files are ~230 MB and are not stored in the repository. They are
# permanently archived at GEO under the accessions below and are downloaded
# here rather than deposited, so that the repository contains only work
# original to this study.
# ---------------------------------------------------------------------------

fetch <- function(url, dest) {
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(dest)) { cat("  have    ", basename(dest), "\n"); return(invisible(TRUE)) }
  cat("  getting ", basename(dest), "\n")
  ok <- tryCatch({ download.file(url, dest, mode = "wb", quiet = TRUE); TRUE },
                 error = function(e) { cat("  FAILED  ", conditionMessage(e), "\n"); FALSE })
  if (!ok && file.exists(dest)) unlink(dest)
  invisible(ok)
}

S <- "https://ftp.ncbi.nlm.nih.gov/geo/series"
P <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms"

cat("GSE262419  discovery, hiPSC-CM TempO-Seq\n")
fetch(paste0(S,"/GSE262nnn/GSE262419/suppl/GSE262419_Plate15.csv.gz"), "data/GSE262419_Plate15.csv.gz")
fetch(paste0(S,"/GSE262nnn/GSE262419/suppl/GSE262419_Plate16.csv.gz"), "data/GSE262419_Plate16.csv.gz")
fetch(paste0(S,"/GSE262nnn/GSE262419/suppl/GSE262419_hash.csv.gz"),    "data/GSE262419_hash.csv.gz")

cat("GSE244740  validation, independent hiPSC-CM\n")
fetch(paste0(S,"/GSE244nnn/GSE244740/matrix/GSE244740_series_matrix.txt.gz"),
      "data/GSE244740/GSE244740_series_matrix.txt.gz")
fetch(paste0(S,"/GSE244nnn/GSE244740/suppl/GSE244740_processed_data_counts.txt.gz"),
      "data/GSE244740/GSE244740_processed_data_counts.txt.gz")

cat("GSE59905  validation, rat heart (GPL5426)\n")
fetch(paste0(S,"/GSE59nnn/GSE59905/matrix/GSE59905-GPL5426_series_matrix.txt.gz"),
      "data/GSE59905/GSE59905-GPL5426_series_matrix.txt.gz")
fetch(paste0(P,"/GPL5nnn/GPL5426/annot/GPL5426.annot.gz"), "data/GPL5426.annot.gz")

cat("\nDone. Any line marked FAILED must be downloaded by hand from the\n",
    "accession page at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<GSE>\n")
