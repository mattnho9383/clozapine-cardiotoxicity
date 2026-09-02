# Clozapine cardiotoxicity: pharmacovigilance and cross-species transcriptomic analysis

Analysis code and derived results for:

> Nho M, Abdel-Latif A, Singh AP. *Clozapine-associated cardiac safety signals are
> not mirrored by contractile gene engagement in a human cardiomyocyte model.*
> [journal, year, DOI once available]

Code archived at https://doi.org/10.5281/zenodo.22255746

The study pairs disproportionality analysis of the FDA Adverse Event Reporting
System with transcriptomic analysis of human iPSC-derived cardiomyocytes exposed
to clozapine, risperidone and haloperidol, plus validation in an independent human
dataset and in rat myocardium.

No new experimental data were generated. Every input is either public or included
here.

---

## What is and is not in this repository

**Included.** All analysis code, all derived result tables, all figures, session
information for each script, and the ten deduplicated FAERS extractions.

**Not included: the GEO source data** (~230 MB). The scripts download it directly
from NCBI at run time. It is permanently archived under the accessions below, so
depositing a copy here would add bulk without adding reproducibility.

**Why the FAERS extractions *are* included.** FAERS is refreshed quarterly. The
extractions here are a snapshot taken 18 December 2024 and cannot be regenerated
identically afterwards — anyone re-running the same query today will get different
counts. They are the one irreplaceable input in the project.

---

## Data sources

| Accession | Description | Platform |
|---|---|---|
| GSE262419 | Discovery. hiPSC-CM, TempO-Seq targeted panel, plates 15 and 16 | TempO-Seq |
| GSE244740 | Validation, human. hiPSC-CM, 42 compounds across 8 concentrations | TempO-Seq |
| GSE59905 | Validation, rat. DrugMatrix cardiac tissue, 95 mg/kg, days 1/3/5 | GPL5426 |

Gene set sizes in MSigDB as used here: GOBP_CARDIAC_MUSCLE_CONTRACTION 156, GOBP_SARCOMERE_ORGANIZATION 51, GOBP_CARDIAC_MUSCLE_TISSUE_DEVELOPMENT 253.

FAERS accessed through OpenVigil 2.1, MedDRA v24, database snapshot 18 December
2024. Case counts: 103,331 clozapine cases within 14,104,742 total distinct cases.

---

## Layout

```
R/          analysis scripts and the reanalysis notebook
data/faers/ ten deduplicated OpenVigil extractions + 2x2 contingency cells
results/    derived tables, one per analysis step, plus sessionInfo
figures/    Figures 1, 2 and S1 as vector PDF and raster PNG
```

### Scripts

| Script | Produces |
|---|---|
| `00_download_geo.R` | Downloads the public GEO inputs. Run once before 01–09 |
| `01_figure1_concentration.R` | Concentration–response for the nine clozapine genes (Fig. 1) |
| `02_figure2_discover.R` | Discovery effect sizes for cross-species comparison |
| `03_figure2.R` | Validation in GSE244740 and GSE59905 |
| `04_figure2_final.R` | Assembles the cross-species comparison table |
| `05_figure2_plot.R` | Cross-species forest plot (Fig. 2) |
| `06_enrichment_screen.R` | Hypergeometric cardiac gene-set enrichment (Table 3) |
| `07_figureS1_pca.R` | Per-plate principal component analysis (Fig. S1) |
| `08_sensitivity_outlier.R` | Vehicle-well sensitivity analysis (plate 16, well B02) |
| `09_faers_tables.R` | Tables 1, 2 and S2 from the FAERS extractions |
| `check_con_direction.R` | Verifies the direction of the GSE244740 concentration series |
| `clozapine_reanalysis.ipynb` | Interactive notebook from which several analyses were originally run |

---

## Reproducing the analysis

```bash
module load R/4.6.0
export R_LIBS_USER=$HOME/Rlibs/R-4.6.0

Rscript R/00_download_geo.R
Rscript R/01_figure1_concentration.R
Rscript R/02_figure2_discover.R
Rscript R/03_figure2.R
Rscript R/04_figure2_final.R
Rscript R/05_figure2_plot.R
Rscript R/06_enrichment_screen.R
Rscript R/07_figureS1_pca.R
Rscript R/08_sensitivity_outlier.R
Rscript R/09_faers_tables.R
```

Run `00_download_geo.R` once to fetch the GEO inputs (~230 MB); scripts 01–08 then
read them from `data/`. Script 09 reads the workbooks in `data/faers/` and requires
`readxl`.

`09_faers_tables.R` ends with assertions against the published values — 1,558
myocarditis cases, 1,097 with a primary-suspect record, 122 deaths, 852
hospitalisations, ROR 35.45, IC 4.79, and all four Evans criteria satisfied for
all ten events. **If an assertion fails, the inputs have changed and the output
should not be used.**

---

## How the FAERS data was prepared

This is the part of the pipeline that was *not* fully scripted, and the
distinction matters for anyone trying to reproduce it.

**1. Extraction — web interface, not scripted.** Queries were run through the
OpenVigil 2.1 web interface, one per MedDRA preferred term, against the 18
December 2024 snapshot. Search terms covered all five marketed names: CLOZAPINE,
CLOZARIL, LEPONEX, FAZACLO and VERSACLOZ. The queries themselves cannot be
reproduced from code here; the exports they produced are deposited instead.

**2. Column-alignment repair.** The OpenVigil tabular export writes the
indication, event and reporter country fields unquoted, so values containing
commas split across cells and displace every subsequent field rightwards. Across
the ten extractions, 996 of 14,176 rows (7.0%) were affected, arising from
fourteen distinct comma-bearing values. Two properties are worth noting: the
defect is invisible to a distinct-value count, because in two files the spurious
insertions and deletions cancelled exactly; and it removes real cases rather than
merely corrupting metadata — tachycardia loses 99 cases (4.8%) unrepaired.

**3. Deduplication to one row per case.** For each `Case_id` a single record was
retained, selected first on suspect-drug attribution — primary suspect where one
existed, falling back to secondary, then to interacting or concomitant — and
thereafter on recency, taking the most recent date of receipt and, where dates
tied, the highest ISR.

**Steps 2 and 3 are documented inside the workbooks themselves.** Each file in
`data/faers/` carries four sheets: the retained rows, every dropped row, a
per-case selection log giving the reason each row was kept, and a README
describing the operation. The deduplication is therefore fully auditable and
reversible from the deposited files, even though the code that performed it is
not included.

**4. Table construction — scripted.** `09_faers_tables.R` takes the deduplicated
workbooks and produces Tables 1, 2 and S2, including all disproportionality
statistics. Every parsing rule in it was validated against the published tables
before being committed.

### One derived input

Cells A and B of each 2 × 2 table are computed from the workbooks. Cells C and D —
counts of each event across all other drugs — come from OpenVigil queries and
cannot be derived from these files, so they are supplied as a documented input in
`data/faers/faers_2x2_cells.csv`. These were recovered by inversion from the
reported proportional reporting ratio and are accurate to within approximately
two cases; resulting PRR values may therefore differ from the published table in
the second decimal place. Information component point estimates are computed by
the script, while the 95% credible interval bounds are taken from OpenVigil 2.1
output, which applies a gamma-posterior approximation, and are recorded in the
same file.

### Known limitation

Table S2's indication grouping is not reproduced by `09_faers_tables.R`. The
mapping from free-text indication strings to the reported categories was applied
manually and is not recoverable from the deposited outputs. All other Table S2
blocks — attribution, sex, age band and reporter country — are generated by the
script and reproduce the published values exactly.

---

## Multi-probe genes in the TempO-Seq panel

The panel carries more than one probe for some genes; *EGLN3*, for example, is
represented by both `EGLN3_88861` and `EGLN3_90114` on plate 16. Where several
probes mapped to one gene symbol, the probe with the highest mean normalized
expression was retained. This rule is applied identically in
`01_figure1_concentration.R`, `06_enrichment_screen.R` and
`08_sensitivity_outlier.R`.

---

## Environment

| Component | Version |
|---|---|
| R | 4.6.0 (2026-04-24) |
| OS | Red Hat Enterprise Linux 8.10 |
| BLAS / LAPACK | FlexiBLAS, NETLIB backend; LAPACK 3.12.1 |
| Bioconductor | 3.23 |
| DESeq2 | 1.52.0 |
| limma | 3.68.5 |
| ggplot2 | 4.0.3 |
| msigdbr | 26.1.1 |
| MSigDB release | as bundled with msigdbr 26.1.1 |
| readxl | 1.4.5 |

Per-script session information is in `results/sessionInfo_*.txt`.

GEO series matrices are parsed directly as text with base R rather than through
GEOquery, which does not build in this environment.

MSigDB gene-set membership changes between releases, and two enrichment p-values
in Table 3 shifted accordingly during reanalysis. The release version is therefore
pinned above and stated in Supplementary Table S1.

---

## Citation

If you use this code, please cite the paper above. The archived release carries
its own DOI: https://doi.org/10.5281/zenodo.22255746

## License

MIT. See `LICENSE`.

## Contact

Correspondence to Anand Prakash Singh. Questions about the code to Matthew Nho.
