#!/bin/bash
#SBATCH --job-name=cloz_fig1
#SBATCH --account=lsa1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/home/%u/clozapine/logs/fig1_%j.out
#SBATCH --error=/home/%u/clozapine/logs/fig1_%j.err

set -euo pipefail
export CLOZ_HOME="${HOME}/clozapine"
export R_LIBS_USER="${HOME}/Rlibs/R-4.6.0"
mkdir -p "${CLOZ_HOME}/logs" "${CLOZ_HOME}/results" "${CLOZ_HOME}/figures"

echo "Job ${SLURM_JOB_ID} on $(hostname), started $(date)"
module purge
module load R/4.6.0
echo "Rscript: $(which Rscript)"

for f in GSE262419_hash.csv.gz GSE262419_Plate15.csv.gz GSE262419_Plate16.csv.gz; do
  [ -s "${CLOZ_HOME}/data/${f}" ] || { echo "MISSING: data/${f}"; exit 1; }
done

cd "${CLOZ_HOME}"
Rscript --vanilla 01_figure1.R
echo "Finished: $(date)"
