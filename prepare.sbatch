#!/bin/bash
#SBATCH --job-name=csrna_prepare
#SBATCH --time=2:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source "${SLURM_SUBMIT_DIR}/homerun/config.env"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CSRNA_CONDA_ENV}"
export PYTHONNOUSERSITE=1

export CSRNA_THREADS="${SLURM_CPUS_PER_TASK:-4}"
# CSRNA_LOG and CSRNA_RUN_LOG_DIR are passed in via --export from submit_array.sh

echo "[$(date)] prepare on $(hostname)" | tee -a "${CSRNA_LOG}"
python -m csrnaseq --only-prepare
echo "[$(date)] prepare done" | tee -a "${CSRNA_LOG}"
