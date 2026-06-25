#!/bin/bash
#SBATCH --job-name=csrna_collect
#SBATCH --time=4:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=32G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source "${SLURM_SUBMIT_DIR}/homerun/config.env"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CSRNA_CONDA_ENV}"
export PYTHONNOUSERSITE=1

export CSRNA_THREADS="${SLURM_CPUS_PER_TASK:-8}"
# CSRNA_LOG and CSRNA_RUN_LOG_DIR passed in via --export from submit_array.sh

for t in homerTools python; do
    command -v "$t" >/dev/null || { echo "MISSING ON PATH: $t"; exit 1; }
done

echo "[$(date)] collect on $(hostname)" | tee -a "${CSRNA_LOG}"
python -m csrnaseq --steps tagdirs bedgraphs tss qc stability report
echo "[$(date)] collect done" | tee -a "${CSRNA_LOG}"
