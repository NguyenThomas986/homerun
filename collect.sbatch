#!/bin/bash
#SBATCH --job-name=csrna_collect
#SBATCH --time=4:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=32G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ./config.env
module load "${CSRNA_CONDA_MODULE}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CSRNA_CONDA_ENV}"
export PYTHONNOUSERSITE=1
mkdir -p "${CSRNA_PROJECT}/logs_slurm"

export CSRNA_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export CSRNA_LOG="${CSRNA_PROJECT}/logs_slurm/collect_${SLURM_JOB_ID}.log"

for t in homerTools python; do
    command -v "$t" >/dev/null || { echo "MISSING ON PATH: $t"; exit 1; }
done

echo "[$(date)] collect on $(hostname)"
python -m csrnaseq --steps tagdirs bedgraphs tss qc stability report
echo "[$(date)] collect done"
