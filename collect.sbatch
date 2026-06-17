#!/bin/bash
#SBATCH --job-name=csrna_collect
#SBATCH --output=logs_slurm/collect_%j.out
#SBATCH --error=logs_slurm/collect_%j.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=32G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ./config.env
module load "${CSRNA_CONDA_MODULE}"
source activate "${CSRNA_CONDA_ENV}"
export PYTHONNOUSERSITE=1
mkdir -p logs_slurm

export CSRNA_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export CSRNA_LOG="${SLURM_SUBMIT_DIR}/logs_slurm/collect_${SLURM_JOB_ID}.log"

for t in homerTools python; do
    command -v "$t" >/dev/null || { echo "MISSING ON PATH: $t"; exit 1; }
done

echo "[$(date)] collect on $(hostname)"
python -m csrnaseq --steps tagdirs bedgraphs tss qc stability report
echo "[$(date)] collect done"
