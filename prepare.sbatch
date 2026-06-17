#!/bin/bash
#SBATCH --job-name=csrna_prepare
#SBATCH --output=logs_slurm/prepare_%j.out
#SBATCH --error=logs_slurm/prepare_%j.err
#SBATCH --time=2:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ./config.env
module load "${CSRNA_CONDA_MODULE}"
source activate "${CSRNA_CONDA_ENV}"
export PYTHONNOUSERSITE=1
mkdir -p logs_slurm

export CSRNA_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export CSRNA_LOG="${SLURM_SUBMIT_DIR}/logs_slurm/prepare_${SLURM_JOB_ID}.log"

echo "[$(date)] prepare on $(hostname)"
python -m csrnaseq --only-prepare
echo "[$(date)] prepare done"

