#!/bin/bash
#SBATCH --job-name=csrna_align
#SBATCH --time=8:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=64G
# Phase 1: one task per sample — trim + align for RawData R1 #SLURM_ARRAY_TASK_ID.
set -euo pipefail

# Args from submit_array.sh: <conda_module> <conda_env> <project> [python flags...]
CONDA_MODULE="$1"; CONDA_ENV="$2"; PROJECT="$3"; PKG_DIR="$4"; shift 4

# Make `python -m csrnaseq` importable without pip install (repo dir on path)
export PYTHONPATH="${PKG_DIR}:${PYTHONPATH:-}"

module load "${CONDA_MODULE}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
export PYTHONNOUSERSITE=1
mkdir -p "${PROJECT}/logs_slurm"
LOG="${PROJECT}/logs_slurm/align_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

# homerTools + python must be present; python validates the chosen aligner/index.
for t in homerTools python; do
    command -v "$t" >/dev/null || { echo "MISSING ON PATH: $t"; exit 1; }
done

echo "[$(date)] align task ${SLURM_ARRAY_TASK_ID} on $(hostname), ${SLURM_CPUS_PER_TASK:-?} threads"
python -m csrnaseq --sample-index "${SLURM_ARRAY_TASK_ID}" --steps trim align \
       --skip-prepare --log-path "${LOG}" "$@"
echo "[$(date)] align task ${SLURM_ARRAY_TASK_ID} done"
