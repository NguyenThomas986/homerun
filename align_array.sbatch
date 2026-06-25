#!/bin/bash
#SBATCH --job-name=csrna_align
#SBATCH --time=8:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=64G
# Phase 1: one task per sample — trim + STAR for RawData R1 file #SLURM_ARRAY_TASK_ID.
# Submitted by submit_array.sh with the right --array range and --dependency.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ./config.env
module load "${CSRNA_CONDA_MODULE}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CSRNA_CONDA_ENV}"
export PYTHONNOUSERSITE=1
mkdir -p "${CSRNA_PROJECT}/logs_slurm"

export CSRNA_THREADS="${SLURM_CPUS_PER_TASK:-20}"
export CSRNA_LOG="${CSRNA_PROJECT}/logs_slurm/align_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

# Required tools for this phase (skewer only if a totalRNA sample lands here):
for t in homerTools STAR python; do
    command -v "$t" >/dev/null || { echo "MISSING ON PATH: $t"; exit 1; }
done

echo "[$(date)] align task ${SLURM_ARRAY_TASK_ID} on $(hostname), ${CSRNA_THREADS} threads"
python -m csrnaseq --sample-index "${SLURM_ARRAY_TASK_ID}" --steps trim align --skip-prepare
echo "[$(date)] align task ${SLURM_ARRAY_TASK_ID} done"
