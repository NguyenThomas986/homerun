#!/bin/bash
#SBATCH --job-name=csrna_map
#SBATCH --output=logs_slurm/map_%A_%a.out
#SBATCH --error=logs_slurm/map_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=64G
# Phase 1: one task per sample — trim + mapping (STAR or HISAT2) for R1 #TASK_ID.
# Aligner + genome index come from config.env (CSRNA_ALIGNER / CSRNA_GENOME_INDEX).
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source "${SLURM_SUBMIT_DIR}/homerun/config.env"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CSRNA_CONDA_ENV}"
mkdir -p logs_slurm

export CSRNA_THREADS="${SLURM_CPUS_PER_TASK:-20}"
export CSRNA_LOG="${SLURM_SUBMIT_DIR}/logs_slurm/map_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

# Aligner tool must be on PATH (STAR or hisat2 depending on CSRNA_ALIGNER)
ALN_TOOL="STAR"; [ "${CSRNA_ALIGNER}" = "hisat2" ] && ALN_TOOL="hisat2"
for t in homerTools "${ALN_TOOL}" python; do
    command -v "$t" >/dev/null || { echo "MISSING ON PATH: $t"; exit 1; }
done

echo "[$(date)] map task ${SLURM_ARRAY_TASK_ID} | aligner=${CSRNA_ALIGNER} | $(hostname)"
python -m csrnaseq --sample-index "${SLURM_ARRAY_TASK_ID}" --steps trim mapping --skip-prepare
