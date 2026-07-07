#!/bin/bash
#SBATCH --job-name=csrna_tagdir
#SBATCH --time=4:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=32G
# Phase 2: one task per sample — build ONE leaf TagDir for RawData R1 #SLURM_ARRAY_TASK_ID.
# Same indexing as align_array.sbatch (one R1 file = one leaf run), so task N
# here builds the TagDir for the SAM file task N of align_array produced.
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
LOG="${PROJECT}/logs_slurm/tagdir_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

for t in makeTagDirectory python; do
    command -v "$t" >/dev/null || { echo "MISSING ON PATH: $t"; exit 1; }
done

echo "[$(date)] tagdir task ${SLURM_ARRAY_TASK_ID} on $(hostname), ${SLURM_CPUS_PER_TASK:-?} threads"
python -m csrnaseq --sample-index "${SLURM_ARRAY_TASK_ID}" --steps tagdirs \
       --skip-prepare --log-path "${LOG}" "$@"
echo "[$(date)] tagdir task ${SLURM_ARRAY_TASK_ID} done"
