#!/bin/bash
#SBATCH --job-name=csrna_tagdircombo
#SBATCH --time=4:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=5GB
# Phase 3: one task per Species/Sample — build that sample's combo TagDir(s)
# (one per assay: csRNA-combo, sRNA-combo, totalRNA-combo). Indexed by
# --group-index into utils.list_samples(cfg), NOT by R1 file — a different
# indexing scheme from align_array/tagdir_array. Only needs the align phase
# to have finished (doesn't depend on tagdir_array's leaf TagDirs).
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
LOG="${PROJECT}/logs_slurm/tagdircombo_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

for t in makeTagDirectory python; do
    command -v "$t" >/dev/null || { echo "MISSING ON PATH: $t"; exit 1; }
done

echo "[$(date)] tagdirs-combo task ${SLURM_ARRAY_TASK_ID} on $(hostname), ${SLURM_CPUS_PER_TASK:-?} threads"
python -m csrnaseq --group-index "${SLURM_ARRAY_TASK_ID}" --steps tagdirs-combo \
       --skip-prepare --log-path "${LOG}" "$@"
echo "[$(date)] tagdirs-combo task ${SLURM_ARRAY_TASK_ID} done"
