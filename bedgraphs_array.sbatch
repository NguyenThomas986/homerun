#!/bin/bash
#SBATCH --job-name=csrna_bedgraphs
#SBATCH --time=2:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=1GB
# Phase 5: one task per Species/Sample — build bedGraphs for that sample's
# leaf and combo TagDirs under Species/TagDirs/. Indexed by --group-index
# into utils.list_samples(cfg). Needs BOTH tagdir_array (leaf TagDirs) and
# tagdirs_combo_array (combo TagDirs) to have finished first.
set -euo pipefail

# Args from submit_array.sh: <conda_module> <conda_env> <project> [python flags...]
CONDA_MODULE="$1"; CONDA_ENV="$2"; PROJECT="$3"; PKG_DIR="$4"; shift 4

# Make `python -m homerun` importable without pip install (repo dir on path)
export PYTHONPATH="${PKG_DIR}/src:${PYTHONPATH:-}"

module load "${CONDA_MODULE}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
export PYTHONNOUSERSITE=1
mkdir -p "${PROJECT}/logs_slurm"
LOG="${PROJECT}/logs_slurm/bedgraphs_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

for t in makeUCSCfile python; do
    command -v "$t" >/dev/null || { echo "MISSING ON PATH: $t"; exit 1; }
done

echo "[$(date)] bedgraphs task ${SLURM_ARRAY_TASK_ID} on $(hostname), ${SLURM_CPUS_PER_TASK:-?} threads"
python -m homerun --group-index "${SLURM_ARRAY_TASK_ID}" --steps bedgraphs \
       --skip-prepare --log-path "${LOG}" "$@"
echo "[$(date)] bedgraphs task ${SLURM_ARRAY_TASK_ID} done"
