#!/bin/bash
#SBATCH --job-name=csrna_tss
#SBATCH --time=4:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=16G
# Phase 6: one task per Species/Sample — call TSRs with findcsRNATSS.pl for
# that one Species/Sample. Indexed by --group-index into
# utils.list_samples(cfg). Needs tagdirs_combo_array (csRNA-combo/sRNA-combo/
# totalRNA-combo TagDirs) to have finished first.
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
LOG="${PROJECT}/logs_slurm/tss_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

for t in findcsRNATSS.pl python; do
    command -v "$t" >/dev/null || { echo "MISSING ON PATH: $t"; exit 1; }
done

echo "[$(date)] tss task ${SLURM_ARRAY_TASK_ID} on $(hostname), ${SLURM_CPUS_PER_TASK:-?} threads"
python -m csrnaseq --group-index "${SLURM_ARRAY_TASK_ID}" --steps tss \
       --skip-prepare --log-path "${LOG}" "$@"
echo "[$(date)] tss task ${SLURM_ARRAY_TASK_ID} done"
