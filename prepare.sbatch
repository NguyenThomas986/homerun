#!/bin/bash
#SBATCH --job-name=csrna_prepare
#SBATCH --time=2:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16G
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
LOG="${PROJECT}/logs_slurm/prepare_${SLURM_JOB_ID}.log"

echo "[$(date)] prepare on $(hostname)"
python -m csrnaseq --only-prepare --log-path "${LOG}" "$@"
echo "[$(date)] prepare done"
