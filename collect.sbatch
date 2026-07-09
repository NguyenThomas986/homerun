#!/bin/bash
#SBATCH --job-name=csrna_collect
#SBATCH --time=4:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=32G
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
LOG="${PROJECT}/logs_slurm/collect_${SLURM_JOB_ID}.log"

for t in homerTools python; do
    command -v "$t" >/dev/null || { echo "MISSING ON PATH: $t"; exit 1; }
done

echo "[$(date)] collect on $(hostname)"
python -m csrnaseq --steps tagdirs-combo bedgraphs tss ritrie qc stability report \
       --log-path "${LOG}" "$@"

echo "[$(date)] cleaning up intermediate files"
find "${PROJECT}" -mindepth 4 -maxdepth 4 -type d \( -name Trimmed -o -name Aligned \) -exec rm -rf {} +
echo "[$(date)] collect done"
