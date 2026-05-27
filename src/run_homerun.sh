#!/bin/bash
#SBATCH --job-name=homerun
#SBATCH --output=homerun_%j.out
#SBATCH --error=homerun_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=normal

# ─────────────────────────────────────────────────────────────────────────────
# Homerun SLURM submission script
#
# Usage:
#   sbatch run_homerun.sh                         # run all samples, all steps
#   sbatch run_homerun.sh --dry-run               # preview only
#   sbatch run_homerun.sh --steps trim star       # only trim + STAR steps
#   sbatch run_homerun.sh --retry-failed          # retry any FAILED samples
#   sbatch run_homerun.sh --samples THP1_rep1     # single sample
# ─────────────────────────────────────────────────────────────────────────────

# ── Environment ───────────────────────────────────────────────────────────────
# Activate your conda/mamba environment here:
# source activate homerun_env
# module load STAR/2.7.10a
# module load hisat2/2.2.1
# module load homer/4.11

# ── Configuration ─────────────────────────────────────────────────────────────
PIPELINE_CSV="homerun_pipeline.csv"   # path to your CSV
HOMERUN_DIR="$(dirname "$0")"         # directory containing Homerun.py
WORKING_PATH="$(pwd)"                 # project root (adjust if needed)

# ── Run ───────────────────────────────────────────────────────────────────────
echo "===== Homerun Pipeline ====="
echo "CSV:          $PIPELINE_CSV"
echo "Working path: $WORKING_PATH"
echo "CPUs:         $SLURM_CPUS_PER_TASK"
echo "Start time:   $(date)"
echo "============================"

python "${HOMERUN_DIR}/Homerun.py" \
    --csv         "${PIPELINE_CSV}" \
    --working-path "${WORKING_PATH}" \
    "$@"   # pass any extra CLI flags (--dry-run, --steps, etc.)

EXIT_CODE=$?
echo "============================"
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
exit $EXIT_CODE
