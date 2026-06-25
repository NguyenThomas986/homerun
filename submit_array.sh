#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# Controller: submits the 3-phase job-array pipeline with dependencies.
#
#   prepare ──afterok──> align_array[0..N-1] ──afterok──> collect
#
# config.env lives NEXT TO this script (in homerun/). Launch from anywhere:
#     /path/to/homerun/submit_array.sh
# RawData / outputs live wherever CSRNA_PROJECT points (separate from homerun/).
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

# Anchor to this script's own directory so config.env + the .sbatch jobs are
# found here no matter where you launch from. cd-ing here (not just sourcing an
# absolute path) also makes SLURM_SUBMIT_DIR = homerun/ for the submitted jobs,
# so their `cd "$SLURM_SUBMIT_DIR"; source ./config.env` resolves to this file too.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

if [ ! -f ./config.env ]; then
    echo "ERROR: no config.env next to submit_array.sh (${SCRIPT_DIR})"
    echo "  Copy config.env.example → config.env in that folder and edit it."
    exit 1
fi
source ./config.env

# Activate env so `python -m csrnaseq --count-samples` works on the login node
module load "${CSRNA_CONDA_MODULE}" 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate "${CSRNA_CONDA_ENV}" 2>/dev/null || true
export PYTHONNOUSERSITE=1

: "${CSRNA_PROJECT:?Set CSRNA_PROJECT in config.env}"

# Guard: refuse to run if CSRNA_PROJECT doesn't exist, so mkdir -p can't
# fabricate a stale cluster path (e.g. /weka/...) on a machine that lacks it.
if [ ! -d "${CSRNA_PROJECT}" ]; then
    echo "ERROR: CSRNA_PROJECT does not exist: ${CSRNA_PROJECT}"
    echo "  Point CSRNA_PROJECT at your real project dir and create it first."
    exit 1
fi
LOG_DIR="${CSRNA_PROJECT}/logs_slurm"
mkdir -p "${LOG_DIR}" "${CSRNA_PROJECT}/RawData"

# Move any loose *_R1*/*_R2* FASTQs from the project root into RawData/ first,
# so --count-samples below sees them.
python -m csrnaseq --stage-raw

# Stage raw data if RawData is empty and a copy source is configured
N=$(python -m csrnaseq --count-samples)
if [ "${N}" -eq 0 ] && [ -n "${CSRNA_COPY_SRC}" ]; then
    echo "RawData empty — copying from CSRNA_COPY_SRC ..."
    cp -r ${CSRNA_COPY_SRC} "${CSRNA_PROJECT}/RawData"/
    N=$(python -m csrnaseq --count-samples)
fi
[ "${N}" -ge 1 ] || { echo "ERROR: no *_R1* FASTQs in ${CSRNA_PROJECT}/RawData"; exit 1; }
echo "Found ${N} sample file(s) → array 0-$((N-1))"

P="--partition=${CSRNA_PARTITION} --mail-user=${CSRNA_EMAIL} --mail-type=ALL"

PREP=$(sbatch --parsable ${P} \
       --output="${LOG_DIR}/prepare-%j.out" --error="${LOG_DIR}/prepare-%j.err" \
       prepare.sbatch)
ARRAY=$(sbatch --parsable ${P} --dependency=afterok:${PREP} \
        --output="${LOG_DIR}/align-%A_%a.out" --error="${LOG_DIR}/align-%A_%a.err" \
        --array=0-$((N-1))%"${CSRNA_ARRAY_THROTTLE}" align_array.sbatch)
COLLECT=$(sbatch --parsable ${P} --dependency=afterok:${ARRAY} \
          --output="${LOG_DIR}/collect-%j.out" --error="${LOG_DIR}/collect-%j.err" \
          collect.sbatch)

echo "Submitted:"
echo "  prepare     = ${PREP}"
echo "  align_array = ${ARRAY}   (tasks 0-$((N-1)), <= ${CSRNA_ARRAY_THROTTLE} concurrent)"
echo "  collect     = ${COLLECT} (runs after all array tasks succeed)"
echo "Watch with: sq   |   logs in ${LOG_DIR}/"
