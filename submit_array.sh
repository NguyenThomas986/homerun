#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# Controller: submits the 3-phase job-array pipeline with dependencies.
#
#   prepare ──afterok──> align_array[0..N-1] ──afterok──> collect
#
# Run on the LOGIN node:   ./submit_array.sh
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd /weka/data/lab/duttke/personal/tnguyen/projects/TestRun4
source ./config.env

# Activate env so `python -m csrnaseq --count-samples` works on the login node
module load "${CSRNA_CONDA_MODULE}" 2>/dev/null || true
source activate "${CSRNA_CONDA_ENV}" 2>/dev/null || true
export PYTHONNOUSERSITE=1

: "${CSRNA_PROJECT:?Set CSRNA_PROJECT in config.env}"
mkdir -p logs_slurm "${CSRNA_PROJECT}/RawData"

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

PREP=$(sbatch --parsable ${P} prepare.sbatch)
ARRAY=$(sbatch --parsable ${P} --dependency=afterok:${PREP} \
               --array=0-$((N-1))%"${CSRNA_ARRAY_THROTTLE}" align_array.sbatch)
COLLECT=$(sbatch --parsable ${P} --dependency=afterok:${ARRAY} collect.sbatch)

echo "Submitted:"
echo "  prepare     = ${PREP}"
echo "  align_array = ${ARRAY}   (tasks 0-$((N-1)), <= ${CSRNA_ARRAY_THROTTLE} concurrent)"
echo "  collect     = ${COLLECT} (runs after all array tasks succeed)"
echo "Watch with: sq   |   logs in ./logs_slurm/"
