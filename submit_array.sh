#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# Controller: submits the 3-phase job-array pipeline with dependencies.
#
#   prepare ──afterok──> align_array[0..N-1] ──afterok──> collect
#
# Run on the LOGIN node:   ./submit_array.sh
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}" 
source "${SLURM_SUBMIT_DIR}/homerun/config.env"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CSRNA_CONDA_ENV}"

# Activate env so `python -m csrnaseq --count-samples` works on the login node
module load "${CSRNA_CONDA_MODULE}" 2>/dev/null || true
source activate "${CSRNA_CONDA_ENV}" 2>/dev/null || true
export PYTHONNOUSERSITE=1

: "${CSRNA_PROJECT:?Set CSRNA_PROJECT in config.env}"

# ── One timestamped log dir for this entire run ───────────────────────────────
RUN_TS=$(date +%Y%m%d_%H%M%S)
LOG_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}/logs_slurm/${RUN_TS}"
mkdir -p "${LOG_DIR}" "${CSRNA_PROJECT}/RawData"
export CSRNA_LOG="${LOG_DIR}/pipeline.log"
echo "[$(date)] run ${RUN_TS} — logs in ${LOG_DIR}" | tee -a "${CSRNA_LOG}"

# Stage raw data if RawData is empty and a copy source is configured
N=$(python -m csrnaseq --count-samples)
if [ "${N}" -eq 0 ] && [ -n "${CSRNA_COPY_SRC}" ]; then
    echo "RawData empty — copying from CSRNA_COPY_SRC ..."
    cp -r ${CSRNA_COPY_SRC} "${CSRNA_PROJECT}/RawData"/
    N=$(python -m csrnaseq --count-samples)
fi
[ "${N}" -ge 1 ] || { echo "ERROR: no *_R1* FASTQs in ${CSRNA_PROJECT}/RawData"; exit 1; }

# Print sample names
echo "Found ${N} sample file(s) → array 0-$((N-1))"
python -m csrnaseq --list-samples | while IFS= read -r line; do
    echo "  [${line%%:*}] ${line#*: }"
done
echo ""

P="--partition=${CSRNA_PARTITION}"
[ -n "${CSRNA_EMAIL}" ] && P="${P} --mail-user=${CSRNA_EMAIL} --mail-type=ALL"

# Pass the shared log dir and log file to every job via --export
E="ALL,CSRNA_RUN_LOG_DIR=${LOG_DIR},CSRNA_LOG=${CSRNA_LOG}"

PREP=$(sbatch --parsable ${P} --export="${E}" \
       --output="${LOG_DIR}/prepare.out" --error="${LOG_DIR}/prepare.err" \
       prepare.sbatch)

ARRAY=$(sbatch --parsable ${P} --export="${E}" \
        --output="${LOG_DIR}/align_%a.out" --error="${LOG_DIR}/align_%a.err" \
        --dependency=afterok:${PREP} \
        --array=0-$((N-1))%"${CSRNA_ARRAY_THROTTLE}" align_array.sbatch)

COLLECT=$(sbatch --parsable ${P} --export="${E}" \
          --output="${LOG_DIR}/collect.out" --error="${LOG_DIR}/collect.err" \
          --dependency=afterok:${ARRAY} collect.sbatch)

echo "Submitted:"
echo "  prepare     = ${PREP}"
echo "  align_array = ${ARRAY}   (tasks 0-$((N-1)), <= ${CSRNA_ARRAY_THROTTLE} concurrent)"
echo "  collect     = ${COLLECT} (runs after all array tasks succeed)"
echo ""
echo "All logs → ${LOG_DIR}/"
echo "  ${LOG_DIR}/pipeline.log   ← combined Python log"
echo "  ${LOG_DIR}/prepare.out"
echo "  ${LOG_DIR}/align_0.out .. align_$((N-1)).out"
echo "  ${LOG_DIR}/collect.out"
echo ""
echo "Watch with: sq"
