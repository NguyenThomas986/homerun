#!/bin/bash
#SBATCH --job-name=csrna_map
#SBATCH --time=8:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=20 --mem=64G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ./config.env
module load "${CSRNA_CONDA_MODULE}"
source activate "${CSRNA_CONDA_ENV}"

export CSRNA_THREADS="${SLURM_CPUS_PER_TASK:-20}"
# CSRNA_LOG and CSRNA_RUN_LOG_DIR passed in via --export from submit_array.sh

ALN_TOOL="STAR"; [ "${CSRNA_ALIGNER}" = "hisat2" ] && ALN_TOOL="hisat2"
for t in homerTools "${ALN_TOOL}" python; do
    command -v "$t" >/dev/null || { echo "MISSING ON PATH: $t"; exit 1; }
done

echo "[$(date)] map task ${SLURM_ARRAY_TASK_ID} | aligner=${CSRNA_ALIGNER} | $(hostname)" | tee -a "${CSRNA_LOG}"
python -m csrnaseq --sample-index "${SLURM_ARRAY_TASK_ID}" --steps trim align --skip-prepare
echo "[$(date)] map task ${SLURM_ARRAY_TASK_ID} done" | tee -a "${CSRNA_LOG}"
