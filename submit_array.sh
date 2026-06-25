#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# Controller: submits the 3-phase job-array pipeline with dependencies.
#
#   prepare ──afterok──> align_array[0..N-1] ──afterok──> collect
#
# No config.env. All settings are command-line flags. Launch from anywhere:
#
#   /path/to/homerun/submit_array.sh \
#       --project /path/to/proj --partition kamiak --conda-env miniComputer \
#       --genome-index /path/to/STARIndex --genome hg38 \
#       [--conda-module anaconda3] [--aligner star|hisat2] [--throttle 16] \
#       [--email you@wsu.edu] [--copy-src '/src/*_R1*'] \
#       [-- <extra python flags, e.g. --trim-min 18 --star-filter-multimap 5000>]
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

usage() { sed -n '2,18p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

# ── Defaults ──────────────────────────────────────────────────────────────────
CONDA_MODULE="anaconda3"
ALIGNER="star"
THROTTLE="16"
EMAIL=""
COPY_SRC=""
PROJECT="" ; PARTITION="" ; CONDA_ENV="" ; GENOME_INDEX="" ; GENOME=""
EXTRA=()

# ── Parse args ────────────────────────────────────────────────────────────────
while [ $# -gt 0 ]; do
    case "$1" in
        --project)       PROJECT="$2"; shift 2 ;;
        --partition)     PARTITION="$2"; shift 2 ;;
        --conda-env)     CONDA_ENV="$2"; shift 2 ;;
        --conda-module)  CONDA_MODULE="$2"; shift 2 ;;
        --genome-index)  GENOME_INDEX="$2"; shift 2 ;;
        --genome)        GENOME="$2"; shift 2 ;;
        --aligner)       ALIGNER="$2"; shift 2 ;;
        --throttle)      THROTTLE="$2"; shift 2 ;;
        --email)         EMAIL="$2"; shift 2 ;;
        --copy-src)      COPY_SRC="$2"; shift 2 ;;
        -h|--help)       usage 0 ;;
        --)              shift; EXTRA=("$@"); break ;;
        *) echo "Unknown option: $1" >&2; usage 1 ;;
    esac
done

# ── Validate required ─────────────────────────────────────────────────────────
miss=""
[ -n "${PROJECT}" ]      || miss="${miss} --project"
[ -n "${PARTITION}" ]    || miss="${miss} --partition"
[ -n "${CONDA_ENV}" ]    || miss="${miss} --conda-env"
[ -n "${GENOME_INDEX}" ] || miss="${miss} --genome-index"
[ -n "${GENOME}" ]       || miss="${miss} --genome"
[ -z "${miss}" ] || { echo "ERROR: missing required:${miss}" >&2; usage 1; }

# Find the .sbatch job files (they sit next to this script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# Make `python -m csrnaseq` importable WITHOUT pip install: SCRIPT_DIR is the
# repo dir that contains the csrnaseq/ package, so it goes on PYTHONPATH here
# (for the login-node calls) and is forwarded to every job below.
export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH:-}"

# Activate env so the login-node python calls below work
module load "${CONDA_MODULE}" 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate "${CONDA_ENV}" 2>/dev/null || true
export PYTHONNOUSERSITE=1

if [ ! -d "${PROJECT}" ]; then
    echo "ERROR: --project does not exist: ${PROJECT}" >&2
    echo "  Create it first (this script will not fabricate the path)." >&2
    exit 1
fi
LOG_DIR="${PROJECT}/logs_slurm"
mkdir -p "${LOG_DIR}" "${PROJECT}/RawData"

# ── Plumbing args (positional) + python flags (forwarded to every phase) ──────
PLUMBING=( "${CONDA_MODULE}" "${CONDA_ENV}" "${PROJECT}" "${SCRIPT_DIR}" )
PY_ARGS=( --project "${PROJECT}" --aligner "${ALIGNER}"
          --genome-index "${GENOME_INDEX}" --genome "${GENOME}" )
[ -n "${COPY_SRC}" ] && PY_ARGS+=( --copy-src "${COPY_SRC}" )
[ ${#EXTRA[@]} -gt 0 ] && PY_ARGS+=( "${EXTRA[@]}" )

# Stage loose FASTQs, then count samples
python -m csrnaseq "${PY_ARGS[@]}" --stage-raw
N=$(python -m csrnaseq "${PY_ARGS[@]}" --count-samples)
if [ "${N}" -eq 0 ] && [ -n "${COPY_SRC}" ]; then
    echo "RawData empty — copying from ${COPY_SRC} ..."
    cp -r ${COPY_SRC} "${PROJECT}/RawData"/
    N=$(python -m csrnaseq "${PY_ARGS[@]}" --count-samples)
fi
[ "${N}" -ge 1 ] || { echo "ERROR: no *_R1* FASTQs in ${PROJECT}/RawData"; exit 1; }
echo "Found ${N} sample file(s) → array 0-$((N-1))"

# ── SLURM options ─────────────────────────────────────────────────────────────
SOPTS="--partition=${PARTITION} --mail-type=ALL"
[ -n "${EMAIL}" ] && SOPTS="${SOPTS} --mail-user=${EMAIL}"

PREP=$(sbatch --parsable ${SOPTS} \
       --output="${LOG_DIR}/prepare-%j.out" --error="${LOG_DIR}/prepare-%j.err" \
       prepare.sbatch "${PLUMBING[@]}" "${PY_ARGS[@]}")
ARRAY=$(sbatch --parsable ${SOPTS} --dependency=afterok:${PREP} \
        --output="${LOG_DIR}/align-%A_%a.out" --error="${LOG_DIR}/align-%A_%a.err" \
        --array=0-$((N-1))%"${THROTTLE}" \
        align_array.sbatch "${PLUMBING[@]}" "${PY_ARGS[@]}")
COLLECT=$(sbatch --parsable ${SOPTS} --dependency=afterok:${ARRAY} \
          --output="${LOG_DIR}/collect-%j.out" --error="${LOG_DIR}/collect-%j.err" \
          collect.sbatch "${PLUMBING[@]}" "${PY_ARGS[@]}")

echo "Submitted:"
echo "  prepare     = ${PREP}"
echo "  align_array = ${ARRAY}   (tasks 0-$((N-1)), <= ${THROTTLE} concurrent)"
echo "  collect     = ${COLLECT} (runs after all array tasks succeed)"
echo "Watch with: sq   |   logs in ${LOG_DIR}/"
