#!/bin/bash

set -euo pipefail

# PYTHONPATH set up front (not just after arg parsing) so --help can shell
# out to `python -m csrnaseq --help` and show every flag valid after `--`.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH:-}"

usage() {
    python -m csrnaseq --help 2>/dev/null || echo "(activate your conda env for the csrnaseq flag list)"
    exit "${1:-0}"
}

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
# (SCRIPT_DIR/PYTHONPATH already set near the top of the script, before arg
# parsing, so --help can shell out to `python -m csrnaseq --help`.)
cd "${SCRIPT_DIR}"

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
mkdir -p "${LOG_DIR}"

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
    echo "RawData empty — running prepare now to copy from ${COPY_SRC} ..."
    # NOTE: --only-prepare also runs ensure_starindex(); if --genome-index
    # already points at an existing, non-empty directory (the normal case)
    # this is a no-op, but if it's missing AND CSRNA_STARINDEX_URL is set,
    # this will download the STARIndex tarball on the login node.
    python -m csrnaseq "${PY_ARGS[@]}" --only-prepare
    N=$(python -m csrnaseq "${PY_ARGS[@]}" --count-samples)
fi
[ "${N}" -ge 1 ] || { echo "ERROR: no *_R1* FASTQs in ${PROJECT}/RawData"; exit 1; }
echo "Found ${N} sample file(s) → array 0-$((N-1))"

S=$(python -m csrnaseq "${PY_ARGS[@]}" --count-groups)
[ "${S}" -ge 1 ] || { echo "ERROR: no Species/Sample groups found under ${PROJECT}"; exit 1; }
echo "Found ${S} Species/Sample group(s) → array 0-$((S-1))"

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
TAGDIR=$(sbatch --parsable ${SOPTS} --dependency=aftercorr:${ARRAY} \
        --output="${LOG_DIR}/tagdir-%A_%a.out" --error="${LOG_DIR}/tagdir-%A_%a.err" \
        --array=0-$((N-1))%"${THROTTLE}" \
        tagdir_array.sbatch "${PLUMBING[@]}" "${PY_ARGS[@]}")
TAGDIR_COMBO=$(sbatch --parsable ${SOPTS} --dependency=afterany:${ARRAY} \
        --output="${LOG_DIR}/tagdircombo-%A_%a.out" --error="${LOG_DIR}/tagdircombo-%A_%a.err" \
        --array=0-$((S-1))%"${THROTTLE}" \
        tagdirs_combo_array.sbatch "${PLUMBING[@]}" "${PY_ARGS[@]}")
TSS=$(sbatch --parsable ${SOPTS} --dependency=afterany:${TAGDIR_COMBO} \
        --output="${LOG_DIR}/tss-%A_%a.out" --error="${LOG_DIR}/tss-%A_%a.err" \
        --array=0-$((S-1))%"${THROTTLE}" \
        tss_array.sbatch "${PLUMBING[@]}" "${PY_ARGS[@]}")
BEDGRAPH=$(sbatch --parsable ${SOPTS} --dependency=afterany:${TAGDIR}:${TAGDIR_COMBO} \
        --output="${LOG_DIR}/bedgraphs-%A_%a.out" --error="${LOG_DIR}/bedgraphs-%A_%a.err" \
        --array=0-$((S-1))%"${THROTTLE}" \
        bedgraphs_array.sbatch "${PLUMBING[@]}" "${PY_ARGS[@]}")
COLLECT=$(sbatch --parsable ${SOPTS} --dependency=afterok:${TSS}:${BEDGRAPH} \
          --output="${LOG_DIR}/collect-%j.out" --error="${LOG_DIR}/collect-%j.err" \
          collect.sbatch "${PLUMBING[@]}" "${PY_ARGS[@]}")

echo "Submitted:"
echo "  prepare             = ${PREP}"
echo "  align_array         = ${ARRAY}         (tasks 0-$((N-1)), <= ${THROTTLE} concurrent)"
echo "  tagdir_array        = ${TAGDIR}        (tasks 0-$((N-1)), <= ${THROTTLE} concurrent)"
echo "  tagdirs_combo_array = ${TAGDIR_COMBO}  (tasks 0-$((S-1)), <= ${THROTTLE} concurrent)"
echo "  tss_array           = ${TSS}           (tasks 0-$((S-1)), <= ${THROTTLE} concurrent)"
echo "  bedgraphs_array     = ${BEDGRAPH}      (tasks 0-$((S-1)), <= ${THROTTLE} concurrent)"
echo "  collect             = ${COLLECT} (runs after tss_array and bedgraphs_array succeed)"
echo "Watch with: sq   |   logs in ${LOG_DIR}/"