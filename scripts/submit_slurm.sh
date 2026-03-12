#!/usr/bin/env bash
set -euo pipefail

MODE="all"
CPUS="16"
MEM="64G"
TIME="24:00:00"
JOB_NAME=""
ACCOUNT=""
PARTITION=""
QOS=""
N_CPUS=""
DRY_RUN="false"
TEST_PARTITION="false"

R_ENV_NAME="${R_ENV_NAME:-spacedeconv-env}"
PY_ENV_NAME="${PY_ENV_NAME:-breast-cancer-rectangle-env}"

usage() {
  cat <<USAGE
Usage: scripts/submit_slurm.sh [options]

Options:
  --mode <all|spacedeconv|rectangle>  Pipeline mode (default: all)
  --cpus <int>                        SLURM cpus-per-task (default: 16)
  --mem <size>                        SLURM memory (default: 64G)
  --time <HH:MM:SS>                   SLURM walltime (default: 24:00:00)
  --job-name <name>                   SLURM job name (default: breast-<mode>)
  --account <name>                    SLURM account (optional)
  --partition <name>                  SLURM partition (optional)
  --qos <name>                        SLURM qos (optional)
  --n-cpus <int>                      Rectangle worker CPUs inside job (optional)
  --test                              Use test partition (-p test) and max 10:00:00 walltime
  --dry-run                           Print sbatch command only
  -h, --help                          Show this help

Env passthrough (optional):
  CONFIG_PATH, REFERENCE_H5AD, CELL_TYPE_COL,
  GENE_SYMBOL_POLICY, CORRECT_MRNA_BIAS, MIN_GENES_OVERLAP,
  R_ENV_NAME, PY_ENV_NAME
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem) MEM="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --job-name) JOB_NAME="$2"; shift 2 ;;
    --account) ACCOUNT="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --qos) QOS="$2"; shift 2 ;;
    --n-cpus) N_CPUS="$2"; shift 2 ;;
    --test) TEST_PARTITION="true"; shift ;;
    --dry-run) DRY_RUN="true"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1"; usage; exit 1 ;;
  esac
done

if [[ "${MODE}" != "all" && "${MODE}" != "spacedeconv" && "${MODE}" != "rectangle" ]]; then
  echo "Invalid --mode: ${MODE}"
  exit 1
fi

if [[ -z "${JOB_NAME}" ]]; then
  JOB_NAME="breast-${MODE}"
fi

hms_to_seconds() {
  local hms="$1"
  IFS=: read -r hh mm ss <<< "${hms}"
  echo $((10#${hh} * 3600 + 10#${mm} * 60 + 10#${ss}))
}

if [[ "${TEST_PARTITION}" == "true" ]]; then
  PARTITION="test"
  max_test_seconds=$((10 * 3600))
  time_seconds=$(hms_to_seconds "${TIME}")
  if (( time_seconds > max_test_seconds )); then
    echo "For --test, --time must be <= 10:00:00. You provided ${TIME}."
    exit 1
  fi
fi

mkdir -p logs/slurm

exports=(
  "ALL"
  "MODE=${MODE}"
  "R_ENV_NAME=${R_ENV_NAME}"
  "PY_ENV_NAME=${PY_ENV_NAME}"
)

for var in CONFIG_PATH REFERENCE_H5AD CELL_TYPE_COL GENE_SYMBOL_POLICY CORRECT_MRNA_BIAS MIN_GENES_OVERLAP; do
  if [[ -n "${!var:-}" ]]; then
    exports+=("${var}=${!var}")
  fi
done

if [[ -n "${N_CPUS}" ]]; then
  exports+=("N_CPUS=${N_CPUS}")
fi

cmd=(
  sbatch
  --job-name "${JOB_NAME}"
  --time "${TIME}"
  --ntasks 1
  --cpus-per-task "${CPUS}"
  --mem "${MEM}"
  --output "logs/slurm/%x_%j.out"
  --error "logs/slurm/%x_%j.err"
  --export "$(IFS=,; echo "${exports[*]}")"
)

[[ -n "${ACCOUNT}" ]] && cmd+=(--account "${ACCOUNT}")
[[ -n "${PARTITION}" ]] && cmd+=(--partition "${PARTITION}")
[[ -n "${QOS}" ]] && cmd+=(--qos "${QOS}")

cmd+=(scripts/run_pipeline.slurm)

echo "Submitting with:"
printf ' %q' "${cmd[@]}"
echo

if [[ "${DRY_RUN}" == "true" ]]; then
  exit 0
fi

"${cmd[@]}"
