#!/usr/bin/env bash
set -euo pipefail

MODE="${MODE:-all}"
R_ENV_NAME="${R_ENV_NAME:-spacedeconv-env}"
PY_ENV_NAME="${PY_ENV_NAME:-breast-cancer-rectangle-env}"
IMMUNEDECONV_ENV_NAME="${IMMUNEDECONV_ENV_NAME:-immunedeconv-env}"

usage() {
  cat <<USAGE
Usage: scripts/run_bulk_tcga_brca.sh [--mode all|prepare|immunedeconv|rectangle|merge|report]
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1"
      usage
      exit 1
      ;;
  esac
done

run_rscript() {
  local env_name="$1"
  local script="$2"
  shift 2
  echo ">>> [R:${env_name}] ${script} $*"
  if [[ -n "${env_name}" ]]; then
    conda run --no-capture-output -n "${env_name}" Rscript "${script}" "$@"
  else
    Rscript "${script}" "$@"
  fi
}

case "${MODE}" in
  all)
    run_rscript "${R_ENV_NAME}" scripts/prepare_bulk_tcga_brca.R
    run_rscript "${IMMUNEDECONV_ENV_NAME}" scripts/run_bulk_immunedeconv.R
    PY_ENV_NAME="${PY_ENV_NAME}" scripts/run_bulk_rectangle.sh
    run_rscript "${R_ENV_NAME}" scripts/merge_bulk_results.R
    run_rscript "${R_ENV_NAME}" scripts/render_bulk_report.R
    ;;
  prepare)
    run_rscript "${R_ENV_NAME}" scripts/prepare_bulk_tcga_brca.R
    ;;
  immunedeconv)
    run_rscript "${IMMUNEDECONV_ENV_NAME}" scripts/run_bulk_immunedeconv.R
    ;;
  rectangle)
    PY_ENV_NAME="${PY_ENV_NAME}" scripts/run_bulk_rectangle.sh
    ;;
  merge)
    run_rscript "${R_ENV_NAME}" scripts/merge_bulk_results.R
    ;;
  report)
    run_rscript "${R_ENV_NAME}" scripts/render_bulk_report.R
    ;;
  *)
    echo "Invalid mode: ${MODE}"
    exit 1
    ;;
esac

echo "Bulk pipeline complete."
