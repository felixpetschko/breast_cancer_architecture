#!/usr/bin/env bash
set -euo pipefail

MODE="${MODE:-all}"
R_ENV_NAME="${R_ENV_NAME:-/scratch/c9881013/.conda_envs/spacedeconv-env}"
PY_ENV_NAME="${PY_ENV_NAME:-/scratch/c9881013/.conda_envs/breast-cancer-rectangle-env}"
IMMUNEDECONV_ENV_NAME="${IMMUNEDECONV_ENV_NAME:-/scratch/c9881013/.conda_envs/immunedeconv-env}"
N_CPUS="${N_CPUS:-32}"

usage() {
  cat <<USAGE
Usage: scripts/run_bulk_cptac_brca.sh [--mode all|prepare|immunedeconv|rectangle|merge|export]
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="$2"
      shift 2
      ;;
    --n-cpus)
      N_CPUS="$2"
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
    if [[ -d "${env_name}" || "${env_name}" == */* ]]; then
      conda run --no-capture-output -p "${env_name}" Rscript "${script}" "$@"
    else
      conda run --no-capture-output -n "${env_name}" Rscript "${script}" "$@"
    fi
  else
    Rscript "${script}" "$@"
  fi
}

run_prepare() {
  run_rscript "${R_ENV_NAME}" scripts/prepare_bulk_cptac_brca.R
}

run_immunedeconv() {
  run_rscript "${IMMUNEDECONV_ENV_NAME}" scripts/run_bulk_immunedeconv.R \
    --input-csv results/bulk_rnaseq/cptac_brca/intermediate/cptac_brca_tpm_samples_by_genes.csv \
    --out-dir results/bulk_rnaseq/cptac_brca/objects \
    --diagnostics-json results/bulk_rnaseq/cptac_brca/objects/deconv_immunedeconv_cptac_brca_diagnostics.json
  mv results/bulk_rnaseq/cptac_brca/objects/deconv_estimate_tcga_brca.csv results/bulk_rnaseq/cptac_brca/objects/deconv_estimate_cptac_brca.csv
  mv results/bulk_rnaseq/cptac_brca/objects/deconv_quantiseq_tcga_brca.csv results/bulk_rnaseq/cptac_brca/objects/deconv_quantiseq_cptac_brca.csv
  mv results/bulk_rnaseq/cptac_brca/objects/deconv_epic_tcga_brca.csv results/bulk_rnaseq/cptac_brca/objects/deconv_epic_cptac_brca.csv
}

run_rectangle() {
  PY_ENV_NAME="${PY_ENV_NAME}" N_CPUS="${N_CPUS}" scripts/run_bulk_cptac_rectangle.sh
}

run_merge() {
  run_rscript "${R_ENV_NAME}" scripts/merge_bulk_cptac_results.R
}

run_export() {
  run_rscript "${R_ENV_NAME}" scripts/export_bulk_cptac_results.R
}

case "${MODE}" in
  all)
    run_prepare
    run_immunedeconv
    run_rectangle
    run_merge
    run_export
    ;;
  prepare)
    run_prepare
    ;;
  immunedeconv)
    run_immunedeconv
    ;;
  rectangle)
    run_rectangle
    ;;
  merge)
    run_merge
    ;;
  export)
    run_export
    ;;
  *)
    echo "Invalid mode: ${MODE}"
    usage
    exit 1
    ;;
esac

echo "CPTAC bulk pipeline complete."
