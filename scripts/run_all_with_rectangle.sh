#!/usr/bin/env bash
set -euo pipefail

slides=(andersson 1142243F 1160920F 4465 44971 4290 4535)

R_ENV_NAME="${R_ENV_NAME:-spacedeconv-env}"
PY_ENV_NAME="${PY_ENV_NAME:-breast-cancer-rectangle-env}"
MODE="${MODE:-all}"
CONFIG_PATH="${CONFIG_PATH:-config/rectangle.yaml}"
RECTANGLE_PY="${RECTANGLE_PY:-scripts/run_rectangle.py}"
REFERENCE_H5AD="${REFERENCE_H5AD:-}"
CELL_TYPE_COL="${CELL_TYPE_COL:-}"
GENE_SYMBOL_POLICY="${GENE_SYMBOL_POLICY:-upper}"
CORRECT_MRNA_BIAS="${CORRECT_MRNA_BIAS:-true}"
MIN_GENES_OVERLAP="${MIN_GENES_OVERLAP:-500}"
N_CPUS="${N_CPUS:-}"

usage() {
  cat <<USAGE
Usage: scripts/run_all_with_rectangle.sh [--mode all|spacedeconv|rectangle]

Modes:
  all          Run baseline + rectangle + merge + report (default)
  spacedeconv  Run baseline + report only
  rectangle    Run rectangle export/deconvolution; merge+report only if baseline outputs exist

Environment overrides:
  R_ENV_NAME, PY_ENV_NAME, CONFIG_PATH, REFERENCE_H5AD, CELL_TYPE_COL,
  GENE_SYMBOL_POLICY, CORRECT_MRNA_BIAS, MIN_GENES_OVERLAP, N_CPUS
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="$2"
      shift 2
      ;;
    --help|-h)
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

if [[ "${MODE}" != "all" && "${MODE}" != "spacedeconv" && "${MODE}" != "rectangle" ]]; then
  echo "Invalid mode: ${MODE}. Use all|spacedeconv|rectangle"
  exit 1
fi

run_rscript() {
  local script="$1"
  shift
  echo ">>> [R] ${script} $*"
  if [[ -n "${R_ENV_NAME}" ]]; then
    conda run --no-capture-output -n "${R_ENV_NAME}" Rscript "${script}" "$@"
  else
    Rscript "${script}" "$@"
  fi
}

run_python() {
  local script="$1"
  shift
  echo ">>> [PY] ${script} $*"
  if [[ -n "${PY_ENV_NAME}" ]]; then
    conda run --no-capture-output -n "${PY_ENV_NAME}" python "${script}" "$@"
  else
    python3 "${script}" "$@"
  fi
}

run_python_quiet() {
  if [[ -n "${PY_ENV_NAME}" ]]; then
    conda run --no-capture-output -n "${PY_ENV_NAME}" python "$@"
  else
    python3 "$@"
  fi
}

read_cfg_value() {
  local key="$1"
  run_python_quiet -c 'import sys,yaml; from pathlib import Path; p=Path(sys.argv[1]); key=sys.argv[2]; cfg=yaml.safe_load(p.read_text()) if p.exists() else {}; print((cfg or {}).get(key,""))' "${CONFIG_PATH}" "${key}"
}

have_baseline_outputs() {
  [[ -f "results/objects/deconv_estimate_andersson.rds" && -f "results/tables/targets_andersson.csv" ]]
}

resolve_rectangle_config() {
  if [[ -z "${REFERENCE_H5AD}" ]]; then
    REFERENCE_H5AD="$(read_cfg_value reference_path)"
  fi

  if [[ -z "${CELL_TYPE_COL}" ]]; then
    CELL_TYPE_COL="$(read_cfg_value cell_type_column)"
    if [[ -z "${CELL_TYPE_COL}" ]]; then
      CELL_TYPE_COL="cell_type"
    fi
  fi

  if [[ -z "${REFERENCE_H5AD}" ]]; then
    echo "Rectangle reference path is empty. Set config/rectangle.yaml reference_path or REFERENCE_H5AD env."
    exit 1
  fi
}

run_baseline() {
  echo "[Baseline] Running spacedeconv"
  run_rscript scripts/run_spacedeconv.R
}

run_rectangle_stage() {
  echo "[Rectangle] Exporting spatial CPM matrices"
  run_rscript scripts/export_spatial_for_rectangle.R --out-dir results/intermediate/rectangle --gene-symbol-policy "${GENE_SYMBOL_POLICY}"

  echo "[Rectangle] Running deconvolution per slide"
  for slide in "${slides[@]}"; do
    bulk_csv="results/intermediate/rectangle/${slide}_bulks_cpm.csv"
    out_csv="results/objects/deconv_rectangle_${slide}.csv"
    diag_json="results/objects/deconv_rectangle_${slide}_diagnostics.json"

    args=(
      --config "${CONFIG_PATH}"
      --reference-h5ad "${REFERENCE_H5AD}"
      --bulks-csv "${bulk_csv}"
      --cell-type-col "${CELL_TYPE_COL}"
      --gene-symbol-policy "${GENE_SYMBOL_POLICY}"
      --min-genes-overlap "${MIN_GENES_OVERLAP}"
      --out-csv "${out_csv}"
      --diagnostics-json "${diag_json}"
    )

    if [[ "${CORRECT_MRNA_BIAS}" == "true" ]]; then
      args+=(--correct-mrna-bias)
    else
      args+=(--no-correct-mrna-bias)
    fi

    if [[ -n "${N_CPUS}" ]]; then
      args+=(--n-cpus "${N_CPUS}")
    fi

    run_python "${RECTANGLE_PY}" "${args[@]}"
  done
}

run_merge_and_report_if_possible() {
  if have_baseline_outputs; then
    echo "[Post] Merging Rectangle output into target tables"
    run_rscript scripts/merge_rectangle_results.R --rectangle-dir results/objects --targets-dir results/tables

    echo "[Post] Rendering report"
    run_rscript scripts/render_report.R
  else
    echo "[Post] Baseline outputs not found; skipping merge and report."
    echo "       Run --mode spacedeconv first (or --mode all) to create baseline outputs."
  fi
}

echo "Mode: ${MODE}"

case "${MODE}" in
  all)
    resolve_rectangle_config
    run_baseline
    run_rectangle_stage
    echo "[Post] Merging Rectangle output into target tables"
    run_rscript scripts/merge_rectangle_results.R --rectangle-dir results/objects --targets-dir results/tables
    echo "[Post] Rendering report"
    run_rscript scripts/render_report.R
    ;;
  spacedeconv)
    run_baseline
    echo "[Post] Rendering report"
    run_rscript scripts/render_report.R
    ;;
  rectangle)
    resolve_rectangle_config
    run_rectangle_stage
    run_merge_and_report_if_possible
    ;;
esac

echo "Pipeline complete."
