#!/usr/bin/env bash
set -euo pipefail

PY_ENV_NAME="${PY_ENV_NAME:-/scratch/c9881013/.conda_envs/breast-cancer-rectangle-env}"
CONFIG_PATH="${CONFIG_PATH:-config/rectangle.yaml}"
RECTANGLE_PY="${RECTANGLE_PY:-scripts/run_rectangle.py}"
REFERENCE_H5AD="${REFERENCE_H5AD:-data/spatial/wu_reference/wu_breast_atlas_rectangle.h5ad}"
CELL_TYPE_COL="${CELL_TYPE_COL:-celltype_minor}"
GENE_SYMBOL_POLICY="${GENE_SYMBOL_POLICY:-upper}"
MIN_GENES_OVERLAP="${MIN_GENES_OVERLAP:-500}"
CORRECT_MRNA_BIAS="${CORRECT_MRNA_BIAS:-true}"
N_CPUS="${N_CPUS:-32}"

PREPARED_CSV="${PREPARED_CSV:-results/bulk_rnaseq/cptac_brca/intermediate/cptac_brca_tpm_samples_by_genes.csv}"
RECTANGLE_INPUT_CSV="${RECTANGLE_INPUT_CSV:-results/bulk_rnaseq/cptac_brca/intermediate/cptac_brca_bulks_cpm.csv}"
OUT_CSV="${OUT_CSV:-results/bulk_rnaseq/cptac_brca/objects/deconv_rectangle_cptac_brca.csv}"
DIAG_JSON="${DIAG_JSON:-results/bulk_rnaseq/cptac_brca/objects/deconv_rectangle_cptac_brca_diagnostics.json}"

if [[ ! -f "${PREPARED_CSV}" ]]; then
  echo "Missing prepared bulk matrix: ${PREPARED_CSV}"
  exit 1
fi

mkdir -p results/bulk_rnaseq/cptac_brca/intermediate results/bulk_rnaseq/cptac_brca/objects
export PREPARED_CSV RECTANGLE_INPUT_CSV

conda_run_python() {
  if [[ -z "${PY_ENV_NAME}" ]]; then
    python3 "$@"
  elif [[ -d "${PY_ENV_NAME}" || "${PY_ENV_NAME}" == */* ]]; then
    conda run --no-capture-output -p "${PY_ENV_NAME}" python "$@"
  else
    conda run --no-capture-output -n "${PY_ENV_NAME}" python "$@"
  fi
}

if [[ -n "${PY_ENV_NAME}" ]]; then
  conda_run_python - <<'PY'
import pandas as pd
import os
inp = os.environ["PREPARED_CSV"]
out = os.environ["RECTANGLE_INPUT_CSV"]
df = pd.read_csv(inp)
if "sample_id" not in df.columns:
    raise ValueError("Prepared bulk matrix must contain sample_id")
df = df.rename(columns={"sample_id": "spot_id"})
df.to_csv(out, index=False)
print(f"Wrote Rectangle CPTAC bulk input: {out}")
PY
else
  python3 - <<'PY'
import pandas as pd
import os
inp = os.environ["PREPARED_CSV"]
out = os.environ["RECTANGLE_INPUT_CSV"]
df = pd.read_csv(inp)
if "sample_id" not in df.columns:
    raise ValueError("Prepared bulk matrix must contain sample_id")
df = df.rename(columns={"sample_id": "spot_id"})
df.to_csv(out, index=False)
print(f"Wrote Rectangle CPTAC bulk input: {out}")
PY
fi

args=(
  --config "${CONFIG_PATH}"
  --reference-h5ad "${REFERENCE_H5AD}"
  --bulks-csv "${RECTANGLE_INPUT_CSV}"
  --cell-type-col "${CELL_TYPE_COL}"
  --gene-symbol-policy "${GENE_SYMBOL_POLICY}"
  --min-genes-overlap "${MIN_GENES_OVERLAP}"
  --out-csv "${OUT_CSV}"
  --diagnostics-json "${DIAG_JSON}"
)

if [[ "${CORRECT_MRNA_BIAS}" == "true" ]]; then
  args+=(--correct-mrna-bias)
else
  args+=(--no-correct-mrna-bias)
fi

if [[ -n "${N_CPUS}" ]]; then
  args+=(--n-cpus "${N_CPUS}")
fi

echo ">>> [PY] ${RECTANGLE_PY} ${args[*]}"
conda_run_python "${RECTANGLE_PY}" "${args[@]}"
