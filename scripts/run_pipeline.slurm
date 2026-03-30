#!/usr/bin/env bash
set -euo pipefail

# This script is intended to be submitted via scripts/submit_slurm.sh
# Required exported vars from sbatch wrapper:
#   MODE, R_ENV_NAME, PY_ENV_NAME, IMMUNEDECONV_ENV_NAME
# Optional exported vars:
#   N_CPUS, CONFIG_PATH, REFERENCE_H5AD, CELL_TYPE_COL,
#   GENE_SYMBOL_POLICY, CORRECT_MRNA_BIAS, MIN_GENES_OVERLAP

module load Anaconda3
source "$(conda info --base)/etc/profile.d/conda.sh"

REPO_DIR="/scratch/c9881013/felix_breast_cancer/breast_cancer"
cd "${REPO_DIR}"

mkdir -p logs/slurm

echo "[$(date)] Starting breast_cancer pipeline on host $(hostname)"
echo "MODE=${MODE:-all}"
echo "R_ENV_NAME=${R_ENV_NAME:-spacedeconv-env}"
echo "PY_ENV_NAME=${PY_ENV_NAME:-breast-cancer-rectangle-env}"
echo "IMMUNEDECONV_ENV_NAME=${IMMUNEDECONV_ENV_NAME:-immunedeconv-env}"
[[ -n "${N_CPUS:-}" ]] && echo "N_CPUS=${N_CPUS}"

case "${MODE:-all}" in
  all|spacedeconv|rectangle)
    scripts/run_all_with_rectangle.sh --mode "${MODE:-all}"
    ;;
  bulk-all)
    scripts/run_bulk_tcga_brca.sh --mode all
    ;;
  bulk-prepare)
    scripts/run_bulk_tcga_brca.sh --mode prepare
    ;;
  bulk-immunedeconv)
    scripts/run_bulk_tcga_brca.sh --mode immunedeconv
    ;;
  bulk-rectangle)
    scripts/run_bulk_tcga_brca.sh --mode rectangle
    ;;
  bulk-merge)
    scripts/run_bulk_tcga_brca.sh --mode merge
    ;;
  bulk-report)
    scripts/run_bulk_tcga_brca.sh --mode report
    ;;
  *)
    echo "Unsupported MODE=${MODE:-all}"
    exit 1
    ;;
esac

echo "[$(date)] Pipeline finished"
