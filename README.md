# breast_cancer

This repository runs a breast-cancer spatial deconvolution analysis for the 7 slides used in the spacedeconv paper (Fig 2 + Suppl Fig 2), focused on:

- ESTIMATE tumor purity
- quanTIseq CD8+ T cells
- EPIC CAFs
- Rectangle cell fractions (single-cell informed, using a Wu et al. reference)

## Data

Required slide objects are stored in `data/`:

- `allresults_minor_andersson.rds`
- `allresults_minor_1142243F.rds`
- `allresults_minor_1160920F.rds`
- `allresults_minor_4465.rds`
- `allresults_minor_44971.rds`
- `allresults_minor_4290.rds`
- `allresults_minor_4535.rds`

## Run

```bash
# Full workflow: baseline + Rectangle + merge + report
scripts/run_all_with_rectangle.sh --mode all

# Baseline only (ESTIMATE/quanTIseq/EPIC + report)
scripts/run_all_with_rectangle.sh --mode spacedeconv

# Rectangle only (uses existing baseline outputs for merge/report if available)
scripts/run_all_with_rectangle.sh --mode rectangle
```

The runner prints live progress from R/Python commands via `conda run --no-capture-output`.

If needed, override env names:

```bash
R_ENV_NAME=spacedeconv-env \
PY_ENV_NAME=breast-cancer-rectangle-env \
scripts/run_all_with_rectangle.sh --mode all
```

## Rectangle Setup

1. Create/activate Rectangle Python environment:
```bash
conda env create -f environment/rectangle.yml
conda activate breast-cancer-rectangle-env
```
2. Ensure the Wu reference file exists at:
`data/reference/wu_breast_atlas_rectangle.h5ad`

This repo already contains scripts to create it from GEO files:
```bash
python3 scripts/convert_wu_matrix_to_h5ad.py \
  --input-dir data/reference/Wu_etal_2021_BRCA_scRNASeq \
  --output-h5ad data/reference/wu_raw.h5ad

python3 scripts/prepare_wu_reference.py \
  --input-h5ad data/reference/wu_raw.h5ad \
  --output-h5ad data/reference/wu_breast_atlas_rectangle.h5ad \
  --cell-type-col celltype_minor \
  --summary-json results/objects/wu_reference_summary.json
```
3. Verify `config/rectangle.yaml` (reference path + label column).

## Outputs

- Per-slide method objects: `results/objects/`
- Per-slide target tables: `results/tables/targets_<slide>.csv`
- Combined table: `results/tables/targets_all_slides.csv`
- Plot panels: `results/plots/slide_<slide>_targets.png`
- Analysis markdown: `report/breast_cancer_analysis.md`
- Rectangle intermediates: `results/intermediate/rectangle/`
- Rectangle diagnostics: `results/objects/deconv_rectangle_<slide>_diagnostics.json`
