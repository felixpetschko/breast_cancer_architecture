# Breast Cancer Architecture

This repository runs a breast-cancer spatial deconvolution analysis for the 7 slides used in the spacedeconv paper (Fig 2 + Suppl Fig 2), focused on:

- ESTIMATE tumor purity
- quanTIseq CD8+ T cells
- EPIC CAFs
- EPIC B cells
- Rectangle cell fractions (single-cell informed, using a Wu et al. reference)

## Data

Required spatial slide objects are stored in `data/spatial/slides/`:

- `allresults_minor_andersson.rds`
- `allresults_minor_1142243F.rds`
- `allresults_minor_1160920F.rds`
- `allresults_minor_4465.rds`
- `allresults_minor_44971.rds`
- `allresults_minor_4290.rds`
- `allresults_minor_4535.rds`

Data layout:

- Spatial Wu reference files in `data/spatial/wu_reference/`:
  - `count_matrix_sparse.mtx`
  - `count_matrix_genes.tsv`
  - `count_matrix_barcodes.tsv`
  - `metadata.csv`
- Bulk TCGA-BRCA files in `data/bulk_rnaseq/tcga_brca/`:
  - `brca_tpm_all_symbols.csv`
  - `brca_counts_all.csv`
  - `meta_data_samples.csv`
  - `clinical.tsv`
  - `tcga_survival_curated_2018.csv`

Source links:

- Spatial slide inputs:
  - https://github.com/ComputationalBiomedicineGroup/spacedeconv_paper/tree/main/data
- Wu scRNA-seq reference files:
  - https://drive.google.com/drive/folders/1WJpUG134anGH1dq12n_f0gt5Kz8vKA37?usp=sharing
- TCGA-BRCA bulk RNA-seq dataset:
  - https://drive.google.com/file/d/1kBY2jAqpKjk6Ptmqq8FUAOGdtBDzQfhN/view

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

## SLURM (easy)

Submit with one command:

```bash
# Full pipeline
scripts/submit_slurm.sh --mode all

# Only spacedeconv
scripts/submit_slurm.sh --mode spacedeconv

# Only rectangle (and use 8 CPUs inside Rectangle deconvolution)
scripts/submit_slurm.sh --mode rectangle --n-cpus 8

# Quick test queue run (partition=test, max walltime 10:00:00)
scripts/submit_slurm.sh --mode rectangle --test --time 10:00:00
```

Common resource flags:

```bash
scripts/submit_slurm.sh \
  --mode all \
  --cpus 32 \
  --mem 128G \
  --time 48:00:00 \
  --account biology
```

Logs are written to `logs/slurm/%x_%j.out` and `logs/slurm/%x_%j.err`.

## Environments

The R environment `spacedeconv-env` can be created based on the `environment.yml` from the upstream spacedeconv repository:

- https://github.com/omnideconv/spacedeconv

The Rectangle Python environment used here is `breast-cancer-rectangle-env`.

```bash
conda env create -f environment/rectangle.yml
conda activate breast-cancer-rectangle-env
```

## Rectangle Setup

1. Activate the Rectangle Python environment:

```bash
conda activate breast-cancer-rectangle-env
```

2. Keep the 4 Wu raw files in `data/spatial/wu_reference/`:

```bash
ls data/spatial/wu_reference
# count_matrix_sparse.mtx  count_matrix_genes.tsv  count_matrix_barcodes.tsv  metadata.csv
```

3. The pipeline auto-builds these generated files if missing:
- `data/spatial/wu_reference/wu_raw.h5ad`
- `data/spatial/wu_reference/wu_breast_atlas_rectangle.h5ad`

Manual build remains possible:

```bash
python3 scripts/convert_wu_matrix_to_h5ad.py \
  --input-dir data/spatial/wu_reference \
  --output-h5ad data/spatial/wu_reference/wu_raw.h5ad

python3 scripts/prepare_wu_reference.py \
  --input-h5ad data/spatial/wu_reference/wu_raw.h5ad \
  --output-h5ad data/spatial/wu_reference/wu_breast_atlas_rectangle.h5ad \
  --cell-type-col celltype_minor \
  --summary-json results/objects/wu_reference_summary.json
```
4. Verify `config/rectangle.yaml` (reference path + label column).
   The `excluded_cell_types` list in that config is applied during reference preparation.
   If you change this list, force rebuild by removing generated files once:
```bash
rm -f data/spatial/wu_reference/wu_raw.h5ad data/spatial/wu_reference/wu_breast_atlas_rectangle.h5ad
```

## Outputs

- Per-slide method objects: `results/objects/`
- Per-slide target tables: `results/tables/targets_<slide>.csv`
- Combined table: `results/tables/targets_all_slides.csv`
- Plot panels: `results/plots/slide_<slide>_targets.png`
- Analysis markdown: `report/breast_cancer_analysis.md`
- Rectangle intermediates: `results/intermediate/rectangle/`
- Rectangle diagnostics: `results/objects/deconv_rectangle_<slide>_diagnostics.json`

Rectangle downstream aggregates used for plotting/tables:

- `rectangle_cd8 = T cells CD8+`
- `rectangle_caf = CAFs MSC iCAF-like + CAFs myCAF-like`
- `rectangle_tumor = Cancer Her2 SC + Cancer LumB SC + Cancer Basal SC + Cancer LumA SC + Unknown`
- `rectangle_bcells = B cells Memory + B cells Naive + Plasmablasts`

## Report Only

To re-render the report from existing results:

```bash
conda run --no-capture-output -n spacedeconv-env Rscript scripts/render_report.R
```
