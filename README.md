# breast_cancer

This repository runs a breast-cancer spacedeconv analysis for the 7 slides used in the spacedeconv paper (Fig 2 + Suppl Fig 2), focused on:

- ESTIMATE tumor purity
- quanTIseq CD8+ T cells
- EPIC CAFs

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
scripts/run_breast_spacedeconv_methods.R
scripts/render_report.R
```

## Outputs

- Per-slide method objects: `results/objects/`
- Per-slide target tables: `results/tables/targets_<slide>.csv`
- Combined table: `results/tables/targets_all_slides.csv`
- Plot panels: `results/plots/slide_<slide>_targets.png`
- Analysis markdown: `report/breast_spacedeconv_analysis.md`
