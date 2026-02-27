#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(spacedeconv)
  library(SpatialExperiment)
})

dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

slides <- c(
  "andersson",
  "1142243F",
  "1160920F",
  "4465",
  "44971",
  "4290",
  "4535"
)

input_files <- file.path("data", paste0("allresults_minor_", slides, ".rds"))
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0) {
  stop(
    "Missing input files:\n",
    paste(missing_inputs, collapse = "\n")
  )
}

pick_col <- function(df, candidates, label) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    stop("Could not find column for ", label, ". Checked: ", paste(candidates, collapse = ", "))
  }
  hit[[1]]
}

extract_table <- function(spe, slide_id) {
  cd <- as.data.frame(colData(spe))
  coords <- as.data.frame(spatialCoords(spe))
  coords$spot_id <- colnames(spe)
  coords <- coords[, c("spot_id", colnames(coords)[colnames(coords) != "spot_id"]), drop = FALSE]

  tumor_col <- pick_col(cd, c("estimate_tumor.purity"), "ESTIMATE tumor purity")
  cd8_col <- pick_col(cd, c("quantiseq_T.cell.CD8."), "quanTIseq CD8")
  caf_col <- pick_col(cd, c("epic_Cancer.associated.fibroblast"), "EPIC CAF")

  out <- data.frame(
    slide_id = slide_id,
    spot_id = colnames(spe),
    tumor_purity_estimate = cd[[tumor_col]],
    cd8_quantiseq = cd[[cd8_col]],
    caf_epic = cd[[caf_col]],
    stringsAsFactors = FALSE
  )
  merge(out, coords, by = "spot_id", all.x = TRUE, sort = FALSE)
}

all_tables <- vector("list", length(slides))
names(all_tables) <- slides

for (slide in slides) {
  message("Processing slide: ", slide)
  in_file <- file.path("data", paste0("allresults_minor_", slide, ".rds"))
  spe <- readRDS(in_file)

  if ("symbol" %in% colnames(as.data.frame(rowData(spe)))) {
    rownames(spe) <- make.names(rowData(spe)$symbol, unique = TRUE)
  } else {
    rownames(spe) <- make.names(rownames(spe), unique = TRUE)
  }

  # Follow the deconvolution setup used in the spacedeconv paper for ESTIMATE/EPIC/quanTIseq.
  spe <- preprocess(spe, min_umi = 87)
  spe <- spacedeconv::normalize(spe)

  deconv_estimate <- deconvolute(spe, method = "estimate", assay_sp = "cpm")
  deconv_epic <- deconvolute(spe, method = "epic", assay_sp = "cpm", tumor = TRUE)
  deconv_quantiseq <- deconvolute(spe, method = "quantiseq", assay_sp = "cpm", tumor = TRUE)

  saveRDS(
    deconv_estimate,
    file = file.path("results/objects", paste0("deconv_estimate_", slide, ".rds"))
  )
  saveRDS(
    deconv_epic,
    file = file.path("results/objects", paste0("deconv_epic_", slide, ".rds"))
  )
  saveRDS(
    deconv_quantiseq,
    file = file.path("results/objects", paste0("deconv_quantiseq_", slide, ".rds"))
  )

  merged <- deconv_estimate
  colData(merged)[, colnames(colData(deconv_epic))] <- colData(deconv_epic)
  colData(merged)[, colnames(colData(deconv_quantiseq))] <- colData(deconv_quantiseq)

  slide_table <- extract_table(merged, slide)
  all_tables[[slide]] <- slide_table

  write.csv(
    slide_table,
    file = file.path("results/tables", paste0("targets_", slide, ".csv")),
    row.names = FALSE
  )
}

combined <- do.call(rbind, all_tables)
write.csv(combined, file = "results/tables/targets_all_slides.csv", row.names = FALSE)

writeLines(
  c(
    "Analysis run complete.",
    paste("Slides processed:", paste(slides, collapse = ", "))
  ),
  con = "results/tables/run_summary.txt"
)
