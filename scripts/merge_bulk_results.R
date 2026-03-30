#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx[[1]] == length(args)) {
    return(default)
  }
  args[[idx[[1]] + 1]]
}

objects_dir <- get_arg("--objects-dir", "results/bulk_rnaseq/objects")
tables_dir <- get_arg("--tables-dir", "results/bulk_rnaseq/tables")
metadata_csv <- get_arg("--metadata-csv", "data/bulk_rnaseq/tcga_brca/meta_data_samples.csv")
diagnostics_json <- get_arg("--diagnostics-json", file.path(tables_dir, "tcga_brca_merge_diagnostics.json"))

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

estimate_path <- file.path(objects_dir, "deconv_estimate_tcga_brca.csv")
quantiseq_path <- file.path(objects_dir, "deconv_quantiseq_tcga_brca.csv")
epic_path <- file.path(objects_dir, "deconv_epic_tcga_brca.csv")
rectangle_path <- file.path(objects_dir, "deconv_rectangle_tcga_brca.csv")

required_files <- c(estimate_path, quantiseq_path, epic_path, rectangle_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing bulk deconvolution outputs:\n", paste(missing_files, collapse = "\n"))
}

tidy_immunedeconv <- function(path, prefix) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!("cell_type" %in% colnames(df))) {
    stop("Missing cell_type column in ", path)
  }
  sample_cols <- setdiff(colnames(df), "cell_type")
  out <- data.frame(sample_id = sample_cols, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(df))) {
    col_name <- paste0(prefix, make.names(df$cell_type[[i]]))
    out[[col_name]] <- as.numeric(df[i, sample_cols, drop = TRUE])
  }
  out
}

estimate <- tidy_immunedeconv(estimate_path, "estimate_")
quantiseq <- tidy_immunedeconv(quantiseq_path, "quantiseq_")
epic <- tidy_immunedeconv(epic_path, "epic_")

rectangle <- read.csv(rectangle_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!("spot_id" %in% colnames(rectangle))) {
  stop("Missing spot_id column in Rectangle output")
}
colnames(rectangle)[colnames(rectangle) == "spot_id"] <- "sample_id"
rect_cols <- setdiff(colnames(rectangle), "sample_id")
colnames(rectangle)[match(rect_cols, colnames(rectangle))] <- paste0("rectangle_", make.names(rect_cols))

raw_merged <- Reduce(
  function(x, y) merge(x, y, by = "sample_id", all = TRUE, sort = FALSE),
  list(estimate, quantiseq, epic, rectangle)
)

metadata_joined <- FALSE
matched_metadata_rows <- 0L
if (file.exists(metadata_csv)) {
  metadata <- read.csv(metadata_csv, check.names = FALSE, stringsAsFactors = FALSE)
  if ("sample_barcode" %in% colnames(metadata)) {
    raw_merged <- merge(raw_merged, metadata, by.x = "sample_id", by.y = "sample_barcode", all.x = TRUE, sort = FALSE)
    metadata_joined <- TRUE
    if ("project_id" %in% colnames(raw_merged)) {
      matched_metadata_rows <- sum(!is.na(raw_merged$project_id))
    }
  }
}

aggregated <- raw_merged

require_cols <- function(df, cols, label) {
  missing <- cols[!cols %in% colnames(df)]
  if (length(missing) > 0) {
    stop("Missing columns for ", label, ": ", paste(missing, collapse = ", "))
  }
}

require_cols(aggregated, "estimate_tumor.purity", "estimate_tumor")
require_cols(aggregated, "quantiseq_T.cell.CD8.", "quantiseq_cd8")
require_cols(aggregated, "epic_Cancer.associated.fibroblast", "epic_caf")
require_cols(aggregated, "epic_B.cell", "epic_bcells")
require_cols(
  aggregated,
  c(
    "rectangle_T.cells.CD8.",
    "rectangle_CAFs.MSC.iCAF.like",
    "rectangle_CAFs.myCAF.like",
    "rectangle_Cancer.Her2.SC",
    "rectangle_Cancer.LumB.SC",
    "rectangle_Cancer.Basal.SC",
    "rectangle_Cancer.LumA.SC",
    "rectangle_Unknown",
    "rectangle_B.cells.Memory",
    "rectangle_B.cells.Naive",
    "rectangle_Plasmablasts"
  ),
  "rectangle aggregates"
)

aggregated$estimate_tumor <- aggregated$estimate_tumor.purity
aggregated$quantiseq_cd8 <- aggregated$quantiseq_T.cell.CD8.
aggregated$epic_caf <- aggregated$epic_Cancer.associated.fibroblast
aggregated$epic_bcells <- aggregated$epic_B.cell
aggregated$rectangle_cd8 <- aggregated$rectangle_T.cells.CD8.
aggregated$rectangle_caf <- rowSums(aggregated[, c("rectangle_CAFs.MSC.iCAF.like", "rectangle_CAFs.myCAF.like")], na.rm = TRUE)
aggregated$rectangle_tumor <- rowSums(
  aggregated[, c("rectangle_Cancer.Her2.SC", "rectangle_Cancer.LumB.SC", "rectangle_Cancer.Basal.SC", "rectangle_Cancer.LumA.SC", "rectangle_Unknown")],
  na.rm = TRUE
)
aggregated$rectangle_bcells <- rowSums(
  aggregated[, c("rectangle_B.cells.Memory", "rectangle_B.cells.Naive", "rectangle_Plasmablasts")],
  na.rm = TRUE
)

raw_out <- file.path(tables_dir, "tcga_brca_all_methods_raw.csv")
agg_out <- file.path(tables_dir, "tcga_brca_all_methods_aggregated.csv")
write.csv(raw_merged, raw_out, row.names = FALSE)
write.csv(aggregated, agg_out, row.names = FALSE)

diagnostics <- list(
  metadata_joined = metadata_joined,
  matched_metadata_rows = matched_metadata_rows,
  n_samples = nrow(aggregated),
  raw_out = raw_out,
  aggregated_out = agg_out
)
writeLines(toJSON(diagnostics, pretty = TRUE, auto_unbox = TRUE), diagnostics_json)

message("Bulk merge complete.")
