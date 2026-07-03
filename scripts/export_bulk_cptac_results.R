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

input_csv <- get_arg("--input-csv", "results/bulk_rnaseq/cptac_brca/tables/cptac_brca_all_methods_aggregated.csv")
mapping_csv <- get_arg("--mapping-csv", "celltype_mapping.csv")
out_csv <- get_arg("--out-csv", "tum_export/deconvolution_results/CPTAC_BRCA_bulk_deconvolution_results_aggregated_rectangle_cell_types.csv")
diagnostics_json <- get_arg("--diagnostics-json", "results/bulk_rnaseq/cptac_brca/tables/cptac_brca_export_diagnostics.json")

if (!file.exists(input_csv)) {
  stop("Missing input CSV: ", input_csv)
}
if (!file.exists(mapping_csv)) {
  stop("Missing mapping CSV: ", mapping_csv)
}

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(diagnostics_json), recursive = TRUE, showWarnings = FALSE)

tab <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
mapping <- read.csv(mapping_csv, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("Fine Annotation", "Coarse Annotation") %in% colnames(mapping))) {
  stop("Mapping file must contain 'Fine Annotation' and 'Coarse Annotation' columns.")
}

fine_to_rectangle <- c(
  rectangle_fine_B_cells_Memory = "rectangle_B.cells.Memory",
  rectangle_fine_B_cells_Naive = "rectangle_B.cells.Naive",
  rectangle_fine_CAFs_MSC_iCAF_like = "rectangle_CAFs.MSC.iCAF.like",
  rectangle_fine_CAFs_myCAF_like = "rectangle_CAFs.myCAF.like",
  rectangle_fine_Cancer_Basal_SC = "rectangle_Cancer.Basal.SC",
  rectangle_fine_Cancer_Her2_SC = "rectangle_Cancer.Her2.SC",
  rectangle_fine_Cancer_LumA_SC = "rectangle_Cancer.LumA.SC",
  rectangle_fine_Cancer_LumB_SC = "rectangle_Cancer.LumB.SC",
  rectangle_fine_DCs = "rectangle_DCs",
  rectangle_fine_Endothelial_ACKR1 = "rectangle_Endothelial.ACKR1",
  rectangle_fine_Endothelial_CXCL12 = "rectangle_Endothelial.CXCL12",
  rectangle_fine_Endothelial_Lymphatic_LYVE1 = "rectangle_Endothelial.Lymphatic.LYVE1",
  rectangle_fine_Endothelial_RGS5 = "rectangle_Endothelial.RGS5",
  rectangle_fine_Macrophage = "rectangle_Macrophage",
  rectangle_fine_Monocyte = "rectangle_Monocyte",
  rectangle_fine_NK_cells = "rectangle_NK.cells",
  rectangle_fine_NKT_cells = "rectangle_NKT.cells",
  rectangle_fine_PVL_Differentiated = "rectangle_PVL.Differentiated",
  rectangle_fine_PVL_Immature = "rectangle_PVL.Immature",
  rectangle_fine_Plasmablasts = "rectangle_Plasmablasts",
  rectangle_fine_T_cells_CD4 = "rectangle_T.cells.CD4.",
  rectangle_fine_T_cells_CD8 = "rectangle_T.cells.CD8.",
  rectangle_fine_Unknown = "rectangle_Unknown"
)

clean_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

base_cols <- c("sample_id", "estimate_tumor", "quantiseq_cd8", "epic_caf", "epic_bcells")
missing_base <- setdiff(base_cols, colnames(tab))
if (length(missing_base) > 0) {
  stop("Missing base export columns: ", paste(missing_base, collapse = ", "))
}

out <- tab[, base_cols, drop = FALSE]
mapping <- mapping[mapping[["Coarse Annotation"]] != "REMOVE", , drop = FALSE]
coarse_levels <- unique(mapping[["Coarse Annotation"]])
aggregated_cols <- character()
for (coarse in coarse_levels) {
  fine_names <- mapping[["Fine Annotation"]][mapping[["Coarse Annotation"]] == coarse]
  rect_cols <- unname(fine_to_rectangle[fine_names])
  if (any(is.na(rect_cols))) {
    stop("Missing fine-to-rectangle mapping for: ", paste(fine_names[is.na(rect_cols)], collapse = ", "))
  }
  missing_rect <- setdiff(rect_cols, colnames(tab))
  if (length(missing_rect) > 0) {
    stop("Missing Rectangle columns for ", coarse, ": ", paste(missing_rect, collapse = ", "))
  }
  out_col <- paste0("rectangle_aggregated_", clean_name(coarse))
  out[[out_col]] <- rowSums(tab[, rect_cols, drop = FALSE], na.rm = TRUE)
  aggregated_cols <- c(aggregated_cols, out_col)
}

write.csv(out, out_csv, row.names = FALSE)

removed_col <- "rectangle_NKT.cells"
sum_with_removed <- NA
if (removed_col %in% colnames(tab)) {
  sum_with_removed <- rowSums(out[, aggregated_cols, drop = FALSE], na.rm = TRUE) + as.numeric(tab[[removed_col]])
}

diagnostics <- list(
  input_csv = input_csv,
  mapping_csv = mapping_csv,
  out_csv = out_csv,
  n_samples = nrow(out),
  n_columns = ncol(out),
  aggregated_rectangle_columns = aggregated_cols,
  removed_rectangle_column = removed_col,
  sum_with_removed_min = if (all(is.na(sum_with_removed))) NA else min(sum_with_removed, na.rm = TRUE),
  sum_with_removed_max = if (all(is.na(sum_with_removed))) NA else max(sum_with_removed, na.rm = TRUE)
)
writeLines(toJSON(diagnostics, pretty = TRUE, auto_unbox = TRUE), con = diagnostics_json)

cat("Wrote CPTAC export:", out_csv, "\n")
cat("Wrote diagnostics:", diagnostics_json, "\n")
