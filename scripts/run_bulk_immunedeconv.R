#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(immunedeconv)
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

input_csv <- get_arg("--input-csv", "results/bulk_rnaseq/intermediate/tcga_brca_tpm_samples_by_genes.csv")
out_dir <- get_arg("--out-dir", "results/bulk_rnaseq/objects")
diagnostics_json <- get_arg("--diagnostics-json", file.path(out_dir, "deconv_immunedeconv_tcga_brca_diagnostics.json"))

if (!file.exists(input_csv)) {
  stop("Missing prepared bulk matrix: ", input_csv)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prepared <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
if (!("sample_id" %in% colnames(prepared))) {
  stop("Prepared bulk matrix must contain a sample_id column.")
}

sample_ids <- prepared$sample_id
expr <- prepared[, setdiff(colnames(prepared), "sample_id"), drop = FALSE]
expr[] <- lapply(expr, as.numeric)
expr <- as.data.frame(t(as.matrix(expr)), check.names = FALSE, stringsAsFactors = FALSE)
colnames(expr) <- sample_ids
rownames(expr) <- colnames(prepared)[colnames(prepared) != "sample_id"]

run_method <- function(method, tumor = TRUE) {
  deconvolute(expr, method = method, tumor = tumor, arrays = FALSE)
}

estimate_raw <- run_method("estimate", tumor = TRUE)
quantiseq_raw <- run_method("quantiseq", tumor = TRUE)
epic_raw <- run_method("epic", tumor = TRUE)

write.csv(estimate_raw, file.path(out_dir, "deconv_estimate_tcga_brca.csv"), row.names = FALSE)
write.csv(quantiseq_raw, file.path(out_dir, "deconv_quantiseq_tcga_brca.csv"), row.names = FALSE)
write.csv(epic_raw, file.path(out_dir, "deconv_epic_tcga_brca.csv"), row.names = FALSE)

diagnostics <- list(
  input_csv = input_csv,
  n_samples = ncol(expr),
  n_genes = nrow(expr),
  methods = c("estimate", "quantiseq", "epic"),
  outputs = list(
    estimate = file.path(out_dir, "deconv_estimate_tcga_brca.csv"),
    quantiseq = file.path(out_dir, "deconv_quantiseq_tcga_brca.csv"),
    epic = file.path(out_dir, "deconv_epic_tcga_brca.csv")
  )
)
writeLines(toJSON(diagnostics, pretty = TRUE, auto_unbox = TRUE), con = diagnostics_json)

cat("Wrote immunedeconv outputs to:", out_dir, "\n")
cat("Wrote diagnostics:", diagnostics_json, "\n")
