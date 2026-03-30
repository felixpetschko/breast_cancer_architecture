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

input_csv <- get_arg("--input-csv", "data/bulk_rnaseq/tcga_brca/brca_tpm_all_symbols.csv")
out_csv <- get_arg("--out-csv", "results/bulk_rnaseq/intermediate/tcga_brca_tpm_samples_by_genes.csv")
summary_json <- get_arg("--summary-json", "results/bulk_rnaseq/intermediate/tcga_brca_input_summary.json")
gene_symbol_policy <- get_arg("--gene-symbol-policy", "upper")

if (!file.exists(input_csv)) {
  stop("Missing input CSV: ", input_csv)
}

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_json), recursive = TRUE, showWarnings = FALSE)

harmonize_genes <- function(values, policy) {
  genes <- as.character(values)
  if (policy == "upper") {
    genes <- toupper(genes)
  } else if (policy == "lower") {
    genes <- tolower(genes)
  }
  make.unique(genes)
}

raw <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(raw) < 2) {
  stop("Input CSV must contain one gene column and at least one sample column.")
}

genes_before <- as.character(raw[[1]])

had_mean_tpm <- "..meanTPM.." %in% colnames(raw)
if (had_mean_tpm) {
  raw <- raw[, colnames(raw) != "..meanTPM..", drop = FALSE]
}

sample_ids <- colnames(raw)[-1]
genes_after <- harmonize_genes(raw[[1]], gene_symbol_policy)
expr <- raw[, -1, drop = FALSE]
expr[] <- lapply(expr, as.numeric)
rownames(expr) <- genes_after

prepared <- data.frame(
  sample_id = sample_ids,
  t(as.matrix(expr)),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.csv(prepared, out_csv, row.names = FALSE)

summary <- list(
  input_csv = input_csv,
  out_csv = out_csv,
  n_samples = length(sample_ids),
  n_genes_before = length(genes_before),
  n_genes_after = nrow(expr),
  had_mean_tpm = had_mean_tpm,
  gene_symbol_policy = gene_symbol_policy
)
writeLines(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), con = summary_json)

cat("Wrote prepared bulk matrix:", out_csv, "\n")
cat("Wrote diagnostics:", summary_json, "\n")
