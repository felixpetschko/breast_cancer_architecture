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

input_csv <- get_arg("--input-csv", "data/bulk_rnaseq/cptac_brca/tpm_matrix_gene_id_name.csv")
cleaned_csv <- get_arg("--cleaned-csv", "data/bulk_rnaseq/cptac_brca/cptac_brca_tpm_gene_symbols.csv")
out_csv <- get_arg("--out-csv", "results/bulk_rnaseq/cptac_brca/intermediate/cptac_brca_tpm_samples_by_genes.csv")
summary_json <- get_arg("--summary-json", "results/bulk_rnaseq/cptac_brca/intermediate/cptac_brca_input_summary.json")
gene_symbol_policy <- get_arg("--gene-symbol-policy", "upper")

if (!file.exists(input_csv)) {
  stop("Missing input CSV: ", input_csv)
}

dir.create(dirname(cleaned_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_json), recursive = TRUE, showWarnings = FALSE)

harmonize_genes <- function(values, policy) {
  genes <- as.character(values)
  if (policy == "upper") {
    genes <- toupper(genes)
  } else if (policy == "lower") {
    genes <- tolower(genes)
  }
  genes
}

raw <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
required_cols <- c("gene_id", "gene_name")
missing_cols <- setdiff(required_cols, colnames(raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns in CPTAC input: ", paste(missing_cols, collapse = ", "))
}

sample_cols <- setdiff(colnames(raw), required_cols)
if (length(sample_cols) == 0) {
  stop("CPTAC input must contain sample columns after gene_id/gene_name.")
}

genes_raw <- as.character(raw$gene_name)
keep_gene <- !is.na(genes_raw) & genes_raw != ""
raw <- raw[keep_gene, , drop = FALSE]
genes <- harmonize_genes(raw$gene_name, gene_symbol_policy)

expr <- raw[, sample_cols, drop = FALSE]
expr[] <- lapply(expr, as.numeric)
mean_tpm <- rowMeans(expr, na.rm = TRUE)
mean_tpm[is.nan(mean_tpm)] <- NA_real_

dedup_frame <- data.frame(
  row_index = seq_len(nrow(raw)),
  gene_name = genes,
  mean_tpm = mean_tpm,
  stringsAsFactors = FALSE
)
dedup_frame <- dedup_frame[order(dedup_frame$gene_name, -dedup_frame$mean_tpm, dedup_frame$row_index), , drop = FALSE]
keep_idx <- dedup_frame$row_index[!duplicated(dedup_frame$gene_name)]

cleaned_expr <- expr[keep_idx, , drop = FALSE]
cleaned_genes <- genes[keep_idx]
rownames(cleaned_expr) <- cleaned_genes

cleaned <- data.frame(
  gene_name = cleaned_genes,
  cleaned_expr,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.csv(cleaned, cleaned_csv, row.names = FALSE)

prepared <- data.frame(
  sample_id = colnames(cleaned_expr),
  t(as.matrix(cleaned_expr)),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.csv(prepared, out_csv, row.names = FALSE)

summary <- list(
  input_csv = input_csv,
  cleaned_csv = cleaned_csv,
  out_csv = out_csv,
  n_samples = length(sample_cols),
  n_rows_input = nrow(raw),
  n_unique_gene_symbols = length(unique(genes)),
  n_rows_after_duplicate_resolution = nrow(cleaned_expr),
  n_removed_empty_gene_names = sum(!keep_gene),
  n_duplicate_rows_removed = nrow(raw) - nrow(cleaned_expr),
  duplicate_resolution = "For duplicated gene_name values, kept the row with the highest mean TPM across all samples.",
  gene_symbol_policy = gene_symbol_policy
)
writeLines(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), con = summary_json)

cat("Wrote cleaned CPTAC gene-symbol matrix:", cleaned_csv, "\n")
cat("Wrote prepared CPTAC bulk matrix:", out_csv, "\n")
cat("Wrote diagnostics:", summary_json, "\n")
