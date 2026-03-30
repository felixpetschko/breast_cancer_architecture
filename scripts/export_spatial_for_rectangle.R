#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
})

slides <- c(
  "andersson",
  "1142243F",
  "1160920F",
  "4465",
  "44971",
  "4290",
  "4535"
)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx[[1]] == length(args)) {
    return(default)
  }
  args[[idx[[1]] + 1]]
}

out_dir <- get_arg("--out-dir", "results/intermediate/rectangle")
gene_symbol_policy <- get_arg("--gene-symbol-policy", "upper")
assay_name <- get_arg("--assay", "counts")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

harmonize_genes <- function(x, policy) {
  genes <- as.character(x)
  if (policy == "upper") {
    genes <- toupper(genes)
  } else if (policy == "lower") {
    genes <- tolower(genes)
  }
  make.unique(genes)
}

for (slide in slides) {
  message("Exporting slide: ", slide)
  in_file <- file.path("data", "spatial", "slides", paste0("allresults_minor_", slide, ".rds"))
  if (!file.exists(in_file)) {
    stop("Missing input file: ", in_file)
  }

  spe <- readRDS(in_file)

  current_assay <- assay_name
  if (!(current_assay %in% assayNames(spe))) {
    current_assay <- assayNames(spe)[[1]]
    message("Requested assay not found, using first assay: ", current_assay)
  }

  counts <- as.matrix(assay(spe, current_assay))

  if ("symbol" %in% colnames(as.data.frame(rowData(spe)))) {
    genes <- rowData(spe)$symbol
  } else {
    genes <- rownames(spe)
  }

  genes <- harmonize_genes(genes, gene_symbol_policy)
  rownames(counts) <- genes

  lib_sizes <- colSums(counts)
  lib_sizes[lib_sizes == 0] <- 1
  cpm <- t(t(counts) / lib_sizes) * 1e6

  out <- as.data.frame(t(cpm), stringsAsFactors = FALSE)
  out$spot_id <- rownames(out)
  out <- out[, c("spot_id", setdiff(colnames(out), "spot_id")), drop = FALSE]

  out_file <- file.path(out_dir, paste0(slide, "_bulks_cpm.csv"))
  write.csv(out, out_file, row.names = FALSE)
}

message("Spatial exports complete: ", out_dir)
