#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rmarkdown)
})

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)

setwd("report")

rmarkdown::render(
  input = "cptac_brca_bulk_rnaseq_analysis.Rmd",
  output_file = "cptac_brca_bulk_rnaseq_analysis.md"
)
