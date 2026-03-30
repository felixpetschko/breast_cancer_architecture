#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rmarkdown)
})

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)

setwd("report")

rmarkdown::render(
  input = "bulk_rnaseq_analysis.Rmd",
  output_file = "bulk_rnaseq_analysis.md"
)
