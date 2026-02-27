#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rmarkdown)
})

rmarkdown::render(
  input = "report/breast_spacedeconv_analysis.Rmd",
  output_file = "breast_spacedeconv_analysis.md",
  output_dir = "report"
)
