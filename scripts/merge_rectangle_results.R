#!/usr/bin/env Rscript

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

rectangle_dir <- get_arg("--rectangle-dir", "results/objects")
targets_dir <- get_arg("--targets-dir", "results/tables")

dir.create(targets_dir, recursive = TRUE, showWarnings = FALSE)

merge_one_slide <- function(slide) {
  base_path <- file.path(targets_dir, paste0("targets_", slide, ".csv"))
  rect_path <- file.path(rectangle_dir, paste0("deconv_rectangle_", slide, ".csv"))

  if (!file.exists(base_path)) {
    stop("Missing base targets table: ", base_path)
  }
  if (!file.exists(rect_path)) {
    stop("Missing Rectangle table: ", rect_path)
  }

  base <- read.csv(base_path, stringsAsFactors = FALSE, check.names = FALSE)
  rect <- read.csv(rect_path, stringsAsFactors = FALSE, check.names = FALSE)

  if (!("spot_id" %in% colnames(base)) || !("spot_id" %in% colnames(rect))) {
    stop("Both tables must contain a spot_id column for slide: ", slide)
  }

  rect_cols <- setdiff(colnames(rect), "spot_id")
  prefixed <- paste0("rectangle_", make.names(rect_cols))
  colnames(rect)[match(rect_cols, colnames(rect))] <- prefixed

  merged <- base
  row_idx <- match(merged$spot_id, rect$spot_id)
  for (col in setdiff(colnames(rect), "spot_id")) {
    merged[[col]] <- rect[[col]][row_idx]
  }

  if (!("slide_id" %in% colnames(merged))) {
    merged$slide_id <- slide
  }

  write.csv(merged, base_path, row.names = FALSE)
  merged
}

all_tables <- vector("list", length(slides))
names(all_tables) <- slides

for (slide in slides) {
  message("Merging Rectangle results for slide: ", slide)
  all_tables[[slide]] <- merge_one_slide(slide)
}

combined <- do.call(rbind, all_tables)
write.csv(combined, file.path(targets_dir, "targets_all_slides.csv"), row.names = FALSE)

message("Rectangle merge complete.")
