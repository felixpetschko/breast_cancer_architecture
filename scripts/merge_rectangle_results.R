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

  # Aggregate Rectangle outputs for downstream plotting/summary:
  # CD8    = T cells CD8+
  # CAF    = CAFs MSC iCAF-like + CAFs myCAF-like
  # Tumor  = Cancer Her2 SC + Cancer LumB SC + Cancer Basal SC + Cancer LumA SC + Unknown
  # B cell = B cells Memory + B cells Naive + Plasmablasts
  req_cd8 <- "rectangle_T.cells.CD8."
  req_caf <- c("rectangle_CAFs.MSC.iCAF.like", "rectangle_CAFs.myCAF.like")
  req_tumor <- c(
    "rectangle_Cancer.Her2.SC",
    "rectangle_Cancer.LumB.SC",
    "rectangle_Cancer.Basal.SC",
    "rectangle_Cancer.LumA.SC",
    "rectangle_Unknown"
  )
  req_bcells <- c(
    "rectangle_B.cells.Memory",
    "rectangle_B.cells.Naive",
    "rectangle_Plasmablasts"
  )

  if (req_cd8 %in% colnames(merged)) {
    merged$rectangle_cd8 <- merged[[req_cd8]]
  }

  if (all(req_caf %in% colnames(merged))) {
    merged$rectangle_caf <- rowSums(merged[, req_caf, drop = FALSE], na.rm = TRUE)
  }

  if (all(req_tumor %in% colnames(merged))) {
    merged$rectangle_tumor <- rowSums(merged[, req_tumor, drop = FALSE], na.rm = TRUE)
  }

  if (all(req_bcells %in% colnames(merged))) {
    merged$rectangle_bcells <- rowSums(merged[, req_bcells, drop = FALSE], na.rm = TRUE)
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
