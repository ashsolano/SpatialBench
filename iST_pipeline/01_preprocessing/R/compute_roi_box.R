# Purpose:  Compute a fixed-size physical region-of-interest (ROI) box,
#           centred on one animal's VisiumHD tissue extent, and save its
#           definition (center + bounds, in microns) for reuse — both by
#           extract_roi.R (to subset VisiumHD/MERSCOPE/Xenium objects to the
#           same physical region) and for later plotting the ROI box on any
#           image of the tissue, without reloading the full VisiumHD object.
# Inputs:   VisiumHD reference Seurat object RDS (config$visiumhd); this
#           animal's VisiumHD sample ID and microns_per_pixel value
#           (config$stalign$visiumhd / microns_per_pixel); ROI size
#           (config$roi$size_um)
# Outputs:  List with animal, visium_sample, mpp, roi_size_um, cx, cy, xmin,
#           xmax, ymin, ymax, saved to --out_rds
# Usage:    Rscript 01_preprocessing/R/compute_roi_box.R \
#             --sample <animal ID, e.g. wt709> --config config/config.yaml \
#             --out_rds <path>

library(optparse)
library(Seurat)
library(yaml)

source("01_preprocessing/R/roi_utils.R")  # must be run from the project root

option_list <- list(
  make_option(c("--sample"),  type = "character", default = NULL,
              help = "Animal ID, e.g. 'wt709', 'ko167' — looked up in config$stalign$visiumhd"),
  make_option(c("--config"),  type = "character", default = "config/config.yaml",
              help = "Path to config.yaml [default: %default]"),
  make_option(c("--out_rds"), type = "character", default = NULL,
              help = "Output path for the ROI box definition RDS")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$sample))  stop("--sample is required")
if (is.null(opt$out_rds)) stop("--out_rds is required")

config <- yaml::read_yaml(opt$config)

visium_sample <- config$stalign$visiumhd[[opt$sample]]
if (is.null(visium_sample)) {
  stop("No config$stalign$visiumhd entry for sample '", opt$sample, "'")
}

mpp <- config$stalign$microns_per_pixel[[visium_sample]]
if (is.null(mpp)) {
  stop("No config$stalign$microns_per_pixel entry for VisiumHD sample '", visium_sample, "'")
}

visium_rds_name <- config$visiumhd$samples[[visium_sample]]
if (is.null(visium_rds_name)) {
  stop("No config$visiumhd$samples entry for VisiumHD sample '", visium_sample, "'")
}
visium_rds <- file.path(config$visiumhd$data_dir, visium_rds_name)
if (!file.exists(visium_rds)) {
  stop("VisiumHD reference RDS not found: ", visium_rds)
}

visium_obj <- readRDS(visium_rds)

roi_center_um <- get_roi_center(visium_obj, mpp)
roi_size_um   <- config$roi$size_um
half          <- roi_size_um / 2

roi_box <- list(
  animal        = opt$sample,
  visium_sample = visium_sample,
  mpp           = mpp,
  roi_size_um   = roi_size_um,
  cx            = roi_center_um$cx,
  cy            = roi_center_um$cy,
  xmin          = roi_center_um$cx - half,
  xmax          = roi_center_um$cx + half,
  ymin          = roi_center_um$cy - half,
  ymax          = roi_center_um$cy + half
)

dir.create(dirname(opt$out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(roi_box, file = opt$out_rds)
