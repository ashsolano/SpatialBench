# Purpose:  Extract a fixed-size physical region of interest (ROI) from a
#           VisiumHD, aligned MERSCOPE, or aligned Xenium Seurat object,
#           using a previously-computed ROI box definition (see
#           compute_roi_box.R) so the same physical region is subset from
#           every platform for a given animal.
# Inputs:   Seurat object RDS (VisiumHD reference, or aligned MERSCOPE/Xenium
#           output of align_binned.R); ROI box definition RDS (output of
#           compute_roi_box.R)
# Outputs:  Subsetted Seurat object saved to --out_rds
# Usage:    Rscript 01_preprocessing/R/extract_roi.R \
#             --input_rds <path> --platform <visium|xenium|merscope> \
#             --roi_box_rds <path> --out_rds <path>

library(optparse)
library(Seurat)

source("01_preprocessing/R/roi_utils.R")  # must be run from the project root

option_list <- list(
  make_option(c("--input_rds"),    type = "character", default = NULL,
              help = "Path to VisiumHD, aligned MERSCOPE, or aligned Xenium Seurat object RDS"),
  make_option(c("--platform"),     type = "character", default = NULL,
              help = "Platform: 'visium', 'xenium', or 'merscope' — selects the assay to subset"),
  make_option(c("--roi_box_rds"),  type = "character", default = NULL,
              help = "Path to this animal's ROI box definition RDS (output of compute_roi_box.R)"),
  make_option(c("--out_rds"),      type = "character", default = NULL,
              help = "Output path for the ROI-subsetted Seurat object RDS")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds))   stop("--input_rds is required")
if (is.null(opt$platform))    stop("--platform is required")
if (is.null(opt$roi_box_rds)) stop("--roi_box_rds is required")
if (is.null(opt$out_rds))     stop("--out_rds is required")

if (!opt$platform %in% c("visium", "xenium", "merscope")) {
  stop("--platform must be 'visium', 'xenium', or 'merscope', got: ", opt$platform)
}

roi_box <- readRDS(opt$roi_box_rds)
obj     <- readRDS(opt$input_rds)

assay_name <- if (opt$platform == "visium") "Spatial.008um" else NULL

obj_roi <- extract_roi(
  obj,
  mpp           = roi_box$mpp,
  roi_center_um = list(cx = roi_box$cx, cy = roi_box$cy),
  roi_size_um   = roi_box$roi_size_um,
  assay_name    = assay_name
)

dir.create(dirname(opt$out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj_roi, file = opt$out_rds)
