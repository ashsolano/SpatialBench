# Purpose:  Create a binned Seurat object from a single MERSCOPE sample.
#           Resolution is a parameter so the same script handles 8um and 16um.
# Inputs:   MERSCOPE data directory (containing detected_transcripts.csv)
# Outputs:  Seurat object saved as <out_dir>/<sample_name>_<resolution>um.rds
# Usage:    Rscript 01_preprocessing/R/create_seurat_binned_merscope.R \
#             --data_dir <path> --sample_name <name> --resolution <8|16> --out_dir <path>

library(optparse)
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)

source("utils/binning_utils.R")  # must be run from the project root

option_list <- list(
  make_option(c("--data_dir"),    type = "character", default = NULL,
              help = "Path to MERSCOPE data directory (containing detected_transcripts.csv)"),
  make_option(c("--sample_name"), type = "character", default = NULL,
              help = "Sample identifier; used as the FOV name and in the output filename"),
  make_option(c("--resolution"),  type = "integer",   default = 8L,
              help = "Bin size in microns [default: %default]"),
  make_option(c("--out_dir"),     type = "character", default = NULL,
              help = "Output directory for the RDS file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$data_dir))    stop("--data_dir is required")
if (is.null(opt$sample_name)) stop("--sample_name is required")
if (is.null(opt$out_dir))     stop("--out_dir is required")

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

message("Sample:     ", opt$sample_name)
message("Resolution: ", opt$resolution, "um")
message("Data dir:   ", opt$data_dir)
message("Output dir: ", opt$out_dir)

start_time <- Sys.time()

obj <- LoadVizgen_binned(
  data.dir   = opt$data_dir,
  resolution = opt$resolution,
  fov        = opt$sample_name,
  assay      = "Vizgen"
)

elapsed <- round(Sys.time() - start_time, digits = 2)
message("Loaded in ", elapsed, " ", attr(elapsed, "units"))

out_file <- file.path(opt$out_dir, paste0(opt$sample_name, "_", opt$resolution, "um.rds"))
saveRDS(obj, file = out_file)
message("Saved: ", out_file)
