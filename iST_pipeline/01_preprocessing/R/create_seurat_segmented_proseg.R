# Purpose:  Create a cell-segmented Seurat object from a single Proseg-segmented sample
#           (Xenium or MERSCOPE). Reads Proseg output files and passes them to
#           myLoadProseg() in utils/segmentation_utils.R.
# Inputs:   Proseg output directory (expected-counts.csv.gz, cell-metadata.csv.gz,
#           transcript-metadata.csv.gz, cell-polygons.geojson.gz)
# Outputs:  Seurat object saved as <out_dir>/<sample_name>_proseg.rds
# Usage:    Rscript 01_preprocessing/R/create_seurat_segmented_proseg.R \
#             --data_dir <path> --sample_name <name> --assay <Vizgen|Xenium> \
#             --out_dir <path>
# Author:   Ashleigh Solano
# Date:     2026-04-20

library(optparse)
library(Seurat)
library(Matrix)

source("utils/segmentation_utils.R")  # must be run from the project root

option_list <- list(
  make_option(c("--data_dir"),    type = "character", default = NULL,
              help = "Path to Proseg output directory (containing expected-counts.csv.gz etc.)"),
  make_option(c("--sample_name"), type = "character", default = NULL,
              help = "Sample identifier; used as the FOV name and in the output filename"),
  make_option(c("--assay"),       type = "character", default = NULL,
              help = "Seurat assay name: 'Vizgen' for MERSCOPE, 'Xenium' for Xenium"),
  make_option(c("--out_dir"),     type = "character", default = NULL,
              help = "Output directory for the RDS file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$data_dir))    stop("--data_dir is required")
if (is.null(opt$sample_name)) stop("--sample_name is required")
if (is.null(opt$assay))       stop("--assay is required (Vizgen or Xenium)")
if (is.null(opt$out_dir))     stop("--out_dir is required")

if (!opt$assay %in% c("Vizgen", "Xenium")) {
  stop("--assay must be 'Vizgen' (MERSCOPE) or 'Xenium'")
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

message("Sample:    ", opt$sample_name)
message("Assay:     ", opt$assay)
message("Data dir:  ", opt$data_dir)
message("Output dir:", opt$out_dir)

start_time <- Sys.time()

obj <- myLoadProseg(
  data.dir = opt$data_dir,
  fov      = opt$sample_name,
  assay    = opt$assay
)

elapsed <- round(Sys.time() - start_time, digits = 2)
message("Loaded in ", elapsed, " ", attr(elapsed, "units"))

out_file <- file.path(opt$out_dir, paste0(opt$sample_name, "_proseg.rds"))
saveRDS(obj, file = out_file)
message("Saved: ", out_file)
