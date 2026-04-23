# Purpose:  Create a cell-segmented Seurat object from a single MERSCOPE sample.
#           Loads cell boundaries and transcript locations via Seurat's LoadVizgen().
#           Blank- negative control genes are separated into a dedicated assay.
# Inputs:   Vizgen output directory (cell_by_gene.csv, detected_transcripts.csv,
#           cell_metadata.csv, cell boundary HDF5 files)
# Outputs:  Seurat object saved as <out_dir>/<sample_name>_<method>.rds
# Usage:    Rscript 01_preprocessing/R/create_seurat_segmented_merscope.R \
#             --data_dir <path> --sample_name <name> --method <default|cellpose> \
#             --out_dir <path>

library(optparse)
library(Seurat)

source("utils/segmentation_utils.R")  # must be run from the project root

option_list <- list(
  make_option(c("--data_dir"),    type = "character", default = NULL,
              help = "Path to Vizgen output directory (containing cell_by_gene.csv etc.)"),
  make_option(c("--sample_name"), type = "character", default = NULL,
              help = "Sample identifier; used as the FOV name and in the output filename"),
  make_option(c("--method"),      type = "character", default = NULL,
              help = "Segmentation method label (e.g. default, cellpose)"),
  make_option(c("--out_dir"),     type = "character", default = NULL,
              help = "Output directory for the RDS file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$data_dir))    stop("--data_dir is required")
if (is.null(opt$sample_name)) stop("--sample_name is required")
if (is.null(opt$method))      stop("--method is required")
if (is.null(opt$out_dir))     stop("--out_dir is required")

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

message("Sample:    ", opt$sample_name)
message("Method:    ", opt$method)
message("Data dir:  ", opt$data_dir)
message("Output dir:", opt$out_dir)

start_time <- Sys.time()

obj <- myLoadVizgen(
  data.dir = opt$data_dir,
  fov      = opt$sample_name,
  assay    = "Vizgen"
)

elapsed <- round(Sys.time() - start_time, digits = 2)
message("Loaded in ", elapsed, " ", attr(elapsed, "units"))

out_file <- file.path(opt$out_dir, paste0(opt$sample_name, "_", opt$method, ".rds"))
saveRDS(obj, file = out_file)
message("Saved: ", out_file)
