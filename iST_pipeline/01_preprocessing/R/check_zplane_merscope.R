# Purpose:  Diagnostic: quantify what a single-plane (z = 3) filter loses versus
#           all z-planes in LoadVizgen_binned (utils/binning_utils.R). Bins one
#           sample twice in memory (z = 3 and z = "all") and compares them.
# Inputs:   --data_dir     MERSCOPE data directory (containing detected_transcripts.csv)
#           --sample_name  Sample identifier (used as FOV name)
# Outputs:  Console report: total transcripts / median nCount / median nFeature /
#           number of bins for z = 3 vs all-z, plus raw per-z-plane transcript counts.
# Usage:    Rscript 01_preprocessing/R/check_zplane_merscope.R \
#             --data_dir <path> --sample_name ctrl172_batch10
#           Run from the project root (sources utils/binning_utils.R).

library(optparse)
library(Seurat)
library(Matrix)
library(data.table)

source("utils/binning_utils.R")

option_list <- list(
  make_option(c("--data_dir"),    type = "character", default = NULL,
              help = "Path to MERSCOPE data directory (containing detected_transcripts.csv)"),
  make_option(c("--sample_name"), type = "character", default = NULL,
              help = "Sample identifier; used as the FOV name"),
  make_option(c("--resolution"),  type = "integer",   default = 8L,
              help = "Bin size in microns [default: %default]"),
  make_option(c("--assay"),       type = "character", default = "Vizgen",
              help = "Assay name to summarise [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$data_dir))    stop("--data_dir is required")
if (is.null(opt$sample_name)) stop("--sample_name is required")

message("Sample:     ", opt$sample_name)
message("Data dir:   ", opt$data_dir)
message("Resolution: ", opt$resolution, "um")

# --- Summarise one binned object --------------------------------------------
summarise_obj <- function(obj, assay) {
  counts <- GetAssayData(obj, assay = assay, layer = "counts")

  ncount_col   <- paste0("nCount_", assay)
  nfeature_col <- paste0("nFeature_", assay)

  ncount <- if (ncount_col %in% colnames(obj@meta.data)) {
    obj@meta.data[[ncount_col]]
  } else {
    Matrix::colSums(counts)
  }

  nfeature <- if (nfeature_col %in% colnames(obj@meta.data)) {
    obj@meta.data[[nfeature_col]]
  } else {
    Matrix::colSums(counts > 0)
  }

  list(
    total_transcripts = sum(counts),
    median_ncount      = median(ncount),
    median_nfeature    = median(nfeature),
    n_cells            = ncol(obj)
  )
}

# --- Bin with z = 3 and z = "all" -------------------------------------------
message("\nRunning LoadVizgen_binned with z = 3 (single plane, for comparison)...")
obj_z3 <- LoadVizgen_binned(
  data.dir   = opt$data_dir,
  resolution = opt$resolution,
  fov        = opt$sample_name,
  assay      = opt$assay,
  z          = 3L
)

message("Running LoadVizgen_binned with z = \"all\" (current default)...")
obj_all <- LoadVizgen_binned(
  data.dir   = opt$data_dir,
  resolution = opt$resolution,
  fov        = opt$sample_name,
  assay      = opt$assay,
  z          = "all"
)

summary_z3  <- summarise_obj(obj_z3,  opt$assay)
summary_all <- summarise_obj(obj_all, opt$assay)

comparison <- data.frame(
  metric = c("total_transcripts", "median_nCount", "median_nFeature", "n_cells_bins"),
  z3     = c(summary_z3$total_transcripts,  summary_z3$median_ncount,  summary_z3$median_nfeature,  summary_z3$n_cells),
  all_z  = c(summary_all$total_transcripts, summary_all$median_ncount, summary_all$median_nfeature, summary_all$n_cells)
)
comparison$diff       <- comparison$all_z - comparison$z3
comparison$pct_change <- round(100 * comparison$diff / comparison$z3, 1)

cat("\n===== Binning comparison (z=3 vs. all-z):", opt$sample_name, "=====\n\n")
print(comparison, row.names = FALSE)

# --- Raw transcripts per z-plane -------------------------------------------
transcripts_file <- file.path(opt$data_dir, "detected_transcripts.csv")
if (!file.exists(transcripts_file)) stop("Transcripts file not found: ", transcripts_file)
mx <- fread(transcripts_file, sep = ",", verbose = FALSE)

if (!"global_z" %in% colnames(mx)) stop("'global_z' column not found in ", transcripts_file)

z_counts <- mx[, .N, by = global_z][order(global_z)]

cat("\n===== Raw transcripts per z-plane:", opt$sample_name, "=====\n\n")
cat("Unique z-planes in raw transcripts file:", nrow(z_counts), "\n")
cat("Total transcripts in raw file:          ", nrow(mx), "\n\n")
print(z_counts)

cat("\nConclusion: compare the raw per-z-plane counts above against the z3 vs.\n")
cat("all-z totals in the comparison table — this shows how many transcripts\n")
cat("(and how much of the resulting nCount/nFeature/cell yield) is discarded\n")
cat("by hardcoding z = 3 instead of using all z-planes.\n")
