# Purpose:  Merge all per-sample preprocessed Seurat objects for one segmentation
#           method into a single combined object ready for embedding and clustering.
# Inputs:   Directory of per-sample RDS files produced by preprocess_sample.R
# Outputs:  Merged Seurat object saved as <out_file>
# Usage:    Rscript 02_spatial_analysis/R/merge_samples.R \
#             --input_dir <path> --method <xenium_default|...> \
#             --out_file <path>


library(optparse)
library(Seurat)
library(future)

# Allow large global objects when merging many samples
options(future.globals.maxSize = 8000 * 1024^3)

option_list <- list(
  make_option(c("--input_dir"), type = "character", default = NULL,
              help = "Directory containing per-sample preprocessed RDS files"),
  make_option(c("--method"),    type = "character", default = NULL,
              help = "Method name (e.g. xenium_default); used to label the merged project"),
  make_option(c("--out_file"),  type = "character", default = NULL,
              help = "Full path for the merged output RDS file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_dir)) stop("--input_dir is required")
if (is.null(opt$method))    stop("--method is required")
if (is.null(opt$out_file))  stop("--out_file is required")

message("Method:      ", opt$method)
message("Input dir:   ", opt$input_dir)
message("Output file: ", opt$out_file)

# --- Discover per-sample RDS files ---
rds_files <- sort(list.files(opt$input_dir, pattern = "\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) {
  stop("No .rds files found in: ", opt$input_dir)
}
message("\nFound ", length(rds_files), " sample file(s):")
for (f in rds_files) message("  ", f)

# --- Load each sample ---
so_list <- lapply(rds_files, function(f) {
  message("Loading: ", basename(f))
  readRDS(f)
})

# --- Validate required metadata columns ---
# All three columns are expected to be set by preprocess_sample.R
required_cols <- c("sample_name", "platform", "seg")
for (i in seq_along(so_list)) {
  missing_cols <- setdiff(required_cols, colnames(so_list[[i]]@meta.data))
  if (length(missing_cols) > 0) {
    stop(
      "Sample file '", basename(rds_files[i]), "' is missing metadata columns: ",
      paste(missing_cols, collapse = ", "), ". ",
      "Ensure preprocess_sample.R ran successfully for this sample."
    )
  }
}

# --- Get sample names for cell ID prefixes ---
# Using the sample_name metadata column ensures prefixes match the biological sample
sample_names <- sapply(so_list, function(so) unique(so$sample_name))
names(so_list) <- sample_names
message("\nMerging ", length(sample_names), " samples:")
for (s in sample_names) message("  ", s)

# --- Merge all objects into one ---
# add.cell.ids prepends sample_name to each cell barcode to ensure uniqueness
# merge.data = TRUE preserves the normalised data slot across all samples
merged <- merge(
  so_list[[1]],
  y            = so_list[-1],
  add.cell.ids = sample_names,
  project      = paste0(opt$method, "_banksy"),
  merge.data   = TRUE
)

message("\nMerged object: ", ncol(merged), " cells, ", nrow(merged), " features")

# --- Save ---
dir.create(dirname(opt$out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(merged, file = opt$out_file)
message("Saved: ", opt$out_file)
message("Done.")
