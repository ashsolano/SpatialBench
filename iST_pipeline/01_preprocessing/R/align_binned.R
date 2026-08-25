# Purpose:  Fit and apply an STalign registration (affine + LDDMM), mapping a
#           filtered, binned Seurat object's bin centroids into VisiumHD
#           target coordinate space. The source (platform) and target
#           (VisiumHD H&E) rasterised images are computed fresh from the
#           filtered Seurat object and the VisiumHD reference object on every
#           run — no static intermediate image files are read. Landmark
#           points were picked manually (see STalign_reticulate.Rmd) and are
#           hardcoded in config.yaml (stalign.landmarks) — this script
#           performs no interactive landmark picking. The actual STalign
#           computation runs in align_binned.py (called via system2()),
#           since STalign is only available in a Python virtualenv.
# Inputs:   Filtered binned Seurat object RDS (output of filter_binned.R, any
#           bin resolution); config.yaml's `stalign` section, which supplies
#           the Python interpreter, the animal-ID -> sample-key mapping per
#           platform, per-sample landmark points, and the VisiumHD low-res
#           scale factor per sample; config.yaml's `visiumhd` section, which
#           supplies the VisiumHD reference RDS directory and per-sample
#           filenames
# Outputs:  Aligned Seurat object saved to --out_rds
# Usage:    Rscript 01_preprocessing/R/align_binned.R \
#             --input_rds <path> --platform <xenium|merscope> \
#             --sample <animal ID, e.g. ko167> --resolution <8|16> \
#             --config config/config.yaml --out_rds <path>

library(optparse)
library(Seurat)
library(yaml)
library(png)

option_list <- list(
  make_option(c("--input_rds"),  type = "character", default = NULL,
              help = "Path to filtered binned Seurat object RDS"),
  make_option(c("--platform"),   type = "character", default = NULL,
              help = "Platform: 'xenium' or 'merscope' — selects the source raster and assay"),
  make_option(c("--sample"),     type = "character", default = NULL,
              help = "Animal ID, e.g. 'wt709', 'ko167' — looked up in config$stalign$matched_samples/visiumhd"),
  make_option(c("--resolution"), type = "integer",   default = NULL,
              help = "Bin size in microns (8 or 16) — used only for logging/temp-file naming"),
  make_option(c("--config"),     type = "character", default = "config/config.yaml",
              help = "Path to config.yaml [default: %default]"),
  make_option(c("--out_rds"),    type = "character", default = NULL,
              help = "Output path for the aligned Seurat object RDS")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds))  stop("--input_rds is required")
if (is.null(opt$platform))   stop("--platform is required")
if (is.null(opt$sample))     stop("--sample is required")
if (is.null(opt$resolution)) stop("--resolution is required")
if (is.null(opt$out_rds))    stop("--out_rds is required")

if (!opt$platform %in% c("xenium", "merscope")) {
  stop("--platform must be 'xenium' or 'merscope', got: ", opt$platform)
}

config <- yaml::read_yaml(opt$config)

python_bin    <- config$stalign$python_bin
assay         <- config$spatial_analysis$platforms[[opt$platform]]$assay

full_sample <- config$stalign$matched_samples[[opt$platform]][[opt$sample]]
if (is.null(full_sample)) {
  stop("No config$stalign$matched_samples entry for platform '", opt$platform,
       "', sample '", opt$sample, "'")
}

visium_sample <- config$stalign$visiumhd[[opt$sample]]
if (is.null(visium_sample)) {
  stop("No config$stalign$visiumhd entry for sample '", opt$sample, "'")
}

lowres_factor <- config$stalign$visiumhd_scalefactors[[visium_sample]]
if (is.null(lowres_factor)) {
  stop("No config$stalign$visiumhd_scalefactors entry for VisiumHD sample '", visium_sample, "'")
}

dir.create(dirname(opt$out_rds), recursive = TRUE, showWarnings = FALSE)

message("Sample:            ", opt$sample, " (", full_sample, ")")
message("Platform:          ", opt$platform)
message("Resolution:        ", opt$resolution, "um")
message("Assay:             ", assay)
message("Input RDS:         ", opt$input_rds)
message("VisiumHD sample:    ", visium_sample)
message("tissue_lowres_scalef applied for ", visium_sample, ": ", lowres_factor)

obj <- readRDS(opt$input_rds)

# --- Load the VisiumHD reference object and extract its H&E image -------------
# Rasterised fresh on every run (rather than reusing a static .npz) so a
# stale target image can never silently diverge from the current VisiumHD
# reference object.
visium_rds_name <- config$visiumhd$samples[[visium_sample]]
if (is.null(visium_rds_name)) {
  stop("No config$visiumhd$samples entry for VisiumHD sample '", visium_sample, "'")
}
visium_rds <- file.path(config$visiumhd$data_dir, visium_rds_name)
if (!file.exists(visium_rds)) {
  stop("VisiumHD reference RDS not found: ", visium_rds)
}

message("Loading VisiumHD reference: ", visium_rds)
visium_obj <- readRDS(visium_rds)

he_img <- visium_obj@images[["slice1.008um"]]@image
message("VisiumHD H&E image dimensions: ", paste(dim(he_img), collapse = " x "))

# png::writePNG() requires values in [0, 1]; Seurat's @image slot is expected
# to already be scaled this way. Fail loudly rather than silently
# clamping/normalising if that assumption ever breaks.
stopifnot(all(he_img >= 0 & he_img <= 1, na.rm = TRUE))

# --- Extract this object's bin centroids (any resolution) ---------------------
fov_name <- names(obj@images)[1]
fov      <- obj@images[[fov_name]]
cents    <- fov@boundaries$centroids@coords
rownames(cents) <- colnames(obj)

ncount_col <- paste0("nCount_", assay)
if (!ncount_col %in% colnames(obj@meta.data)) {
  stop("Metadata column '", ncount_col, "' not found; check --platform")
}

coords_df <- data.frame(
  cell_id = rownames(cents),
  x       = cents[, 1],
  y       = cents[, 2],
  nCount  = obj@meta.data[rownames(cents), ncount_col]
)

# --- Landmarks from config -> CSV for the Python step --------------------------
# Landmark points are hardcoded in config$stalign$landmarks (one [y, x] pair
# per row, per sample/platform) — column order is (y, x) throughout, do not
# swap to (x, y) here or in align_binned.py.
landmarks_entry <- config$stalign$landmarks[[opt$sample]][[opt$platform]]
if (is.null(landmarks_entry)) {
  stop("No config$stalign$landmarks entry for sample '", opt$sample,
       "', platform '", opt$platform, "'")
}

points_source <- do.call(rbind, lapply(landmarks_entry$source, as.numeric))
points_target <- do.call(rbind, lapply(landmarks_entry$target, as.numeric))
colnames(points_source) <- c("y", "x")
colnames(points_target) <- c("y", "x")

# --- Write temp CSVs/PNG, run align_binned.py, read the result back -----------
tmp_dir              <- tempdir()
tag                   <- paste0(full_sample, "_", opt$resolution, "um")
coords_csv            <- file.path(tmp_dir, paste0(tag, "_coords.csv"))
source_landmarks_csv  <- file.path(tmp_dir, paste0(tag, "_source_landmarks.csv"))
target_landmarks_csv  <- file.path(tmp_dir, paste0(tag, "_target_landmarks.csv"))
target_image_png      <- file.path(tmp_dir, paste0(tag, "_target_image.png"))
output_csv             <- file.path(tmp_dir, paste0(tag, "_aligned_coords.csv"))

write.csv(coords_df, coords_csv, row.names = FALSE)
write.csv(as.data.frame(points_source), source_landmarks_csv, row.names = FALSE)
write.csv(as.data.frame(points_target), target_landmarks_csv, row.names = FALSE)
png::writePNG(he_img, target_image_png)

align_script <- "01_preprocessing/Python/align_binned.py"  # must be run from the project root

message("Running STalign registration (affine + LDDMM)...")
status <- system2(
  python_bin,
  args = c(
    align_script,
    "--target_image",     target_image_png,
    "--source_landmarks", source_landmarks_csv,
    "--target_landmarks", target_landmarks_csv,
    "--coords_csv",       coords_csv,
    "--output_csv",       output_csv
  )
)

if (status != 0) {
  stop("align_binned.py failed with exit status ", status)
}

aligned_df <- read.csv(output_csv)
aligned_df <- aligned_df[match(coords_df$cell_id, aligned_df$cell_id), ]

if (anyNA(aligned_df$cell_id)) {
  stop("Some bins were not returned by align_binned.py; alignment output is incomplete")
}

# --- Rescale into the VisiumHD low-res pixel space -----------------------------
aligned_df$x <- aligned_df$x / lowres_factor
aligned_df$y <- aligned_df$y / lowres_factor

# --- Write the aligned coordinates back into the object's FOV ------------------
new_coords <- as.matrix(aligned_df[, c("y", "x")])
colnames(new_coords) <- c("x", "y")
rownames(new_coords) <- aligned_df$cell_id

obj@images[[fov_name]]@boundaries$centroids@coords <- new_coords

unlink(c(coords_csv, source_landmarks_csv, target_landmarks_csv, target_image_png, output_csv))

saveRDS(obj, file = opt$out_rds)
message("Saved: ", opt$out_rds)
