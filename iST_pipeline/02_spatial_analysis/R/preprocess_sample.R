# Purpose:  Per-sample preprocessing for spatial analysis: QC filtering,
#           normalisation, scaling, and BANKSY spatial feature extraction.
#           Takes a Seurat object from 01_preprocessing as input.
# Inputs:   RDS file produced by 01_preprocessing (cell-segmented Seurat object);
#           platform, segmentation method, and analysis parameters passed as CLI args
# Outputs:  Seurat object with BANKSY assay saved as <out_dir>/<sample_name>.rds
# Usage:    Rscript 02_spatial_analysis/R/preprocess_sample.R \
#             --input_rds <path> --sample_name <name> \
#             --platform <xenium|merscope> --seg <default|cellpose|proseg> \
#             --assay <Xenium|Vizgen> \
#             --qc_min_counts <int> --scale_factor <int> \
#             --banksy_lambda <float> --banksy_k_geom <int> \
#             --out_dir <path>


library(optparse)
library(Seurat)
library(SeuratWrappers)
library(Banksy)
library(future)

# Allow large global objects when using future-backed parallelism (e.g. RunBanksy)
options(future.globals.maxSize = 8000 * 1024^3)

option_list <- list(
  make_option(c("--input_rds"),     type = "character", default = NULL,
              help = "Path to input RDS file from 01_preprocessing"),
  make_option(c("--sample_name"),   type = "character", default = NULL,
              help = "Sample identifier (stored as metadata and used in output filename)"),
  make_option(c("--platform"),      type = "character", default = NULL,
              help = "Sequencing platform: xenium or merscope"),
  make_option(c("--seg"),           type = "character", default = NULL,
              help = "Segmentation method: default, cellpose, or proseg"),
  make_option(c("--assay"),         type = "character", default = NULL,
              help = "Seurat assay name: Xenium (Xenium) or Vizgen (MERSCOPE)"),
  make_option(c("--qc_min_counts"), type = "integer",   default = NULL,
              help = "Minimum transcript count per cell (strict greater-than threshold)"),
  make_option(c("--scale_factor"),  type = "numeric",   default = NULL,
              help = "Scale factor for LogNormalize"),
  make_option(c("--banksy_lambda"), type = "numeric",   default = NULL,
              help = "BANKSY lambda: spatial weight (0 = non-spatial, 1 = fully spatial)"),
  make_option(c("--banksy_k_geom"), type = "integer",   default = NULL,
              help = "BANKSY k_geom: number of geometric neighbours"),
  make_option(c("--out_dir"),       type = "character", default = NULL,
              help = "Output directory for the processed RDS file")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Validate required arguments ---
if (is.null(opt$input_rds))     stop("--input_rds is required")
if (is.null(opt$sample_name))   stop("--sample_name is required")
if (is.null(opt$platform))      stop("--platform is required")
if (is.null(opt$seg))           stop("--seg is required")
if (is.null(opt$assay))         stop("--assay is required")
if (is.null(opt$qc_min_counts)) stop("--qc_min_counts is required")
if (is.null(opt$scale_factor))  stop("--scale_factor is required")
if (is.null(opt$banksy_lambda)) stop("--banksy_lambda is required")
if (is.null(opt$banksy_k_geom)) stop("--banksy_k_geom is required")
if (is.null(opt$out_dir))       stop("--out_dir is required")

if (!opt$assay %in% c("Xenium", "Vizgen")) {
  stop("--assay must be 'Xenium' or 'Vizgen'")
}
if (!opt$platform %in% c("xenium", "merscope")) {
  stop("--platform must be 'xenium' or 'merscope'")
}
if (!opt$seg %in% c("default", "cellpose", "proseg")) {
  stop("--seg must be 'default', 'cellpose', or 'proseg'")
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

message("Sample:        ", opt$sample_name)
message("Platform:      ", opt$platform)
message("Seg method:    ", opt$seg)
message("Assay:         ", opt$assay)
message("QC min counts: ", opt$qc_min_counts)
message("Scale factor:  ", opt$scale_factor)
message("BANKSY lambda: ", opt$banksy_lambda)
message("BANKSY k_geom: ", opt$banksy_k_geom)
message("Input RDS:     ", opt$input_rds)
message("Output dir:    ", opt$out_dir)

# --- Load Seurat object ---
message("\nLoading RDS...")
obj <- readRDS(opt$input_rds)
message("  Loaded: ", ncol(obj), " cells, ", nrow(obj), " features")

# --- Attach sample metadata ---
# Store platform, seg method, and sample name as separate columns so that
# downstream filtering and ggplot faceting can use each dimension independently
obj$sample_name <- opt$sample_name
obj$platform    <- opt$platform
obj$seg         <- opt$seg

# Derive condition from the sample_name prefix: wt / ko / ctrl
obj$condition <- ifelse(startsWith(opt$sample_name, "wt"), "wt",
                 ifelse(startsWith(opt$sample_name, "ko"), "ko", "ctrl"))

# Derive animal_id by stripping the batch and region suffix
# Handles both single-underscore (wt709_batch27) and double-underscore (wt709__batch34__0032118)
obj$animal_id <- sub("_{1,2}batch.*", "", opt$sample_name)

# --- Set default assay to the platform's native assay ---
DefaultAssay(obj) <- opt$assay

# --- QC filter: minimum transcript count ---
# The count column name is nCount_<assay> (e.g. nCount_Xenium, nCount_Vizgen)
ncount_col <- paste0("nCount_", opt$assay)
if (!ncount_col %in% colnames(obj@meta.data)) {
  stop(
    "QC column '", ncount_col, "' not found in metadata. ",
    "Check that --assay matches the assay name in the input RDS."
  )
}

n_before <- ncol(obj)
keep     <- colnames(obj)[obj@meta.data[[ncount_col]] > opt$qc_min_counts]
obj      <- subset(obj, cells = keep)
n_after  <- ncol(obj)
message(
  "\nQC filter (", ncount_col, " > ", opt$qc_min_counts, "): ",
  n_before, " -> ", n_after, " cells (removed ", n_before - n_after, ")"
)

# --- Normalise ---
message("\nNormalising (LogNormalize, scale.factor = ", opt$scale_factor, ")...")
obj <- NormalizeData(
  obj,
  normalization.method = "LogNormalize",
  scale.factor         = opt$scale_factor,
  verbose              = FALSE
)

# --- Scale all features in the native assay ---
# Populates scale.data for the native assay; RunBanksy reads the data slot
message("Scaling data...")
obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)

# --- Run BANKSY spatial feature extraction ---
# Creates a new BANKSY assay combining gene expression with neighbourhood context
message(
  "\nRunning BANKSY (lambda = ", opt$banksy_lambda,
  ", k_geom = ", opt$banksy_k_geom, ")..."
)
obj <- RunBanksy(
  obj,
  assay    = opt$assay,
  slot     = "data",
  features = rownames(obj[[opt$assay]]),
  lambda   = opt$banksy_lambda,
  k_geom   = opt$banksy_k_geom,
  verbose  = FALSE
)

# --- Save ---
out_file <- file.path(opt$out_dir, paste0(opt$sample_name, ".rds"))
saveRDS(obj, file = out_file)
message("\nSaved: ", out_file)
message("Done.")
