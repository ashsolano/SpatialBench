# Purpose:  Dimensionality reduction, Harmony batch correction by sample, and
#           Louvain clustering for a merged spatial dataset. Operates on the
#           BANKSY assay produced by preprocess_sample.R.
# Inputs:   Merged Seurat object from merge_samples.R (must contain BANKSY assay
#           and sample_name metadata column)
# Outputs:  Embedded and clustered Seurat object saved as <out_file>, containing:
#             pca_banksy  — PCA on BANKSY features
#             h_sampleid           — Harmony-corrected embedding (grouped by sample_name)
#             snn_harmony          — shared nearest-neighbour graph
#             umap_harmony         — UMAP on Harmony embedding
#             cluster_harmony_sampleid — Louvain cluster assignments
# Usage:    Rscript 02_spatial_analysis/R/embed_harmony.R \
#             --input_rds <path> --method <xenium_default|...> \
#             --npcs <int> --dims <int> --resolution <float> \
#             --umap_seed <int> --out_file <path>


library(optparse)
library(Seurat)
library(SeuratWrappers)
library(harmony)
library(future)

# Allow large global objects when running Harmony and UMAP
options(future.globals.maxSize = 8000 * 1024^3)

option_list <- list(
  make_option(c("--input_rds"),  type = "character", default = NULL,
              help = "Path to merged RDS file from merge_samples.R"),
  make_option(c("--method"),     type = "character", default = NULL,
              help = "Method name (e.g. xenium_default); used for logging"),
  make_option(c("--npcs"),       type = "integer",   default = NULL,
              help = "Number of principal components for PCA"),
  make_option(c("--dims"),       type = "integer",   default = NULL,
              help = "Dimensions used for Harmony, FindNeighbors, and UMAP (1:dims)"),
  make_option(c("--resolution"), type = "numeric",   default = NULL,
              help = "Louvain clustering resolution"),
  make_option(c("--umap_seed"),  type = "integer",   default = NULL,
              help = "Random seed for UMAP reproducibility"),
  make_option(c("--out_file"),   type = "character", default = NULL,
              help = "Full path for the embedded output RDS file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds))  stop("--input_rds is required")
if (is.null(opt$method))     stop("--method is required")
if (is.null(opt$npcs))       stop("--npcs is required")
if (is.null(opt$dims))       stop("--dims is required")
if (is.null(opt$resolution)) stop("--resolution is required")
if (is.null(opt$umap_seed))  stop("--umap_seed is required")
if (is.null(opt$out_file))   stop("--out_file is required")

message("Method:      ", opt$method)
message("npcs:        ", opt$npcs)
message("dims:        ", opt$dims)
message("resolution:  ", opt$resolution)
message("umap_seed:   ", opt$umap_seed)
message("Input RDS:   ", opt$input_rds)
message("Output file: ", opt$out_file)

# --- Load merged object ---
message("\nLoading merged object...")
obj <- readRDS(opt$input_rds)
message("  Loaded: ", ncol(obj), " cells, ", nrow(obj), " features")

# --- Validate prerequisites ---
if (!"BANKSY" %in% Assays(obj)) {
  stop(
    "BANKSY assay not found in merged object. ",
    "Ensure preprocess_sample.R ran RunBanksy() for all samples."
  )
}
if (!"sample_name" %in% colnames(obj@meta.data)) {
  stop(
    "'sample_name' column not found in metadata. ",
    "Ensure preprocess_sample.R attached sample metadata correctly."
  )
}

# --- Set BANKSY as default assay for all downstream steps ---
DefaultAssay(obj) <- "BANKSY"

# --- Scale BANKSY features across the full merged object ---
# Must be done after merging because per-sample scaling in preprocess_sample.R
# only scaled within each sample; here we scale across all cells jointly
message("\nScaling BANKSY features...")
obj <- ScaleData(obj, assay = "BANKSY", verbose = FALSE)

# --- PCA on BANKSY features ---
message("Running PCA (npcs = ", opt$npcs, ")...")
feats <- rownames(obj[["BANKSY"]])
obj <- RunPCA(
  obj,
  assay          = "BANKSY",
  features       = feats,
  npcs           = opt$npcs,
  reduction.name = "pca_banksy",
  verbose        = FALSE
)

# --- Harmony batch correction by sample ---
# Corrects for sample-to-sample technical variation while preserving
# biological signal; grouped by sample_name (one correction per sample)
message("Running Harmony (group by sample_name, dims = 1:", opt$dims, ")...")
obj <- RunHarmony(
  obj,
  group.by.vars  = "sample_name",
  reduction      = "pca_banksy",
  reduction.save = "h_sampleid",
  assay.use      = "BANKSY",
  project.dim    = FALSE,
  verbose        = FALSE
)

# --- Shared nearest-neighbour graph on Harmony embedding ---
message("Finding neighbours (dims = 1:", opt$dims, ")...")
obj <- FindNeighbors(
  obj,
  reduction  = "h_sampleid",
  dims       = 1:opt$dims,
  graph.name = "snn_harmony",
  verbose    = FALSE
)

# --- Louvain clustering (algorithm 4 = igraph implementation) ---
message("Clustering (resolution = ", opt$resolution, ")...")
obj <- FindClusters(
  obj,
  graph.name = "snn_harmony",
  resolution = opt$resolution,
  algorithm  = 4,
  method     = "igraph",
  verbose    = FALSE
)

# Store clusters under a descriptive name; seurat_clusters is overwritten
# by subsequent FindClusters calls, so the named column is safer for downstream use
obj$cluster_harmony_sampleid <- obj$seurat_clusters
message("  Clusters found: ", nlevels(obj$seurat_clusters))

# --- UMAP on Harmony embedding ---
message("Running UMAP (seed = ", opt$umap_seed, ")...")
obj <- RunUMAP(
  obj,
  reduction      = "h_sampleid",
  dims           = 1:opt$dims,
  reduction.name = "umap_harmony",
  seed.use       = opt$umap_seed,
  verbose        = FALSE
)

# --- Save ---
dir.create(dirname(opt$out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, file = opt$out_file)
message("\nSaved: ", opt$out_file)
message("Done.")
