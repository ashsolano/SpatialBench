# Purpose: GC B cell subclustering. Subsets cells annotated as "GC B cells"
#          from an annotated Seurat object and runs an independent embedding:
#          ScaleData -> PCA (pca.gc) -> Harmony (harmony.gc, by sample_name)
#          -> SNN graph -> Leiden clustering -> UMAP (umap.gc).
#          Cluster labels are stored in gc_subcluster.
# Inputs:  Annotated Seurat object from annotate_clusters.R (must contain
#          BANKSY assay, sample_name metadata, and cell_type == "GC B cells")
# Outputs: Subclustered Seurat object (GC B cells only) saved as <out_file>,
#          containing: pca.gc, harmony.gc, snn.gc, umap.gc, gc_subcluster
# Usage:   Rscript 02_spatial_analysis/R/subcluster_gc.R \
#            --input_rds <path> --method <name> \
#            --dims <int> --resolution <float> \
#            --umap_seed <int> --out_file <path>

library(optparse)
library(Seurat)
library(SeuratWrappers)
library(harmony)
library(future)

# Allow large global objects for Harmony and UMAP
options(future.globals.maxSize = 8000 * 1024^3)

option_list <- list(
  make_option(c("--input_rds"),  type = "character", default = NULL,
              help = "Path to annotated RDS from annotate_clusters.R"),
  make_option(c("--method"),     type = "character", default = NULL,
              help = "Method name (e.g. xenium_default); used for logging"),
  make_option(c("--dims"),       type = "integer",   default = NULL,
              help = "Dimensions for FindNeighbors and RunUMAP (1:dims)"),
  make_option(c("--resolution"), type = "numeric",   default = NULL,
              help = "Louvain clustering resolution"),
  make_option(c("--umap_seed"),  type = "integer",   default = NULL,
              help = "Random seed for UMAP reproducibility"),
  make_option(c("--out_file"),   type = "character", default = NULL,
              help = "Full path for the subclustered output RDS file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds))  stop("--input_rds is required")
if (is.null(opt$method))     stop("--method is required")
if (is.null(opt$dims))       stop("--dims is required")
if (is.null(opt$resolution)) stop("--resolution is required")
if (is.null(opt$umap_seed))  stop("--umap_seed is required")
if (is.null(opt$out_file))   stop("--out_file is required")

message("Method:      ", opt$method)
message("dims:        ", opt$dims)
message("resolution:  ", opt$resolution)
message("umap_seed:   ", opt$umap_seed)
message("Input RDS:   ", opt$input_rds)
message("Output file: ", opt$out_file)

# --- Load annotated object ---
message("\nLoading annotated object...")
obj <- readRDS(opt$input_rds)
message("  Loaded: ", ncol(obj), " cells total")

# --- Validate prerequisites ---
if (!"BANKSY" %in% Assays(obj)) {
  stop("BANKSY assay not found. Ensure preprocess_sample.R ran RunBanksy().")
}
if (!"sample_name" %in% colnames(obj@meta.data)) {
  stop("'sample_name' column not found. Ensure merge_samples.R attached sample metadata.")
}
if (!"cell_type" %in% colnames(obj@meta.data)) {
  stop("'cell_type' column not found. Ensure annotate_clusters.R ran successfully.")
}
if (!"GC B cells" %in% obj$cell_type) {
  stop(
    "No cells labelled 'GC B cells' found in cell_type. ",
    "Check annotations in config['annotations']['", opt$method, "']."
  )
}

# --- Subset to GC B cells ---
message("\nSubsetting to GC B cells...")
gc_subset <- subset(obj, subset = cell_type == "GC B cells")
message("  GC B cells: ", ncol(gc_subset))

# --- Set BANKSY as default assay ---
DefaultAssay(gc_subset) <- "BANKSY"

# --- Scale BANKSY features on the GC subset ---
# Rescale within the subset so that the PCA reflects GC B cell variation,
# not the full-object scale that was computed in embed_harmony.R
message("\nScaling BANKSY features...")
gc_subset <- ScaleData(gc_subset, assay = "BANKSY", verbose = FALSE)

# --- PCA on BANKSY features ---
# npcs fixed at 30 (same as full-object PCA) to retain sufficient variance;
# gc_dims (default 15) controls how many of those PCs are used downstream
message("Running PCA (npcs = 30)...")
feats <- rownames(gc_subset[["BANKSY"]])
gc_subset <- RunPCA(
  gc_subset,
  assay          = "BANKSY",
  features       = feats,
  npcs           = 30,
  reduction.name = "pca.gc",
  verbose        = FALSE
)

# --- Harmony batch correction by sample ---
message("Running Harmony (group by sample_name, dims = 1:", opt$dims, ")...")
gc_subset <- RunHarmony(
  gc_subset,
  group.by.vars  = "sample_name",
  reduction      = "pca.gc",
  reduction.save = "harmony.gc",
  assay.use      = "BANKSY",
  project.dim    = FALSE,
  verbose        = FALSE
)

# --- Shared nearest-neighbour graph on Harmony embedding ---
message("Finding neighbours (dims = 1:", opt$dims, ")...")
gc_subset <- FindNeighbors(
  gc_subset,
  reduction  = "harmony.gc",
  dims       = 1:opt$dims,
  graph.name = "snn.gc",
  verbose    = FALSE
)

# --- Leiden clustering ---
message("Clustering (resolution = ", opt$resolution, ")...")
gc_subset <- FindClusters(
  gc_subset,
  graph.name = "snn.gc",
  resolution = opt$resolution,
  algorithm  = 4,
  method     = "igraph",
  verbose    = FALSE
)

# Store under a named column; seurat_clusters is overwritten by later
# FindClusters calls so the named column is safer for downstream use
gc_subset$gc_subcluster <- gc_subset$seurat_clusters
message("  GC subclusters found: ", nlevels(gc_subset$seurat_clusters))

# --- UMAP on Harmony embedding ---
message("Running UMAP (seed = ", opt$umap_seed, ")...")
gc_subset <- RunUMAP(
  gc_subset,
  reduction      = "harmony.gc",
  dims           = 1:opt$dims,
  reduction.name = "umap.gc",
  seed.use       = opt$umap_seed,
  verbose        = FALSE
)

# --- Save ---
dir.create(dirname(opt$out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(gc_subset, file = opt$out_file)
message("\nSaved: ", opt$out_file)
message("Done.")
