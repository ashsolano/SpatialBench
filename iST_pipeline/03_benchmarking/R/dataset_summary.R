# Purpose:  Compute dataset-level summary metrics across VisiumHD, MERSCOPE,
#           and Xenium binning objects. Extracts per-sample bins, transcripts,
#           sparsity, and common-gene equivalents for each resolution.
# Inputs:   config/config.yaml  (visiumhd, spatial_analysis, bin_resolutions)
#           results/01_preprocessing/merscope_{res}um/{sample}_{res}um.rds
#           results/01_preprocessing/xenium_{res}um/{sample}_{res}um.rds
# Outputs:  results/03_benchmarking/dataset_summary/metrics.rds

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(yaml)
  library(optparse)
})

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--config"),  type = "character",
              default = "config/config.yaml",
              help    = "Path to config.yaml [default: %default]"),
  make_option(c("--out_dir"), type = "character",
              default = "results/03_benchmarking/dataset_summary",
              help    = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Proportion of zero entries in a count matrix
compute_sparsity <- function(count_matrix) {
  sum(count_matrix == 0) / length(count_matrix)
}

# Extract summary metrics from a VisiumHD Seurat object.
# The object carries both Spatial.008um and Spatial.016um assays, so both
# resolutions are computed in a single call.
extract_visium_hd_metrics <- function(seurat_obj, sample_name, common_genes) {
  meta     <- seurat_obj@meta.data

  bins_8   <- sum(!is.na(meta$nCount_Spatial.008um))
  bins_16  <- sum(!is.na(meta$nCount_Spatial.016um))
  genes_8  <- length(rownames(seurat_obj@assays[["Spatial.008um"]]))
  genes_16 <- length(rownames(seurat_obj@assays[["Spatial.016um"]]))
  tx_8     <- sum(meta$nCount_Spatial.008um,  na.rm = TRUE)
  tx_16    <- sum(meta$nCount_Spatial.016um, na.rm = TRUE)

  mat_8  <- GetAssayData(seurat_obj, assay = "Spatial.008um", slot = "counts")
  mat_16 <- GetAssayData(seurat_obj, assay = "Spatial.016um", slot = "counts")
  sp_8   <- compute_sparsity(mat_8)
  sp_16  <- compute_sparsity(mat_16)

  # Subset to common genes for cross-platform comparisons
  cg_8  <- intersect(common_genes, rownames(seurat_obj@assays[["Spatial.008um"]]))
  cg_16 <- intersect(common_genes, rownames(seurat_obj@assays[["Spatial.016um"]]))
  tx_cg_8   <- sum(mat_8[cg_8,  , drop = FALSE], na.rm = TRUE)
  tx_cg_16  <- sum(mat_16[cg_16, , drop = FALSE], na.rm = TRUE)
  sp_cg_8   <- compute_sparsity(mat_8[cg_8,  , drop = FALSE])
  sp_cg_16  <- compute_sparsity(mat_16[cg_16, , drop = FALSE])

  data.frame(
    Sample     = sample_name,
    Technology = "VisiumHD",
    Metric     = rep(c("Genes", "Bins", "Transcripts", "Sparsity",
                       "Transcripts (Common)", "Sparsity (Common)"), each = 2),
    Binning    = rep(c("16µm", "8µm"), times = 6),
    Value      = c(genes_16, genes_8,
                   bins_16,  bins_8,
                   tx_16,    tx_8,
                   sp_16,    sp_8,
                   tx_cg_16, tx_cg_8,
                   sp_cg_16, sp_cg_8)
  )
}

# Extract summary metrics from a MERSCOPE (Vizgen assay) binning Seurat object.
# `binning` is a display label, e.g. "8µm", passed from the caller.
extract_vizgen_metrics <- function(seurat_obj, sample_name, common_genes, binning) {
  meta  <- seurat_obj@meta.data
  bins  <- nrow(meta)
  genes <- length(rownames(seurat_obj@assays[["Vizgen"]]))
  tx    <- sum(meta$nCount_Vizgen, na.rm = TRUE)

  mat    <- GetAssayData(seurat_obj, assay = "Vizgen", slot = "counts")
  sp_all <- compute_sparsity(mat)

  cg     <- intersect(common_genes, rownames(seurat_obj@assays[["Vizgen"]]))
  mat_cg <- mat[cg, , drop = FALSE]
  tx_cg  <- sum(mat_cg, na.rm = TRUE)
  sp_cg  <- compute_sparsity(mat_cg)

  data.frame(
    Sample     = sample_name,
    Technology = "MERSCOPE",
    Metric     = c("Genes", "Bins", "Transcripts", "Sparsity",
                   "Transcripts (Common)", "Sparsity (Common)"),
    Binning    = rep(binning, times = 6),
    Value      = c(genes, bins, tx, sp_all, tx_cg, sp_cg)
  )
}

# Extract summary metrics from a Xenium binning Seurat object.
# `binning` is a display label, e.g. "16µm", passed from the caller.
extract_xenium_metrics <- function(seurat_obj, sample_name, common_genes, binning) {
  meta  <- seurat_obj@meta.data
  bins  <- nrow(meta)
  genes <- length(rownames(seurat_obj@assays[["Xenium"]]))
  tx    <- sum(meta$nCount_Xenium, na.rm = TRUE)

  mat    <- GetAssayData(seurat_obj, assay = "Xenium", slot = "counts")
  sp_all <- compute_sparsity(mat)

  cg     <- intersect(common_genes, rownames(seurat_obj@assays[["Xenium"]]))
  mat_cg <- mat[cg, , drop = FALSE]
  tx_cg  <- sum(mat_cg, na.rm = TRUE)
  sp_cg  <- compute_sparsity(mat_cg)

  data.frame(
    Sample     = sample_name,
    Technology = "Xenium",
    Metric     = c("Genes", "Bins", "Transcripts", "Sparsity",
                   "Transcripts (Common)", "Sparsity (Common)"),
    Binning    = rep(binning, times = 6),
    Value      = c(genes, bins, tx, sp_all, tx_cg, sp_cg)
  )
}

# ---------------------------------------------------------------------------
# Load VisiumHD samples
# ---------------------------------------------------------------------------
message("Loading VisiumHD samples...")

visiumhd_dir     <- cfg$visiumhd$data_dir
visiumhd_samples <- cfg$visiumhd$samples

visiumhd_objs <- lapply(names(visiumhd_samples), function(samp) {
  path <- file.path(visiumhd_dir, visiumhd_samples[[samp]])
  message("  ", samp, ": ", path)
  obj <- readRDS(path)
  DefaultAssay(obj) <- "Spatial.008um"
  obj
})
names(visiumhd_objs) <- names(visiumhd_samples)

# ---------------------------------------------------------------------------
# Load MERSCOPE binning objects
# ---------------------------------------------------------------------------
message("Loading MERSCOPE binning objects...")

merscope_samples <- cfg$spatial_analysis$merscope_samples
bin_resolutions  <- cfg$bin_resolutions

merscope_objs <- lapply(bin_resolutions, function(res) {
  bin_dir <- file.path(cfg$output_dir, "01_preprocessing", paste0("merscope_", res, "um"))
  lapply(merscope_samples, function(samp) {
    path <- file.path(bin_dir, paste0(samp, "_", res, "um.rds"))
    message("  ", samp, " @ ", res, "um: ", path)
    readRDS(path)
  }) |> setNames(merscope_samples)
}) |> setNames(as.character(bin_resolutions))

# ---------------------------------------------------------------------------
# Load Xenium binning objects
# ---------------------------------------------------------------------------
message("Loading Xenium binning objects...")

xenium_samples <- cfg$spatial_analysis$xenium_default_samples

xenium_objs <- lapply(bin_resolutions, function(res) {
  bin_dir <- file.path(cfg$output_dir, "01_preprocessing", paste0("xenium_", res, "um"))
  lapply(xenium_samples, function(samp) {
    path <- file.path(bin_dir, paste0(samp, "_", res, "um.rds"))
    message("  ", samp, " @ ", res, "um: ", path)
    readRDS(path)
  }) |> setNames(xenium_samples)
}) |> setNames(as.character(bin_resolutions))

# ---------------------------------------------------------------------------
# Determine common genes across all three platforms
# ---------------------------------------------------------------------------
message("Computing common gene set...")

# Use one representative object from each platform
ref_visiumhd <- visiumhd_objs[[1]]
ref_merscope <- merscope_objs[[as.character(bin_resolutions[1])]][[1]]
ref_xenium   <- xenium_objs[[as.character(bin_resolutions[1])]][[1]]

genes_visiumhd <- rownames(ref_visiumhd@assays[["Spatial.008um"]])
genes_merscope <- rownames(ref_merscope@assays[["Vizgen"]])
genes_xenium   <- rownames(ref_xenium@assays[["Xenium"]])
common_genes   <- Reduce(intersect, list(genes_visiumhd, genes_merscope, genes_xenium))
message("  Common genes: ", length(common_genes))

# ---------------------------------------------------------------------------
# Extract metrics
# ---------------------------------------------------------------------------
message("Extracting VisiumHD metrics...")
visiumhd_data <- dplyr::bind_rows(
  lapply(names(visiumhd_objs), function(samp) {
    extract_visium_hd_metrics(visiumhd_objs[[samp]], samp, common_genes)
  })
)

message("Extracting MERSCOPE metrics...")
merscope_data <- dplyr::bind_rows(
  lapply(bin_resolutions, function(res) {
    label <- paste0(res, "µm")
    lapply(merscope_samples, function(samp) {
      extract_vizgen_metrics(
        merscope_objs[[as.character(res)]][[samp]], samp, common_genes, label
      )
    })
  })
)

message("Extracting Xenium metrics...")
xenium_data <- dplyr::bind_rows(
  lapply(bin_resolutions, function(res) {
    label <- paste0(res, "µm")
    lapply(xenium_samples, function(samp) {
      extract_xenium_metrics(
        xenium_objs[[as.character(res)]][[samp]], samp, common_genes, label
      )
    })
  })
)

# ---------------------------------------------------------------------------
# Combine and save
# ---------------------------------------------------------------------------
combined_data <- dplyr::bind_rows(visiumhd_data, merscope_data, xenium_data)

message("Total rows: ", nrow(combined_data))

out_file <- file.path(opt$out_dir, "metrics.rds")
saveRDS(combined_data, out_file)
message("Saved: ", out_file)

# Save platform gene lists so downstream figure scripts (fig1.R) can build the
# Venn diagram without reloading the full Seurat objects.
gene_lists <- list(
  VisiumHD = genes_visiumhd,
  Xenium   = genes_xenium,
  MERSCOPE = genes_merscope
)
gene_lists_file <- file.path(opt$out_dir, "gene_lists.rds")
saveRDS(gene_lists, gene_lists_file)
message("Saved: ", gene_lists_file)
