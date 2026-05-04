# Purpose:  Compute per-bin QC metadata (nCount, nFeature) for VisiumHD,
#           MERSCOPE, and Xenium at all configured bin resolutions, for both
#           the full platform gene set and the three-platform common gene subset.
#           Also computes per-cell intersect-gene metadata for VisiumHD, FLEX
#           snRNA-seq, and scRNA-seq (scGEM), used for Figure 2 FLEX comparison.
# Inputs:   config/config.yaml  (visiumhd, spatial_analysis, bin_resolutions)
#           results/01_preprocessing/merscope_{res}um/{sample}_{res}um.rds
#           results/01_preprocessing/xenium_{res}um/{sample}_{res}um.rds
#           results/03_benchmarking/dataset_summary/gene_lists.rds  (common genes)
#           cfg$scrna$path  (scFlex_seu.rds — FLEX snRNA-seq reference)
#           --sc_rds        (scGEM_seu.rds  — scRNA-seq; suggest adding to config.yaml
#                            as scrna$sc_path for consistency)
# Outputs:  results/03_benchmarking/qc_metrics/metadata_combined.rds
#               (per-bin rows: Sample, nCount, nFeature, platform, Subset, bin_size)
#           results/03_benchmarking/qc_metrics/metadata_flex_scrna.rds
#               (per-cell rows: barcode, nCount_intersect, nFeature_intersect,
#                Sample, platform)
#           results/03_benchmarking/qc_metrics/genes_intersect_flex.rds

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(yaml)
  library(optparse)
})

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--config"),     type = "character",
              default = "config/config.yaml",
              help    = "Path to config.yaml [default: %default]"),
  make_option(c("--out_dir"),    type = "character",
              default = "results/03_benchmarking/qc_metrics",
              help    = "Output directory [default: %default]"),
  make_option(c("--gene_lists"), type = "character",
              default = "results/03_benchmarking/dataset_summary/gene_lists.rds",
              help    = "Path to gene_lists.rds from dataset_summary.R [default: %default]"),
  make_option(c("--sc_rds"),     type = "character",
              default = NULL,
              help    = "Path to scGEM_seu.rds (10x scRNA-seq). Flag: add as scrna$sc_path in config.yaml")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opt$gene_lists)) stop("--gene_lists not found: ", opt$gene_lists)

cfg <- yaml::read_yaml(opt$config)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Seurat v4/v5 compatible raw counts extraction
get_counts_mat <- function(so, assay) {
  tryCatch(
    GetAssayData(so, assay = assay, layer  = "counts"),
    error = function(e) GetAssayData(so, assay = assay, slot = "counts")
  )
}

# Build per-bin metadata for one platform / bin-size / gene-subset combination.
# gene_subset = NULL uses all genes; otherwise restricts to the supplied vector.
# Returns a data frame with one row per bin/cell.
build_meta_platform <- function(obj_list, sample_names, assay_name,
                                platform_tag, subset_tag, bin_size_label,
                                gene_subset = NULL) {
  do.call(rbind, lapply(sample_names, function(s) {
    obj     <- obj_list[[s]]
    mat     <- get_counts_mat(obj, assay_name)
    keep_g  <- if (is.null(gene_subset)) rownames(mat) else intersect(gene_subset, rownames(mat))
    sub_mat <- mat[keep_g, , drop = FALSE]
    data.frame(
      Sample   = s,
      nCount   = Matrix::colSums(sub_mat),
      nFeature = Matrix::colSums(sub_mat > 0),
      platform = platform_tag,
      Subset   = subset_tag,
      bin_size = bin_size_label,
      stringsAsFactors = FALSE
    )
  }))
}

# Compute per-barcode nCount and nFeature restricted to the intersect gene set.
# Returns a data frame with columns: barcode, nCount_intersect, nFeature_intersect.
calc_intersect_qc <- function(so, genes, assay) {
  mat <- get_counts_mat(so, assay)
  g   <- intersect(genes, rownames(mat))
  mat <- mat[g, , drop = FALSE]
  data.frame(
    barcode            = colnames(mat),
    nCount_intersect   = Matrix::colSums(mat),
    nFeature_intersect = Matrix::colSums(mat > 0)
  )
}

# ---------------------------------------------------------------------------
# Load VisiumHD samples
# ---------------------------------------------------------------------------
message("Loading VisiumHD samples...")

visiumhd_dir     <- cfg$visiumhd$data_dir
visiumhd_samples <- cfg$visiumhd$samples   # named list: sample_key -> filename

visiumhd_objs <- lapply(names(visiumhd_samples), function(samp) {
  path <- file.path(visiumhd_dir, visiumhd_samples[[samp]])
  message("  ", samp, ": ", path)
  readRDS(path)
})
names(visiumhd_objs) <- names(visiumhd_samples)

# ---------------------------------------------------------------------------
# Load MERSCOPE binning objects
# ---------------------------------------------------------------------------
message("Loading MERSCOPE binning objects...")

merscope_samples <- cfg$spatial_analysis$merscope_samples
bin_resolutions  <- cfg$bin_resolutions

merscope_objs <- setNames(
  lapply(bin_resolutions, function(res) {
    bin_dir <- file.path(cfg$output_dir, "01_preprocessing", paste0("merscope_", res, "um"))
    setNames(
      lapply(merscope_samples, function(samp) {
        path <- file.path(bin_dir, paste0(samp, "_", res, "um.rds"))
        message("  ", samp, " @ ", res, "um: ", path)
        readRDS(path)
      }),
      merscope_samples
    )
  }),
  as.character(bin_resolutions)
)

# ---------------------------------------------------------------------------
# Load Xenium binning objects
# ---------------------------------------------------------------------------
message("Loading Xenium binning objects...")

xenium_samples <- cfg$spatial_analysis$xenium_default_samples

xenium_objs <- setNames(
  lapply(bin_resolutions, function(res) {
    bin_dir <- file.path(cfg$output_dir, "01_preprocessing", paste0("xenium_", res, "um"))
    setNames(
      lapply(xenium_samples, function(samp) {
        path <- file.path(bin_dir, paste0(samp, "_", res, "um.rds"))
        message("  ", samp, " @ ", res, "um: ", path)
        readRDS(path)
      }),
      xenium_samples
    )
  }),
  as.character(bin_resolutions)
)

# ---------------------------------------------------------------------------
# Load common genes (pre-computed by dataset_summary.R)
# ---------------------------------------------------------------------------
gene_lists   <- readRDS(opt$gene_lists)
common_genes <- Reduce(intersect, list(gene_lists$VisiumHD,
                                       gene_lists$MERSCOPE,
                                       gene_lists$Xenium))
message("Three-platform common genes: ", length(common_genes))

# ---------------------------------------------------------------------------
# Build spatial platform metadata
# ---------------------------------------------------------------------------
# For each bin resolution, build metadata for All genes and the common-gene
# subset ("90") for VisiumHD, MERSCOPE, and Xenium.
# VisiumHD assay name encodes the bin size: Spatial.008um, Spatial.016um, etc.
message("Building spatial platform metadata across bin resolutions...")

meta_list <- list()

for (res in bin_resolutions) {
  bin_label  <- paste0(res, "um")
  vis_assay  <- sprintf("Spatial.%03dum", res)   # e.g. "Spatial.008um"
  mer_assay  <- "Vizgen"
  xen_assay  <- "Xenium"
  mer_objs_r <- merscope_objs[[as.character(res)]]
  xen_objs_r <- xenium_objs[[as.character(res)]]
  vis_samps  <- names(visiumhd_objs)
  mer_samps  <- names(mer_objs_r)
  xen_samps  <- names(xen_objs_r)

  message("  Building VisiumHD @ ", bin_label)
  meta_list[[paste0("vis_all_",  res)]] <- build_meta_platform(
    visiumhd_objs, vis_samps, vis_assay, "VisiumHD", "All", bin_label)
  meta_list[[paste0("vis_cg_",   res)]] <- build_meta_platform(
    visiumhd_objs, vis_samps, vis_assay, "VisiumHD", "90",  bin_label, common_genes)

  message("  Building MERSCOPE @ ", bin_label)
  meta_list[[paste0("mer_all_",  res)]] <- build_meta_platform(
    mer_objs_r, mer_samps, mer_assay, "MERSCOPE", "All", bin_label)
  meta_list[[paste0("mer_cg_",   res)]] <- build_meta_platform(
    mer_objs_r, mer_samps, mer_assay, "MERSCOPE", "90",  bin_label, common_genes)

  message("  Building Xenium @ ", bin_label)
  meta_list[[paste0("xen_all_",  res)]] <- build_meta_platform(
    xen_objs_r, xen_samps, xen_assay, "Xenium",   "All", bin_label)
  meta_list[[paste0("xen_cg_",   res)]] <- build_meta_platform(
    xen_objs_r, xen_samps, xen_assay, "Xenium",   "90",  bin_label, common_genes)
}

# Combine and apply ordered factor levels for consistent plotting
metadata_combined <- dplyr::bind_rows(meta_list) %>%
  dplyr::mutate(
    platform = factor(platform, levels = c("VisiumHD", "MERSCOPE", "Xenium")),
    Subset   = factor(Subset,   levels = c("All", "90"))
  )

message("Saving metadata_combined.rds (", nrow(metadata_combined), " rows)...")
saveRDS(metadata_combined, file.path(opt$out_dir, "metadata_combined.rds"))

# ---------------------------------------------------------------------------
# Build FLEX / scRNA-seq comparison metadata
# ---------------------------------------------------------------------------
# Intersect-gene QC across VisiumHD 8um, FLEX snRNA-seq, and 10x scRNA-seq.
# Uses raw counts ("RNA" assay) for QC; VisiumHD uses "Spatial.008um".
# NOTE: HTO_demuxmix is the hashtag demultiplexing column in the scGEM object.
#       If this column name changes, update the transmute() call below.
if (!is.null(opt$sc_rds)) {

  message("Loading FLEX snRNA-seq: ", cfg$scrna$path)
  flex_seu <- readRDS(cfg$scrna$path)

  message("Loading scRNA-seq: ", opt$sc_rds)
  sc_seu <- readRDS(opt$sc_rds)

  # Intersect gene set across the three data types
  vis_genes  <- rownames(get_counts_mat(visiumhd_objs[[1]], "Spatial.008um"))
  flex_genes <- rownames(get_counts_mat(flex_seu, "RNA"))
  sc_genes   <- rownames(get_counts_mat(sc_seu,   "RNA"))
  genes_intersect <- Reduce(intersect, list(vis_genes, flex_genes, sc_genes))
  message("Intersect genes (VisiumHD × FLEX × scRNA): ", length(genes_intersect))

  # VisiumHD: one row per bin per sample
  meta_vis_flex <- dplyr::bind_rows(
    lapply(names(visiumhd_objs), function(samp) {
      calc_intersect_qc(visiumhd_objs[[samp]], genes_intersect, "Spatial.008um") %>%
        dplyr::mutate(Sample = samp, platform = "VisiumHD")
    })
  )

  # FLEX snRNA-seq: sample label from cfg$scrna$sample_col
  flex_col  <- cfg$scrna$sample_col
  meta_flex <- calc_intersect_qc(flex_seu, genes_intersect, "RNA") %>%
    dplyr::left_join(
      flex_seu@meta.data %>%
        tibble::rownames_to_column("barcode") %>%
        dplyr::transmute(barcode, Sample = paste0("FLEX_", .data[[flex_col]])),
      by = "barcode"
    ) %>%
    dplyr::mutate(platform = "FLEX")

  # 10x scRNA-seq: sample label from HTO demultiplexing column
  meta_sc <- calc_intersect_qc(sc_seu, genes_intersect, "RNA") %>%
    dplyr::left_join(
      sc_seu@meta.data %>%
        tibble::rownames_to_column("barcode") %>%
        dplyr::transmute(barcode,
                         Sample = paste0("scRNA_", as.character(.data$HTO_demuxmix))),
      by = "barcode"
    ) %>%
    dplyr::filter(!is.na(Sample)) %>%
    dplyr::mutate(platform = "scRNAseq")

  metadata_flex_scrna <- dplyr::bind_rows(meta_vis_flex, meta_flex, meta_sc) %>%
    dplyr::mutate(
      platform = factor(platform, levels = c("scRNAseq", "FLEX", "VisiumHD"))
    )

  message("Saving metadata_flex_scrna.rds (", nrow(metadata_flex_scrna), " rows)...")
  saveRDS(metadata_flex_scrna,
          file.path(opt$out_dir, "metadata_flex_scrna.rds"))
  saveRDS(genes_intersect,
          file.path(opt$out_dir, "genes_intersect_flex.rds"))

} else {
  message("--sc_rds not provided: skipping FLEX/scRNA comparison metadata.")
  message("  Re-run with --sc_rds /path/to/scGEM_seu.rds to generate metadata_flex_scrna.rds")
}

message("Done. Outputs written to: ", opt$out_dir)
