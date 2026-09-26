# Purpose:  Shared helpers for per-bin QC metadata (nCount / nFeature),
#           sourced by qc_metrics.R (full tissue) and qc_metrics_roi.R
#           (ROI-extracted objects).
# Inputs:   None (function definitions only)
# Outputs:  None
# Author:   Ashleigh Solano
# Date:     2026-09-26

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

# Three-platform common genes from dataset_summary.R's gene_lists.rds
load_common_genes <- function(gene_lists_path) {
  gene_lists <- readRDS(gene_lists_path)
  Reduce(intersect, list(gene_lists$VisiumHD, gene_lists$MERSCOPE, gene_lists$Xenium))
}

# Apply ordered factor levels used by all Figure 2 QC plots
order_qc_factors <- function(df) {
  dplyr::mutate(
    df,
    platform = factor(platform, levels = c("VisiumHD", "MERSCOPE", "Xenium")),
    Subset   = factor(Subset,   levels = c("All", "90"))
  )
}
