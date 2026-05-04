# Purpose:  Build pseudobulk count matrices across VisiumHD, MERSCOPE, and Xenium
#           using 8µm bins restricted to the common gene set. Finds genes present in
#           every sample of every platform, sums counts per sample, and exports a
#           DGEList (edgeR) with sample metadata for downstream gene-comparison panels.
# Inputs:   config/config.yaml
#           results/01_preprocessing/merscope_8um/{sample}_8um.rds
#           results/01_preprocessing/xenium_8um/{sample}_8um.rds
#           cfg$visiumhd$data_dir / cfg$visiumhd$samples  (VisiumHD Seurat objects)
# Outputs:  results/03_benchmarking/gene_comparison/dge.rds
#           results/03_benchmarking/gene_comparison/counts_mat.rds


suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
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
              default = "results/03_benchmarking/gene_comparison",
              help    = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Retrieve the counts layer/slot, supporting both Seurat v4 and v5 APIs.
get_counts_compat <- function(sobj, assay) {
  tryCatch(
    GetAssayData(sobj, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(sobj, assay = assay, slot = "counts")
  )
}

# Sum counts across all non-empty bins for a single sample, restricted to
# the supplied gene set. Returns a one-row-per-gene tibble suitable for
# bind_rows() into the full pseudobulk data frame.
get_pseudobulk_common <- function(sobj, sample_name, platform, genes) {
  assay <- DefaultAssay(sobj)
  mat   <- get_counts_compat(sobj, assay)

  # Drop bins with zero total counts — they contribute nothing to the sum
  keep <- Matrix::colSums(mat) > 0
  if (any(!keep))
    message("    ", sample_name, ": removing ", sum(!keep), " empty bins")

  mat   <- mat[genes, keep, drop = FALSE]
  cnts  <- Matrix::rowSums(mat)

  tibble(feature  = names(cnts),
         count    = as.numeric(cnts),
         Platform = platform,
         SampleID = sample_name)
}

# ---------------------------------------------------------------------------
# Load VisiumHD samples
# ---------------------------------------------------------------------------
message("Loading VisiumHD samples...")

visiumhd_dir     <- cfg$visiumhd$data_dir
visiumhd_samples <- cfg$visiumhd$samples

visiumhd_8um <- lapply(names(visiumhd_samples), function(samp) {
  path <- file.path(visiumhd_dir, visiumhd_samples[[samp]])
  message("  ", samp, ": ", path)
  obj <- readRDS(path)
  DefaultAssay(obj) <- "Spatial.008um"
  obj
}) |> setNames(names(visiumhd_samples))

# ---------------------------------------------------------------------------
# Load MERSCOPE 8µm binning objects
# ---------------------------------------------------------------------------
message("Loading MERSCOPE 8µm binning objects...")

merscope_samples <- cfg$spatial_analysis$merscope_samples
bin_res          <- cfg$bin_resolutions[[1]]   # 8

merscope_8um <- lapply(merscope_samples, function(samp) {
  path <- file.path(cfg$output_dir, "01_preprocessing",
                    paste0("merscope_", bin_res, "um"),
                    paste0(samp, "_", bin_res, "um.rds"))
  message("  ", samp, ": ", path)
  obj <- readRDS(path)
  DefaultAssay(obj) <- "Vizgen"
  obj
}) |> setNames(merscope_samples)

# ---------------------------------------------------------------------------
# Load Xenium 8µm binning objects
# ---------------------------------------------------------------------------
message("Loading Xenium 8µm binning objects...")

xenium_samples <- cfg$spatial_analysis$xenium_default_samples

xenium_8um <- lapply(xenium_samples, function(samp) {
  path <- file.path(cfg$output_dir, "01_preprocessing",
                    paste0("xenium_", bin_res, "um"),
                    paste0(samp, "_", bin_res, "um.rds"))
  message("  ", samp, ": ", path)
  obj <- readRDS(path)
  DefaultAssay(obj) <- "Xenium"
  obj
}) |> setNames(xenium_samples)

# ---------------------------------------------------------------------------
# Identify common genes across all samples and all platforms
# ---------------------------------------------------------------------------
message("Finding common genes across all platforms...")

# Within-platform intersection first, then across platforms.
# This ensures every sample in every platform carries the common gene set.
vis_genes <- Reduce(intersect,
                    lapply(visiumhd_8um,  function(s) rownames(get_counts_compat(s, DefaultAssay(s)))))
mer_genes <- Reduce(intersect,
                    lapply(merscope_8um, function(s) rownames(get_counts_compat(s, DefaultAssay(s)))))
xen_genes <- Reduce(intersect,
                    lapply(xenium_8um,  function(s) rownames(get_counts_compat(s, DefaultAssay(s)))))

common_genes <- Reduce(intersect, list(vis_genes, mer_genes, xen_genes))
message("  VisiumHD genes: ", length(vis_genes))
message("  MERSCOPE genes: ", length(mer_genes))
message("  Xenium genes:   ", length(xen_genes))
message("  Common genes:   ", length(common_genes))

# ---------------------------------------------------------------------------
# Build pseudobulk data frame
# ---------------------------------------------------------------------------
message("Building pseudobulk count matrix...")

pb_df <- bind_rows(
  map_dfr(names(visiumhd_8um),  ~ get_pseudobulk_common(visiumhd_8um[[.]],  ., "VisiumHD", common_genes)),
  map_dfr(names(merscope_8um),  ~ get_pseudobulk_common(merscope_8um[[.]],  ., "MERSCOPE", common_genes)),
  map_dfr(names(xenium_8um),    ~ get_pseudobulk_common(xenium_8um[[.]],    ., "Xenium",   common_genes))
)

# Pivot to genes × samples count matrix, filling any gaps with 0
counts_mat <- pb_df %>%
  pivot_wider(names_from  = c(Platform, SampleID),
              names_sep   = "_",
              values_from = count,
              values_fill = list(count = 0)) %>%
  column_to_rownames("feature") %>%
  as.matrix()

# ---------------------------------------------------------------------------
# Build sample metadata (col_info) and derive condition + animal ID
# ---------------------------------------------------------------------------
# Column names of counts_mat are "{Platform}_{SampleID}".
# Type is derived from the sample name prefix (MERSCOPE/Xenium: wt/ko/ctrl)
# or from the numeric animal ID suffix (VisiumHD: batch33_NNN).
# FLAG: animal ID sets below are project-specific; update if new samples are added.

col_info <- tibble(Combined = colnames(counts_mat)) %>%
  separate(col    = Combined,
           into   = c("Platform", "SampleID"),
           sep    = "_",
           extra  = "merge",
           remove = FALSE) %>%
  mutate(
    IDnum = case_when(
      # VisiumHD: SampleID like "batch33_167" — take the trailing number
      Platform == "VisiumHD" ~ sub(".*_(\\d+)$", "\\1", SampleID),
      # MERSCOPE/Xenium: SampleID like "wt709_batch13" — take the leading number
      TRUE                   ~ sub("^[a-z]+(\\d+)_.*", "\\1", SampleID)
    ),
    Type = case_when(
      # MERSCOPE/Xenium: uppercase the alpha prefix (wt -> WT, ko -> KO, ctrl -> CTRL)
      Platform != "VisiumHD" ~ toupper(sub("^([a-z]+)\\d+.*", "\\1", SampleID)),
      # VisiumHD: infer from animal ID
      IDnum %in% c("709", "710", "713") ~ "WT",
      IDnum %in% c("166", "167", "168") ~ "KO",
      IDnum %in% c("172", "173", "174") ~ "CTRL",
      TRUE ~ "Unknown"
    ),
    Type = factor(Type, levels = c("WT", "KO", "CTRL"))
  ) %>%
  column_to_rownames("Combined") %>%
  select(SampleID, Platform, Type, IDnum)

# ---------------------------------------------------------------------------
# Create DGEList
# ---------------------------------------------------------------------------
dge <- DGEList(counts = counts_mat, samples = col_info)
message("DGEList: ", nrow(dge), " genes x ", ncol(dge), " samples")

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
dge_path  <- file.path(opt$out_dir, "dge.rds")
mat_path  <- file.path(opt$out_dir, "counts_mat.rds")

saveRDS(dge,        dge_path)
saveRDS(counts_mat, mat_path)

message("Saved: ", dge_path)
message("Saved: ", mat_path)
message("Done.")
