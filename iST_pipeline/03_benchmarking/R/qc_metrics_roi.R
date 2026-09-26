# Purpose:  Compute per-bin QC metadata (nCount, nFeature) for the fixed-size
#           ROI-extracted 8µm objects of the matched animals (VisiumHD,
#           MERSCOPE, Xenium), for both the full platform gene set and the
#           three-platform common gene subset. ROI-restricted counterpart of
#           qc_metrics.R, used for the Figure 2 ROI QC boxplots.
# Inputs:   config/config.yaml  (stalign$visiumhd, stalign$matched_samples)
#           results/01_preprocessing/visium_8um_roi/{visium_sample}_8um_roi_{roi_label}.rds
#           results/01_preprocessing/merscope_8um_roi/{sample}_8um_roi_{roi_label}.rds
#           results/01_preprocessing/xenium_8um_roi/{sample}_8um_roi_{roi_label}.rds
#           results/03_benchmarking/dataset_summary/gene_lists.rds  (common genes)
# Outputs:  results/03_benchmarking/qc_metrics/metadata_roi{roi_label}.rds
#               (per-bin rows: Sample, nCount, nFeature, platform, Subset,
#                bin_size, animal — same columns as metadata_combined.rds
#                plus the matched animal ID)
# Author:   Ashleigh Solano
# Date:     2026-09-26

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(yaml)
  library(optparse)
})

source("03_benchmarking/R/utils/qc_utils.R")  # must be run from the project root

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
  make_option(c("--roi_label"),  type = "character", default = "2mm",
              help    = "ROI size label used in the ROI filenames, e.g. '2mm' [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opt$gene_lists)) stop("--gene_lists not found: ", opt$gene_lists)

cfg <- yaml::read_yaml(opt$config)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Matched animals and their per-platform sample keys
# ---------------------------------------------------------------------------
# Same animal -> sample mapping used by align_binned.R and extract_roi.R
animals <- names(cfg$stalign$visiumhd)

platforms <- list(
  VisiumHD = list(keys  = unlist(cfg$stalign$visiumhd[animals]),
                  dir   = "visium_8um_roi",   assay = "Spatial.008um"),
  MERSCOPE = list(keys  = unlist(cfg$stalign$matched_samples$merscope[animals]),
                  dir   = "merscope_8um_roi", assay = "Vizgen"),
  Xenium   = list(keys  = unlist(cfg$stalign$matched_samples$xenium[animals]),
                  dir   = "xenium_8um_roi",   assay = "Xenium")
)

# ---------------------------------------------------------------------------
# Load common genes (pre-computed by dataset_summary.R)
# ---------------------------------------------------------------------------
common_genes <- load_common_genes(opt$gene_lists)
message("Three-platform common genes: ", length(common_genes))

# ---------------------------------------------------------------------------
# Build ROI metadata per platform (All genes + common-gene subset "90")
# ---------------------------------------------------------------------------
meta_list <- list()

for (plat in names(platforms)) {
  p <- platforms[[plat]]
  if (length(p$keys) != length(animals)) {
    stop("Missing config$stalign sample key(s) for ", plat)
  }

  # Load this platform's ROI objects, one per matched animal
  message("Loading ", plat, " ROI objects...")
  objs <- setNames(lapply(p$keys, function(key) {
    path <- file.path(cfg$output_dir, "01_preprocessing", p$dir,
                      paste0(key, "_8um_roi_", opt$roi_label, ".rds"))
    message("  ", key, ": ", path)
    readRDS(path)
  }), p$keys)

  meta_list[[paste0(plat, "_all")]] <- build_meta_platform(
    objs, p$keys, p$assay, plat, "All", "8um")
  meta_list[[paste0(plat, "_cg")]]  <- build_meta_platform(
    objs, p$keys, p$assay, plat, "90",  "8um", common_genes)

  rm(objs); gc()
}

# ---------------------------------------------------------------------------
# Combine, tag the matched animal, and save
# ---------------------------------------------------------------------------
# Lookup: platform sample key (e.g. "wt709_batch13") -> animal ID ("wt709")
key_to_animal <- setNames(rep(animals, length(platforms)),
                          unlist(lapply(platforms, function(p) unname(p$keys))))

metadata_roi <- order_qc_factors(dplyr::bind_rows(meta_list)) %>%
  dplyr::mutate(animal = unname(key_to_animal[Sample]))

out_path <- file.path(opt$out_dir, paste0("metadata_roi", opt$roi_label, ".rds"))
message("Saving ", basename(out_path), " (", nrow(metadata_roi), " rows)...")
saveRDS(metadata_roi, out_path)

message("Done. Output written to: ", opt$out_dir)
