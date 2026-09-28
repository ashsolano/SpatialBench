# Purpose:  Background vs target count summaries and false discovery rate (FDR)
#           for MERSCOPE and Xenium post-QC 8µm bins (empty-bin + DBSCAN filtered
#           by filter_binned_{merscope,xenium}); also combines the pre-computed
#           Moran's I results into one RDS.
# Inputs:   config/config.yaml  (spatial_analysis, qc_backgrounds, output_dir)
#           results/01_preprocessing/merscope_8um_filtered/{sample}_8um_filtered.rds
#           results/01_preprocessing/xenium_8um_filtered/{sample}_8um_filtered.rds
#           --moransi_mer  path to Moran's I RDS for MERSCOPE (from moransi.R)
#           --moransi_xen  path to Moran's I RDS for Xenium   (from moransi.R)
# Outputs:  results/03_benchmarking/qc_backgrounds/background_per_sample.rds
#               (per-sample × assay total counts: platform, Sample, Assay, total_calls)
#           results/03_benchmarking/qc_backgrounds/background_summary.rds
#               (median and IQR per platform × assay)
#           results/03_benchmarking/qc_backgrounds/fdr_results.rds
#               (per-sample FDR: Sample, Platform, FDR)
#           results/03_benchmarking/qc_backgrounds/moransi_combined.rds
#               (combined Moran's I with platform column; only if both --moransi_* provided)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(yaml)
  library(optparse)
})

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--config"), type = "character",
              default = "config/config.yaml",
              help    = "Path to config.yaml [default: %default]"),
  make_option(c("--out_dir"), type = "character",
              default = "results/03_benchmarking/qc_backgrounds",
              help    = "Output directory [default: %default]"),
  make_option(c("--moransi_mer"), type = "character",
              default = NULL,
              help    = "Path to pre-computed Moran's I RDS for MERSCOPE (optional)"),
  make_option(c("--moransi_xen"), type = "character",
              default = NULL,
              help    = "Path to pre-computed Moran's I RDS for Xenium (optional)")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load 8um binning objects
# ---------------------------------------------------------------------------
# Samples without background probe signal are excluded via
# config["qc_backgrounds"]["exclude_samples"] (see config.yaml for reasons)
bg_cfg           <- cfg$qc_backgrounds
merscope_samples <- setdiff(cfg$spatial_analysis$merscope_samples,      unlist(bg_cfg$exclude_samples$merscope))
xenium_samples   <- setdiff(cfg$spatial_analysis$xenium_default_samples, unlist(bg_cfg$exclude_samples$xenium))
walk(c("merscope", "xenium"), function(p) {
  excl <- unlist(bg_cfg$exclude_samples[[p]])
  if (length(excl) > 0) message("Excluding ", p, " samples (no background signal): ", paste(excl, collapse = ", "))
})

# Read the filtered 8um object for each sample of one platform; stops if any
# sample is missing so that all samples are guaranteed to be retained
load_filtered_8um <- function(platform, samples) {
  in_dir <- file.path(cfg$output_dir, "01_preprocessing", paste0(platform, "_8um_filtered"))
  paths  <- file.path(in_dir, paste0(samples, "_8um_filtered.rds"))
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) stop("Missing filtered 8um objects:\n  ", paste(missing, collapse = "\n  "))
  walk2(samples, paths, ~ message("  ", .x, ": ", .y))
  setNames(map(paths, readRDS), samples)
}

# Stop if any retained sample lacks a required metadata column (e.g. a sample
# with no background assay that has not been listed in exclude_samples)
check_features <- function(obj_list, features) {
  iwalk(obj_list, function(obj, samp) {
    absent <- setdiff(features, colnames(obj@meta.data))
    if (length(absent) > 0) {
      stop(samp, " is missing: ", paste(absent, collapse = ", "),
           ". Add it to qc_backgrounds$exclude_samples in config.yaml if intended.")
    }
  })
}

message("Loading MERSCOPE filtered 8um objects...")
merscope_8um <- load_filtered_8um("merscope", merscope_samples)

message("Loading Xenium filtered 8um objects...")
xenium_8um   <- load_filtered_8um("xenium", xenium_samples)

message("Samples loaded: MERSCOPE = ", length(merscope_8um), ", Xenium = ", length(xenium_8um))

# ---------------------------------------------------------------------------
# Background count summaries
# ---------------------------------------------------------------------------
# Metadata column names are derived from the Seurat assay names in the binned
# objects. MERSCOPE blank probes are in the "Blanks" assay. Xenium
# "Unassigned Codeword" features are loaded by Seurat into the "BlankCodeword"
# assay. If assay naming changes upstream, update these vectors.
mer_features <- c("nCount_Vizgen", "nCount_Blanks")
xen_features <- c(
  "nCount_Xenium",
  "nCount_ControlCodeword",
  "nCount_ControlProbe",
  "nCount_BlankCodeword"
)
check_features(merscope_8um, mer_features)
check_features(xenium_8um,   xen_features)

summarise_counts <- function(obj_list, features, platform_name) {
  tibble(
    Sample   = names(obj_list),
    platform = platform_name,
    Counts   = map(obj_list, ~ colSums(.x@meta.data[, features, drop = FALSE], na.rm = TRUE))
  ) %>%
    unnest_wider(Counts) %>%
    pivot_longer(
      cols      = all_of(features),
      names_to  = "Assay",
      values_to = "total_calls"
    )
}

mer_df <- summarise_counts(merscope_8um, mer_features, "MERSCOPE")
xen_df <- summarise_counts(xenium_8um,   xen_features, "Xenium")

plot_df <- bind_rows(mer_df, xen_df) %>%
  mutate(
    Assay = recode(
      Assay,
      nCount_Vizgen          = "Gene",
      nCount_Xenium          = "Gene",
      nCount_Blanks          = "Blanks",
      nCount_ControlCodeword = "Control Codeword",
      nCount_ControlProbe    = "Control Probe",
      nCount_BlankCodeword   = "Unassigned"
    ),
    Assay = factor(
      Assay,
      levels = c("Gene", "Blanks", "Control Codeword", "Control Probe", "Unassigned")
    )
  )

background_per_sample <- plot_df %>%
  group_by(platform, Assay, Sample) %>%
  summarise(total_calls = sum(total_calls, na.rm = TRUE), .groups = "drop") %>%
  mutate(total_calls = round(total_calls)) %>%
  arrange(platform, Assay, Sample)

background_summary <- background_per_sample %>%
  group_by(platform, Assay) %>%
  summarise(
    n_samples    = dplyr::n(),
    median_total = median(total_calls, na.rm = TRUE),
    IQR_total    = IQR(total_calls, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    median_total = round(median_total),
    IQR_total    = round(IQR_total)
  ) %>%
  arrange(platform, Assay)

message("Background summary:")
print(background_summary)

# ---------------------------------------------------------------------------
# FDR computation
# ---------------------------------------------------------------------------
# Probe-panel constants from config["qc_backgrounds"]:
#   MERSCOPE: 89 blank probes, 91 target genes.
#   Xenium:   380 unassigned codeword probes (background), 100 target genes.
n_bg_mer <- bg_cfg$n_bg$merscope
n_tg_mer <- bg_cfg$n_tg$merscope
n_bg_xen <- bg_cfg$n_bg$xenium
n_tg_xen <- bg_cfg$n_tg$xenium

# FDR = (background calls / n background probes) / (target calls / n target genes) * 100
compute_fdr <- function(obj_list, bg_feat, tg_feat, n_bg, n_tg, platform_name) {
  tibble(
    Sample   = names(obj_list),
    Platform = platform_name,
    bg_calls = map_dbl(obj_list, ~ sum(.x@meta.data[, bg_feat], na.rm = TRUE)),
    tg_calls = map_dbl(obj_list, ~ sum(.x@meta.data[, tg_feat], na.rm = TRUE))
  ) %>%
    mutate(FDR = (bg_calls / n_bg) * (n_tg / tg_calls) * 100) %>%
    select(Sample, Platform, FDR)
}

mer_fdr <- compute_fdr(
  merscope_8um, "nCount_Blanks", "nCount_Vizgen",
  n_bg_mer, n_tg_mer, "MERSCOPE"
)
xen_fdr <- compute_fdr(
  xenium_8um, "nCount_BlankCodeword", "nCount_Xenium",
  n_bg_xen, n_tg_xen, "Xenium"
)

fdr_results <- bind_rows(mer_fdr, xen_fdr)

message("FDR results:")
print(fdr_results, n = Inf)

# ---------------------------------------------------------------------------
# Moran's I (pass-through: load moransi.R outputs and combine)
# ---------------------------------------------------------------------------
if (!is.null(opt$moransi_mer) && !is.null(opt$moransi_xen)) {
  if (!file.exists(opt$moransi_mer)) stop("--moransi_mer not found: ", opt$moransi_mer)
  if (!file.exists(opt$moransi_xen)) stop("--moransi_xen not found: ", opt$moransi_xen)

  message("Loading Moran's I results...")
  mrs <- readRDS(opt$moransi_mer) %>% mutate(platform = "MERSCOPE")
  xen <- readRDS(opt$moransi_xen) %>% mutate(platform = "Xenium")
  moransi_combined <- bind_rows(mrs, xen)

  saveRDS(moransi_combined, file.path(opt$out_dir, "moransi_combined.rds"))
  message("Saved: moransi_combined.rds")
} else {
  message("--moransi_mer / --moransi_xen not provided: skipping Moran's I.")
  message("  Re-run with both flags to generate moransi_combined.rds.")
}

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
saveRDS(background_per_sample, file.path(opt$out_dir, "background_per_sample.rds"))
message("Saved: background_per_sample.rds")

saveRDS(background_summary, file.path(opt$out_dir, "background_summary.rds"))
message("Saved: background_summary.rds")

saveRDS(fdr_results, file.path(opt$out_dir, "fdr_results.rds"))
message("Saved: fdr_results.rds")

message("Done. Outputs written to: ", opt$out_dir)
