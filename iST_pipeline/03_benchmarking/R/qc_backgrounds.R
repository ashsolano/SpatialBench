# Purpose:  Compute background vs target count summaries and false discovery rate
#           (FDR) for MERSCOPE and Xenium 8µm binning objects. Also loads
#           pre-computed Moran's I results and combines them into a single RDS
#           for downstream figure generation.
# Inputs:   config/config.yaml  (spatial_analysis, output_dir)
#           results/01_preprocessing/merscope_8um/{sample}_8um.rds
#           results/01_preprocessing/xenium_8um/{sample}_8um.rds
#           --moransi_mer  path to pre-computed Moran's I RDS for MERSCOPE
#           --moransi_xen  path to pre-computed Moran's I RDS for Xenium
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
merscope_samples <- cfg$spatial_analysis$merscope_samples
xenium_samples   <- cfg$spatial_analysis$xenium_default_samples

message("Loading MERSCOPE 8um objects...")
mer_dir      <- file.path(cfg$output_dir, "01_preprocessing", "merscope_8um")
merscope_8um <- setNames(
  lapply(merscope_samples, function(samp) {
    path <- file.path(mer_dir, paste0(samp, "_8um.rds"))
    message("  ", samp, ": ", path)
    readRDS(path)
  }),
  merscope_samples
)

message("Loading Xenium 8um objects...")
xen_dir     <- file.path(cfg$output_dir, "01_preprocessing", "xenium_8um")
xenium_8um  <- setNames(
  lapply(xenium_samples, function(samp) {
    path <- file.path(xen_dir, paste0(samp, "_8um.rds"))
    message("  ", samp, ": ", path)
    readRDS(path)
  }),
  xenium_samples
)

# ---------------------------------------------------------------------------
# Background count summaries
# ---------------------------------------------------------------------------
# FLAG: nCount_BlankProbe and nCount_Vizgen are MERSCOPE-specific metadata
#       column names derived from the Seurat object assay names. If assay
#       naming changes upstream, update mer_features here.
mer_features <- c("nCount_Vizgen", "nCount_BlankProbe")
xen_features <- c(
  "nCount_Xenium",
  "nCount_ControlCodeword",
  "nCount_ControlProbe",
  "nCount_Unassigned"
)

# Sum each metadata column across cells to get per-sample totals per assay type
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

# Recode raw metadata column names to human-readable assay labels
plot_df <- bind_rows(mer_df, xen_df) %>%
  mutate(
    Assay = recode(
      Assay,
      nCount_Vizgen          = "Gene",
      nCount_Xenium          = "Gene",
      nCount_BlankProbe      = "Blanks",
      nCount_ControlCodeword = "Control Codeword",
      nCount_ControlProbe    = "Control Probe",
      nCount_Unassigned      = "Unassigned"
    ),
    Assay = factor(
      Assay,
      levels = c("Gene", "Blanks", "Control Codeword", "Control Probe", "Unassigned")
    )
  )

# One summed total per sample × assay for boxplot-level data
background_per_sample <- plot_df %>%
  group_by(platform, Assay, Sample) %>%
  summarise(total_calls = sum(total_calls, na.rm = TRUE), .groups = "drop") %>%
  mutate(total_calls = round(total_calls)) %>%
  arrange(platform, Assay, Sample)

# Median/IQR summary across samples per platform × assay
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
# FLAG: n_bg and n_tg are probe-panel constants specific to these assays.
#       MERSCOPE: 89 blank probes, 91 target genes.
#       Xenium:   380 unassigned codeword probes (background), 100 target genes.
#       Consider moving these values into config.yaml if the probe panel changes
#       between experiments.
n_bg_mer <- 89
n_tg_mer <- 91
n_bg_xen <- 380
n_tg_xen <- 100

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
  merscope_8um, "nCount_BlankProbe", "nCount_Vizgen",
  n_bg_mer, n_tg_mer, "MERSCOPE"
)
xen_fdr <- compute_fdr(
  xenium_8um, "nCount_Unassigned", "nCount_Xenium",
  n_bg_xen, n_tg_xen, "Xenium"
)

fdr_results <- bind_rows(mer_fdr, xen_fdr)

message("FDR results:")
print(fdr_results)

# ---------------------------------------------------------------------------
# Moran's I (pass-through: load pre-computed RDS files and combine)
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
