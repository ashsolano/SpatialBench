# Purpose:  Extended figure for Figure 2: per-sample background panels, one
#           packed row per platform (MERSCOPE top, Xenium bottom), one facet per
#           animal, shared y limits across rows.
#           Panels: background_counts_persample (total counts per assay type);
#           moransi_persample (Moran's I, gene vs background features).
#           Animals missing from a platform must match config
#           qc_backgrounds.exclude_samples (checked).
# Supports: fig2_background.R — background signal is low and spatially
#           unstructured relative to genes, per sample (not only pooled).
# Legend notes:
#           - MERSCOPE ctrl174/ko166: no blank probes in batch 6 run (not shown).
#           - Moran's I "Background" for Xenium = unassigned codewords only; the
#             counts panel also shows Control Codeword and Control Probe.
#           - MERSCOPE Moran's I values reflect z-plane-corrected data, so they
#             are higher than in earlier drafts.
# Inputs:   results/03_benchmarking/qc_backgrounds/background_per_sample.rds
#           results/03_benchmarking/qc_backgrounds/moransi_combined.rds
#           config/config.yaml (sample lists and exclude_samples)
# Outputs:  extended_figures/fig2ext_background/background_counts_persample.pdf
#           extended_figures/fig2ext_background/moransi_persample.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(optparse)
})

source("04_manuscript/R/utils/theme.R")
source("04_manuscript/R/utils/palettes.R")
source("04_manuscript/R/utils/extended_helpers.R")

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--qc_backgrounds_dir"), type = "character",
              default = "results/03_benchmarking/qc_backgrounds",
              help    = "Directory containing qc_backgrounds.R outputs [default: %default]"),
  make_option(c("--config"), type = "character",
              default = "config/config.yaml",
              help    = "Pipeline config (sample lists, exclude_samples) [default: %default]"),
  make_option(c("--out_dir"), type = "character",
              default = "extended_figures/fig2ext_background",
              help    = "Output directory for panel PDFs [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# First level is the top row
platform_levels <- c("MERSCOPE", "Xenium")

# Animal column order (ko, wt, ctrl) and upstream exclusions from config
animal_levels    <- animal_levels_from_config(cfg)
excluded_samples <- excluded_samples_from_config(cfg)

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
bg_per_sample_path <- file.path(opt$qc_backgrounds_dir, "background_per_sample.rds")
moransi_path       <- file.path(opt$qc_backgrounds_dir, "moransi_combined.rds")
if (!file.exists(bg_per_sample_path)) stop("Not found: ", bg_per_sample_path)
if (!file.exists(moransi_path))       stop("Not found: ", moransi_path)

message("Loading background per-sample data...")
counts_df <- readRDS(bg_per_sample_path) %>%
  mutate(
    animal     = factor(sample_to_animal(Sample), levels = animal_levels),
    platform   = factor(platform, levels = platform_levels),
    # fill_group drives scale_fill_platform(): platform colour for gene counts, grey for background
    fill_group = ifelse(Assay == "Gene", as.character(platform), "Background")
  )

message("Loading Moran's I combined data...")
moransi_raw <- readRDS(moransi_path)

# Drop features whose Moran's I is undefined (reported, not silent)
na_rows <- moransi_raw %>% filter(is.na(observed))
message("Dropping ", nrow(na_rows), " features with NA Moran's I: ",
        paste0(na_rows$platform, " ", na_rows$sample, " ", na_rows$feature, collapse = "; "))

moransi_df <- moransi_raw %>%
  filter(!is.na(observed)) %>%
  mutate(
    type       = recode(type, Target = "Gene", Background = "Background"),
    type       = factor(type, levels = c("Gene", "Background")),
    animal     = factor(sample_to_animal(sample), levels = animal_levels),
    platform   = factor(platform, levels = platform_levels),
    fill_group = ifelse(type == "Gene", as.character(platform), "Background")
  )

# ---------------------------------------------------------------------------
# Missing animals must match config exclude_samples (see extended_helpers.R)
# ---------------------------------------------------------------------------
check_missing_animals(counts_df,  "Sample", "Counts",
                      animal_levels, excluded_samples, platform_levels)
check_missing_animals(moransi_df, "sample", "Moran's I",
                      animal_levels, excluded_samples, platform_levels)

# ---------------------------------------------------------------------------
# Shared y limits across both platform rows
# ---------------------------------------------------------------------------
# Counts: from 10^0 up to the next even power of ten above the largest total
counts_ymax   <- 10^(2 * ceiling(log10(max(counts_df$total_calls)) / 2))
counts_limits <- c(1, counts_ymax)
counts_breaks <- 10^seq(0, log10(counts_ymax), by = 2)

# Moran's I: span of boxplot whiskers over every sample x type box, so hidden
# outliers do not stretch the axis (applied with coord_cartesian, which zooms
# without dropping data from the box statistics)
whisker_range <- moransi_df %>%
  group_by(platform, animal, type) %>%
  summarise(
    lo = boxplot.stats(observed)$stats[1],
    hi = boxplot.stats(observed)$stats[5],
    .groups = "drop"
  )
moransi_limits <- c(min(0, whisker_range$lo), max(whisker_range$hi))

# Same fixed panel height in both PDFs so boxes match; page heights follow
# from each layout
panel_h_mm <- 20

# ---------------------------------------------------------------------------
# Panel 1: per-sample background counts
# ---------------------------------------------------------------------------
message("Building background_counts_persample panel...")

# Assay levels dropped per row so each row shows only its own assay types
plot_counts_row <- function(plat) {
  row_df <- counts_df %>% filter(platform == plat) %>% mutate(Assay = droplevels(Assay))

  ggplot(row_df, aes(x = Assay, y = total_calls, fill = fill_group)) +
    geom_col(width = 0.7) +
    scale_fill_platform() +
    facet_animals() +
    coord_cartesian(clip = "off") +
    scale_y_log10(
      "Total counts/sample",
      limits = counts_limits,
      breaks = counts_breaks,
      labels = trans_format("log10", math_format(10^.x)),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(x = NULL, title = plat) +
    theme_sb() +
    theme_bg_panels +
    theme_ext
}

p_counts <- wrap_plots(map(platform_levels, plot_counts_row), ncol = 1)

save_fixed_panels(p_counts, file.path(opt$out_dir, "background_counts_persample.pdf"),
                  panel_h_mm = panel_h_mm, n_rows = length(platform_levels))

# ---------------------------------------------------------------------------
# Panel 2: per-sample Moran's I
# ---------------------------------------------------------------------------
message("Building moransi_persample panel...")

plot_moransi_row <- function(plat) {
  row_df <- moransi_df %>% filter(platform == plat)

  ggplot(row_df, aes(x = type, y = observed, fill = fill_group)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
    # thinner outline than the default so small facets keep a visible fill
    geom_boxplot(outlier.shape = NA, width = 0.7, linewidth = 0.3) +
    scale_fill_platform() +
    facet_animals() +
    coord_cartesian(ylim = moransi_limits, clip = "off") +
    labs(x = NULL, y = "Observed Moran's I", title = plat) +
    theme_sb() +
    theme_bg_panels +
    theme_ext
}

p_moransi <- wrap_plots(map(platform_levels, plot_moransi_row), ncol = 1)

save_fixed_panels(p_moransi, file.path(opt$out_dir, "moransi_persample.pdf"),
                  panel_h_mm = panel_h_mm, n_rows = length(platform_levels))

message("Done. All panels written to: ", opt$out_dir)
