# Purpose:  Figure 2 — QC background panels.
#           Panel 1 (background_vs_target): total counts by assay type (gene vs
#                    background signals) per platform, boxplot + per-sample jitter.
#           Panel 2 (moransi): Moran's I by assay type, faceted by platform.
#                    Only produced if moransi_combined.rds is present.
#           Panel 3 (fdr): false discovery rate by platform, bar chart with
#                    raw sample points.
#           Panel 4 (probe_scurves): probe rank S-curves coloured by mean count,
#                    with overlap probes highlighted and labelled.
#           Each panel saved as a separate PDF to figures/fig2/.
#           Adapted from fig2_qcbackgrounds_v2.R, fig2_moransI_fdr.R, and
#           probe_rank_platform_v2.R.
# Inputs:   results/03_benchmarking/qc_backgrounds/background_per_sample.rds
#           results/03_benchmarking/qc_backgrounds/fdr_results.rds
#           results/03_benchmarking/qc_backgrounds/moransi_combined.rds  (optional)
#           results/03_benchmarking/probe_rank/ranked_plat.rds
#           results/03_benchmarking/probe_rank/label_pool_genes.csv
# Outputs:  figures/fig2/background_vs_target.pdf
#           figures/fig2/moransi.pdf  (only if moransi_combined.rds present)
#           figures/fig2/fdr.pdf
#           figures/fig2/probe_scurves.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(optparse)
})

source("04_manuscript/R/utils/theme.R")
source("04_manuscript/R/utils/palettes.R")

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--qc_backgrounds_dir"), type = "character",
              default = "results/03_benchmarking/qc_backgrounds",
              help    = "Directory containing qc_backgrounds.R outputs [default: %default]"),
  make_option(c("--probe_rank_dir"), type = "character",
              default = "results/03_benchmarking/probe_rank",
              help    = "Directory containing probe_rank.R outputs [default: %default]"),
  make_option(c("--out_dir"), type = "character",
              default = "figures/fig2",
              help    = "Output directory for panel PDFs [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

bg_dir   <- opt$qc_backgrounds_dir
rank_dir <- opt$probe_rank_dir
out_dir  <- opt$out_dir

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
bg_per_sample_path <- file.path(bg_dir, "background_per_sample.rds")
fdr_path           <- file.path(bg_dir, "fdr_results.rds")
moransi_path       <- file.path(bg_dir, "moransi_combined.rds")
ranked_plat_path   <- file.path(rank_dir, "ranked_plat.rds")
label_pool_path    <- file.path(rank_dir, "label_pool_genes.csv")

if (!file.exists(bg_per_sample_path)) stop("Not found: ", bg_per_sample_path)
if (!file.exists(fdr_path))           stop("Not found: ", fdr_path)
if (!file.exists(ranked_plat_path))   stop("Not found: ", ranked_plat_path)
if (!file.exists(label_pool_path))    stop("Not found: ", label_pool_path)

message("Loading background per-sample data...")
background_per_sample <- readRDS(bg_per_sample_path) %>%
  # fill_group drives scale_fill_platform(): platform colour for gene counts, grey for background
  mutate(fill_group = ifelse(Assay == "Gene", platform, "Background"))

message("Loading FDR results...")
fdr_results <- readRDS(fdr_path) %>%
  mutate(platform = factor(Platform, levels = c("MERSCOPE", "Xenium")))

message("Loading probe rank table...")
ranked_plat <- readRDS(ranked_plat_path)
label_pool  <- utils::read.csv(label_pool_path)$feature

# ---------------------------------------------------------------------------
# Shared theme additions for all panels in this script
# ---------------------------------------------------------------------------
theme_bg_panels <- theme(
  axis.text.x      = element_text(angle = 45, hjust = 1),
  legend.position  = "none",
  axis.line        = element_line(colour = "black"),
  strip.background = element_rect(colour = "black", fill = "white"),
  strip.text       = element_text(face = "italic"),
  panel.background = element_blank(),
  panel.border     = element_blank(),
  plot.background  = element_blank(),
  panel.spacing    = unit(2, "mm")
)

# ---------------------------------------------------------------------------
# Panel 1: Background vs target counts
# ---------------------------------------------------------------------------
message("Building background_vs_target panel...")

total_df <- background_per_sample %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

p_background <- ggplot(total_df, aes(x = Assay, y = total_calls, fill = fill_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_point(
    aes(group = Sample),
    position = position_jitter(width = 0.10, height = 0),
    shape    = 21,
    colour   = "grey20",
    stroke   = 0.35,
    size     = 2.0,
    alpha    = 0.80
  ) +
  scale_fill_platform() +
  scale_y_log10(
    "Total counts/sample",
    labels = trans_format("log10", math_format(10^.x)),
    breaks = trans_breaks("log10", function(x) 10^x)
  ) +
  facet_wrap(~ platform, scales = "free_x", nrow = 1) +
  labs(x = "Assay") +
  theme_sb() +
  theme_bg_panels

ggsave(
  file.path(out_dir, "background_vs_target.pdf"),
  p_background,
  width  = dims$full_w,
  height = dims$half_w * 1.4,
  units  = "mm",
  device = cairo_pdf,
  bg     = "white"
)
message("Saved: background_vs_target.pdf")

# ---------------------------------------------------------------------------
# Panel 2: Moran's I
# ---------------------------------------------------------------------------
if (file.exists(moransi_path)) {
  message("Loading Moran's I combined data...")
  moransi_combined <- readRDS(moransi_path)

  df_moransi <- moransi_combined %>%
    mutate(
      type       = recode(type, Target = "Gene", Background = "Background"),
      type       = factor(type, levels = c("Gene", "Background")),
      fill_group = ifelse(type == "Gene", platform, "Background")
    )

  p_moransi <- ggplot(df_moransi, aes(x = type, y = observed, fill = fill_group)) +
    geom_boxplot(
      outlier.shape = NA,
      position      = position_dodge(width = 0.8),
      width         = 0.7
    ) +
    facet_wrap(~ platform, scales = "free_x", nrow = 1) +
    scale_fill_platform() +
    labs(x = "Assay", y = "Observed Moran's I") +
    theme_sb() +
    theme_bg_panels

  ggsave(
    file.path(out_dir, "moransi.pdf"),
    p_moransi,
    width  = dims$full_w,
    height = dims$half_w * 1.4,
    units  = "mm",
    device = cairo_pdf,
    bg     = "white"
  )
  message("Saved: moransi.pdf")
} else {
  message("Skipping Moran's I panel: ", moransi_path, " not found.")
  message("  Re-run qc_backgrounds.R with --moransi_mer and --moransi_xen to generate it.")
}

# ---------------------------------------------------------------------------
# Panel 3: FDR bar chart
# ---------------------------------------------------------------------------
message("Building FDR panel...")

fdr_stats <- fdr_results %>%
  group_by(platform) %>%
  summarise(
    meanFDR = mean(FDR),
    sdFDR   = sd(FDR),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("%.2f±%.2f%%", meanFDR, sdFDR))

p_fdr <- ggplot(fdr_stats, aes(x = platform, y = meanFDR, fill = platform)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = meanFDR - sdFDR, ymax = meanFDR + sdFDR),
    width = 0.2
  ) +
  geom_jitter(
    data   = fdr_results,
    aes(x = platform, y = FDR),
    width  = 0.15,
    shape  = 21,
    fill   = "white",
    colour = "black",
    size   = 1,
    alpha  = 0.8
  ) +
  geom_text(
    data  = fdr_stats,
    aes(label = label, y = meanFDR + sdFDR + 1),
    vjust = 0,
    size  = 3
  ) +
  scale_fill_platform() +
  labs(x = "Platform", y = "FDR (%)") +
  theme_sb() +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.line        = element_line(colour = "black"),
    panel.background = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA),
    plot.background  = element_blank()
  )

ggsave(
  file.path(out_dir, "fdr.pdf"),
  p_fdr,
  width  = dims$half_w * 0.6,
  height = dims$half_w * 1.2,
  units  = "mm",
  device = cairo_pdf,
  bg     = "white"
)
message("Saved: fdr.pdf")

# ---------------------------------------------------------------------------
# Panel 4: Probe S-curves
# ---------------------------------------------------------------------------
message("Building probe S-curves panel...")

# Label points for the subset of probes in the pre-computed label pool
labels_df <- ranked_plat %>%
  filter(Type == "Target", feature %in% label_pool)

# Shared y-limits across both platforms with small padding
limits_y <- range(ranked_plat$y, finite = TRUE)
pad      <- diff(limits_y) * 0.04
limits_y <- c(max(0, limits_y[1] - pad), limits_y[2] + pad)

# Shaded background region covering mean counts up to the bg95 threshold
thr_df   <- ranked_plat %>% distinct(platform, bg95) %>% mutate(yline = log10(bg95 + 1))
shade_df <- thr_df %>% transmute(platform, xmin = 0, xmax = 1, ymin = -Inf, ymax = yline)

# Burgundy used for overlap targets that fall within the background threshold
overlap_color <- "#800020"
point_size    <- 1.9

p_probe <- ggplot() +
  geom_rect(
    data        = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill        = "grey90",
    alpha       = 0.5
  ) +
  geom_point(
    data  = ranked_plat %>% filter(Type == "Blank"),
    aes(x = rank_frac, y = y),
    shape = 17, colour = "grey50", size = point_size, alpha = 0.7
  ) +
  geom_point(
    data  = ranked_plat %>% filter(Type == "Target", !is_overlap),
    aes(x = rank_frac, y = y, colour = y),
    shape = 16, size = point_size, alpha = 0.85
  ) +
  scale_colour_viridis_c(
    name   = "log10(mean count + 1)",
    option = "D",
    guide  = guide_colourbar(
      title.position = "top",
      barwidth       = unit(3, "mm"),
      barheight      = unit(20, "mm")
    )
  ) +
  geom_point(
    data   = ranked_plat %>% filter(Type == "Target", is_overlap),
    aes(x = rank_frac, y = y),
    shape  = 21,
    fill   = overlap_color,
    colour = "white",
    stroke = 0.25,
    size   = point_size + 0.2,
    alpha  = 0.95
  ) +
  ggrepel::geom_text_repel(
    data               = labels_df,
    aes(x = rank_frac, y = y, label = feature),
    colour             = "black",
    size               = 2.2,
    segment.colour     = "grey50",
    segment.curvature  = -0.2,
    segment.ncp        = 5,
    segment.angle      = 90,
    box.padding        = 0.25,
    point.padding      = 0.12,
    min.segment.length = 0,
    force              = 1.2,
    max.overlaps       = Inf,
    seed               = 123
  ) +
  facet_wrap(~ platform, nrow = 1, scales = "free_x", drop = TRUE) +
  scale_x_continuous(
    "Probe rank (percentile)",
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = limits_y, clip = "off") +
  labs(y = expression(log[10] ~ "(mean count + 1)")) +
  theme_sb() +
  theme(
    legend.position  = "right",
    legend.direction = "vertical",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.line        = element_line(colour = "black"),
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text.x     = element_text(size = 8, face = "italic"),
    panel.background = element_blank(),
    panel.border     = element_blank(),
    plot.background  = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(2, "mm"),
    plot.margin      = margin(5.5, 18, 5.5, 5.5, "pt")
  )

ggsave(
  file.path(out_dir, "probe_scurves.pdf"),
  p_probe,
  width  = dims$full_w,
  height = dims$half_w,
  units  = "mm",
  device = cairo_pdf,
  bg     = "white"
)
message("Saved: probe_scurves.pdf")

message("Done. All panels written to: ", out_dir)
