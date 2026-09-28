# Purpose:  Fig 2 extended: per-sample probe rank vs background, and Moran's I of the low-rank genes.
# Supports: fig2_background.R — low-rank MERSCOPE genes sit in background per sample; structure differs by platform.
# Inputs:   results/03_benchmarking/probe_rank/{ranked_sample.rds,label_pool_genes.csv}, results/03_benchmarking/qc_backgrounds/moransi_combined.rds, config/config.yaml
# Outputs:  extended_figures/fig2ext_probe_rank/{probe_rank_persample,moransi_low_rank_genes}.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(optparse)
})

source("04_manuscript/R/utils/theme.R")
source("04_manuscript/R/utils/palettes.R")
source("04_manuscript/R/utils/extended_helpers.R")

set.seed(42)

# ---- CLI arguments ----
option_list <- list(
  make_option(c("--probe_rank_dir"), type = "character",
              default = "results/03_benchmarking/probe_rank",
              help    = "Directory containing probe_rank.R outputs [default: %default]"),
  make_option(c("--qc_backgrounds_dir"), type = "character",
              default = "results/03_benchmarking/qc_backgrounds",
              help    = "Directory containing qc_backgrounds.R outputs [default: %default]"),
  make_option(c("--config"), type = "character",
              default = "config/config.yaml",
              help    = "Pipeline config (sample lists, exclude_samples) [default: %default]"),
  make_option(c("--out_dir"), type = "character",
              default = "extended_figures/fig2ext_probe_rank",
              help    = "Output directory for panel PDFs [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

platform_levels <- c("MERSCOPE", "Xenium")

animal_levels    <- animal_levels_from_config(cfg)
excluded_samples <- excluded_samples_from_config(cfg)

# ---- Load inputs ----
ranked_sample_path <- file.path(opt$probe_rank_dir, "ranked_sample.rds")
label_pool_path    <- file.path(opt$probe_rank_dir, "label_pool_genes.csv")
moransi_path       <- file.path(opt$qc_backgrounds_dir, "moransi_combined.rds")
walk(c(ranked_sample_path, label_pool_path, moransi_path),
     function(p) if (!file.exists(p)) stop("Not found: ", p))

message("Loading per-sample probe rank table...")
ranked_df <- readRDS(ranked_sample_path) %>%
  mutate(
    animal   = factor(sample_to_animal(Sample), levels = animal_levels),
    platform = factor(platform, levels = platform_levels),
    class    = case_when(
      Type == "Blank" ~ "background",
      is_overlap      ~ "Gene in background",
      TRUE            ~ "Target"
    )
  )

low_rank_genes <- sort(utils::read.csv(label_pool_path)$feature)
message("Low-rank gene set (", length(low_rank_genes), "): ", paste(low_rank_genes, collapse = ", "))

message("Loading Moran's I combined data...")
moransi_df <- readRDS(moransi_path) %>%
  filter(type == "Target", feature %in% low_rank_genes)

na_rows <- moransi_df %>% filter(is.na(observed))
message("Dropping ", nrow(na_rows), " low-rank gene values with NA Moran's I",
        if (nrow(na_rows) > 0) paste0(": ", paste(na_rows$sample, na_rows$feature, collapse = "; ")) else "")

moransi_df <- moransi_df %>%
  filter(!is.na(observed)) %>%
  mutate(
    platform = factor(platform, levels = platform_levels),
    feature  = factor(feature, levels = low_rank_genes)
  )

missing_genes <- setdiff(low_rank_genes, unique(as.character(moransi_df$feature)))
if (length(missing_genes) > 0) stop("No Moran's I for: ", paste(missing_genes, collapse = ", "))

# ---- Missing-animal checks ----
check_missing_animals(ranked_df,  "Sample", "Probe rank",
                      animal_levels, excluded_samples, platform_levels)
check_missing_animals(moransi_df, "sample", "Moran's I",
                      animal_levels, excluded_samples, platform_levels)

# ---- Panel 1: per-sample probe rank curves ----
message("Building probe_rank_persample panel...")

probe_page_h_mm <- 68

probe_base_size <- 6
small_text_pt   <- probe_base_size - 1

point_size       <- 0.5
blank_point_size <- 0.42
overlap_size     <- 0.9
label_size_pt    <- 5
max_labels       <- 8
overlap_colour   <- "#800020"

probe_ylim <- c(0, max(ranked_df$y, na.rm = TRUE) * 1.04)

threshold_df <- ranked_df %>%
  distinct(platform, animal, bg95) %>%
  mutate(yline = log10(bg95 + 1))

label_choice <- ranked_df %>%
  filter(class == "Gene in background") %>%
  mutate(priority = if_else(feature %in% low_rank_genes, 1L, 2L)) %>%
  group_by(platform, animal) %>%
  arrange(priority, desc(rank), .by_group = TRUE) %>%
  slice_head(n = max_labels) %>%
  ungroup() %>%
  select(platform, animal, feature) %>%
  mutate(is_label = TRUE)

message("Labels per sample: ",
        paste0(label_choice %>% count(platform, animal) %>%
                 transmute(x = paste0(platform, " ", animal, " ", n)) %>% pull(x),
               collapse = ", "))

label_df <- ranked_df %>%
  left_join(label_choice, by = c("platform", "animal", "feature")) %>%
  mutate(is_label = coalesce(is_label, FALSE)) %>%
  arrange(platform, animal, rank_frac) %>%
  group_by(platform, animal) %>%
  mutate(
    side     = if_else(cumsum(is_label) %% 3 == 1, "above", "below"),
    label    = if_else(is_label, feature, ""),
    nudge_y  = case_when(!is_label ~ 0, side == "above" ~ 1.2, TRUE ~ -1.8)
  ) %>%
  ungroup()

class_shapes <- c("background" = 17, "Gene in background" = 21)

plot_probe_row <- function(plat) {
  row_df   <- ranked_df   %>% filter(platform == plat)
  row_lab  <- label_df    %>% filter(platform == plat)
  row_thr  <- threshold_df %>% filter(platform == plat)
  show_key <- any(row_df$class == "Gene in background")

  ggplot(row_df, aes(x = rank_frac, y = y)) +
    geom_hline(data = row_thr, aes(yintercept = yline),
               linetype = "dashed", colour = "grey30", linewidth = 0.25) +
    geom_point(data = filter(row_df, class == "background"),
               aes(shape = "background"),
               colour = "grey50", size = blank_point_size, alpha = 0.7,
               show.legend = show_key) +
    geom_point(data = filter(row_df, class == "Target"),
               aes(colour = y), shape = 16, size = point_size, alpha = 0.85) +
    geom_point(data = filter(row_df, class == "Gene in background"),
               aes(shape = "Gene in background"),
               fill = overlap_colour, colour = "white", stroke = 0.1,
               size = overlap_size, show.legend = show_key) +
    # Labels are pre-selected; max.overlaps = Inf so none are dropped
    geom_text_repel(
      data               = row_lab,
      aes(label = label),
      position           = position_nudge_repel(x = 0, y = row_lab$nudge_y),
      size               = label_size_pt / .pt,
      fontface           = "italic",
      colour             = "black",
      point.size         = point_size,
      segment.colour     = "grey65",
      segment.size       = 0.1,
      segment.curvature  = -0.15,
      segment.ncp        = 3,
      segment.angle      = 20,
      box.padding        = 0.12,
      point.padding      = 0.05,
      min.segment.length = unit(1.5, "mm"),
      force              = 3,
      force_pull         = 0.5,
      max.overlaps       = Inf,
      max.iter           = 1e5,
      max.time           = 5,
      seed               = 123
    ) +
    scale_colour_viridis_c(
      name   = expression(log[10] ~ "(total count + 1)"),
      option = "D",
      guide  = guide_colourbar(
        order          = 1,
        title.position = "left",
        title.theme    = element_text(size = small_text_pt, angle = 90, hjust = 0.5),
        barwidth       = unit(1.5, "mm"),
        barheight      = unit(10, "mm")
      )
    ) +
    scale_shape_manual(
      name   = NULL,
      values = class_shapes,
      labels = c("background" = "background", "Gene in background" = "Gene in\nbackground"),
      guide  = guide_legend(
        order        = 2,
        override.aes = list(colour = c("grey50", "white"),
                            fill   = c(NA, overlap_colour),
                            size   = c(0.8, 1.1))
      )
    ) +
    facet_animals() +
    scale_x_continuous(
      "Probe rank (percentile)",
      breaks = seq(0, 1, by = 0.25),
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    # Limits on the coord, not the scale, so nudged label obstacles are kept
    coord_cartesian(xlim = c(0, 1), ylim = probe_ylim, clip = "off") +
    labs(y = expression(log[10] ~ "(total count + 1)"), title = plat) +
    theme_sb(base_size = probe_base_size) +
    theme_ext +
    theme(
      plot.title       = element_text(size = probe_base_size),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = small_text_pt),
      axis.text.y      = element_text(size = small_text_pt),
      axis.ticks       = element_line(linewidth = 0.2),
      axis.ticks.length = unit(0.6, "mm"),
      strip.text       = element_text(face = "italic", size = small_text_pt,
                                      margin = margin(0.8, 0, 0.8, 0, "mm")),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.15),
      panel.spacing    = unit(2.5, "mm"),
      legend.position  = "right",
      legend.text      = element_text(size = small_text_pt),
      legend.key.size  = unit(2, "mm"),
      legend.spacing.y = unit(1, "mm"),
      legend.margin    = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1, "mm")
    )
}

p_probe <- wrap_plots(map(platform_levels, plot_probe_row), ncol = 1)

panel_h_probe_mm <- panel_h_for_page(p_probe, probe_page_h_mm, n_rows = length(platform_levels))

save_fixed_panels(p_probe, file.path(opt$out_dir, "probe_rank_persample.pdf"),
                  panel_h_mm = panel_h_probe_mm, n_rows = length(platform_levels))

# ---- Panel 2: Moran's I of the low-rank genes ----
message("Building moransi_low_rank_genes panel...")

p_moransi <- ggplot(moransi_df, aes(x = platform, y = observed)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  geom_boxplot(aes(fill = platform), outlier.shape = NA, width = 0.6,
               linewidth = 0.3, alpha = 0.5) +
  geom_point(aes(colour = platform),
             position = position_jitter(width = 0.15, height = 0, seed = 42),
             size = 0.9, alpha = 0.8) +
  scale_fill_platform() +
  scale_colour_platform() +
  facet_wrap(~ feature, ncol = 4) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Observed Moran's I") +
  theme_sb() +
  theme_bg_panels +
  theme_ext +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(
  file.path(opt$out_dir, "moransi_low_rank_genes.pdf"),
  p_moransi,
  width  = dims$full_w,
  height = dims$half_w * 1.05,
  units  = "mm",
  device = cairo_pdf,
  bg     = "white"
)
message("Saved: moransi_low_rank_genes.pdf")

message("Done. All panels written to: ", opt$out_dir)
