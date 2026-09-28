# Purpose:  Fig 2 extended: per-bin counts/genes violins per animal, one PDF per bin size.
# Supports: fig2_qc.R — per-bin distributions are consistent across animals.
# Inputs:   results/03_benchmarking/qc_metrics/metadata_combined.rds, config/config.yaml
# Outputs:  extended_figures/fig2ext_qc_violins/qc_violins_{8um,16um}.pdf

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

# ---- CLI arguments ----
option_list <- list(
  make_option(c("--metadata"), type = "character",
              default = "results/03_benchmarking/qc_metrics/metadata_combined.rds",
              help    = "Per-bin QC metadata from qc_metrics.R [default: %default]"),
  make_option(c("--config"), type = "character",
              default = "config/config.yaml",
              help    = "Pipeline config (sample lists, gene_comparison.animals) [default: %default]"),
  make_option(c("--out_dir"), type = "character",
              default = "extended_figures/fig2ext_qc_violins",
              help    = "Output directory for panel PDFs [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opt$metadata)) stop("Not found: ", opt$metadata)
cfg <- yaml::read_yaml(opt$config)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

columns <- tibble::tribble(
  ~column,   ~platform,  ~Subset, ~title,
  "vis_all", "VisiumHD", "All",   "Visium HD\nall genes",
  "vis_90",  "VisiumHD", "90",    "Visium HD\n90 common genes",
  "mer_90",  "MERSCOPE", "90",    "MERSCOPE\n90 common genes",
  "xen_90",  "Xenium",   "90",    "Xenium\n90 common genes"
)
metrics <- c(nCount = "Counts\nper bin", nFeature = "Genes\nper bin")

trim_quantile  <- 0.99
label_band_mm  <- 2.4
violin_adjust  <- 2
genotype_gap   <- 1.25
panel_h_mm     <- 12

# ---- Load per-bin metadata and attach animal IDs ----
message("Loading per-bin QC metadata...")
meta <- readRDS(opt$metadata) %>%
  mutate(platform = as.character(platform), Subset = as.character(Subset)) %>%
  inner_join(select(columns, column, platform, Subset), by = c("platform", "Subset"))

vis_map        <- visiumhd_animal_map(cfg)
spatial_levels <- animal_levels_from_config(cfg)
animal_levels  <- order_animals(union(unname(vis_map), spatial_levels))

sample_tbl <- meta %>%
  distinct(platform, Sample) %>%
  mutate(
    animal   = if_else(platform == "VisiumHD", unname(vis_map[Sample]), sample_to_animal(Sample)),
    genotype = factor(toupper(sub("[0-9]+$", "", animal)), levels = c("KO", "WT", "CTRL"))
  )
unknown <- sample_tbl %>% filter(is.na(animal) | !animal %in% animal_levels | is.na(genotype))
if (nrow(unknown) > 0) {
  stop("Samples with no config animal: ", paste(unknown$platform, unknown$Sample, collapse = ", "))
}

meta <- meta %>%
  left_join(sample_tbl, by = c("platform", "Sample")) %>%
  mutate(animal = factor(animal, levels = animal_levels))

no_exclusions <- list(VisiumHD = character(0), MERSCOPE = character(0), Xenium = character(0))
check_missing_animals(meta, "animal", "QC violins", order_animals(unname(vis_map)),
                      no_exclusions, "VisiumHD")
check_missing_animals(meta, "animal", "QC violins", spatial_levels,
                      no_exclusions, c("MERSCOPE", "Xenium"))

# ---- Medians ----
animal_medians <- meta %>%
  group_by(bin_size, column, animal) %>%
  summarise(across(all_of(names(metrics)), median), .groups = "drop")

column_medians <- animal_medians %>%
  group_by(bin_size, column) %>%
  summarise(across(all_of(names(metrics)), median), .groups = "drop")

message("Per-animal medians:")
print(as.data.frame(animal_medians))
message("Column medians (median of per-animal medians):")
print(as.data.frame(column_medians))

# Half rounds up (round() would round half to even)
fmt_median <- function(x) comma(floor(x + 0.5), accuracy = 1)

# ---- Y axes (coord_cartesian, so boxplots use all bins) ----
axis_group <- function(col) if_else(col == "vis_all", "vis_all", "shared_90")

animal_trim <- meta %>%
  group_by(bin_size, column, animal) %>%
  summarise(across(all_of(names(metrics)), ~ unname(quantile(.x, trim_quantile))), .groups = "drop")

three_breaks <- function(top) {
  steps <- as.vector(outer(c(1, 1.5, 2, 2.5, 3, 4, 5), 10^(0:5)))
  step  <- max(steps[2 * steps <= top])
  c(0, step, 2 * step)
}

axis_spec <- function(bin, metric, group) {
  trims <- filter(animal_trim, bin_size == bin, axis_group(column) == group)[[metric]]
  genes_90 <- metric == "nFeature" && group == "shared_90"
  top <- if (genes_90) 90 else max(trims)
  list(
    ylim   = c(-0.02 * top, top),
    breaks = if (genes_90) c(0, 45, 90) else three_breaks(top)
  )
}

# ---- Panels ----
median_label <- function(label) {
  txt <- grid::textGrob(label, x = unit(1, "npc"), y = unit(1, "npc") + unit(0.3, "mm"),
                        just = c("right", "bottom"),
                        gp = grid::gpar(fontsize = 6, fontfamily = "Arial"))
  # Blank background, else the inset hides the panels
  wrap_elements(full = txt, clip = FALSE) +
    theme(plot.background = element_blank(), plot.margin = margin(0, 0, 0, 0))
}

plot_panel <- function(panel_df, trim_df, metric, col, spec, median_val,
                       first_col, top_row, bottom_row) {
  col_info <- filter(columns, column == col)

  trims     <- setNames(trim_df[[metric]], as.character(trim_df$animal))
  violin_df <- filter(panel_df, .data[[metric]] <= trims[as.character(animal)])

  p <- ggplot(panel_df, aes(x = animal, y = .data[[metric]])) +
    geom_violin(data = violin_df, fill = pal_muted[[col_info$platform]], alpha = 0.7, colour = NA,
                scale = "width", width = 0.85, adjust = violin_adjust, trim = TRUE) +
    geom_hline(yintercept = median_val, linetype = "dashed", linewidth = 0.3, colour = "grey20") +
    geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.25) +
    facet_grid(. ~ genotype, scales = "free_x", space = "free_x", switch = "x") +
    # 0.075 each side = one unit per animal; col_widths() relies on this
    scale_x_discrete(labels = function(x) sub("^[a-z]+", "", x),
                     expand = expansion(add = 0.075)) +
    scale_y_continuous(breaks = spec$breaks, labels = label_comma(), expand = expansion(0)) +
    coord_cartesian(ylim = spec$ylim) +
    labs(x = NULL, y = if (first_col) metrics[[metric]] else NULL,
         title = if (top_row) col_info$title else NULL) +
    theme_sb() +
    theme_ext +
    theme(
      panel.border     = element_blank(),
      axis.line        = element_line(colour = "black", linewidth = 0.3),
      strip.background = element_blank(),
      strip.placement  = "outside",
      axis.ticks.length = unit(0.6, "mm"),
      axis.text.x      = element_text(margin = margin(0.3, 0, 0, 0, "mm")),
      axis.text.y      = element_text(margin = margin(0, 0.3, 0, 0, "mm")),
      strip.text.x     = element_text(size = 6, margin = margin(0.2, 0, 0, 0, "mm")),
      strip.switch.pad.grid = unit(0, "mm"),
      panel.spacing.x  = unit(genotype_gap, "mm"),
      plot.title       = element_text(lineheight = 0.9, margin = margin(0, 0, label_band_mm, 0, "mm")),
      plot.margin      = margin(if (top_row) 0.5 else 0.5 + label_band_mm, 1.5, 0.3, 1.5, "mm"),
      legend.position  = "none"
    )

  if (!bottom_row) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())
  }
  p + inset_element(median_label(paste("median =", fmt_median(median_val))),
                    left = 0, bottom = 0, right = 1, top = 1,
                    align_to = "panel", clip = FALSE, on_top = TRUE, ignore_tag = TRUE)
}

plot_bin_size <- function(bin) {
  bin_df <- filter(meta, bin_size == bin)
  grid   <- tidyr::expand_grid(metric = names(metrics), column = columns$column)

  above <- map_dfr(names(metrics), function(metric) {
    map_dfr(columns$column, function(col) {
      top  <- axis_spec(bin, metric, axis_group(col))$ylim[2]
      vals <- bin_df[[metric]][bin_df$column == col]
      tibble::tibble(bin_size = bin, metric = metric, column = col, axis_top = top,
                     pct_above = 100 * mean(vals > top))
    })
  })
  message("% of bins above the axis top")
  print(as.data.frame(above))

  panels <- pmap(grid, function(metric, column) {
    median_val <- column_medians %>%
      filter(bin_size == bin, .data$column == !!column) %>%
      pull(all_of(metric))
    plot_panel(
      panel_df   = filter(bin_df, .data$column == !!column) %>% mutate(animal = droplevels(animal)),
      trim_df    = filter(animal_trim, bin_size == bin, .data$column == !!column),
      metric     = metric,
      col        = column,
      spec       = axis_spec(bin, metric, axis_group(column)),
      median_val = median_val,
      first_col  = column == columns$column[1],
      top_row    = metric == names(metrics)[1],
      bottom_row = metric == tail(names(metrics), 1)
    )
  })

  n_animals  <- map_int(columns$column, ~ n_distinct(bin_df$animal[bin_df$column == .x]))
  n_gaps     <- map_int(columns$column, ~ n_distinct(bin_df$genotype[bin_df$column == .x])) - 1L
  col_widths <- function(slot_mm) unit(n_animals * slot_mm + n_gaps * genotype_gap, "mm")
  layout_for <- function(slot_mm) {
    wrap_plots(panels, nrow = length(metrics), byrow = TRUE, widths = col_widths(slot_mm))
  }

  probe_slot <- 5
  overhead   <- layout_width_mm(
    layout_for(probe_slot) + plot_layout(heights = unit(rep(panel_h_mm, length(metrics)), "mm"))
  ) - sum(n_animals) * probe_slot
  slot_mm <- (dims$full_w - overhead) / sum(n_animals)
  message(sprintf("Violin slot width: %.2f mm", slot_mm))
  layout_for(slot_mm)
}

# ---- Save one PDF per bin size ----
walk(unique(meta$bin_size), function(bin) {
  message("Building QC violins @ ", bin, "...")
  save_fixed_panels(plot_bin_size(bin), file.path(opt$out_dir, paste0("qc_violins_", bin, ".pdf")),
                    panel_h_mm = panel_h_mm, n_rows = length(metrics))
})

message("Done. All panels written to: ", opt$out_dir)
