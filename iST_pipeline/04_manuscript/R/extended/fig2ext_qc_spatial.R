# Purpose:  Extended figure for Figure 2: per-bin spatial maps of counts and genes
#           per bin (all genes, 8 µm), one PDF per platform, one tile per animal
#           in genotype groups KO, WT, CTRL; each tile is a separate image in the
#           PDF and is also written as a transparent PNG.
# Supports: fig2_qc.R — per-bin counts and genes are spatially consistent within
#           and across animals (same bins as fig2ext_qc_violins.R). ROI-level QC
#           for the matched animals: 03_benchmarking/R/qc_metrics_roi.R.
# Inputs:   results/03_benchmarking/qc_spatial/bins_8um.rds (qc_spatial.R)
#           config/config.yaml, config/qc_spatial.yaml
# Outputs:  extended_figures/fig2ext_qc_spatial/qc_spatial_{xenium,merscope,visiumhd}_8um.pdf
#           extended_figures/fig2ext_qc_spatial/tiles/{platform}_{animal}_{metric}.png

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

option_list <- list(
  make_option("--bins",       type = "character", default = "results/03_benchmarking/qc_spatial/bins_8um.rds"),
  make_option("--config",     type = "character", default = "config/config.yaml"),
  make_option("--qc_spatial", type = "character", default = "config/qc_spatial.yaml"),
  make_option("--out_dir",    type = "character", default = "extended_figures/fig2ext_qc_spatial"),
  make_option("--bin_size",   type = "integer",   default = 8L)
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg  <- yaml::read_yaml(opt$config)
disp <- yaml::read_yaml(opt$qc_spatial)$qc_spatial
res  <- opt$bin_size
tile_dir <- file.path(opt$out_dir, "tiles")
dir.create(tile_dir, recursive = TRUE, showWarnings = FALSE)

platforms <- c(Xenium = "Xenium, all genes", MERSCOPE = "MERSCOPE, all genes",
               VisiumHD = "Visium HD, all genes")
metrics   <- c(nCount = "Counts\nper bin", nFeature = "Genes\nper bin")

tile_pad_um  <- 200
tile_gap_mm  <- 0.6
gap_mm       <- 1.5      # between genotype groups
header_mm    <- 5
row_gap_mm   <- 1
title_w_mm   <- 7
legend_w_mm  <- 11
margin_mm    <- 1
page_h_mm    <- 50.8     # the scale-bar row takes up the remainder
scale_mm     <- 4.5      # probe value, solved below
scale_bar_um <- 1000

# ---------------------------------------------------------------------------
# Bins and animals
# ---------------------------------------------------------------------------
bins <- readRDS(opt$bins) %>%
  mutate(platform = as.character(platform)) %>%
  filter(main_section)

vis_map    <- visiumhd_animal_map(cfg)
sample_tbl <- bins %>%
  distinct(platform, Sample) %>%
  mutate(animal   = if_else(platform == "VisiumHD", unname(vis_map[Sample]), sample_to_animal(Sample)),
         genotype = factor(toupper(sub("[0-9]+$", "", animal)), levels = c("KO", "WT", "CTRL")))
if (anyNA(sample_tbl$animal) || anyNA(sample_tbl$genotype)) stop("Samples with no config animal")

no_exclusions <- list(VisiumHD = character(0), MERSCOPE = character(0), Xenium = character(0))
check_missing_animals(mutate(sample_tbl, Sample = animal), "Sample", "QC spatial",
                      order_animals(unname(vis_map)), no_exclusions, "VisiumHD")
check_missing_animals(mutate(sample_tbl, Sample = animal), "Sample", "QC spatial",
                      animal_levels_from_config(cfg), no_exclusions, c("MERSCOPE", "Xenium"))

# ---------------------------------------------------------------------------
# Orientation and centring
# ---------------------------------------------------------------------------
# Display frame is y down: mirror left-right, then rotate clockwise
orient <- function(x, y, rotate, flip) {
  if (!rotate %in% c(0, 90, 180, 270)) stop("rotate must be 0/90/180/270, got ", rotate)
  if (flip) x <- -x
  for (i in seq_len(rotate / 90)) {
    x_new <- -y
    y     <- x
    x     <- x_new
  }
  list(x = x, y = y)
}

sample_orientation <- function(platform, sample, animal, aligned) {
  if (aligned) {
    extra <- disp$matched_rotate[[animal]] %||% 0
    return(list(rotate = (disp$matched_base$rotate + extra) %% 360, flip = disp$matched_base$flip))
  }
  o <- disp$orientation[[tolower(platform)]][[sample]]
  if (is.null(o)) stop("No orientation for ", platform, " ", sample, " in ", opt$qc_spatial)
  list(rotate = o$rotate, flip = isTRUE(o$flip))
}

bins <- bins %>%
  left_join(select(sample_tbl, platform, Sample, animal), by = c("platform", "Sample")) %>%
  group_by(platform, Sample) %>%
  group_modify(function(df, key) {
    o  <- sample_orientation(key$platform, key$Sample, df$animal[1], df$aligned[1])
    xy <- orient(df$x_um, df$y_um, o$rotate, o$flip)
    mutate(df, x_disp = xy$x, y_disp = xy$y)
  }) %>%
  ungroup()

# Aligned animals share their VisiumHD section centre (registered across
# platforms); others are centred on their own extent
centres <- bins %>%
  group_by(platform, Sample, animal, aligned) %>%
  summarise(cx = mean(range(x_disp)), cy = mean(range(y_disp)), .groups = "drop")
vis_centres <- centres %>% filter(platform == "VisiumHD") %>% select(animal, vis_cx = cx, vis_cy = cy)
centres <- centres %>%
  left_join(vis_centres, by = "animal") %>%
  mutate(cx = if_else(aligned, vis_cx, cx), cy = if_else(aligned, vis_cy, cy))

bins <- bins %>%
  left_join(select(centres, platform, Sample, cx, cy), by = c("platform", "Sample")) %>%
  mutate(x_disp = x_disp - cx, y_disp = y_disp - cy)

half_span <- max(abs(c(bins$x_disp, bins$y_disp))) + res / 2
tile_um   <- ceiling((2 * half_span + 2 * tile_pad_um) / 200) * 200

# Tile coordinates: origin at the tile's top-left, y down
bins <- mutate(bins, x_tile = x_disp + tile_um / 2, y_tile = y_disp + tile_um / 2)

# ---------------------------------------------------------------------------
# Colour caps
# ---------------------------------------------------------------------------
caps <- bins %>%
  group_by(platform) %>%
  summarise(across(all_of(names(metrics)), ~ list(unname(quantile(.x, disp$caps_quantiles)))),
            .groups = "drop")
cap_for <- function(plat, metric) caps[[metric]][[which(caps$platform == plat)]]
message("Colour caps:")
walk(caps$platform, function(p) walk(names(metrics), function(m)
  message(sprintf("  %-8s %-8s %s", p, m, paste(round(cap_for(p, m), 1), collapse = " - ")))))

fill_scale <- function(limits, guide = "none") {
  scale_fill_viridis_c(option = "D", limits = limits, oob = squish, na.value = NA,
                       breaks = limits, labels = label_comma(accuracy = 1), guide = guide)
}

# ---------------------------------------------------------------------------
# Tile images
# ---------------------------------------------------------------------------
# One colour matrix per tile (row 1 = top, NA = no bin): each image pixel is
# the mean of the bins whose centres fall in it. Built directly rather than by
# a raster device so the image covers exactly its tile
tile_image <- function(df, metric, limits) {
  n_img <- round(tile_mm / 25.4 * disp$dpi)
  px_um <- tile_um / n_img
  idx   <- (floor(df$x_tile / px_um)) * n_img + floor(df$y_tile / px_um) + 1
  sums  <- rowsum(as.numeric(df[[metric]]), idx)
  n     <- rowsum(rep(1, nrow(df)), idx)

  sc <- fill_scale(limits)
  sc$train(limits)
  img <- matrix(NA_character_, n_img, n_img)
  img[as.integer(rownames(sums))] <- sc$map(sums[, 1] / n[, 1])
  img
}

tile_png <- function(img, file) {
  rgba <- grDevices::col2rgb(ifelse(is.na(img), "#00000000", img), alpha = TRUE) / 255
  png::writePNG(aperm(array(rgba, c(4, nrow(img), ncol(img))), c(2, 3, 1)), file, dpi = disp$dpi)
}

sample_bins <- function(plat, sample) filter(bins, platform == plat, Sample == sample)

# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------
widest   <- names(which.max(table(sample_tbl$platform)))
n_tiles  <- sum(sample_tbl$platform == widest)
set_tile_mm <- function(mm) {
  tile_mm  <<- mm
  um_to_mm <<- mm / tile_um
}
set_tile_mm(15)   # probe; solved from the measured layout below

tile_layout <- function(plat) {
  sample_tbl %>%
    filter(platform == plat) %>%
    mutate(animal = factor(animal, levels = order_animals(animal))) %>%
    arrange(animal) %>%
    mutate(label  = sub("^[a-z]+", "", animal),
           gap_mm = case_when(row_number() == 1         ~ 0,
                              genotype != lag(genotype) ~ gap_mm,
                              TRUE                      ~ tile_gap_mm),
           x0_mm  = (row_number() - 1) * tile_mm + cumsum(gap_mm))
}

tile_widths <- function(layout) as.vector(rbind(layout$gap_mm, tile_mm))[-1]

# coord_cartesian, not coord_fixed(): fixed-aspect panels make patchwork drop
# the absolute mm sizes, and each cell is already an exact square
map_tile <- function(img, plat) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = tile_um, ymin = 0, ymax = tile_um,
             fill = pal_muted_light[[plat]], colour = NA) +
    annotation_raster(img, xmin = 0, xmax = tile_um, ymin = 0, ymax = tile_um, interpolate = FALSE) +
    coord_cartesian(xlim = c(0, tile_um), ylim = c(0, tile_um), expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# Vector gradient colourbar, taken from a ggplot legend
colourbar <- function(limits) {
  p <- ggplot(tibble(x = 0, y = 0, v = limits), aes(x, y, fill = v)) +
    geom_raster() +
    fill_scale(limits, guide = guide_colourbar(
      display = "gradient", barwidth = unit(1.5, "mm"), barheight = unit(0.7 * tile_mm, "mm"))) +
    labs(fill = NULL) +
    theme_sb() +
    theme(legend.ticks = element_blank(), legend.margin = margin(0, 0, 0, 0),
          legend.justification = c(0, 0.5))
  gt <- ggplotGrob(p)
  wrap_elements(full = gt$grobs[[which(gt$layout$name == "guide-box-right")]]) +
    theme(plot.margin = margin(0, 0, 0, 1.5, "mm"))
}

text_cell <- function(label, rot = 0) {
  wrap_elements(full = grid::textGrob(label, rot = rot, gp = grid::gpar(
    fontsize = 7, fontfamily = "Arial", lineheight = 0.9))) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# Platform title, genotype labels over a rule per group, animal numbers (x in mm)
header_plot <- function(layout, width_mm, title) {
  groups <- layout %>%
    group_by(genotype) %>%
    summarise(x_start = min(x0_mm), x_end = max(x0_mm) + tile_mm, .groups = "drop")
  inset <- 0.03 * tile_mm
  ggplot() +
    annotate("text", x = (groups$x_start + groups$x_end) / 2, y = 0.8,
             label = as.character(groups$genotype), size = 6 / .pt, family = "Arial") +
    annotate("segment", x = groups$x_start + inset, xend = groups$x_end - inset,
             y = 0.5, yend = 0.5, linewidth = 0.3) +
    annotate("text", x = layout$x0_mm + tile_mm / 2, y = 0.25,
             label = layout$label, size = 6 / .pt, family = "Arial") +
    scale_x_continuous(limits = c(0, width_mm), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(title = title) +
    theme_void(base_family = "Arial") +
    theme(plot.title  = element_text(face = "italic", size = 7, hjust = 0,
                                     margin = margin(0, 0, 0.5, 0, "mm")),
          plot.margin = margin(0, 0, 0, 0))
}

# Scale bar under the last tile; y in mm below the tile row
scale_bar_plot <- function() {
  bar_x2 <- 0.95 * tile_um
  ggplot() +
    annotate("segment", x = bar_x2 - scale_bar_um, xend = bar_x2, y = -1, yend = -1,
             linewidth = 0.4, lineend = "butt") +
    annotate("text", x = bar_x2 - scale_bar_um / 2, y = -1.6, label = paste(scale_bar_um / 1000, "mm"),
             vjust = 1, size = 6 / .pt, family = "Arial") +
    scale_x_continuous(limits = c(0, tile_um), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-scale_mm, 0), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# One flat patchwork grid with zero-margin cells so the mm sizes are exact.
# Columns: row title, tiles and (empty) gap columns, colourbar. Rows: header,
# metric rows with a gap row between, scale bar
build_platform <- function(plat) {
  layout     <- tile_layout(plat)
  widths     <- c(title_w_mm, tile_widths(layout), legend_w_mm)
  n_col      <- length(widths)
  tile_col   <- 1 + seq(1, 2 * nrow(layout) - 1, by = 2)
  metric_row <- 2 * seq_along(metrics)
  n_row      <- max(metric_row) + 1

  cells <- list(header_plot(layout, sum(tile_widths(layout)), platforms[[plat]]))
  areas <- list(area(1, 2, 1, n_col - 1))
  add   <- function(cell, a) {
    cells[[length(cells) + 1]] <<- cell
    areas[[length(areas) + 1]] <<- a
  }
  iwalk(names(metrics), function(metric, i) {
    limits <- cap_for(plat, metric)
    add(text_cell(metrics[[metric]], rot = 90), area(metric_row[i], 1))
    walk2(layout$Sample, tile_col, function(sample, col) {
      add(map_tile(tile_image(sample_bins(plat, sample), metric, limits), plat), area(metric_row[i], col))
    })
    add(colourbar(limits), area(metric_row[i], n_col))
  })
  add(scale_bar_plot(), area(n_row, max(tile_col)))

  heights <- c(header_mm, head(rep(c(tile_mm, row_gap_mm), length(metrics)), -1), scale_mm)
  wrap_plots(cells) +
    plot_layout(design = do.call(c, areas), widths = unit(widths, "mm"), heights = unit(heights, "mm")) +
    plot_annotation(theme = theme(plot.margin = margin(margin_mm, margin_mm, margin_mm, margin_mm, "mm")))
}

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
# Measure text on a cairo device (Arial), not the default pdf() device
cairo_pdf(tempfile(fileext = ".pdf"))

# Widths and heights are linear in the tile size and scale-bar row, so one
# measurement of the fixed overhead solves both
probe_mm <- tile_mm
set_tile_mm((dims$full_w - (layout_width_mm(build_platform(widest)) - n_tiles * probe_mm)) / n_tiles)
scale_mm <- scale_mm + page_h_mm - layout_height_mm(build_platform(widest))
if (scale_mm < 3) stop("page_h_mm = ", page_h_mm, " leaves ", round(scale_mm, 1), " mm for the scale bar")
message(sprintf("Tile %.2f mm for %.1f mm of tissue; scale-bar row %.1f mm", tile_mm, tile_um / 1000, scale_mm))

walk(names(platforms), function(plat) {
  p      <- build_platform(plat)
  page_w <- layout_width_mm(p)
  file   <- file.path(opt$out_dir, sprintf("qc_spatial_%s_%dum.pdf", tolower(plat), res))
  ggsave(file, p, width = page_w, height = layout_height_mm(p), units = "mm", device = cairo_pdf, bg = "white")
  message(sprintf("Saved: %s (%.1f mm wide, %.1f MB)", basename(file), page_w, file.size(file) / 1e6))

  layout <- tile_layout(plat)
  walk2(layout$Sample, as.character(layout$animal), function(sample, animal) {
    walk(names(metrics), function(metric) {
      tile_png(tile_image(sample_bins(plat, sample), metric, cap_for(plat, metric)),
               file.path(tile_dir, sprintf("%s_%s_%s.png", tolower(plat), animal, metric)))
    })
  })
})
message("Done: ", opt$out_dir)
