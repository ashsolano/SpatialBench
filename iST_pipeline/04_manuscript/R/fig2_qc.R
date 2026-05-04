# Purpose:  Figure 2 QC panels — platform QC panels.
#           Spatial scatter montage coloured by nCount and nFeature across
#           VisiumHD, MERSCOPE, and Xenium (3 platforms × 4 samples, 8µm bins).
#           Per-platform boxplots of median counts/bin and genes/bin (8µm bins,
#           90-gene common subset). Per-cell intersect-gene boxplots comparing
#           VisiumHD, FLEX snRNA-seq, and scRNA-seq (optional; requires
#           metadata_flex_scrna.rds from qc_metrics.R --sc_rds).
#           Each panel saved as a separate PDF.
#           Adapted from montage_v2.R, fig2_qcmetrics.R, flex_visiumhd_qc_v2.R.
# Inputs:   config/config.yaml
#           results/01_preprocessing/merscope_8um/{sample}_8um.rds
#           results/01_preprocessing/xenium_8um/{sample}_8um.rds
#           results/03_benchmarking/qc_metrics/metadata_combined.rds
#           results/03_benchmarking/qc_metrics/metadata_flex_scrna.rds  (optional)
#           results/03_benchmarking/qc_metrics/genes_intersect_flex.rds (optional)
# Outputs:  figures/fig2/spatial_ncount.pdf
#           figures/fig2/spatial_nfeature.pdf
#           figures/fig2/qc_counts_spatial.pdf
#           figures/fig2/qc_genes_spatial.pdf
#           figures/fig2/qc_counts_flex_all.pdf       (if FLEX metadata present)
#           figures/fig2/qc_genes_flex_all.pdf        (if FLEX metadata present)
#           figures/fig2/qc_counts_flex_visiumhd.pdf  (if FLEX metadata present)
#           figures/fig2/qc_genes_flex_visiumhd.pdf   (if FLEX metadata present)


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(dbscan)
  library(scattermore)
  library(gridExtra)
  library(grid)
  library(cowplot)
  library(colorspace)
  library(scales)
  library(yaml)
  library(optparse)
})

source("04_manuscript/R/utils/theme.R")
source("04_manuscript/R/utils/palettes.R")

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--config"),    type = "character",
              default = "config/config.yaml",
              help    = "Path to config.yaml [default: %default]"),
  make_option(c("--input_dir"), type = "character",
              default = "results/03_benchmarking/qc_metrics",
              help    = "Directory containing qc_metrics.R outputs [default: %default]"),
  make_option(c("--out_dir"),   type = "character",
              default = "figures/fig2",
              help    = "Output directory for panel PDFs [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

meta_combined_path   <- file.path(opt$input_dir, "metadata_combined.rds")
meta_flex_path       <- file.path(opt$input_dir, "metadata_flex_scrna.rds")
genes_intersect_path <- file.path(opt$input_dir, "genes_intersect_flex.rds")

if (!file.exists(meta_combined_path)) stop("Not found: ", meta_combined_path)

cfg <- yaml::read_yaml(opt$config)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# SPATIAL SCATTER MONTAGE (nCount / nFeature per bin)
# ===========================================================================

# ---------------------------------------------------------------------------
# Sample mapping: display name -> config sample key
# ---------------------------------------------------------------------------
# These four samples are selected for Fig 2 from each platform.
# FLAG: consider adding a fig2_samples section to config.yaml to avoid
#       hardcoding these mappings here.
#
# VisiumHD keys match cfg$visiumhd$samples.
# MERSCOPE/Xenium keys match cfg$spatial_analysis$merscope_samples /
#   xenium_default_samples; bin objects are at results/01_preprocessing/.
vis_sample_map <- c(
  WT709 = "batch33_709",
  WT713 = "batch33_713",
  KO167 = "batch33_167",
  KO168 = "batch33_168"
)
mer_sample_map <- c(
  WT709 = "wt709_batch13",
  WT713 = "wt713_batch13",
  KO167 = "ko167_batch10",
  KO168 = "ko168_batch9"
)
xen_sample_map <- c(
  WT709 = "wt709_batch27",
  WT713 = "wt713_batch24",
  KO167 = "ko167_batch24",
  KO168 = "ko168_batch27"
)

display_samples <- names(vis_sample_map)   # c("WT709","WT713","KO167","KO168")

# ---------------------------------------------------------------------------
# Load Seurat objects for spatial montage
# ---------------------------------------------------------------------------
message("Loading VisiumHD samples...")
visiumhd_dir <- cfg$visiumhd$data_dir

visiumhd_objs <- setNames(
  lapply(vis_sample_map, function(key) {
    path <- file.path(visiumhd_dir, cfg$visiumhd$samples[[key]])
    message("  ", key, ": ", path)
    readRDS(path)
  }),
  display_samples
)

message("Loading MERSCOPE 8um samples...")
mer_8um_dir <- file.path(cfg$output_dir, "01_preprocessing", "merscope_8um")

merscope_objs <- setNames(
  lapply(mer_sample_map, function(key) {
    path <- file.path(mer_8um_dir, paste0(key, "_8um.rds"))
    message("  ", key, ": ", path)
    readRDS(path)
  }),
  display_samples
)

message("Loading Xenium 8um samples...")
xen_8um_dir <- file.path(cfg$output_dir, "01_preprocessing", "xenium_8um")

xenium_objs <- setNames(
  lapply(xen_sample_map, function(key) {
    path <- file.path(xen_8um_dir, paste0(key, "_8um.rds"))
    message("  ", key, ": ", path)
    readRDS(path)
  }),
  display_samples
)

# ---------------------------------------------------------------------------
# Platform colour palette (muted, lightened for tile backgrounds)
# ---------------------------------------------------------------------------
# pal_muted is defined by theme.R; lighten() from colorspace
pal_bg <- colorspace::lighten(c(
  VisiumHD = pal_muted[["VisiumHD"]],
  MERSCOPE = pal_muted[["MERSCOPE"]],
  Xenium   = pal_muted[["Xenium"]]
), amount = 0.2)

# ---------------------------------------------------------------------------
# Spatial data extraction helpers
# ---------------------------------------------------------------------------

# Remove DBSCAN noise and return only the dominant cluster.
# FLAG: eps=8 and minPts=10 are calibrated for 8um-bin coordinates;
#       verify if coordinate units change across datasets.
filter_noise <- function(df, xcol = "x", ycol = "y", eps = 8, minPts = 10) {
  cl   <- dbscan::dbscan(as.matrix(df[, c(xcol, ycol)]), eps = eps, minPts = minPts)$cluster
  keep <- as.integer(names(which.max(table(cl[cl > 0]))))
  df[cl == keep, , drop = FALSE]
}

centre_window      <- function(win) c(x = mean(win$x), y = mean(win$y))
recentre_to_window <- function(df, win) {
  cen_df  <- c(x = mean(df$x), y = mean(df$y))
  cen_win <- centre_window(win)
  df %>% mutate(x = x + (cen_win["x"] - cen_df["x"]),
                y = y + (cen_win["y"] - cen_df["y"]))
}

# Rotate coordinates 90 degrees counter-clockwise (used to align KO168)
rotate90_df <- function(df) {
  df %>% transmute(x = df$y, y = -df$x, across(-c(x, y), identity))
}

make_square <- function(win) {
  xr <- win$x; yr <- win$y
  dx <- diff(xr); dy <- diff(yr)
  if (dx > dy) {
    mid <- mean(yr); yr <- mid + c(-dx / 2, dx / 2)
  } else {
    mid <- mean(xr); xr <- mid + c(-dy / 2, dy / 2)
  }
  list(x = xr, y = yr)
}

pad_window <- function(win, frac = 0.10) {
  dx <- diff(win$x); dy <- diff(win$y)
  list(x = win$x + c(-dx, dx) * frac,
       y = win$y + c(-dy, dy) * frac)
}

# Label grobs for the row (platform) and column (sample) axes
row_label <- function(label) {
  ggplot() +
    annotate("text", x = 0.05, y = 0.5, label = label,
             family = "Arial", fontface = "bold", size = 5,
             angle = 90, hjust = 0.5) +
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    theme_void()
}

col_label <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.05, label = label,
             family = "Arial", fontface = "bold", size = 5, vjust = 0) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme_void()
}

# ---------------------------------------------------------------------------
# Extract coordinate + QC data frames from Seurat objects
# ---------------------------------------------------------------------------

# VisiumHD: tissue coordinates scaled by the low-res image scale factor.
# Filters bins with fewer than 5 counts to remove empty spots.
make_df_visium <- function(seu, bin_suffix = "008um") {
  assay_nm <- paste0("Spatial.", bin_suffix)
  px       <- GetTissueCoordinates(seu, assay = assay_nm)
  coords   <- as.data.frame(px)[, 1:2]
  colnames(coords) <- c("x_px", "y_px")

  # Locate the matching image and retrieve its scale factor
  img_idx <- grep(paste0("\\.", bin_suffix, "$"), names(seu@images))
  sf      <- seu@images[[img_idx]]@scale.factors[["lowres"]]

  coords <- coords %>%
    transmute(x = x_px * sf, y = y_px * sf)

  meta <- FetchData(seu, vars = c(paste0("nCount_",   assay_nm),
                                   paste0("nFeature_", assay_nm)))
  colnames(meta) <- c("nCount", "nFeature")

  bind_cols(coords, meta) %>%
    mutate(spot_id = rownames(.)) %>%
    filter(nCount >= 5)
}

# MERSCOPE / Xenium: centroid coordinates converted from pixels to µm.
# FLAG: scale_um = 0.02721088 is instrument-specific (MERSCOPE pixel pitch).
#       Verify against your microscope calibration if data changes.
make_df_image <- function(seu, prefix, scale_um = 0.02721088) {
  raw    <- seu@images[[1]]@boundaries[["centroids"]]@coords
  coords <- as.data.frame(raw)[, 1:2] %>%
    setNames(c("x_px", "y_px")) %>%
    transmute(x = x_px * scale_um, y = y_px * scale_um)

  coords$nCount   <- seu@meta.data[[paste0("nCount_",   prefix)]]
  coords$nFeature <- seu@meta.data[[paste0("nFeature_", prefix)]]
  coords$spot_id  <- rownames(coords)
  coords
}

# ---------------------------------------------------------------------------
# Extract all data frames
# ---------------------------------------------------------------------------
message("Extracting VisiumHD coordinates and QC metrics...")
vis_dfs       <- setNames(lapply(visiumhd_objs, make_df_visium), display_samples)
vis_dfs_clean <- lapply(vis_dfs, filter_noise)

message("Extracting MERSCOPE coordinates and QC metrics...")
mer_dfs <- setNames(
  lapply(merscope_objs, make_df_image, prefix = "Vizgen"),
  display_samples
)

message("Extracting Xenium coordinates and QC metrics...")
xen_dfs <- setNames(
  lapply(xenium_objs, make_df_image, prefix = "Xenium"),
  display_samples
)

# ---------------------------------------------------------------------------
# Compute shared plot window and per-platform colour limits
# ---------------------------------------------------------------------------
all_df <- bind_rows(
  bind_rows(vis_dfs_clean, .id = "id") %>% mutate(platform = "VisiumHD"),
  bind_rows(mer_dfs,       .id = "id") %>% mutate(platform = "MERSCOPE"),
  bind_rows(xen_dfs,       .id = "id") %>% mutate(platform = "Xenium")
)

# Square, padded window that encompasses all three platforms
global_sq  <- make_square(list(x = range(all_df$x), y = range(all_df$y)))
global_pad <- pad_window(global_sq, frac = 0.10)
vis_sq <- mer_sq <- xen_sq <- global_pad

# Per-platform colour limits: 1st–95th percentile per metric
qr_from <- function(lst, col) {
  v <- unlist(lapply(lst, `[[`, col))
  c(quantile(v, 0.01, na.rm = TRUE), quantile(v, 0.95, na.rm = TRUE))
}

quant_ranges <- list(
  VisiumHD = list(nCount   = qr_from(vis_dfs_clean, "nCount"),
                  nFeature = qr_from(vis_dfs_clean, "nFeature")),
  MERSCOPE = list(nCount   = qr_from(mer_dfs, "nCount"),
                  nFeature = qr_from(mer_dfs, "nFeature")),
  Xenium   = list(nCount   = qr_from(xen_dfs, "nCount"),
                  nFeature = qr_from(xen_dfs, "nFeature"))
)

# ---------------------------------------------------------------------------
# Core tile plot function
# ---------------------------------------------------------------------------
# Draws one spatial scatter tile with a 1 mm scale bar.
# limits clamps the viridis colour scale for cross-sample comparability.
plot_tile <- function(df, platform, win, value_col, limits, flip_y = FALSE) {
  sb     <- round(1000 / 8)   # ~125 units ≈ 1 mm at 8 µm/unit
  dx     <- diff(win$x); dy <- diff(win$y)
  xgap   <- 0.05 * dx;  ygap <- 0.05 * dy
  bar_y  <- win$y[1] + ygap
  bar_x2 <- win$x[2] - xgap; bar_x1 <- bar_x2 - sb
  tick_h <- 0.03 * dy
  y0 <- bar_y - tick_h / 2; y1 <- bar_y + tick_h / 2

  p <- ggplot(df, aes(x, y, colour = .data[[value_col]])) +
    geom_scattermore(pointsize = 2) +
    scale_colour_viridis_c(
      option    = "D",
      direction = 1,
      limits    = limits,
      oob       = scales::squish,
      guide     = "none"
    ) +
    coord_fixed(xlim = win$x, ylim = win$y, expand = FALSE) +
    theme_void(base_family = "Arial", base_size = 5) +
    theme(
      panel.background = element_rect(fill   = pal_bg[[platform]],
                                      colour = pal_bg[[platform]]),
      plot.margin      = margin(0, 0, 0, 0)
    )

  if (flip_y) p <- p + scale_y_reverse()

  p +
    annotate("segment", x = bar_x1, xend = bar_x2, y = bar_y,  yend = bar_y,  linewidth = 0.3) +
    annotate("segment", x = bar_x1, xend = bar_x1, y = y0,     yend = y1,     linewidth = 0.3) +
    annotate("segment", x = bar_x2, xend = bar_x2, y = y0,     yend = y1,     linewidth = 0.3)
}

# ---------------------------------------------------------------------------
# Colourbar legend builder
# ---------------------------------------------------------------------------
make_legend <- function(vals, title_txt) {
  lims  <- quantile(vals, c(0.05, 0.95), na.rm = TRUE)
  ticks <- ceiling(lims)
  p <- ggplot() +
    geom_point(aes(1, 1, colour = mean(lims)), size = 0) +
    scale_colour_viridis_c(
      option = "D", direction = 1,
      limits = lims, breaks = ticks, labels = ticks,
      guide  = guide_colorbar(
        title          = title_txt,
        title.position = "top",
        barheight      = unit(20, "mm"),
        barwidth       = unit(2,  "mm"),
        ticks          = TRUE
      )
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title    = element_text(size = 6),
      legend.text     = element_text(size = 5),
      legend.margin   = margin(0, 0, 0, 0)
    )
  suppressWarnings(cowplot::get_legend(p))
}

# ---------------------------------------------------------------------------
# Panel builder: returns a gtable for one metric (nCount or nFeature)
# ---------------------------------------------------------------------------
# KO168 is rotated 90° and re-centred to match the other samples' orientation.
build_montage_panel <- function(metric = c("nCount", "nFeature")) {
  metric    <- match.arg(metric)
  leg_title <- if (metric == "nCount") "Transcripts" else "Genes"

  # Build tile grobs for each platform
  vis_grobs <- lapply(display_samples, function(s) {
    df <- vis_dfs_clean[[s]]
    if (s == "KO168") df <- rotate90_df(df) %>% recentre_to_window(vis_sq)
    ggplotGrob(plot_tile(df, "VisiumHD", vis_sq, metric,
                         limits = quant_ranges[["VisiumHD"]][[metric]]))
  })

  mer_grobs <- lapply(display_samples, function(s) {
    df <- mer_dfs[[s]]
    if (s == "KO168") df <- rotate90_df(df) %>% recentre_to_window(mer_sq)
    ggplotGrob(plot_tile(df, "MERSCOPE", mer_sq, metric,
                         limits = quant_ranges[["MERSCOPE"]][[metric]],
                         flip_y = TRUE))
  })

  xen_grobs <- lapply(display_samples, function(s) {
    df <- xen_dfs[[s]]
    if (s == "KO168") df <- rotate90_df(df) %>% recentre_to_window(xen_sq)
    ggplotGrob(plot_tile(df, "Xenium", xen_sq, metric,
                         limits = quant_ranges[["Xenium"]][[metric]],
                         flip_y = TRUE))
  })

  # Label and spacer grobs
  blank_grob      <- ggplotGrob(ggplot() + theme_void())
  col_label_grobs <- lapply(display_samples, function(s) ggplotGrob(col_label(s)))
  row_label_grobs <- lapply(c("VisiumHD", "MERSCOPE", "Xenium"),
                            function(lbl) ggplotGrob(row_label(lbl)))

  # Per-platform colourbar legends
  leg_vis <- make_legend(unlist(lapply(vis_dfs_clean, `[[`, metric)), leg_title)
  leg_mer <- make_legend(unlist(lapply(mer_dfs,       `[[`, metric)), leg_title)
  leg_xen <- make_legend(unlist(lapply(xen_dfs,       `[[`, metric)), leg_title)

  spacer_row <- replicate(5, nullGrob(), simplify = FALSE)

  # Assemble: row of column labels, spacer, then one platform row per platform
  all_grobs <- c(
    list(blank_grob), col_label_grobs,                            #  1– 5
    spacer_row,                                                    #  6–10
    list(row_label_grobs[[1]]), vis_grobs, list(leg_vis),         # 11–16
    spacer_row,                                                    # 17–21
    list(row_label_grobs[[2]]), mer_grobs, list(leg_mer),         # 22–27
    spacer_row,                                                    # 28–32
    list(row_label_grobs[[3]]), xen_grobs, list(leg_xen)          # 33–38
  )

  layout_mat <- rbind(
    c( 1,  2,  3,  4,  5, NA),
    c( 6,  7,  8,  9, 10, NA),
    c(11, 12, 13, 14, 15, 16),
    c(17, 18, 19, 20, 21, NA),
    c(22, 23, 24, 25, 26, 27),
    c(28, 29, 30, 31, 32, NA),
    c(33, 34, 35, 36, 37, 38)
  )

  widths_mm  <- c(12, 30, 30, 30, 30, 10)
  heights_mm <- c( 6,  2, 30,  2, 30,  2, 30)

  arrangeGrob(
    grobs         = all_grobs,
    layout_matrix = layout_mat,
    widths        = unit(widths_mm,  "mm"),
    heights       = unit(heights_mm, "mm"),
    padding       = unit(0, "mm")
  )
}

# Wrap a gtable in a ggdraw canvas for ggsave compatibility
save_montage_panel <- function(grob, filename, width_mm = 144, height_mm = 108) {
  p <- cowplot::ggdraw() +
    cowplot::draw_grob(grob, x = 0, y = 0, width = 1, height = 1) +
    theme(plot.margin = grid::unit(c(0, 0, 0, 0), "pt"))
  ggsave(filename, p,
         device = cairo_pdf,
         width  = width_mm,
         height = height_mm,
         units  = "mm",
         bg     = "white")
  message("Saved: ", filename)
}

message("Building nCount spatial montage...")
p_ncount <- build_montage_panel("nCount")
save_montage_panel(p_ncount, file.path(opt$out_dir, "spatial_ncount.pdf"))

message("Building nFeature spatial montage...")
p_nfeature <- build_montage_panel("nFeature")
save_montage_panel(p_nfeature, file.path(opt$out_dir, "spatial_nfeature.pdf"))

# ===========================================================================
# QC METRIC BOXPLOTS (spatial platforms + FLEX comparison)
# ===========================================================================

# ---------------------------------------------------------------------------
# Helper: compute 4-tick y-axis breaks from the data maximum
# ---------------------------------------------------------------------------
four_ticks <- function(x) {
  step  <- ceiling((x / 3) / 10) * 10
  upper <- step * 3
  list(breaks = seq(0, upper, by = step), upper = upper)
}

# ---------------------------------------------------------------------------
# Spatial platform QC panels (VisiumHD / MERSCOPE / Xenium)
# ---------------------------------------------------------------------------
# Filtered to 8µm bins and the 90-gene common subset, then summarised to one
# median value per sample × platform before plotting.

subtitle_spatial <- "8 μm bins · 90 common genes"

message("Loading metadata_combined...")
metadata_combined <- readRDS(meta_combined_path)

make_spatial_qc_panel <- function(yvar, ylab) {

  df <- metadata_combined %>%
    dplyr::filter(Subset == "90", bin_size == "8um") %>%
    dplyr::group_by(Sample, platform) %>%
    dplyr::summarise(val = median(.data[[yvar]]), .groups = "drop") %>%
    # fix order for consistent jitter position across runs
    dplyr::arrange(platform, Sample) %>%
    dplyr::mutate(
      platform = factor(platform, levels = c("VisiumHD", "MERSCOPE", "Xenium")),
      subtitle = subtitle_spatial
    )

  ax <- four_ticks(max(df$val, na.rm = TRUE))

  ggplot(df, aes(x = platform, y = val)) +
    geom_boxplot(
      aes(fill = platform),
      outlier.shape = NA,
      width         = 0.7
    ) +
    geom_point(
      aes(fill = platform, group = Sample),
      position = position_jitter(width = 0.10, height = 0),
      shape  = 21,
      colour = "grey20",
      stroke = 0.35,
      size   = 2.0,
      alpha  = 0.80
    ) +
    scale_fill_platform() +
    scale_y_continuous(
      limits = c(0, ax$upper),
      breaks = ax$breaks,
      labels = label_number(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = "Platform", y = ylab) +
    facet_wrap(~ subtitle, ncol = 1) +
    theme_sb() +
    theme(
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
}

message("Building spatial QC panels...")

p_spatial_counts <- make_spatial_qc_panel("nCount",   "Median counts/bin")
p_spatial_genes  <- make_spatial_qc_panel("nFeature", "Median genes/bin")

ggsave(file.path(opt$out_dir, "qc_counts_spatial.pdf"),
       p_spatial_counts,
       width = dims$half_w, height = dims$half_w * 1.4,
       units = "mm", device = cairo_pdf, bg = "white")
message("Saved: qc_counts_spatial.pdf")

ggsave(file.path(opt$out_dir, "qc_genes_spatial.pdf"),
       p_spatial_genes,
       width = dims$half_w, height = dims$half_w * 1.4,
       units = "mm", device = cairo_pdf, bg = "white")
message("Saved: qc_genes_spatial.pdf")

# ---------------------------------------------------------------------------
# FLEX / scRNA-seq intersect-gene QC panels
# ---------------------------------------------------------------------------
# Skipped gracefully if metadata_flex_scrna.rds was not produced (--sc_rds
# was omitted from qc_metrics.R).
if (!file.exists(meta_flex_path)) {
  message("Skipping FLEX panels: ", meta_flex_path, " not found.")
  message("  Re-run qc_metrics.R with --sc_rds to generate this file.")
  message("Done. Outputs written to: ", opt$out_dir)
  quit(save = "no", status = 0)
}

message("Loading metadata_flex_scrna and genes_intersect_flex...")
metadata_flex_scrna <- readRDS(meta_flex_path)
genes_intersect     <- readRDS(genes_intersect_path)
subtitle_flex       <- paste0("Intersect genes (n = ", length(genes_intersect), ")")

# Colour palette for the FLEX comparison platforms (not in scale_fill_platform)
pal_flex <- c(
  scRNAseq = "#333333",
  FLEX     = "#888888",
  VisiumHD = pal_muted[["VisiumHD"]]
)

# Build one intersect-gene QC panel, optionally subsetting to specific platforms.
# platforms_keep controls both the data filter and the fill scale.
make_flex_qc_panel <- function(yvar, ylab,
                               platforms_keep = c("scRNAseq", "FLEX", "VisiumHD")) {

  df <- metadata_flex_scrna %>%
    dplyr::filter(platform %in% platforms_keep) %>%
    dplyr::group_by(Sample, platform) %>%
    dplyr::summarise(val = median(.data[[yvar]]), .groups = "drop") %>%
    dplyr::mutate(
      platform = factor(platform, levels = platforms_keep),
      subtitle = subtitle_flex
    )

  ax <- four_ticks(max(df$val, na.rm = TRUE))

  ggplot(df, aes(x = platform, y = val)) +
    geom_boxplot(
      aes(fill = platform),
      outlier.shape = NA,
      width         = 0.7
    ) +
    geom_point(
      aes(fill = platform, group = Sample),
      position = position_jitter(width = 0.08, height = 0),
      shape  = 21,
      colour = "grey20",
      stroke = 0.35,
      size   = 2.0,
      alpha  = 0.80
    ) +
    scale_fill_manual(values = pal_flex[platforms_keep]) +
    scale_y_continuous(
      limits = c(0, ax$upper),
      breaks = ax$breaks,
      labels = label_number(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = "Platform", y = ylab) +
    facet_wrap(~ subtitle, ncol = 1) +
    theme_sb() +
    theme(
      legend.position  = "none",
      axis.text.x      = element_text(angle = 45, hjust = 1),
      axis.line        = element_line(colour = "black"),
      strip.background = element_rect(colour = "black", fill = "white"),
      strip.text       = element_text(face = "italic"),
      panel.background = element_blank(),
      panel.border     = element_blank()
    )
}

message("Building FLEX QC panels (all platforms)...")

p_flex_counts_all <- make_flex_qc_panel(
  "nCount_intersect",
  "Median counts per bin/cell (intersect genes)"
)
p_flex_genes_all  <- make_flex_qc_panel(
  "nFeature_intersect",
  "Median genes per bin/cell (intersect genes)"
)

ggsave(file.path(opt$out_dir, "qc_counts_flex_all.pdf"),
       p_flex_counts_all,
       width = dims$half_w, height = dims$half_w * 1.4,
       units = "mm", device = cairo_pdf, bg = "white")
message("Saved: qc_counts_flex_all.pdf")

ggsave(file.path(opt$out_dir, "qc_genes_flex_all.pdf"),
       p_flex_genes_all,
       width = dims$half_w, height = dims$half_w * 1.4,
       units = "mm", device = cairo_pdf, bg = "white")
message("Saved: qc_genes_flex_all.pdf")

message("Building FLEX QC panels (FLEX + VisiumHD only)...")

p_flex_counts_fxvs <- make_flex_qc_panel(
  "nCount_intersect",
  "Median counts per bin/cell (intersect genes)",
  platforms_keep = c("FLEX", "VisiumHD")
)
p_flex_genes_fxvs  <- make_flex_qc_panel(
  "nFeature_intersect",
  "Median genes per bin/cell (intersect genes)",
  platforms_keep = c("FLEX", "VisiumHD")
)

ggsave(file.path(opt$out_dir, "qc_counts_flex_visiumhd.pdf"),
       p_flex_counts_fxvs,
       width = dims$half_w * 0.6, height = dims$half_w * 1.4,
       units = "mm", device = cairo_pdf, bg = "white")
message("Saved: qc_counts_flex_visiumhd.pdf")

ggsave(file.path(opt$out_dir, "qc_genes_flex_visiumhd.pdf"),
       p_flex_genes_fxvs,
       width = dims$half_w * 0.6, height = dims$half_w * 1.4,
       units = "mm", device = cairo_pdf, bg = "white")
message("Saved: qc_genes_flex_visiumhd.pdf")

message("Done. All panels written to: ", opt$out_dir)
