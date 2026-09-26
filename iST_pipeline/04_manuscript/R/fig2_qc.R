# Purpose:  Figure 2 QC panels — platform QC panels.
#           Spatial scatter montage coloured by nCount and nFeature across
#           VisiumHD, MERSCOPE, and Xenium (3 platforms × 4 samples, 8µm bins).
#           Per-platform dot plots (one shaped point per animal + median
#           line) of median counts/bin and genes/bin (8µm bins,
#           90-gene common subset) for the matched animals, for full tissue and
#           for the fixed-size ROI (shared y-axis). Per-cell intersect-gene boxplots comparing
#           VisiumHD, FLEX snRNA-seq, and scRNA-seq (optional; requires
#           metadata_flex_scrna.rds from qc_metrics.R --sc_rds).
#           Each panel saved as a separate PDF.
#           Adapted from montage_v2.R, fig2_qcmetrics.R, flex_visiumhd_qc_v2.R.
#           Montage coordinates: MERSCOPE/Xenium bins are read from the
#           STalign-aligned objects, whose centroids are already in the
#           matched VisiumHD sample's full-resolution pixel space. All three
#           platforms are converted to µm with that VisiumHD sample's
#           microns_per_pixel (config stalign$microns_per_pixel).
# Inputs:   config/config.yaml
#           results/01_preprocessing/merscope_8um_aligned/{sample}_8um_aligned.rds
#           results/01_preprocessing/xenium_8um_aligned/{sample}_8um_aligned.rds
#           results/03_benchmarking/qc_metrics/metadata_combined.rds
#           results/03_benchmarking/qc_metrics/metadata_roi{roi_label}.rds
#           results/03_benchmarking/qc_metrics/metadata_flex_scrna.rds  (optional)
#           results/03_benchmarking/qc_metrics/genes_intersect_flex.rds (optional)
# Outputs:  figures/fig2/spatial_ncount.pdf
#           figures/fig2/spatial_nfeature.pdf
#           figures/fig2/qc_counts_spatial.pdf
#           figures/fig2/qc_genes_spatial.pdf
#           figures/fig2/qc_counts_spatial_roi{roi_label}.pdf
#           figures/fig2/qc_genes_spatial_roi{roi_label}.pdf
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
              help    = "Output directory for panel PDFs [default: %default]"),
  make_option(c("--samples"),   type = "character",
              default = "WT709,WT713,KO167,KO168",
              help    = "Comma-separated montage samples (display names) [default: %default]"),
  make_option(c("--roi_label"), type = "character", default = "2mm",
              help    = "ROI size label of metadata_roi{label}.rds and ROI panel filenames [default: %default]"),
  make_option(c("--montage_only"), action = "store_true", default = FALSE,
              help    = "Only build the spatial montages (skip boxplots); for checking")
)
opt <- parse_args(OptionParser(option_list = option_list))

meta_combined_path   <- file.path(opt$input_dir, "metadata_combined.rds")
meta_roi_path        <- file.path(opt$input_dir, paste0("metadata_roi", opt$roi_label, ".rds"))
meta_flex_path       <- file.path(opt$input_dir, "metadata_flex_scrna.rds")
genes_intersect_path <- file.path(opt$input_dir, "genes_intersect_flex.rds")

if (!opt$montage_only) {
  for (f in c(meta_combined_path, meta_roi_path)) {
    if (!file.exists(f)) stop("Not found: ", f)
  }
}

cfg <- yaml::read_yaml(opt$config)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# SPATIAL SCATTER MONTAGE (nCount / nFeature per bin)
# ===========================================================================

# ---------------------------------------------------------------------------
# Sample mapping: display name -> config sample key
# ---------------------------------------------------------------------------
# Display names are upper-case animal IDs (e.g. "KO168"); the matching
# VisiumHD / MERSCOPE / Xenium sample keys come from the same config$stalign
# mapping that align_binned.R used, so montage and alignment cannot diverge.
display_samples <- trimws(strsplit(opt$samples, ",")[[1]])
animal_ids      <- setNames(tolower(display_samples), display_samples)

lookup_keys <- function(mapping) {
  keys <- vapply(animal_ids, function(a) {
    k <- mapping[[a]]
    if (is.null(k)) stop("No config$stalign entry for animal '", a, "'")
    k
  }, character(1))
  setNames(keys, display_samples)
}

vis_sample_map <- lookup_keys(cfg$stalign$visiumhd)
mer_sample_map <- lookup_keys(cfg$stalign$matched_samples$merscope)
xen_sample_map <- lookup_keys(cfg$stalign$matched_samples$xenium)

# µm per full-resolution pixel for each display sample's matched VisiumHD
# sample — shared by all three platforms, since aligned MERSCOPE/Xenium
# centroids live in that VisiumHD sample's pixel space.
um_per_px <- setNames(
  vapply(vis_sample_map, function(k) as.numeric(cfg$stalign$microns_per_pixel[[k]]),
         numeric(1)),
  display_samples
)

# ---------------------------------------------------------------------------
# Load Seurat objects for spatial montage (read-only; never modified/saved)
# ---------------------------------------------------------------------------
# Read one RDS per display sample from a directory + filename template
load_objs <- function(sample_map, dir, suffix = "") {
  setNames(
    lapply(sample_map, function(key) {
      path <- file.path(dir, paste0(key, suffix))
      message("  ", key, ": ", path)
      readRDS(path)
    }),
    names(sample_map)
  )
}

message("Loading VisiumHD samples...")
visiumhd_objs <- load_objs(
  setNames(unlist(cfg$visiumhd$samples[vis_sample_map]), display_samples),
  cfg$visiumhd$data_dir
)

message("Loading MERSCOPE 8um aligned samples...")
merscope_objs <- load_objs(
  mer_sample_map,
  file.path(cfg$output_dir, "01_preprocessing", "merscope_8um_aligned"),
  "_8um_aligned.rds"
)

message("Loading Xenium 8um aligned samples...")
xenium_objs <- load_objs(
  xen_sample_map,
  file.path(cfg$output_dir, "01_preprocessing", "xenium_8um_aligned"),
  "_8um_aligned.rds"
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

# Remove DBSCAN noise and return only the dominant cluster (VisiumHD only;
# MERSCOPE/Xenium were already DBSCAN-filtered by filter_binned.R).
# eps is in µm: 80 µm matches the previous eps = 8 low-res VisiumHD pixels
# (8 px / 0.0272 lowres_scalef * 0.274 µm/px ≈ 80 µm).
filter_noise <- function(df, xcol = "x", ycol = "y", eps = 80, minPts = 10) {
  cl   <- dbscan::dbscan(as.matrix(df[, c(xcol, ycol)]), eps = eps, minPts = minPts)$cluster
  keep <- as.integer(names(which.max(table(cl[cl > 0]))))
  df[cl == keep, , drop = FALSE]
}

# Centre of a data frame's x/y bounding box
extent_centre <- function(df) c(x = mean(range(df$x)), y = mean(range(df$y)))

# Display-only 90° rotation, (x, y) -> (y, -x) about `centre`. Used for KO168
# so its tissue matches the other samples' orientation. Applied identically
# to all three platforms about the matched VisiumHD tissue centre, so their
# relative alignment is preserved. Operates on the plotting data frame only.
rotate90_df <- function(df, centre) {
  dx <- df$x - centre[["x"]]
  dy <- df$y - centre[["y"]]
  df %>% mutate(x = centre[["x"]] + dy,
                y = centre[["y"]] - dx)
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

# VisiumHD: native full-resolution pixel coordinates of the 8 µm bins,
# converted to µm. Counts are matched to coordinates by barcode.
# Filters bins with fewer than 5 counts to remove empty spots.
make_df_visium <- function(seu, um_px, bin_suffix = "008um") {
  assay_nm <- paste0("Spatial.", bin_suffix)
  img_nm   <- paste0("slice1.", bin_suffix)
  tc       <- GetTissueCoordinates(seu, image = img_nm)

  meta <- seu@meta.data[tc$cell, c(paste0("nCount_",   assay_nm),
                                   paste0("nFeature_", assay_nm))]
  if (anyNA(meta[[1]])) stop("VisiumHD barcodes missing from meta.data for ", img_nm)

  data.frame(x        = tc$x * um_px,
             y        = tc$y * um_px,
             nCount   = meta[[1]],
             nFeature = meta[[2]],
             spot_id  = tc$cell) %>%
    filter(nCount >= 5)
}

# MERSCOPE / Xenium: STalign-aligned centroids (VisiumHD full-res pixels),
# converted to µm. Counts are matched to centroids by bin ID.
make_df_image <- function(seu, prefix, um_px) {
  raw <- seu@images[[1]]@boundaries[["centroids"]]@coords
  ids <- rownames(raw)
  if (is.null(ids)) stop("Aligned centroids have no bin IDs; cannot match counts")

  meta <- seu@meta.data[ids, c(paste0("nCount_",   prefix),
                               paste0("nFeature_", prefix))]
  if (anyNA(meta[[1]])) stop("Centroid bin IDs missing from meta.data (", prefix, ")")

  data.frame(x        = raw[, "x"] * um_px,
             y        = raw[, "y"] * um_px,
             nCount   = meta[[1]],
             nFeature = meta[[2]],
             spot_id  = ids)
}

# ---------------------------------------------------------------------------
# Extract all data frames
# ---------------------------------------------------------------------------
message("Extracting VisiumHD coordinates and QC metrics...")
vis_dfs <- setNames(
  lapply(display_samples, function(s) make_df_visium(visiumhd_objs[[s]], um_per_px[[s]])),
  display_samples
)

message("DBSCAN noise filtering on VisiumHD (eps = 80 µm)...")
vis_dfs_clean <- lapply(vis_dfs, filter_noise)
for (s in display_samples) {
  message(sprintf("  %s: kept %d / %d bins (%.2f%%)", s, nrow(vis_dfs_clean[[s]]),
                  nrow(vis_dfs[[s]]), 100 * nrow(vis_dfs_clean[[s]]) / nrow(vis_dfs[[s]])))
}

message("Extracting MERSCOPE coordinates and QC metrics...")
mer_dfs <- setNames(
  lapply(display_samples, function(s)
    make_df_image(merscope_objs[[s]], "Vizgen", um_per_px[[s]])),
  display_samples
)

message("Extracting Xenium coordinates and QC metrics...")
xen_dfs <- setNames(
  lapply(display_samples, function(s)
    make_df_image(xenium_objs[[s]], "Xenium", um_per_px[[s]])),
  display_samples
)

# ---------------------------------------------------------------------------
# Display-only KO168 rotation
# ---------------------------------------------------------------------------
# All three platforms are rotated about the same centre (the matched VisiumHD
# tissue extent), so they stay registered to each other. Seurat objects and
# aligned coordinates are untouched; only these plotting data frames change.
rotate_samples <- intersect("KO168", display_samples)

for (s in rotate_samples) {
  rot_centre <- extent_centre(vis_dfs_clean[[s]])
  message("Rotating ", s, " 90° for display about VisiumHD centre (",
          round(rot_centre[["x"]]), ", ", round(rot_centre[["y"]]), ") µm")
  vis_dfs_clean[[s]] <- rotate90_df(vis_dfs_clean[[s]], rot_centre)
  mer_dfs[[s]]       <- rotate90_df(mer_dfs[[s]],       rot_centre)
  xen_dfs[[s]]       <- rotate90_df(xen_dfs[[s]],       rot_centre)
}

# ---------------------------------------------------------------------------
# Compute shared plot window and per-platform colour limits
# ---------------------------------------------------------------------------
# Window is computed after rotation so rotated tissue is never clipped.
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
# Draws one spatial scatter tile with a 1 mm scale bar (coordinates in µm).
# limits clamps the viridis colour scale for cross-sample comparability.
# All platforms share one coordinate space, so no per-platform axis flip.
plot_tile <- function(df, platform, win, value_col, limits) {
  sb     <- 1000              # 1 mm in µm
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
# KO168 data frames were already rotated for display (see above).
build_montage_panel <- function(metric = c("nCount", "nFeature")) {
  metric    <- match.arg(metric)
  leg_title <- if (metric == "nCount") "Transcripts" else "Genes"

  # Build tile grobs for each platform
  platform_tiles <- function(dfs, platform, win) {
    lapply(display_samples, function(s) {
      ggplotGrob(plot_tile(dfs[[s]], platform, win, metric,
                           limits = quant_ranges[[platform]][[metric]]))
    })
  }
  vis_grobs <- platform_tiles(vis_dfs_clean, "VisiumHD", vis_sq)
  mer_grobs <- platform_tiles(mer_dfs,       "MERSCOPE", mer_sq)
  xen_grobs <- platform_tiles(xen_dfs,       "Xenium",   xen_sq)

  # Label and spacer grobs
  blank_grob      <- ggplotGrob(ggplot() + theme_void())
  col_label_grobs <- lapply(display_samples, function(s) ggplotGrob(col_label(s)))
  row_label_grobs <- lapply(c("VisiumHD", "MERSCOPE", "Xenium"),
                            function(lbl) ggplotGrob(row_label(lbl)))

  # Per-platform colourbar legends
  leg_vis <- make_legend(unlist(lapply(vis_dfs_clean, `[[`, metric)), leg_title)
  leg_mer <- make_legend(unlist(lapply(mer_dfs,       `[[`, metric)), leg_title)
  leg_xen <- make_legend(unlist(lapply(xen_dfs,       `[[`, metric)), leg_title)

  # Assemble a grid with one column per sample: header row of column labels,
  # then one row per platform (row label, tiles, legend), separated by spacers.
  n_col      <- length(display_samples) + 1           # row-label column + tiles
  spacer_row <- replicate(n_col, nullGrob(), simplify = FALSE)
  platform_rows <- list(
    list(row_label_grobs[[1]], vis_grobs, leg_vis),
    list(row_label_grobs[[2]], mer_grobs, leg_mer),
    list(row_label_grobs[[3]], xen_grobs, leg_xen)
  )

  all_grobs  <- c(list(blank_grob), col_label_grobs)
  layout_mat <- rbind(c(seq_len(n_col), NA))
  for (pr in platform_rows) {
    # Spacer row (no legend cell)
    idx        <- length(all_grobs) + seq_len(n_col)
    all_grobs  <- c(all_grobs, spacer_row)
    layout_mat <- rbind(layout_mat, c(idx, NA))
    # Platform row: label, one tile per sample, legend
    row_grobs  <- c(list(pr[[1]]), pr[[2]], list(pr[[3]]))
    idx        <- length(all_grobs) + seq_along(row_grobs)
    all_grobs  <- c(all_grobs, row_grobs)
    layout_mat <- rbind(layout_mat, idx)
  }

  widths_mm  <- c(12, rep(30, length(display_samples)), 10)
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
# Width defaults to the 4-sample layout (12 + 4 × 30 + 10 mm, plus margin)
save_montage_panel <- function(grob, filename,
                               width_mm  = 14 + 30 * length(display_samples) + 10,
                               height_mm = 108) {
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

if (opt$montage_only) {
  message("--montage_only: skipping QC boxplots. Outputs in: ", opt$out_dir)
  quit(save = "no", status = 0)
}

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
# Two parallel versions of the same comparison: full tissue
# (metadata_combined.rds) and the fixed-size ROI (metadata_roi{label}.rds).
# Both are restricted to the matched animals, filtered to 8µm bins and the
# 90-gene common subset, then summarised to one median value per
# sample × platform. Full and ROI panels share a y-axis per metric.

message("Loading metadata_combined and ROI metadata...")
metadata_combined <- readRDS(meta_combined_path)
metadata_roi      <- readRDS(meta_roi_path)

# Keep only the matched animals (config$stalign) on every platform, so all
# three platforms show the same animals. Stops if any platform is incomplete.
matched_animals <- names(cfg$stalign$visiumhd)
matched_keys <- list(
  VisiumHD = unlist(cfg$stalign$visiumhd[matched_animals]),
  MERSCOPE = unlist(cfg$stalign$matched_samples$merscope[matched_animals]),
  Xenium   = unlist(cfg$stalign$matched_samples$xenium[matched_animals])
)

filter_matched <- function(meta) {
  out <- meta %>%
    dplyr::filter(
      (platform == "VisiumHD" & Sample %in% matched_keys$VisiumHD) |
      (platform == "MERSCOPE" & Sample %in% matched_keys$MERSCOPE) |
      (platform == "Xenium"   & Sample %in% matched_keys$Xenium)
    )
  n_per_platform <- tapply(out$Sample, as.character(out$platform), dplyr::n_distinct)
  if (length(n_per_platform) != 3 || any(n_per_platform != length(matched_animals))) {
    stop("Matched-animal filter did not give ", length(matched_animals),
         " samples per platform: ",
         paste(names(n_per_platform), n_per_platform, sep = "=", collapse = ", "))
  }
  out
}

# Platform sample key (e.g. "wt709_batch13") -> display animal ID ("WT709")
key_to_animal <- setNames(rep(toupper(matched_animals), length(matched_keys)),
                          unlist(lapply(matched_keys, unname)))

# One point shape per animal, shared across all three platforms
animal_shapes <- c(WT709 = 21, WT713 = 22, KO167 = 24, KO168 = 23)

# One median per sample × platform (8µm bins, 90-gene subset)
summarise_qc <- function(meta, yvar) {
  meta %>%
    dplyr::filter(Subset == "90", bin_size == "8um") %>%
    dplyr::group_by(Sample, platform) %>%
    dplyr::summarise(val = median(.data[[yvar]]), .groups = "drop") %>%
    # fix order for consistent jitter position across runs
    dplyr::arrange(platform, Sample) %>%
    dplyr::mutate(
      platform = factor(platform, levels = c("VisiumHD", "MERSCOPE", "Xenium")),
      animal   = factor(unname(key_to_animal[Sample]), levels = names(animal_shapes))
    )
}

qc_sets <- list(
  full = filter_matched(metadata_combined),
  roi  = filter_matched(metadata_roi)
)

subtitles <- list(
  full = "8 μm bins · 90 common genes",
  roi  = paste0("8 μm bins · 90 common genes · ", opt$roi_label, " ROI")
)

make_spatial_qc_panel <- function(df, ylab, subtitle, ax) {

  df <- df %>% dplyr::mutate(subtitle = subtitle)

  ggplot(df, aes(x = platform, y = val)) +
    # Median across the four animals: horizontal black line per platform
    stat_summary(
      fun       = median,
      fun.min   = median,
      fun.max   = median,
      geom      = "errorbar",
      width     = 0.35,
      linewidth = 0.6,
      colour    = "black"
    ) +
    # One point per animal: platform colour fill, animal-specific shape.
    # Each animal gets a fixed horizontal slot (same order on every platform)
    # so animals with tied medians never overlap.
    geom_point(
      aes(fill = platform, shape = animal, group = animal),
      position = position_dodge(width = 0.30),
      colour = "grey20",
      stroke = 0.35,
      size   = 2.0,
      alpha  = 0.80
    ) +
    scale_fill_platform(guide = "none") +
    scale_shape_manual(
      values = animal_shapes,
      name   = "Animal",
      guide  = guide_legend(override.aes = list(fill = "grey70", alpha = 1))
    ) +
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
      legend.position  = "right",
      axis.line        = element_line(colour = "black"),
      strip.background = element_rect(colour = "black", fill = "white"),
      strip.text       = element_text(face = "italic"),
      panel.background = element_blank(),
      panel.border     = element_blank(),
      plot.background  = element_blank(),
      panel.spacing    = unit(2, "mm")
    )
}

message("Building spatial QC panels (full tissue + ", opt$roi_label, " ROI)...")

qc_metrics_plot <- list(
  counts = list(yvar = "nCount",   ylab = "Median counts/bin"),
  genes  = list(yvar = "nFeature", ylab = "Median genes/bin")
)

for (m in names(qc_metrics_plot)) {
  spec <- qc_metrics_plot[[m]]
  summ <- lapply(qc_sets, summarise_qc, yvar = spec$yvar)

  # Shared y-axis: computed from the maximum across full and ROI versions
  ax <- four_ticks(max(unlist(lapply(summ, `[[`, "val")), na.rm = TRUE))

  for (set in names(summ)) {
    fname <- if (set == "full") {
      paste0("qc_", m, "_spatial.pdf")
    } else {
      paste0("qc_", m, "_spatial_roi", opt$roi_label, ".pdf")
    }
    p <- make_spatial_qc_panel(summ[[set]], spec$ylab, subtitles[[set]], ax)
    ggsave(file.path(opt$out_dir, fname), p,
           width = dims$half_w, height = dims$half_w * 1.4,
           units = "mm", device = cairo_pdf, bg = "white")
    message("Saved: ", fname, "  (", paste(table(summ[[set]]$platform), collapse = "/"),
            " samples per VisiumHD/MERSCOPE/Xenium)")
  }
}

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
