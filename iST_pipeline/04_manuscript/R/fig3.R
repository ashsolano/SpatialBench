# Purpose:  Figure 3 — Segmentation quality comparison across methods and platforms.
#           Produces eight panels from pre-computed segmentation_quality.R outputs:
#           per-platform metric boxplots (cell count, area, transcripts per cell,
#           transcripts assigned), UMAP embeddings coloured by cell type, cell type
#           composition bar charts, MECR boxplots, negative marker purity dotplot,
#           and negative marker purity heatmaps per cell type.
# Inputs:   results/03_benchmarking/segmentation_quality/
#             metrics_long.rds, umap_coords.rds, cell_type_counts.rds,
#             mecr_table.rds, purity_summary.rds, purity_ct.rds
# Outputs:  figures/fig3/
#             metrics_merscope.pdf, metrics_xenium.pdf
#             umap_merscope.pdf,    umap_xenium.pdf
#             cell_counts.pdf,      mecr.pdf
#             purity_dotplot.pdf,   purity_heatmap.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(optparse)
  library(RColorBrewer)
})

source("04_manuscript/R/utils/theme.R")
source("04_manuscript/R/utils/palettes.R")

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--input_dir"), type = "character",
              default = "results/03_benchmarking/segmentation_quality",
              help    = "Directory containing segmentation_quality.R outputs [default: %default]"),
  make_option(c("--out_dir"),   type = "character",
              default = "figures/fig3",
              help    = "Output directory for panel PDFs [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# Check required inputs
required_files <- c("metrics_long.rds", "umap_coords.rds", "cell_type_counts.rds",
                    "mecr_table.rds", "purity_summary.rds", "purity_ct.rds")
for (f in required_files) {
  path <- file.path(opt$input_dir, f)
  if (!file.exists(path)) stop("Required input not found: ", path)
}

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
message("Loading inputs from: ", opt$input_dir)

metrics_long     <- readRDS(file.path(opt$input_dir, "metrics_long.rds"))
umap_coords      <- readRDS(file.path(opt$input_dir, "umap_coords.rds"))
cell_type_counts <- readRDS(file.path(opt$input_dir, "cell_type_counts.rds"))
mecr_table       <- readRDS(file.path(opt$input_dir, "mecr_table.rds"))
purity_summary   <- readRDS(file.path(opt$input_dir, "purity_summary.rds"))
purity_ct        <- readRDS(file.path(opt$input_dir, "purity_ct.rds"))

# ---------------------------------------------------------------------------
# Shared aesthetics
# ---------------------------------------------------------------------------
scale_fill_seg   <- function(...) scale_fill_manual(values = pal_segmentation, ...)
scale_colour_seg <- function(...) scale_colour_manual(values = pal_segmentation, ...)

# ---------------------------------------------------------------------------
# Panels 1–2: metric boxplots (one per platform)
# ---------------------------------------------------------------------------
# Build shared y-axis limits per metric so both platforms are comparable
metric_limits <- metrics_long %>%
  dplyr::group_by(metric) %>%
  dplyr::summarise(lo = min(value, na.rm = TRUE),
                   hi = max(value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    pad = (hi - lo) * 0.05,
    lo  = dplyr::if_else(metric == "Transcripts assigned (%)", 0,   lo - pad),
    hi  = dplyr::if_else(metric == "Transcripts assigned (%)", 100, hi + pad)
  )

metric_lim <- function(m) {
  r <- dplyr::filter(metric_limits, metric == m)
  c(r$lo, r$hi)
}

plot_metric_panel <- function(df_platform) {
  metric_name <- as.character(unique(df_platform$metric))
  lim         <- metric_lim(metric_name)
  brk         <- if (metric_name == "Transcripts assigned (%)") {
    seq(0, 100, 25)
  } else {
    pretty(lim, n = 4)
  }

  lbl <- if (metric_name == "Total cell count") {
    scales::label_number(scale_cut = scales::cut_short_scale())
  } else if (metric_name == "Transcripts assigned (%)") {
    scales::label_number(accuracy = 1)
  } else {
    scales::label_number(big.mark = ",")
  }

  ggplot(df_platform, aes(x = method, y = value)) +
    geom_boxplot(aes(fill = method), width = 0.7, outlier.shape = NA,
                 colour = "black", alpha = 0.40, coef = 2) +
    geom_jitter(aes(fill = method), width = 0.15, shape = 21,
                colour = "black", size = dims$point, alpha = 0.95) +
    scale_fill_seg() +
    scale_y_continuous(limits = lim, breaks = brk, oob = scales::oob_keep,
                       labels = lbl) +
    labs(x = NULL, y = metric_name) +
    theme_sb() +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.line        = element_line(colour = "black"),
      panel.border     = element_rect(colour = "black", fill = NA),
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.position  = "none"
    )
}

build_metrics_figure <- function(platform_name) {
  df <- dplyr::filter(metrics_long, platform == platform_name)
  plots <- lapply(split(df, df$metric), plot_metric_panel)
  patchwork::wrap_plots(plots, nrow = 1)
}

message("Building metrics panels...")

p_mer_metrics <- build_metrics_figure("MERSCOPE")
p_xen_metrics <- build_metrics_figure("Xenium")

ggsave(file.path(opt$out_dir, "metrics_merscope.pdf"),
       p_mer_metrics, width = dims$full_w, height = dims$half_w * 0.9,
       units = "mm", device = cairo_pdf, bg = "white")
message("  Saved: metrics_merscope.pdf")

ggsave(file.path(opt$out_dir, "metrics_xenium.pdf"),
       p_xen_metrics, width = dims$full_w, height = dims$half_w * 0.9,
       units = "mm", device = cairo_pdf, bg = "white")
message("  Saved: metrics_xenium.pdf")

# ---------------------------------------------------------------------------
# Panels 3–4: UMAP coloured by cell type
# ---------------------------------------------------------------------------

build_umap_panel <- function(platform_name) {
  df     <- dplyr::filter(umap_coords, platform == platform_name)
  segs   <- levels(df$segmentation)
  # Subsample for rendering speed (>100k points gets slow)
  set.seed(42)
  if (nrow(df) > 150000) df <- dplyr::slice_sample(df, n = 150000)

  plots <- lapply(segs, function(seg) {
    d <- dplyr::filter(df, segmentation == seg)
    ggplot(d, aes(x = UMAP1, y = UMAP2, colour = cell_type)) +
      geom_point(size = 0.15, alpha = 0.5, stroke = 0) +
      scale_colour_manual(values = pal_cell_type, na.value = "grey80") +
      labs(title = seg, x = "UMAP1", y = "UMAP2", colour = NULL) +
      theme_sb() +
      theme(
        aspect.ratio        = 1,
        legend.position     = "none",
        plot.title.position = "plot",
        plot.title          = element_text(face = "bold", hjust = 0.5),
        axis.line           = element_line(colour = "black"),
        panel.background    = element_blank(),
        panel.border        = element_blank(),
        plot.background     = element_blank()
      )
  })

  patchwork::wrap_plots(plots, nrow = 1)
}

message("Building UMAP panels...")

p_mer_umap <- build_umap_panel("MERSCOPE")
p_xen_umap <- build_umap_panel("Xenium")

ggsave(file.path(opt$out_dir, "umap_merscope.pdf"),
       p_mer_umap, width = dims$full_w * 0.65, height = dims$half_w * 0.7,
       units = "mm", device = cairo_pdf, bg = "white")
message("  Saved: umap_merscope.pdf")

ggsave(file.path(opt$out_dir, "umap_xenium.pdf"),
       p_xen_umap, width = dims$full_w * 0.65, height = dims$half_w * 0.7,
       units = "mm", device = cairo_pdf, bg = "white")
message("  Saved: umap_xenium.pdf")

# ---------------------------------------------------------------------------
# Panel 5: cell type composition bar chart (mean cells per sample ± SEM)
# ---------------------------------------------------------------------------
message("Building cell counts panel...")

cell_summ <- cell_type_counts %>%
  dplyr::group_by(platform, segmentation, cell_type) %>%
  dplyr::summarise(
    mean_cells = mean(n_cells, na.rm = TRUE),
    sem_cells  = sd(n_cells, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups    = "drop"
  )

dodge_w <- 0.78
bar_w   <- 0.85

p_cell_counts <- ggplot(cell_summ,
                        aes(x = cell_type, y = mean_cells, fill = segmentation)) +
  geom_col(position = position_dodge(width = dodge_w), width = bar_w,
           colour = NA, alpha = 0.85) +
  geom_errorbar(
    aes(ymin = mean_cells - sem_cells, ymax = mean_cells + sem_cells),
    position  = position_dodge(width = dodge_w),
    width     = 0.18,
    linewidth = dims$line,
    colour    = "black"
  ) +
  coord_flip() +
  facet_wrap(~ platform, nrow = 1, scales = "free_x") +
  scale_fill_seg() +
  scale_y_continuous(
    labels = scales::label_number(scale_cut = scales::cut_si(""), accuracy = 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(x = NULL, y = "Mean cells per sample ± SEM", fill = NULL) +
  theme_sb() +
  theme(
    legend.position  = "top",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold")
  )

ggsave(file.path(opt$out_dir, "cell_counts.pdf"),
       p_cell_counts, width = dims$full_w, height = dims$half_w * 0.9,
       units = "mm", device = cairo_pdf, bg = "white")
message("  Saved: cell_counts.pdf")

# ---------------------------------------------------------------------------
# Panel 6: MECR boxplot + jitter, faceted by platform
# ---------------------------------------------------------------------------
message("Building MECR panel...")

p_mecr <- ggplot(mecr_table, aes(x = segmentation, y = MECR)) +
  geom_boxplot(aes(fill = segmentation), width = 0.7, outlier.shape = NA,
               colour = "black", alpha = 0.40, coef = 2) +
  geom_jitter(aes(fill = segmentation), width = 0.15, shape = 21,
              colour = "black", size = dims$point, alpha = 0.95) +
  scale_fill_seg() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
                     labels = scales::label_number(accuracy = 0.1)) +
  facet_wrap(~ platform, nrow = 1) +
  labs(x = NULL, y = "MECR") +
  theme_sb() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.line        = element_line(colour = "black"),
    panel.border     = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(colour = "black"),
    legend.position  = "none"
  )

ggsave(file.path(opt$out_dir, "mecr.pdf"),
       p_mecr, width = dims$half_w * 1.1, height = dims$half_w,
       units = "mm", device = cairo_pdf, bg = "white")
message("  Saved: mecr.pdf")

# ---------------------------------------------------------------------------
# Panel 7: purity dotplot (tx_assigned vs negative marker purity)
# ---------------------------------------------------------------------------
# Skipped gracefully if purity_summary is empty (scRNA ref was not available)
if (nrow(purity_summary) == 0) {
  message("purity_summary.rds is empty — skipping purity panels.")
  message("  Re-run segmentation_quality.R with a valid segmentation_comp.scrna_path.")
  message("Done. Outputs written to: ", opt$out_dir)
  quit(save = "no", status = 0)
}

message("Building purity dotplot panel...")

p_purity_dot <- ggplot(
  purity_summary,
  aes(x = tx_assigned, y = purity, fill = segmentation, shape = platform)
) +
  geom_point(size = 3.4, alpha = 0.95, colour = "black", stroke = 0.25) +
  scale_fill_seg() +
  scale_shape_manual(values = c(MERSCOPE = 21, Xenium = 24)) +
  scale_x_continuous(expand = expansion(mult = 0.08),
                     labels = scales::label_percent(accuracy = 1)) +
  scale_y_continuous(expand = expansion(mult = 0.08),
                     limits = c(0, 1)) +
  guides(
    shape = guide_legend(override.aes = list(fill = "white")),
    fill  = guide_legend(override.aes = list(shape = 21, colour = "black"))
  ) +
  labs(x = "Fraction assigned transcripts", y = "Negative-marker purity",
       shape = NULL, fill = NULL) +
  theme_sb() +
  theme(
    legend.position = "right",
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )

ggsave(file.path(opt$out_dir, "purity_dotplot.pdf"),
       p_purity_dot, width = dims$half_w * 1.4, height = dims$half_w,
       units = "mm", device = cairo_pdf, bg = "white")
message("  Saved: purity_dotplot.pdf")

# ---------------------------------------------------------------------------
# Panel 8: negative marker purity heatmap per cell type (ComplexHeatmap)
# ---------------------------------------------------------------------------
message("Building purity heatmap panel...")

if (nrow(purity_ct) == 0) {
  message("purity_ct.rds is empty — skipping purity heatmap.")
  message("Done. Outputs written to: ", opt$out_dir)
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

blues <- RColorBrewer::brewer.pal(9, "Blues")[1:7]
purity_col_fun <- circlize::colorRamp2(seq(0, 1, length.out = length(blues)), blues)

# Row order: cell types sorted by mean purity across all methods (ascending)
ct_order <- purity_ct %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(mu = mean(purity, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(mu) %>%
  dplyr::pull(cell_type)

# Helper: build ComplexHeatmap for one platform
make_purity_heatmap <- function(purity_ct_one, title, row_levels,
                                base_family = "Arial") {
  purity_ct_one <- dplyr::mutate(purity_ct_one,
    segmentation = factor(segmentation, levels = c("Default", "Proseg", "Cellpose")),
    cell_type    = as.character(cell_type)
  )

  mat_df <- purity_ct_one %>%
    dplyr::select(cell_type, segmentation, purity) %>%
    tidyr::pivot_wider(names_from = segmentation, values_from = purity) %>%
    dplyr::right_join(tibble::tibble(cell_type = row_levels), by = "cell_type") %>%
    dplyr::mutate(cell_type = factor(cell_type, levels = row_levels)) %>%
    dplyr::arrange(cell_type)

  mat_df <- as.data.frame(mat_df)
  rownames(mat_df) <- mat_df$cell_type
  mat_df$cell_type <- NULL
  for (col in c("Default", "Proseg", "Cellpose")) {
    if (!col %in% colnames(mat_df)) mat_df[[col]] <- NA_real_
  }
  mat <- as.matrix(mat_df)[, c("Default", "Proseg", "Cellpose"), drop = FALSE]

  # Row colour annotation (cell type)
  ct_pal     <- pal_cell_type
  missing_ct <- setdiff(rownames(mat), names(ct_pal))
  if (length(missing_ct) > 0) {
    ct_pal <- c(ct_pal, setNames(rep("grey80", length(missing_ct)), missing_ct))
  }

  ha_row <- ComplexHeatmap::rowAnnotation(
    `Cell type` = rownames(mat),
    col = list(`Cell type` = ct_pal),
    show_annotation_name = FALSE,
    width = unit(4, "mm"), show_legend = FALSE
  )

  ha_col <- ComplexHeatmap::HeatmapAnnotation(
    Segmentation = colnames(mat),
    col = list(Segmentation = pal_segmentation),
    show_annotation_name = FALSE,
    annotation_height = unit(4, "mm"), show_legend = TRUE
  )

  ComplexHeatmap::Heatmap(
    mat,
    name          = "Purity",
    col           = purity_col_fun,
    na_col        = "grey90",
    top_annotation = ha_col,
    left_annotation = ha_row,
    show_row_names       = TRUE,
    row_names_side       = "left",
    row_names_gp         = gpar(fontsize = 7, fontfamily = base_family),
    row_names_max_width  = unit(6, "cm"),
    show_column_names    = TRUE,
    column_names_gp      = gpar(fontsize = 7, fontfamily = base_family),
    cluster_rows         = FALSE,
    cluster_columns      = FALSE,
    rect_gp              = gpar(col = "white", lwd = 1.2),
    column_title         = title,
    column_title_gp      = gpar(fontsize = 8, fontface = "bold", fontfamily = base_family),
    heatmap_legend_param = list(
      title     = "Purity",
      at        = c(0, 0.25, 0.5, 0.75, 1),
      labels    = c("0", "0.25", "0.5", "0.75", "1"),
      title_gp  = gpar(fontsize = 6.5, fontfamily = base_family),
      labels_gp = gpar(fontsize = 6,   fontfamily = base_family)
    ),
    cell_fun = function(j, i, x, y, w, h, fill) {
      v       <- mat[i, j]
      txt_col <- if (!is.na(v) && v >= 0.80) "white" else "black"
      label   <- if (is.na(v)) "NA" else sprintf("%.2f", v)
      grid.text(label, x, y, gp = gpar(fontsize = 6, col = txt_col,
                                        fontfamily = base_family))
    }
  )
}

ht_mer <- make_purity_heatmap(
  dplyr::filter(purity_ct, platform == "MERSCOPE"),
  title = "MERSCOPE", row_levels = ct_order
)
ht_xen <- make_purity_heatmap(
  dplyr::filter(purity_ct, platform == "Xenium"),
  title = "Xenium", row_levels = ct_order
)

heatmap_path <- file.path(opt$out_dir, "purity_heatmap.pdf")
cairo_pdf(heatmap_path, width = (154 / 25.4), height = (110 / 25.4))
ComplexHeatmap::draw(
  ht_mer + ht_xen,
  heatmap_legend_side  = "right",
  annotation_legend_side = "right",
  merge_legends        = TRUE,
  padding              = unit(c(2, 2, 2, 2), "mm")
)
dev.off()
message("  Saved: purity_heatmap.pdf")

message("Done. All panels written to: ", opt$out_dir)
