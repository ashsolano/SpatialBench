# Purpose:  Manuscript Figure 1 — cross-platform dataset summary.
#           Panel A: horizontal bar chart of bins, transcripts, and sparsity for
#           VisiumHD, MERSCOPE, and Xenium at 8µm and 16µm binning resolutions.
#           Panel B: Venn diagram of gene panel overlap across the three platforms.
#           Panel C: density scatter plots of pseudobulk log10(CPM+1) correlation
#           between 10X FLEX scRNA-seq and each ST platform (WT samples, 8µm).
# Inputs:   results/03_benchmarking/dataset_summary/metrics.rds
#           results/03_benchmarking/dataset_summary/gene_lists.rds
#           results/03_benchmarking/scrna_correlation/avg_expr.rds
# Outputs:  figures/fig1/fig1_barplot.pdf
#           figures/fig1/fig1_venn.pdf
#           figures/fig1/fig1_scrna_correlation.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(ggtext)
  library(ggVennDiagram)
  library(patchwork)
  library(optparse)
})

source("04_manuscript/R/utils/theme.R")
source("04_manuscript/R/utils/palettes.R")

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--input_rds"),     type = "character", default = NULL,
              help = "Path to metrics.rds produced by dataset_summary.R"),
  make_option(c("--gene_lists"),    type = "character", default = NULL,
              help = "Path to gene_lists.rds produced by dataset_summary.R"),
  make_option(c("--scrna_cor_rds"), type = "character",
              default = "results/03_benchmarking/scrna_correlation/avg_expr.rds",
              help    = "Path to avg_expr.rds produced by scrna_correlation.R [default: %default]"),
  make_option(c("--out_dir"),       type = "character",
              default = "figures/fig1",
              help    = "Output directory for figure PDFs [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds))               stop("--input_rds is required")
if (is.null(opt$gene_lists))              stop("--gene_lists is required")
if (!file.exists(opt$scrna_cor_rds)) stop("--scrna_cor_rds not found: ", opt$scrna_cor_rds)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load pre-computed data
# ---------------------------------------------------------------------------
combined_data <- readRDS(opt$input_rds)
gene_lists    <- readRDS(opt$gene_lists)

# ---------------------------------------------------------------------------
# Factor ordering for bar plot
# ---------------------------------------------------------------------------
metric_levels  <- c("Bins", "Transcripts", "Transcripts (Common)",
                    "Sparsity", "Sparsity (Common)")
binning_levels <- c("8µm", "16µm")

# Drop the Genes metric and any rows with NA binning
combined_data <- combined_data |>
  dplyr::filter(Metric != "Genes", !is.na(Binning)) |>
  dplyr::mutate(
    Metric  = factor(Metric,  levels = metric_levels),
    Binning = factor(Binning, levels = binning_levels)
  )

# ---------------------------------------------------------------------------
# Summary statistics — mean ± SEM per Technology × Metric × Binning
# ---------------------------------------------------------------------------
summary_data <- combined_data |>
  dplyr::group_by(Technology, Metric, Binning) |>
  dplyr::summarise(
    MeanValue = mean(Value, na.rm = TRUE),
    SD        = sd(Value,   na.rm = TRUE),
    SEM       = SD / sqrt(dplyr::n()),
    .groups   = "drop"
  ) |>
  dplyr::mutate(
    Metric  = factor(Metric,  levels = metric_levels),
    Binning = factor(Binning, levels = binning_levels)
  )

# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

# Density scatter plot of pseudobulk log10(CPM+1) for one platform vs scRNA-seq.
# color_low/color_high are the gradient endpoints for the 2D density contours.
generate_density_plot <- function(expr_data, cor_value, n_genes,
                                  platform_name, color_low, color_high,
                                  axis_limits) {
  x_pos <- axis_limits[2] - 0.05 * diff(axis_limits)
  y_pos <- axis_limits[1] + 0.05 * diff(axis_limits)

  ggplot(expr_data, aes(x = scRNA, y = ST)) +
    geom_point(color = "grey80", size = 0.4, alpha = 0.35, na.rm = TRUE) +
    stat_density_2d(
      aes(fill = after_stat(level), alpha = after_stat(level)),
      geom = "polygon", color = "black", linewidth = 0.3,
      contour = TRUE, bins = 8, adjust = 1.5, na.rm = TRUE
    ) +
    scale_fill_gradient(low = color_low, high = color_high) +
    scale_alpha(range = c(0.2, 0.75), guide = "none") +
    geom_abline(slope = 1, intercept = 0,
                color = "black", linewidth = 0.4, linetype = "dashed") +
    scale_x_continuous(limits = axis_limits, expand = expansion(0)) +
    scale_y_continuous(limits = axis_limits, expand = expansion(0)) +
    coord_fixed() +
    annotate("text", x = x_pos, y = y_pos,
             label = paste0("R = ", round(cor_value, 2), "\nn = ", n_genes),
             size = 3, hjust = 1, vjust = 0, color = "black") +
    labs(x = "10X FLEX WT log10(CPM + 1)",
         y = paste0(platform_name, " WT log10(CPM + 1)")) +
    theme_sb() +
    theme(
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
      aspect.ratio    = 1,
      legend.position = "none"
    )
}

# Format numbers with SI suffix and one decimal place (e.g. 1.2M, 300K)
fmt_si_1dp <- function(x) {
  out <- character(length(x))
  for (i in seq_along(x)) {
    v  <- x[i]
    if (is.na(v)) { out[i] <- NA_character_; next }
    av <- abs(v)
    if      (av >= 1e9) out[i] <- paste0(round(v / 1e9, 1), "B")
    else if (av >= 1e6) out[i] <- paste0(round(v / 1e6, 1), "M")
    else if (av >= 1e3) out[i] <- paste0(round(v / 1e3, 1), "K")
    else                out[i] <- as.character(round(v, 1))
  }
  out
}

# Horizontal bar chart with mean bars, SEM error bars, and individual points.
# Uses scale_fill_platform() from theme.R so colours stay consistent with the
# rest of the manuscript.
plot_fun <- function(df, summary_df, title) {
  pos <- position_dodge(width = 0.8)

  ggplot() +
    geom_col(
      data     = summary_df,
      aes(y = Technology, x = MeanValue, fill = Technology, group = Technology),
      position = pos, orientation = "y", width = 0.7, alpha = 0.8
    ) +
    geom_errorbar(
      data     = summary_df,
      aes(y    = Technology, x = MeanValue,
          xmin = MeanValue - SEM, xmax = MeanValue + SEM, group = Technology),
      position = pos, orientation = "y", width = 0.2, colour = "black"
    ) +
    geom_point(
      data     = df,
      aes(y = Technology, x = Value, group = Technology),
      position = pos, orientation = "y",
      size = 0.3, shape = 21, colour = "black", fill = "black"
    ) +
    geom_text(
      data     = summary_df,
      aes(
        y     = Technology, x = MeanValue,
        label = ifelse(grepl("Sparsity", Metric),
                       sprintf("%.2f", MeanValue),
                       fmt_si_1dp(MeanValue)),
        group = Technology
      ),
      position = pos, orientation = "y", hjust = -0.3,
      size = 6 / .pt, family = "Arial"
    ) +
    scale_fill_platform() +
    facet_grid(Binning ~ Metric, scales = "free_x", space = "fixed") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.50))) +
    theme_sb() +
    theme(
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.title.x    = element_blank(),
      axis.title.y    = element_blank(),
      strip.text.x    = element_text(size = 7, family = "Arial", face = "bold"),
      strip.text.y    = element_text(size = 7, family = "Arial", face = "bold"),
      panel.spacing   = unit(0.8, "lines"),
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      legend.position = "right",
      plot.title      = element_text(size = 7, family = "Arial", face = "bold")
    ) +
    labs(title = title, y = "", x = "")
}

# ---------------------------------------------------------------------------
# Panel A — bar plot
# ---------------------------------------------------------------------------
p_barplot <- plot_fun(
  combined_data,
  summary_data,
  "Cross-platform dataset summary (8µm vs 16µm binning)"
)

# ---------------------------------------------------------------------------
# Panel B — Venn diagram of gene panel overlap
# ---------------------------------------------------------------------------
genes_visium   <- gene_lists$VisiumHD
genes_xenium   <- gene_lists$Xenium
genes_merscope <- gene_lists$MERSCOPE

venn_list <- list(
  "Visium HD" = genes_visium,
  "Xenium"    = genes_xenium,
  "MERSCOPE"  = genes_merscope
)

set_labels <- c(
  "Visium HD" = paste0("Visium HD\n(", length(genes_visium), " genes)"),
  "Xenium"    = paste0("Xenium\n(",    length(genes_xenium), " genes)"),
  "MERSCOPE"  = paste0("MERSCOPE\n(", length(genes_merscope), " genes)")
)

p_venn <- ggVennDiagram(
  venn_list,
  set_labels  = set_labels,
  label_alpha = 0,
  label       = "count",
  edge_size   = 1,
  set_color   = c("#6DAA2C", "#4477AA", "#AA3377")
) +
  scale_fill_gradient(low = "white", high = "white") +
  theme_void() +
  theme(legend.position = "none")

# ---------------------------------------------------------------------------
# Save panels separately
# ---------------------------------------------------------------------------
ggsave(
  filename = file.path(opt$out_dir, "fig1_barplot.pdf"),
  plot     = p_barplot,
  width    = dims$full_w,
  height   = dims$height / 3,
  units    = "mm",
  device   = cairo_pdf
)
message("Saved: ", file.path(opt$out_dir, "fig1_barplot.pdf"))

ggsave(
  filename = file.path(opt$out_dir, "fig1_venn.pdf"),
  plot     = p_venn,
  width    = dims$half_w,
  height   = dims$half_w,
  units    = "mm",
  device   = cairo_pdf
)
message("Saved: ", file.path(opt$out_dir, "fig1_venn.pdf"))

# ---------------------------------------------------------------------------
# Panel C — scRNA-seq vs ST pseudobulk correlation density plots
# ---------------------------------------------------------------------------

# Light variants of pal_muted for the density gradient low ends
pal_muted_light <- c(
  VisiumHD = "#e5f5d6",   # near-white olive-green
  MERSCOPE = "#f5d9ec",   # near-white plum
  Xenium   = "#d9e8f5"    # near-white indigo
)

avg_expr <- readRDS(opt$scrna_cor_rds)

# Shared axis limits across all three platforms for visual comparability
all_values <- c(
  avg_expr$VisiumHD$data$scRNA, avg_expr$VisiumHD$data$ST,
  avg_expr$MERSCOPE$data$scRNA, avg_expr$MERSCOPE$data$ST,
  avg_expr$Xenium$data$scRNA,   avg_expr$Xenium$data$ST
)
axis_limits <- c(floor(min(all_values, na.rm = TRUE)),
                 ceiling(max(all_values, na.rm = TRUE)))

p_visiumhd_cor <- generate_density_plot(
  avg_expr$VisiumHD$data, avg_expr$VisiumHD$correlation, avg_expr$VisiumHD$n_genes,
  platform_name = "Visium HD",
  color_low  = pal_muted_light["VisiumHD"],
  color_high = pal_muted["VisiumHD"],
  axis_limits = axis_limits
)

p_merscope_cor <- generate_density_plot(
  avg_expr$MERSCOPE$data, avg_expr$MERSCOPE$correlation, avg_expr$MERSCOPE$n_genes,
  platform_name = "MERSCOPE",
  color_low  = pal_muted_light["MERSCOPE"],
  color_high = pal_muted["MERSCOPE"],
  axis_limits = axis_limits
)

p_xenium_cor <- generate_density_plot(
  avg_expr$Xenium$data, avg_expr$Xenium$correlation, avg_expr$Xenium$n_genes,
  platform_name = "Xenium",
  color_low  = pal_muted_light["Xenium"],
  color_high = pal_muted["Xenium"],
  axis_limits = axis_limits
)

p_fig1c <- p_visiumhd_cor + p_merscope_cor + p_xenium_cor +
  plot_layout(ncol = 3)

ggsave(
  filename = file.path(opt$out_dir, "fig1_scrna_correlation.pdf"),
  plot     = p_fig1c,
  width    = dims$full_w,
  height   = dims$full_w / 3,   # square panels side by side
  units    = "mm",
  device   = cairo_pdf,
  bg       = "white"
)
message("Saved: ", file.path(opt$out_dir, "fig1_scrna_correlation.pdf"))
