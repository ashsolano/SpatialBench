# Purpose:  Figure 2 gene comparison panels: pseudobulk MDS plot and
#           per-gene average-expression scatter plots across platforms.
#           Reads the DGEList produced by gene_comparison.R and saves each
#           panel as a separate PDF.
#           Adapted from fig2_pseudobulk_mds.R and fig2_avgexpr_scatter_v2.R.
# Inputs:   results/03_benchmarking/gene_comparison/dge.rds
# Outputs:  figures/fig2/pseudobulk_mds.pdf
#           figures/fig2/avgexpr_scatter_visiumhd_merscope.pdf
#           figures/fig2/avgexpr_scatter_visiumhd_xenium.pdf
#           figures/fig2/avgexpr_scatter_merscope_xenium.pdf

suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(vegan)
  library(optparse)
})

source("04_manuscript/R/utils/theme.R")
source("04_manuscript/R/utils/palettes.R")

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--input_dir"), type = "character",
              default = "results/03_benchmarking/gene_comparison",
              help    = "Directory containing gene_comparison.R outputs [default: %default]"),
  make_option(c("--out_dir"),   type = "character",
              default = "figures/fig2",
              help    = "Output directory for panel PDFs [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dge_path <- file.path(opt$input_dir, "dge.rds")
if (!file.exists(dge_path)) stop("Not found: ", dge_path)

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load and normalise DGEList
# ---------------------------------------------------------------------------
message("Loading DGEList...")
dge    <- readRDS(dge_path) |> calcNormFactors()

# Metadata: rownames are Platform_SampleID (column names of counts_mat)
meta <- as.data.frame(dge$samples) |>
  rownames_to_column("ColName")    # ColName = join key matching expr_mat cols

# ===========================================================================
# PANEL: Pseudobulk MDS
# ===========================================================================

# ---------------------------------------------------------------------------
# Build MDS plot
# Computes euclidean distances on log-CPM, runs classical MDS, and plots
# each sample as a filled shape coloured by platform and shaped by condition.
# ---------------------------------------------------------------------------
build_pseudobulk_mds <- function(dge, meta) {

  logcpm   <- cpm(dge, log = TRUE, prior.count = 2)
  dist_mat <- dist(t(logcpm), method = "euclidean")
  mds_res  <- cmdscale(dist_mat, k = 2, eig = TRUE)

  # Percent variance explained from positive eigenvalues only
  eigs_pos <- mds_res$eig[mds_res$eig > 0]
  var_ex   <- eigs_pos / sum(eigs_pos)
  pct1     <- round(var_ex[1] * 100, 1)
  pct2     <- round(var_ex[2] * 100, 1)

  # Combine MDS coordinates with sample metadata
  mds_df <- as.data.frame(mds_res$points) |>
    setNames(c("Dim1", "Dim2")) |>
    rownames_to_column("ColName") |>
    left_join(meta, by = "ColName") |>
    mutate(strip_label = paste0(length(rownames(dge)), " Common genes"))

  ggplot(mds_df, aes(x = Dim1, y = Dim2)) +

    geom_point(aes(fill = Platform, shape = Type),
               size   = 5,
               colour = "black",
               stroke = 0.5) +

    geom_label_repel(aes(label = IDnum),
                     fill              = "white",
                     colour            = "black",
                     size              = 2,
                     box.padding       = unit(0.3, "lines"),
                     point.padding     = unit(0.2, "lines"),
                     segment.color     = "grey70",
                     segment.curvature = 0.2,
                     segment.ncp       = 5,
                     segment.angle     = 45,
                     max.overlaps      = Inf) +

    scale_fill_platform(name = "Platform") +

    scale_shape_manual(
      name   = "Sample type",
      values = c(WT   = 21,   # filled circle
                 KO   = 24,   # filled triangle up
                 CTRL = 22)   # filled square
    ) +

    # Strip label showing the gene count
    facet_grid(~ strip_label) +

    labs(
      x = paste0("Dim1 (", pct1, "%)"),
      y = paste0("Dim2 (", pct2, "%)")
    ) +

    theme_sb() +
    theme(
      legend.position   = "right",
      legend.key.size   = unit(3, "mm"),
      legend.key.width  = unit(3, "mm"),
      legend.key.height = unit(3, "mm"),
      legend.title      = element_text(size = 8),
      legend.text       = element_text(size = 7),
      legend.spacing.y  = unit(0.5, "mm"),
      legend.spacing.x  = unit(0.5, "mm"),
      legend.margin     = margin(2, 2, 2, 2),
      panel.background  = element_blank(),
      panel.border      = element_rect(colour = "black", fill = NA),
      axis.line         = element_line(colour = "black"),
      strip.background  = element_rect(colour = "black", fill = "white"),
      strip.text        = element_text(face = "italic"),
      plot.background   = element_blank(),
      panel.spacing     = unit(2, "mm")
    ) +

    guides(
      fill = guide_legend(
        override.aes = list(shape = 21, colour = "black", size = 4),
        keywidth     = unit(3, "mm"),
        keyheight    = unit(3, "mm")
      ),
      shape = guide_legend(
        override.aes = list(fill = "white", colour = "black", size = 4),
        keywidth     = unit(3, "mm"),
        keyheight    = unit(3, "mm")
      )
    )
}

message("Building pseudobulk MDS panel...")
p_mds <- build_pseudobulk_mds(dge, meta)

ggsave(
  filename = file.path(opt$out_dir, "pseudobulk_mds.pdf"),
  plot     = p_mds,
  width    = dims$half_w * 1.3,
  height   = dims$half_w,
  units    = "mm",
  device   = cairo_pdf,
  bg       = "white"
)
message("Saved: pseudobulk_mds.pdf")

# ===========================================================================
# PANELS: Average expression scatter plots (one per platform pair)
# ===========================================================================

# ---------------------------------------------------------------------------
# Compute log10(CPM+1) expression matrix and global axis limits once,
# so all three scatter panels use the same scale for fair comparison.
# ---------------------------------------------------------------------------
message("Computing expression matrix and global axis limits...")

cpm_mat  <- cpm(dge, log = FALSE)
expr_mat <- log10(cpm_mat + 1)
expr_lab <- "log10(CPM+1)"

# Derive per-sample platform membership for pivoting
avg_all <- as.data.frame(expr_mat) |>
  rownames_to_column("Gene") |>
  pivot_longer(-Gene, names_to = "ColName", values_to = "Expr") |>
  inner_join(meta, by = "ColName") |>
  group_by(Gene, Platform) |>
  summarise(meanExpr = mean(Expr), .groups = "drop") |>
  pivot_wider(names_from = Platform, values_from = meanExpr)

# Global min/max across all platform-average columns
platform_cols <- setdiff(colnames(avg_all), "Gene")
axis_limits <- c(
  min(as.matrix(avg_all[, platform_cols]), na.rm = TRUE),
  max(as.matrix(avg_all[, platform_cols]), na.rm = TRUE)
)

# ---------------------------------------------------------------------------
# Build one scatter panel for a pair of platforms.
# Points are per-gene average log10(CPM+1); density contours show the
# joint distribution; the 5 most divergent genes are labelled.
# ---------------------------------------------------------------------------
build_avgexpr_scatter <- function(dge, meta, plat_x, plat_y,
                                  axis_limits, panel_title = NULL,
                                  prior_count = 2, ncontours = 5) {

  cpm_mat  <- cpm(dge, log = FALSE)
  expr_mat <- log10(cpm_mat + 1)

  avg_df <- as.data.frame(expr_mat) |>
    rownames_to_column("Gene") |>
    pivot_longer(-Gene, names_to = "ColName", values_to = "Expr") |>
    inner_join(meta, by = "ColName") |>
    filter(Platform %in% c(plat_x, plat_y)) |>
    group_by(Gene, Platform) |>
    summarise(meanExpr = mean(Expr), .groups = "drop") |>
    pivot_wider(names_from = Platform, values_from = meanExpr) |>
    rename(meanX = all_of(plat_x), meanY = all_of(plat_y))

  # Annotation position: R value in upper-left
  rng   <- diff(axis_limits)
  x_lbl <- axis_limits[1] + 0.03 * rng
  y_lbl <- axis_limits[2] - 0.03 * rng

  Rval <- cor(avg_df$meanX, avg_df$meanY, use = "pairwise.complete.obs")
  Rtxt <- paste0("R = ", round(Rval, 2))

  # Top 5 genes with the largest absolute difference between platforms
  top5 <- avg_df |>
    mutate(diff = abs(meanX - meanY)) |>
    slice_max(diff, n = 5)

  if (!is.null(panel_title)) {
    avg_df$Panel <- panel_title
    facet_call   <- facet_grid(~ Panel)
  } else {
    facet_call <- NULL
  }

  p <- ggplot(avg_df, aes(meanX, meanY)) +
    geom_point(size = 1, colour = "grey80") +
    stat_density_2d(colour = "black", bins = ncontours, linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    annotate("text",
             x = x_lbl, y = y_lbl, label = Rtxt,
             hjust = 0, vjust = 1, size = 3) +
    geom_text_repel(
      data          = top5,
      aes(label     = Gene),
      size          = 2.5,
      colour        = "blue",
      box.padding   = unit(0.2, "lines"),
      point.padding = unit(0.1, "lines")
    ) +
    labs(
      x = paste0(expr_lab, " (", plat_x, ")"),
      y = paste0(expr_lab, " (", plat_y, ")")
    ) +
    coord_equal() +
    scale_x_continuous(limits = axis_limits, expand = expansion(0)) +
    scale_y_continuous(limits = axis_limits, expand = expansion(0)) +
    theme_sb() +
    theme(
      aspect.ratio     = 1,
      strip.background = element_rect(colour = "black", fill = "white"),
      panel.background = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA),
      axis.line        = element_line(colour = "black"),
      strip.text       = element_text(face = "italic")
    )

  if (!is.null(facet_call)) p <- p + facet_call
  p
}

# ---------------------------------------------------------------------------
# Generate and save one panel per platform pair
# ---------------------------------------------------------------------------
pairs <- list(
  list(x = "VisiumHD", y = "MERSCOPE",
       file = "avgexpr_scatter_visiumhd_merscope.pdf",
       title = "VisiumHD vs MERSCOPE"),
  list(x = "VisiumHD", y = "Xenium",
       file = "avgexpr_scatter_visiumhd_xenium.pdf",
       title = "VisiumHD vs Xenium"),
  list(x = "MERSCOPE", y = "Xenium",
       file = "avgexpr_scatter_merscope_xenium.pdf",
       title = "MERSCOPE vs Xenium")
)

scatter_w <- dims$half_w * 0.65   # ~55 mm — square panels

for (pr in pairs) {
  message("Building scatter: ", pr$x, " vs ", pr$y, "...")
  p <- build_avgexpr_scatter(
    dge         = dge,
    meta        = meta,
    plat_x      = pr$x,
    plat_y      = pr$y,
    axis_limits = axis_limits,
    panel_title = pr$title
  )
  ggsave(
    filename = file.path(opt$out_dir, pr$file),
    plot     = p,
    width    = scatter_w,
    height   = scatter_w,
    units    = "mm",
    device   = cairo_pdf,
    bg       = "white"
  )
  message("Saved: ", pr$file)
}

message("Done. All panels written to: ", opt$out_dir)
