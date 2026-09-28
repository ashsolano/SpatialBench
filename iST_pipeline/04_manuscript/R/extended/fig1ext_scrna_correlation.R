# Purpose:  Extended figure for Figure 1: pseudobulk logCPM of each spatial
#           platform (y) vs each single-cell reference (x) as one 2 x 6 grid
#           (rows 10X FLEX / 10X 3' GEX; columns Visium HD, MERSCOPE, Xenium on
#           pairwise gene sets, then on the 90-gene common set).
# Supports: fig1.R — spatial vs single-cell agreement holds for a second
#           reference (3' GEX) and on the gene set common to all three platforms.
# Legend notes:
#           - Matched WT animals (709, 713); axis titles omit "WT". 3' GEX WT708
#             and MERSCOPE/Xenium WT710 excluded.
#           - logCPM = log10(CPM + 1). Value = mean over animals of per-animal
#             pseudobulk log10(CPM + 1); library size = the pairwise shared genes.
#             R = Pearson on gene-level means; n = genes.
#           - Density contours use genes with logCPM > 0.1 in both datasets; all
#             genes are shown as points and included in R and n.
#           - Common set = Visium HD ∩ MERSCOPE ∩ Xenium panels (90 genes); values
#             subset from the pairwise sets, not renormalised.
#           - Axes 0–6 as in Fig 1c. MERSCOPE/Xenium values ~1 log10 unit higher
#             than in the earlier draft. Visium HD vs 3' GEX n = 19059 (earlier
#             draft 19053).
# Inputs:   results/03_benchmarking/scrna_correlation/correlation_by_reference.rds
# Outputs:  extended_figures/fig1ext_scrna_correlation/scrna_correlation_grid.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(optparse)
})

source("04_manuscript/R/utils/theme.R")
source("04_manuscript/R/utils/palettes.R")
source("04_manuscript/R/utils/plot_helpers.R")

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--in_file"), type = "character",
              default = "results/03_benchmarking/scrna_correlation/correlation_by_reference.rds",
              help    = "correlation_by_reference.rds from scrna_correlation.R [default: %default]"),
  make_option(c("--out_dir"), type = "character",
              default = "extended_figures/fig1ext_scrna_correlation",
              help    = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opt$in_file)) stop("--in_file not found: ", opt$in_file)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Figure constants
# ---------------------------------------------------------------------------
reference_labels <- c(FLEX = "10X FLEX",               GEX = "10X 3' GEX")
axis_ref_labels  <- c(FLEX = "10X FLEX",               GEX = "10X 3'GEX")
row_strip_labels <- c(FLEX = "snRNA-seq\n10X FLEX",    GEX = "scRNA-seq\n10X 3' GEX")
platform_labels  <- c(VisiumHD = "Visium HD", MERSCOPE = "MERSCOPE", Xenium = "Xenium")
gene_set_levels  <- c("pairwise", "common")

axis_limits <- c(0, 6)
axis_breaks <- c(0, 2, 4, 6)
page_w_mm   <- dims$full_w
gap_w_mm    <- 3                # between the two column groups
strip_w_mm  <- 5.5
header_h_mm <- 4
bar_h_mm    <- 5
bar_size_mm <- c(w = 6, h = 1)

# Text sizes (pt): the reference's sizes scaled from 210 to 170 mm wide
txt <- list(
  base       = 6,    # tick labels = base - 1
  axis_title = 5.5,
  annot      = 5,    # R / n
  label      = 6.5,  # headers and row strips
  bar_title  = 5,
  bar_label  = 4.5
)
contour_lw <- 0.1
density_cutoff <- 0.1   # contours use genes with log10(CPM + 1) > cutoff in both datasets
point_alpha    <- 0.25  # grey points (Fig 1c: 0.35)

# ---------------------------------------------------------------------------
# Load and check data
# ---------------------------------------------------------------------------
corr_by_ref  <- readRDS(opt$in_file)
corr_data    <- corr_by_ref$data
corr_summary <- corr_by_ref$summary
common_genes <- corr_by_ref$common_genes

message("WT animals: ", paste(corr_by_ref$samples$wt_ids, collapse = ", "))
purrr::iwalk(corr_by_ref$samples[names(platform_labels)], function(s, plat) {
  message("  ", plat, " samples: ", paste(s, collapse = ", "))
})

expected <- tidyr::expand_grid(reference = names(reference_labels),
                               platform  = names(platform_labels),
                               gene_set  = gene_set_levels)
missing  <- dplyr::anti_join(expected, corr_summary,
                             by = c("reference", "platform", "gene_set"))
if (nrow(missing) > 0 || nrow(corr_summary) != nrow(expected)) {
  stop("correlation_by_reference.rds does not contain exactly one row per ",
       "reference × platform × gene set")
}

n_common <- corr_summary$n_genes[corr_summary$gene_set == "common"]
if (any(n_common != length(common_genes))) {
  stop("Common-gene panels have n = {", paste(unique(n_common), collapse = ", "),
       "}, expected ", length(common_genes))
}
message("Common gene set: ", length(common_genes), " genes")

# Fail rather than silently crop genes outside the axes
value_range <- range(c(corr_data$scRNA, corr_data$ST), na.rm = TRUE)
message("Expression range: ", paste(round(value_range, 3), collapse = " – "))
if (value_range[1] < axis_limits[1] || value_range[2] > axis_limits[2]) {
  stop("Expression values fall outside the 0–6 axis range: ",
       paste(round(value_range, 3), collapse = " – "))
}

message("R and n per panel:")
print(as.data.frame(dplyr::mutate(corr_summary, r = round(r, 3))))

# ---------------------------------------------------------------------------
# Density scatter panels
# ---------------------------------------------------------------------------
make_panel <- function(ref, gene_set, plat) {
  panel_data <- dplyr::filter(corr_data, reference == ref,
                              gene_set == !!gene_set, platform == plat)
  stats      <- dplyr::filter(corr_summary, reference == ref,
                              gene_set == !!gene_set, platform == plat)
  density_data <- dplyr::filter(panel_data, scRNA > density_cutoff, ST > density_cutoff)
  generate_density_plot(
    panel_data, stats$r, stats$n_genes,
    platform_name     = platform_labels[[plat]],
    color_low         = pal_density[[plat]][["low"]],
    color_high        = pal_density[[plat]][["high"]],
    axis_limits       = axis_limits,
    adjust            = 1.5,
    bins              = density_contour_bins[[plat]],
    x_label           = paste0("logCPM (", axis_ref_labels[[ref]], ")"),
    y_label           = paste0("logCPM (", platform_labels[[plat]], ")"),
    annot_pos         = "top-left",
    annot_size        = txt$annot / .pt,
    contour_linewidth = contour_lw,
    axis_breaks       = axis_breaks,
    base_size         = txt$base,
    density_data      = density_data,
    point_alpha       = point_alpha
  ) +
    theme(axis.title = element_text(size = txt$axis_title))
}

# Row-major: FLEX then GEX; pairwise then common
panel_grid <- tidyr::expand_grid(reference = names(reference_labels),
                                 gene_set  = gene_set_levels,
                                 platform  = names(platform_labels))
panels <- purrr::pmap(panel_grid, function(reference, gene_set, platform) {
  make_panel(reference, gene_set, platform)
})

# ---------------------------------------------------------------------------
# Column-group headers, row strips and colour bars
# ---------------------------------------------------------------------------
# area = "panel" aligns the grob with the panel area of its row / column;
# box_margin replaces wrap_elements()' default 5.5 pt margin
wrap_grob <- function(grob, area = c("full", "panel"), box_margin = margin(0, 0, 0, 0)) {
  area    <- match.arg(area)
  wrapped <- if (area == "full") wrap_elements(full = grob) else wrap_elements(panel = grob)
  wrapped + theme(plot.margin = box_margin)
}

label_box <- function(label, fill, rot = 0) {
  grobTree(
    roundrectGrob(r = unit(0.8, "mm"), gp = gpar(fill = fill, col = NA)),
    textGrob(label, rot = rot,
             gp = gpar(fontfamily = "Arial", fontface = "bold",
                       fontsize = txt$label, lineheight = 0.9))
  )
}

# "Low [gradient] High" over "Point density", as in the reference; levels are
# relative within each panel, so no numbers
colour_bar <- function(color_low, color_high) {
  bar_w  <- unit(bar_size_mm[["w"]], "mm")
  bar_h  <- unit(bar_size_mm[["h"]], "mm")
  bar_y  <- unit(1, "npc") - unit(1.5, "mm")
  label_gp <- gpar(fontfamily = "Arial", fontsize = txt$bar_label)
  grobTree(
    rasterGrob(matrix(scales::seq_gradient_pal(color_low, color_high)(seq(0, 1, length.out = 300)), nrow = 1),
               x = unit(0.5, "npc"), y = bar_y, width = bar_w, height = bar_h,
               interpolate = TRUE),
    textGrob("Low",  x = unit(0.5, "npc") - 0.5 * bar_w - unit(0.5, "mm"), y = bar_y,
             hjust = 1, gp = label_gp),
    textGrob("High", x = unit(0.5, "npc") + 0.5 * bar_w + unit(0.5, "mm"), y = bar_y,
             hjust = 0, gp = label_gp),
    textGrob("Point density", x = unit(0.5, "npc"), y = bar_y - 0.5 * bar_h - unit(0.4, "mm"),
             vjust = 1,
             gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = txt$bar_title))
  )
}

headers <- list(
  wrap_grob(label_box("Pairwise platform-specific gene sets", fill = "grey85"),
            box_margin = margin(0, 0, 1.5, 0, "mm")),
  wrap_grob(label_box(paste0("Platform intersection gene set (", length(common_genes), " common)"),
                      fill = "grey85"),
            box_margin = margin(0, 0, 1.5, 0, "mm"))
)
row_strips <- purrr::map(row_strip_labels, function(label) {
  wrap_grob(label_box(label, fill = pal_muted[["snRNAseq"]], rot = -90),
            area = "panel", box_margin = margin(0, 0, 0, 1, "mm"))
})
bars <- purrr::map(rep(names(platform_labels), times = length(gene_set_levels)), function(plat) {
  wrap_grob(colour_bar(pal_density[[plat]][["low"]], pal_density[[plat]][["high"]]),
            area = "panel")
})

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------
# A/B headers; C–H, I–N panel rows; Q/R row strips; S–X colour bars; # empty.
# Areas are filled in alphabetical order of their letters.
design <- "
AAA#BBB#
CDE#FGHQ
IJK#LMNR
STU#VWX#
"

build_grid <- function(panel_mm) {
  wrap_plots(c(headers, panels, row_strips, bars), design = design) +
    plot_layout(
      widths  = unit(c(rep(panel_mm, 3), gap_w_mm, rep(panel_mm, 3), strip_w_mm), "mm"),
      heights = unit(c(header_h_mm, panel_mm, panel_mm, bar_h_mm), "mm")
    )
}

# Measured on a cairo device so text metrics match the saved PDF
layout_size_mm <- function(p) {
  cairo_pdf(tempfile(fileext = ".pdf"), width = 10, height = 10)
  on.exit(dev.off())
  gt <- patchwork::patchworkGrob(p)
  c(width  = convertWidth(sum(gt$widths),   "mm", valueOnly = TRUE),
    height = convertHeight(sum(gt$heights), "mm", valueOnly = TRUE))
}

# Panel side that makes the page exactly page_w_mm wide
probe_mm <- 20
overhead <- layout_size_mm(build_grid(probe_mm))[["width"]] - 6 * probe_mm
panel_mm <- (page_w_mm - overhead) / 6
if (panel_mm <= 0) stop("Page width leaves no room for panels")

p_grid <- build_grid(panel_mm)
page   <- layout_size_mm(p_grid)

out_file <- file.path(opt$out_dir, "scrna_correlation_grid.pdf")
ggsave(
  out_file, p_grid,
  width  = page[["width"]],
  height = page[["height"]],
  units  = "mm",
  device = cairo_pdf,
  bg     = "white"
)
message("Saved: ", out_file,
        sprintf(" (%.1f x %.1f mm, panels %.1f mm)", page[["width"]], page[["height"]], panel_mm))
