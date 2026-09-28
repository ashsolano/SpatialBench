# Purpose:  Shared plotting helpers for manuscript figure scripts
#           (scRNA-seq vs ST density scatter). Source after theme.R.
# Inputs:   none
# Outputs:  density_contour_bins, generate_density_plot() (in calling environment)

# Contour levels per platform; fewer for Visium HD (~15-19k genes)
density_contour_bins <- c(VisiumHD = 6, MERSCOPE = 8, Xenium = 8)

# Pseudobulk log10(CPM + 1) density scatter, platform (y) vs single-cell
# reference (x). Defaults reproduce Fig 1c. density_data (default: expr_data)
# feeds the contours only; points, R and n always use expr_data.
generate_density_plot <- function(expr_data, cor_value, n_genes,
                                  platform_name, color_low, color_high,
                                  axis_limits, adjust = 1.5, bins = 8,
                                  x_label = "10X FLEX WT log10(CPM + 1)",
                                  y_label = paste0(platform_name, " WT log10(CPM + 1)"),
                                  annot_pos = c("bottom-right", "top-left"),
                                  annot_size  = 3,
                                  contour_linewidth = 0.3,
                                  axis_breaks = NULL,
                                  base_size   = 7,
                                  density_data = NULL,
                                  point_alpha  = 0.35) {
  annot_pos <- match.arg(annot_pos)

  inset <- 0.05 * diff(axis_limits)
  annot <- if (annot_pos == "bottom-right") {
    list(x = axis_limits[2] - inset, y = axis_limits[1] + inset, hjust = 1, vjust = 0)
  } else {
    list(x = axis_limits[1] + inset, y = axis_limits[2] - inset, hjust = 0, vjust = 1)
  }

  breaks_scales <- if (!is.null(axis_breaks)) {
    list(scale_x_continuous(breaks = axis_breaks),
         scale_y_continuous(breaks = axis_breaks))
  }

  ggplot(expr_data, aes(x = scRNA, y = ST)) +
    geom_point(color = "grey80", size = 0.4, alpha = point_alpha, na.rm = TRUE) +
    stat_density_2d(
      data = density_data,
      aes(fill = after_stat(level), alpha = after_stat(level)),
      geom = "polygon", color = "black", linewidth = contour_linewidth,
      contour = TRUE, bins = bins, adjust = adjust, na.rm = TRUE
    ) +
    scale_fill_gradient(low = color_low, high = color_high) +
    scale_alpha(range = c(0.2, 0.75), guide = "none") +
    breaks_scales +
    geom_abline(slope = 1, intercept = 0,
                color = "black", linewidth = 0.4, linetype = "dashed") +
    # Widen the KDE grid to the full axis range without dropping points
    # (scale limits would set out-of-range values to NA)
    expand_limits(x = axis_limits, y = axis_limits) +
    coord_fixed(xlim = axis_limits, ylim = axis_limits, expand = FALSE) +
    annotate("text", x = annot$x, y = annot$y,
             label = paste0("R = ", round(cor_value, 2), "\nn = ", n_genes),
             size = annot_size, hjust = annot$hjust, vjust = annot$vjust,
             color = "black") +
    labs(x = x_label, y = y_label) +
    theme_sb(base_size = base_size) +
    theme(
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
      aspect.ratio    = 1,
      legend.position = "none"
    )
}
