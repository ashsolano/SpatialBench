# Visualisation helper functions.
#
# Key spatial plot functions:
#   plot_spatial_binned() — aggregates 8 µm bins into a coarser grid (default
#     300 bins per axis) using mean (continuous) or majority vote (discrete),
#     then renders with geom_raster().  Fast for full-tissue views.
#   plot_spatial_tiles() — renders every 8 µm bin as an individual tile
#     (hardcoded 29.3 px bin size).  Slower but preserves full resolution;

plot_samples <- function(spes, spot_metric, labs, title,
                         discrete = FALSE, size = 0.25) {
  if (is.list(spes)) {
    p <- sapply(spes,
      function(spe) {
        tibble::tibble(
          counts = spot_metric(spe),
          "x_coor" = spatialCoords(spe)[, 1],
          "y_coor" = spatialCoords(spe)[, 2],
          sample = metadata(spe)$sample
        )
      },
      simplify = FALSE
    ) %>%
      do.call(rbind, .) |>
      ggplot(aes(x = x_coor, y = y_coor, col = counts)) +
      facet_wrap(~sample, ncol = 2)
  } else {
    p <- spatialCoords(spes) |>
      data.frame() |>
      ggplot(aes(
        x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        col = spot_metric(spes)
      )) +
      facet_wrap(~ spes$sample_id, ncol = 2)
  }
  p <- p +
    geom_point(size = size) + coord_fixed() +
    xlab("") + ylab("") + theme_minimal() +
    labs + title +
    if (discrete) {
      scale_colour_discrete()
    } else {
      scale_colour_gradientn(
        colours = colorRampPalette(rev(brewer.pal(11, "Spectral")))(100)
      )
    }
  return(p)
}

plot_pca_density <- function(spe, dimname = "PCA", dims = c(1, 2), col,
                             density_alpha = 0.3, points_alpha = 0.7) {
  data <- data.frame(
    x = reducedDim(spe, dimname)[, dims[1]],
    y = reducedDim(spe, dimname)[, dims[2]],
    col = if (missing(col)) {
      spe$sample_id
    } else {
      col
    }
  )
  p <- ggplot(data, aes(x = x, y = y, col = col)) +
    geom_point(alpha = points_alpha)
  xdens <- axis_canvas(p, axis = "x") +
    geom_density(data = data, aes(x = x, fill = col, col = col), alpha = density_alpha)
  ydens <- axis_canvas(p, axis = "y", coord_flip = TRUE) +
    geom_density(data = data, aes(x = y, fill = col, col = col), alpha = density_alpha) +
    coord_flip()

  return(
    p |>
      cowplot::insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") |>
      cowplot::insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") |>
      ggdraw()
  )
}

plot_complexheatmap <- function(spe, cluster_label, marker_genes_df,
                                counts_fn = SingleCellExperiment::logcounts) {
  m <- spe[marker_genes_df$EnsembleID, order(cluster_label)] |>
    counts_fn() |>
    as.matrix() |>
    scale()

  col_fun <- circlize::colorRamp2(
    seq(min(m, na.rm = TRUE), mean(m, na.rm = TRUE) + 2 * sd(m, na.rm = TRUE), length.out = 100),
    grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100)
  )

  return(ComplexHeatmap::Heatmap(
    m,
    col = col_fun,
    column_split = sort(cluster_label),
    row_split = marker_genes_df$zone,
    cluster_rows = F, cluster_columns = F,
    show_column_names = F,
    use_raster = T, raster_quality = 5,
    raster_resize_mat = mean,
    row_labels = rowData(spe[marker_genes_df$EnsembleID, ])$symbol,
    heatmap_legend_param = list(
      title = "relative expression\n",
      title_gp = grid::gpar(fontsize = 16, fontface = "plain")
    ),
    row_title_rot = 0, column_title_rot = 45
  ))
}

logfc_mtx <- function(spe, cluster_label, marker_genes_df) {
  args.grid <- expand.grid(
    zone = as.character(unique(marker_genes_df$zone)),
    cluster = as.character(unique(cluster_label))
  ) |>
    tibble::as_tibble()

  parallel::mcmapply(function(zone, cluster) {
    rows <- marker_genes_df[marker_genes_df$zone == zone, "EnsembleID"]
    cols <- cluster_label == cluster
    cluster_score <- sum(SingleCellExperiment::counts(spe[rows, cols])) / sum(cols)
    other_score <- sum(SingleCellExperiment::counts(spe[rows, !cols])) / sum(!cols)
    tibble::tibble_row(
      zone = zone,
      cluster = cluster,
      logfc = log2(cluster_score / other_score)
    )
  }, args.grid$zone, args.grid$cluster, mc.cores = 4, SIMPLIFY = F) |>
    dplyr::bind_rows()
}

# deprecated
# too slow
logfc_mtx_gene <- function(spe, cluster_label, marker_genes_df) {
  args.grid <- expand.grid(
    zone = as.character(unique(marker_genes_df$zone)),
    cluster = as.character(unique(cluster_label))
  ) |>
    tibble::as_tibble()

  parallel::mcmapply(function(zone, cluster) {
    rows <- marker_genes_df[marker_genes_df$zone == zone, "gene"]
    cols <- cluster_label == cluster
    cluster_score <- sum(SingleCellExperiment::counts(spe[rows, cols])) / sum(cols)
    other_score <- sum(SingleCellExperiment::counts(spe[rows, !cols])) / sum(!cols)
    tibble::tibble_row(
      zone = zone,
      cluster = cluster,
      logfc = log2(cluster_score / other_score)
    )
  }, args.grid$zone, args.grid$cluster, mc.cores = 4, SIMPLIFY = F) |>
    dplyr::bind_rows()
}


# faster implementation of logfc_mtx_gene with S4Arrays::rowsum / colsum
logfc_tb_zone_rowsum <- function(spe, cluster_label, marker_genes_df) {
  marker_genes_df <- marker_genes_df |>
    as.data.frame() |>
    dplyr::arrange(zone)
  if (!all(marker_genes_df$gene %in% rownames(spe))) {
    warning(
      sprintf(
        "%s not in rownames(spe)",
        paste0(marker_genes_df$gene[!marker_genes_df$gene %in% rownames(spe)], collapse = ", ")
      )
    )
    marker_genes_df <- marker_genes_df[marker_genes_df$gene %in% rownames(spe), ]
  }
  spe <- spe[marker_genes_df$gene, ]
  mtx <- SingleCellExperiment::counts(spe) |>
    as("CsparseMatrix") |>
    SparseArray::rowsum(group = marker_genes_df$zone, reorder = TRUE) |>
    S4Arrays::colsum(group = cluster_label, reorder = TRUE)

  cell_numbers <- table(cluster_label) |>
    as.data.frame() |>
    (\(x) setNames(x$Freq, x$cluster_label))()

  args.grid <- tidyr::expand_grid(
    zone = unique(marker_genes_df$zone),
    cluster = unique(cluster_label),
  ) |>
    tibble::as_tibble()

  mapply(function(zone, cluster) {
    score <- mtx[zone, cluster] / sum(cluster_label == cluster)
    other_score <- sum(mtx[zone, colnames(mtx) != cluster]) / sum(cluster_label != cluster)
    tibble::tibble_row(
      zone = zone,
      cluster = cluster,
      logfc = log2(score / other_score),
      expr = score
    )
  }, args.grid$zone, args.grid$cluster, SIMPLIFY = FALSE) |>
    dplyr::bind_rows() |>
    dplyr::mutate(zone = factor(zone, levels = levels(marker_genes_df$zone)))
}

cowplot_plotVisium <- function(spe, col) {
  spe$col <- col
  cowplot::plot_grid(
    plotlist = unique(spe$sample_id) |>
      lapply(function(i) {
        spe[, spe$sample_id == i] |>
          ggspavis::plotVisium(annotate = "col")
      })
  )
}

# https://github.com/XiaoZhangryy/iSC.MEB/blob/233471ca6a96e2b5ada79f150d3574c1cd906e15/R/Visualization.r
iSC.MEB_palette <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plot_spe_raster <- function(
    spe,
    feature,
    feature_name = "Feature",
    continuous_palette = rev(RColorBrewer::brewer.pal(11, "Spectral")),
    discrete_palette = NULL) {
  stopifnot(length(feature) == ncol(spe))

  # Build tibble with spatial coords + feature
  tb <- spatialCoords(spe) |>
    data.frame() |>
    as_tibble() |>
    mutate(
      sample_id = colData(spe)$sample_id,
      feature = feature
    )

  # Compute bin sizes (fixed 300 bins per axis per sample)
  bin_sizes <- tb |>
    group_by(sample_id) |>
    summarise(
      bin_size_x = (max(x) - min(x)) / 300,
      bin_size_y = (max(y) - min(y)) / 300,
      .groups = "drop"
    )

  # Snap coords to grid & aggregate
  tb_grid <- tb |>
    left_join(bin_sizes, by = "sample_id") |>
    mutate(
      x = round(x / bin_size_x) * bin_size_x,
      y = round(y / bin_size_y) * bin_size_y
    ) |>
    group_by(sample_id, x, y) |>
    summarise(
      feature = if (is.numeric(feature)) {
        mean(feature, na.rm = TRUE)
      } else {
        f <- na.omit(feature)
        if (length(f) == 0) {
          NA
        } else {
          names(sort(table(f), decreasing = TRUE))[1]
        }
      },
      .groups = "drop"
    )

  # Shared theme
  base_theme <- theme_minimal() +
    theme(
      plot.margin = unit(c(1, 0, 0, 0), "cm"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14),
      strip.text = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14)
    )

  # Choose fill scale depending on continuous vs discrete
  if (is.numeric(feature)) {
    fill_scale <- scale_fill_gradientn(
      colours = colorRampPalette(continuous_palette)(100)
    )
    guide <- guides(fill = guide_colorbar(title = feature_name))
  } else {
    if (is.null(discrete_palette)) {
      # default discrete palette: hue
      n_levels <- length(unique(na.omit(tb_grid$feature)))
      discrete_palette <- scales::hue_pal()(n_levels)
    }
    fill_scale <- scale_fill_manual(values = discrete_palette)
    guide <- guides(fill = guide_legend(title = feature_name))
  }

  # Plot
  ggplot(tb_grid, aes(x = x, y = y, fill = feature)) +
    geom_raster() +
    coord_fixed() +
    facet_wrap(~sample_id, ncol = 2) +
    base_theme +
    fill_scale +
    guide +
    xlab("") +
    ylab("") +
    ggtitle(feature_name)
}


plot_spatial_binned <- function(
    spe, value, # column in colData(spe) to plot
    n_bins = 300, # resolution in x and y
    palette # optional custom palette
    ) {
  # stopifnot(value_col %in% colnames(colData(spe)))
  if (length(value) == 1) {
    stopifnot(value %in% colnames(colData(spe)))
    tb <- spatialCoords(spe) |>
      data.frame() |>
      as_tibble() |>
      mutate(
        sample_id = spe$sample_id,
        value = colData(spe)[[value]]
      )
  } else {
    stopifnot(length(value) == ncol(spe))
    tb <- spatialCoords(spe) |>
      data.frame() |>
      as_tibble() |>
      mutate(
        sample_id = spe$sample_id,
        value = value
      )
  }


  # bin sizes per sample
  bin_sizes <- tb |>
    group_by(sample_id) |>
    summarise(
      bin_size_x = (max(x) - min(x)) / n_bins,
      bin_size_y = (max(y) - min(y)) / n_bins,
      .groups = "drop"
    )

  # binning
  tb_grid <- tb |>
    left_join(bin_sizes, by = "sample_id") |>
    mutate(
      x = round(x / bin_size_x) * bin_size_x,
      y = round(y / bin_size_y) * bin_size_y
    )

  # decide categorical vs numerical
  if (is.numeric(tb_grid$value)) {
    tb_grid <- tb_grid |>
      group_by(sample_id, x, y) |>
      summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
    fill_scale <- scale_fill_viridis_c()
  } else {
    tb_grid <- tb_grid |>
      group_by(sample_id, x, y) |>
      summarise(
        value = names(sort(table(value), decreasing = TRUE))[1],
        .groups = "drop"
      )
    if (missing(palette)) {
      fill_scale <- scale_fill_brewer(palette = "Set3")
    } else {
      fill_scale <- scale_fill_manual(values = palette)
    }
  }

  # shared theme
  base_theme <- theme_minimal() +
    theme(
      plot.margin = unit(c(1, 0, 0, 0), "cm"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14),
      strip.text = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14)
    )

  ggplot(tb_grid, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    coord_fixed() +
    facet_wrap(~sample_id, ncol = 2) +
    fill_scale +
    base_theme +
    xlab("") +
    ylab("")
}

plot_spatial_tiles <- function(
    spe, value, # column in colData(spe) to plot
    palette # optional custom palette
    ) {
  # stopifnot(value_col %in% colnames(colData(spe)))
  if (length(value) == 1) {
    stopifnot(value %in% colnames(colData(spe)))
    tb <- spatialCoords(spe) |>
      data.frame() |>
      as_tibble() |>
      mutate(
        sample_id = spe$sample_id,
        value = colData(spe)[[value]]
      )
  } else {
    stopifnot(length(value) == ncol(spe))
    tb <- spatialCoords(spe) |>
      data.frame() |>
      as_tibble() |>
      mutate(
        sample_id = spe$sample_id,
        value = value
      )
  }

  tb <- tb |>
    mutate(
      x = x %/% 29.3,
      y = y %/% 29.3
    )

  # decide categorical vs numerical
  if (is.numeric(tb$value)) {
    fill_scale <- scale_fill_viridis_c()
  } else {
    if (missing(palette)) {
      fill_scale <- scale_fill_brewer(palette = "Set3")
    } else {
      fill_scale <- scale_fill_manual(values = palette)
    }
  }

  # shared theme
  base_theme <- theme_minimal() +
    theme(
      plot.margin = unit(c(1, 0, 0, 0), "cm"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14),
      strip.text = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14)
    )

  p <- ggplot(tb, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    coord_fixed() +
    fill_scale +
    base_theme +
    xlab("") +
    ylab("")
  if (length(unique(tb$sample_id)) > 1) {
    p <- p + facet_wrap(~sample_id, ncol = 2)
  }

  return(p)
}


# simple Venn diagram using eulerr (two sets only)
simple_2sets_eulerr <- function(x, y) {
  # Find the intersection and unique elements
  intersection <- length(intersect(x, y))
  x_only <- length(setdiff(x, y))
  y_only <- length(setdiff(y, x))

  # Create the named vector for the Venn diagram
  venn_data <- c(
    "x only" = x_only,
    "y only" = y_only,
    "Intersection" = intersection
  )

  # Fit the Venn diagram using euler
  fit <- eulerr::euler(venn_data)

  # Plot the Venn diagram
  return(plot(fit))
}

# Wrapper function that takes multiple vectors and returns a logical matrix
strings_to_matrix <- function(...) {
  # Capture all the input vectors
  input_list <- list(...)

  # Combine all vectors into one unique set of values (all unique values across all input vectors)
  all_values <- unique(unlist(input_list))

  # Initialize a logical matrix
  result_matrix <- matrix(FALSE, nrow = length(all_values), ncol = length(input_list))

  # Set column names to the names of the input arguments
  colnames(result_matrix) <- names(input_list)

  # Set row names to the unique values
  rownames(result_matrix) <- all_values

  # Loop through each input vector and fill in the logical matrix
  for (i in seq_along(input_list)) {
    result_matrix[, i] <- all_values %in% input_list[[i]]
  }

  # append set sizes to column names
  colnames(result_matrix) <- paste0(
    colnames(result_matrix),
    " (",
    colSums(result_matrix),
    ")"
  )

  # Return the logical matrix
  return(result_matrix)
}

# helper function to add arial font from my folder
add_arial <- function() {
  sysfonts::font_add(
    family = "Arial",
    regular = "~/fonts/Arial/Arial.ttf",
    bold = "~/fonts/Arial/Arial Bold.ttf",
    italic = "~/fonts/Arial/Arial Italic.ttf",
    bolditalic = "~/fonts/Arial/Arial Bold Italic.ttf"
  )
  showtext::showtext_auto()
}


plot_volcano <- function(
    efit, p.val.cutoff = 0.05, logFC.cutoff = 1, first = "KO",
    hl_gene, box.padding = 0.35, point.padding = 0.5, size = 5, nudge_x = -2, nudge_y = 0.4) {
  tb <- limma::topTable(efit, n = Inf) |>
    mutate(
      DE_Status = case_when(
        adj.P.Val < p.val.cutoff & logFC > logFC.cutoff ~
          sprintf("Significant (Upregulated in %s)", first),
        adj.P.Val < p.val.cutoff & logFC < -logFC.cutoff ~
          sprintf("Significant (Downregulated in %s)", first),
        TRUE ~ "Not Significant"
      )
    )

  p <- tb |>
    ggplot(aes(x = logFC, y = -log10(adj.P.Val))) +
    ggrastr::geom_point_rast(
      data = filter(tb, DE_Status == "Not Significant"),
      aes(color = DE_Status), size = 0.5, alpha = 0.5
    ) +
    geom_point(
      data = filter(tb, DE_Status != "Not Significant"),
      aes(color = DE_Status), size = 2, alpha = 0.8
    )

  if (!missing(hl_gene)) {
    # highlight the specified gene in yellow and add a label
    p <- p +
      geom_point(
        data = filter(tb, symbol == hl_gene),
        aes(x = logFC, y = -log10(adj.P.Val)),
        color = "yellow", size = 3
      ) +
      ggrepel::geom_text_repel(
        data = subset(tb, symbol == hl_gene),
        aes(label = symbol),
        # box.padding = 0.35, point.padding = 0.5,
        # size = 5, nudge_x = -2, nudge_y = 0.4
        box.padding = box.padding, point.padding = point.padding,
        size = size, nudge_x = nudge_x, nudge_y = nudge_y
      )
  }
  p +
    scale_color_manual(
      values =
        setNames(
          c("#b2182b", "#2166ac", "grey"),
          c(
            sprintf("Significant (Upregulated in %s)", first),
            sprintf("Significant (Downregulated in %s)", first),
            "Not Significant"
          )
        ),
      drop = FALSE
    ) +
    # Add reference lines
    geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
    geom_vline(xintercept = -1, linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    # Center the x-axis
    scale_x_continuous(
      limits = c(
        -max(abs(topTable(efit, n = Inf)$logFC)),
        max(abs(topTable(efit, n = Inf)$logFC))
      )
    ) +
    # Add plot title and axis labels
    labs(
      # title = "Germinal Centre KO vs. WT",
      x = "Log2 Fold Change",
      y = "-Log10(Adjusted P-value)",
      color = "DE Status"
    ) +
    # Use a theme with a border
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(hjust = 0.5)
    )
}

plot_volcano_facet <- function(
    efits, id = "Cell-type", p.val.cutoff = 0.05, logFC.cutoff = 1, first = "KO",
    hl_gene, box.padding = 0.35, point.padding = 0.5, size = 5, nudge_x = -2, nudge_y = 5,
    scales = "fixed", ncol = NULL, nrow = NULL) {

  tb <- lapply(efits, function(x) {
    limma::topTable(x, n = Inf) |>
      mutate(
        DE_Status = case_when(
          adj.P.Val < p.val.cutoff & logFC > logFC.cutoff ~
            sprintf("Significant (Upregulated in %s)", first),
          adj.P.Val < p.val.cutoff & logFC < -logFC.cutoff ~
            sprintf("Significant (Downregulated in %s)", first),
          TRUE ~ "Not Significant"
        )
      ) |>
      as_tibble()
  }) |>
    bind_rows(.id = id) |>
    mutate(
      !!id := factor(!!sym(id), levels = names(efits))
    )

  hl_tb <- tb |>
    filter(symbol %in% hl_gene)

  p <- ggplot(tb, aes(x = logFC, y = -log10(adj.P.Val))) +
    # Rasterized layer for non-significant points
    ggrastr::geom_point_rast(
      data = filter(tb, DE_Status == "Not Significant"),
      aes(color = DE_Status), size = 0.5, alpha = 0.5
    ) +
    # add cutoff lines
    geom_vline(xintercept = c(-logFC.cutoff, logFC.cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(p.val.cutoff), linetype = "dashed") +
    # Standard layer for significant points
    geom_point(
      data = filter(tb, DE_Status != "Not Significant"),
      aes(color = DE_Status), size = 1, alpha = 0.8
    ) +
    # Highlight gene (Yellow point)
    geom_point(
      data = hl_tb,
      color = "yellow", size = 2
    ) +
    # Highlight gene label
    ggrepel::geom_text_repel(
      data = hl_tb,
      aes(label = symbol),
      box.padding = box.padding, point.padding = point.padding,
      size = size, nudge_x = nudge_x, nudge_y = nudge_y
    ) +
    facet_wrap(
      as.formula(sprintf("~`%s`", id)),
      scales = scales,
      ncol = ncol,
      nrow = nrow
    ) +
    xlim(-max(abs(tb$logFC)), max(abs(tb$logFC))) +
    scale_color_manual(
      values =
        setNames(
          c("#b2182b", "#2166ac", "grey"),
          c(
            sprintf("Significant (Upregulated in %s)", first),
            sprintf("Significant (Downregulated in %s)", first),
            "Not Significant"
          )
        ),
      drop = FALSE
    ) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      strip.text = element_text(face = "bold", size = 14, hjust = 0), # facet titles
      strip.background = element_blank()
    ) +
    labs(x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)", color = "DE Status")
  p
}
