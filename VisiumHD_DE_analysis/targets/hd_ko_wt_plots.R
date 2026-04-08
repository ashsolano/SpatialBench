# Visualisation targets for the Tbx21 KO vs wild-type Visium HD experiment.
# All plots are saved as PDFs in output/.  Targets with format = "file" return
# the output file path so that targets tracks file freshness automatically.

list(
  # ---------------------------------------------------------------------------
  # Colour palettes
  # ---------------------------------------------------------------------------
  # color palette from Ash
  # /vast/projects/SpatialBench/analysis/sub_workflows/spatialbench_manuscript/R_scripts/panels/fig4_umap_clusters.R
  # converted to rgb instead of previous rgba
  tar_target(
    ash_palette,
    c(
      "Erythrocytes"   = "#00A087",
      "Naive B cells"  = "#E9967A",
      "GC B Cells"     = "#525ecc",
      "Plasma B cells" = "#480607",
      "T cells"        = "#4DBBD5",
      "NK cells"       = "#7E6148",
      "ILC"            = "#802268",
      "Monocytes"      = "#F0E685",
      "Macrophages"    = "#CCEBC5",
      "DC"             = "#386CB0",
      "Granulocytes"   = "#cc52c0",
      "Stem cells"     = "#F0027F"
    )
  ),

  # Rahcel's color palette
  tar_target(rachel_palette,
    deployment = "main",
    c(
      "#9d183299", "#00A087FF", "#525ecc99", "#E9967A",
      "#4DBBD5FF", "#480607", "#386CB0", "#7E6148FF",
      "#ADB6B6FF", "#F0027F", "#F0E685FF", "#cc52c099",
      "#CCEBC5", "#DC143C", "#802268FF", "#FFBB78",
      "#3C548899", "#189b3699", "#cc9b5299", "#1B1919FF",
      "#6acc5299", "#cc527699", "#42857599", "#868686FF",
      "#006400", "#0A47FFFF", "#3B1B53FF", "#749B58FF",
      "#ccc05299", "#cc765299", "#1f0fdf99", "#00D68FFF",
      "#14FFB1FF", "#3C5488FF", "#8F7700FF", "#164194FF",
      "#0094CDFF", "#FFDAB9", "#FF00FF", "#00FFFF", "green",
      "blue", "#003C67FF", "#52a7cc99", "yellow"
    )
  ),

  # ---------------------------------------------------------------------------
  # QC spatial maps
  # plots covering all four tissue sections:
  #   libsize  — total UMI count per bin (Spectral palette).
  #   mt_pct   — mitochondrial gene fraction per bin.
  # ---------------------------------------------------------------------------
  # QC
  # UMI counts
  tar_target(
    spe_hd_ko_wt_libsize_8um,
    format = "file",
    {
      p <- plot_spatial_binned(
        spe = spe_hd_ko_wt_cut_8um_labeled,
        value = "sum",
      ) +
        scale_fill_gradientn(
          colours = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100)
        ) +
        guides(fill = guide_colorbar(title = "UMI counts"))

      file_path <- file.path("output", "spe_hd_ko_wt_libsize_8um.pdf")
      ggsave(file_path, p, width = 20, height = 20, dpi = 300, units = "cm")
      file_path
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c40g")
    )
  ),
  tar_target(
    spe_hd_ko_wt_mt_pct_8um,
    format = "file",
    {
      p <- plot_spatial_binned(
        spe = spe_hd_ko_wt_cut_8um_labeled,
        value = "subsets_mt_percent",
      ) +
        scale_fill_gradientn(
          colours = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100)
        ) +
        guides(fill = guide_colorbar(title = "% mitochondrial genes"))

      file_path <- file.path("output", "spe_hd_ko_wt_mt_pct_8um.pdf")
      ggsave(file_path, p, width = 20, height = 20, dpi = 300, units = "cm")
      file_path
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c40g")
    )
  ),

  # ---------------------------------------------------------------------------
  # Zoomed spatial view of a single germinal centre in sample 709 (WT).
  # Three panels side by side: RCTD cluster labels, light zone signature score,
  # and dark zone signature score (sum of raw counts for zone-specific genes).
  # ---------------------------------------------------------------------------
  tar_target(
    zoomed_spatial_plot,
    format = "file",
    {
      # subset to a gc region in 709
      spe <- spe_hd_ko_wt_cut_8um_labeled
      spe <- spe[, spe$sample_id == "hd_wt_8um_709"]
      spe <- spe[, spatialCoords(spe)[, 1] > 5500 & spatialCoords(spe)[, 1] < 7500 &
                   spatialCoords(spe)[, 2] > 10900 & spatialCoords(spe)[, 2] < 12900]
      spe$cluster <- fct_drop(spe$cluster)

      lz_sum <- colSums(SingleCellExperiment::counts(spe[dz_lz_signatures$gene[dz_lz_signatures$zone == "LZ"], ]))
      dz_sum <- colSums(SingleCellExperiment::counts(spe[dz_lz_signatures$gene[dz_lz_signatures$zone == "DZ"], ]))

      p1 <- plot_spatial_tiles(
        spe,
        "cluster",
        setNames(ash_palette, nm = NULL)
      ) +
        ggtitle("") +
        guides(fill = guide_legend(title = "Cluster"))
      p2 <- plot_spatial_tiles(spe, lz_sum) +
        ggtitle("Light zone signature genes expression") +
        guides(fill = guide_colorbar(title = "Sum"))
      p3 <- plot_spatial_tiles(spe, dz_sum) +
        ggtitle("Dark zone signature genes expression") +
        guides(fill = guide_colorbar(title = "Sum"))
      p <- p1 | p2 | p3

      file_path <- file.path("output", "zoomed_spatial_plot.pdf")
      ggsave(
        filename = file_path,
        plot = p,
        width = 40, height = 10, dpi = 300, units = "cm"
      )
      file_path
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c40g")
    )
  ),

  # ---------------------------------------------------------------------------
  # DE result plots — T-bet target gene enrichment
  # MA plots (mean expression vs log-fold-change) and barcode plots showing
  # the rank positions of T-bet target genes in each cell type's KO vs WT
  # DE analysis, ordered by ROAST p-value significance.
  # Points/bars are coloured by prior directional DE status (up/down in
  # Tbx21 KO from the reference bulk RNA-seq study).
  # ---------------------------------------------------------------------------
  # DE result plots
  # T-bet
  tar_target(HD_KO_WT_Tbet_MAplot_all,
    format = "file",
    deployment = "main",
    efits_hd_ko_wt_cut_8um |>
      # subset(names(efits_hd_ko_wt) != "B cell") |>
      (\(efits) {
        out_file <- file.path("output", "HD_KO_WT_Tbet_MAplot_all.pdf")
        pdf(out_file, width = 10, height = 16)
        par(mfrow = c(3, 2))
        for (efit in efits) {
          status <- efit$genes |>
            dplyr::left_join(prev_de, by = c("symbol" = "Symbol")) |>
            mutate(status = case_when(
              is.na(logFC) ~ "Other",
              symbol == "Tbx21" ~ "T-bet",
              logFC > 0 ~ "T-bet up-regulated",
              logFC < 0 ~ "T-bet down-regulated"
            )) |>
            dplyr::pull(status) |>
            factor()
          limma::plotMA(efit,
            status = status,
            main = paste0("KO vs. WT MA plot of ", efit$cluster, " DE analysis"),
            values = c("T-bet up-regulated", "T-bet down-regulated", "T-bet"),
            hl.col = c("red", "blue", "blue"),
            hl.pch = c(16, 16, 17),
            hl.cex = c(0.7, 0.7, 1)
          )
        }
        dev.off()
        return(out_file)
      })()
  ),
  tar_target(
    HD_KO_WT_Tbet_barcodeplot_all,
    format = "file",
    deployment = "main",

    efits_hd_ko_wt_cut_8um |>
      (\(efits) {
        add_arial()
        out_file <- file.path("output", "HD_KO_WT_Tbet_barcodeplot_all.pdf")
        pdf(out_file, width = 12, height = 16)
        par(
          mfrow = c(5, 3),
          family = "Arial",
          mai = c(0, 0, 0, 0),
          mar = c(3, 1, 1, 0) + 1.5
        )
        # order by significance
        clusters <- roast_table_tbet_hd_ko_wt_cut_8um |>
          arrange(P.Value) |>
          pull('cell-type')
        for (i in seq_along(efits)) {
          efit <- efits[[clusters[i]]]
          status <- efit$genes |>
            dplyr::left_join(prev_de, by = c("symbol" = "Symbol")) |>
            mutate(
              status = case_when(
                is.na(logFC) ~ "Other",
                logFC > 0 ~ "T-bet up-regulated",
                logFC < 0 ~ "T-bet down-regulated"
              )
            ) |>
            dplyr::pull(status) |>
            factor(levels = c("T-bet up-regulated", "T-bet down-regulated", "Other"))
          barcodeplot(efit$t[, "KovsWt"],
            index = status == "T-bet up-regulated", index2 = status == "T-bet down-regulated",
            # main = sprintf("Barcode plot of T-bet signature genes in KO vs. WT %s DE analysis", efit$cluster),
            main = ""
          )
          mtext(sprintf("%s. %s", LETTERS[i], efit$cluster), side = 3, adj = 0, line = 1, font = 2, cex = 1.2)
        }
        dev.off()
        return(out_file)
      })()

  ),
  # ---------------------------------------------------------------------------
  # DE result plots — sex-linked gene enrichment
  # MA and barcode plots for Y-chromosome and X-inactivation-escape genes
  # ---------------------------------------------------------------------------
  # Sex specific genes
  tar_target(HD_KO_WT_sex_MAplot_all,
    format = "file",
    deployment = "main",
    efits_hd_ko_wt_cut_8um |>
      # subset(names(efits_hd_ko_wt) != "B cell") |>
      (\(efits) {
        out_file <- file.path("output", "HD_KO_WT_sex_MAplot_all.pdf")
        pdf(out_file, width = 10, height = 16)
        par(mfrow = c(3, 2))
        for (efit in efits) {
          chr <- efit$gene$chromosome_name
          chr <- ifelse(efit$gene$symbol %in% xie_genes, "XiE", chr)
          chr[!chr %in% c("XiE", "Y")] <- "Other"
          chr <- factor(chr,
            levels = c("XiE", "Y", "Other"),
            labels = c("XiE gene", "chrY gene", "Other")
          )
          limma::plotMA(efit,
            status = chr,
            coef = "KovsWt",
            main = paste0("Male vs. Female\nMA plot of ", efit$cluster, " DE analysis"),
            hl.col = c("blue", "red", "black")
          )
        }
        dev.off()
        return(out_file)
      })()
  ),
  tar_target(
    HD_KO_WT_sex_barcodeplot_all,
    format = "file",
    deployment = "main",
    efits_hd_ko_wt_cut_8um |>
      (\(efits) {
        add_arial()
        out_file <- file.path("output", "HD_KO_WT_sex_barcodeplot_all.pdf")
        pdf(out_file, width = 28, height = 20)
        par(
          mfrow = c(4, 4),
          family = "Arial"
        )
        clusters <- roast_table_sex_hd_ko_wt_cut_8um |>
          arrange(P.Value) |>
          pull('cell-type')
        for (i in seq_along(efits)) {
          efit <- efits[[clusters[i]]]
          chr <- efit$gene$chromosome_name
          chr <- ifelse(efit$gene$symbol %in% xie_genes, "XiE", chr)
          chr[!chr %in% c("XiE", "Y")] <- "Other"
          chr <- factor(chr,
            levels = c("XiE", "Y", "Other")
          )
          barcodeplot(
            efit$t[, "KovsWt"],
            index = chr == "Y", index2 = chr == "XiE", main = ""
          )
          mtext(sprintf("%s. %s", LETTERS[i], efit$cluster), side = 3, adj = 0, line = 1, font = 2, cex = 1.8)
        }
        dev.off()
        return(out_file)
      })()
  ),

  # ---------------------------------------------------------------------------
  # DZ vs LZ barcode plot (WT samples only)
  # ---------------------------------------------------------------------------
  tar_target(
    HD_WT_dzlz_barcodeplot,
    format = "file",
    deployment = "main",
    {
      outfile <- file.path("output", "dz_lz_barcodeplot.pdf")
      add_arial()
      pdf(outfile, width = 8, height = 6)
      par(family = "Arial")
      efit <- efits_dz_lz_cut_8um[["wt"]]
      status <- efit$genes |>
        dplyr::mutate(symbol = toupper(symbol)) |>
        dplyr::left_join(bulk_dz_lz_signatures, by = c("symbol" = "symbol")) |>
        dplyr::mutate(
          status = case_when(
            is.na(z) ~ "Other",
            z > 0 ~ "Up in previous DZ vs. LZ signature",
            z < 0 ~ "Down in previous DZ vs. LZ signature"
          )
        ) |>
        dplyr::pull(status)
      barcodeplot(efit$t[, "DZvsLZ"],
        index = status == "Up in previous DZ vs. LZ signature",
        index2 = status == "Down in previous DZ vs. LZ signature",
        main = ""
      )
      dev.off()
      outfile
    }
  ),

  # ---------------------------------------------------------------------------
  # Spatial cluster / zone maps
  # spe_zone_hd_ko_wt_ggplot_8um_ash     — RCTD cell type labels
  # spe_zone_hd_ko_wt_ggplot_8um — same but different styling...
  # ---------------------------------------------------------------------------
  tar_target(
    spe_zone_hd_ko_wt_ggplot_8um,
    format = "file",
    {
      p <- plot_spatial_binned(
        spe = spe_hd_ko_wt_cut_8um_labeled,
        value = "cluster",
        palette = rachel_palette
      )

      file_path <- file.path("output", "spe_zone_hd_ko_wt_8um.pdf")
      ggsave(file_path, p, width = 20, height = 20, dpi = 300, units = "cm")
      file_path
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c40g")
    )
  ),
  # in ash's palette, with black background
  tar_target(
    spe_zone_hd_ko_wt_ggplot_8um_ash,
    format = "file",
    {
      spe <- spe_hd_ko_wt_cut_8um_labeled
      labels <- spe$cluster |>
        fct_recode(
          "GC B Cells" = "Light zone",
          "GC B Cells" = "Dark zone",
          "GC B Cells" = "Germinal centre (unassigned)"
        )
      p <- plot_spatial_binned(
        spe = spe,
        value = labels,
        palette = ash_palette
      ) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, color = "white"),
          strip.text = element_text(size = 14, color = "white"),
          legend.text = element_text(size = 14, color = "white"),
          legend.title = element_text(size = 14, color = "white"),
          plot.background = element_rect(fill = "black", color = NA),
          panel.background = element_rect(fill = "black", color = NA),
          panel.grid = element_blank()
        ) +
        annotation_custom(
          grob = grid::linesGrob(gp = grid::gpar(col = "white", lwd = 2)),
          xmin = 5, xmax = 3768, ymin = -5005, ymax = -5000
        )
      file_path <- file.path("output", "spe_zone_hd_ko_wt_8um_ash.pdf")
      ggsave(
        filename = file_path,
        plot = p,
        width = 20,
        height = 20,
        dpi = 300,
        units = "cm"
      )
      file_path
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c40g")
    )
  ),
  # ---------------------------------------------------------------------------
  # UMAP visualisations
  # banksy_umap — the joint Banksy UMAP.
  # umap_hd_ko_wt_8um — an expression-only UMAP (PCA on log-normalised
  #   counts) as an alternative view.
  # Both are coloured by RCTD cell type label using ash_palette.
  # ---------------------------------------------------------------------------
  # UMAP
  # Using the banksy UMAP
  tar_target(
    banksy_umap_hd_ko_wt_8um,
    format = "file",
    {
      spe <- spe_hd_ko_wt_cut_8um_labeled
      df <- cbind(
        as.data.frame(reducedDim(spe, "UMAP")),
        as.data.frame(colData(spe))
      ) |>
        mutate(
          cluster = fct_recode(cluster,
            "GC B Cells" = "Light zone",
            "GC B Cells" = "Dark zone",
            "GC B Cells" = "Germinal centre (unassigned)"
          )
        )

      p <- df |>
        ggplot(aes(x = V1, y = V2, color = cluster)) +
        ggrastr::geom_point_rast(size = 0.01, alpha = 0.5) +
        scale_color_manual(values = ash_palette) +
        theme_classic() +
        labs(
          title = "Banksy UMAP of VisiumHD samples",
          x = "UMAP1",
          y = "UMAP2",
          color = "Cell Type"
        ) +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12)
        ) +
        # Large points in the legend
        guides(color = guide_legend(override.aes = list(
          size = 4,
          alpha = 1
        )))
      file_path <- file.path("output", "banksy_umap_hd_ko_wt_8um.pdf")
      ggsave(file_path, p, width = 12, height = 9, dpi = 300)
      file_path
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c120g")
    )
  ),

  # Computing UMAP without spatial information
  tar_target(
    umap_df_hd_ko_wt_8um,
    {
      spe <- spe_hd_ko_wt_cut_8um_labeled[, spe_hd_ko_wt_cut_8um_labeled$filter]
      spe <- computeLibraryFactors(spe) |>
        logNormCounts()
      spe <- spe |>
        scater::runPCA() |>
        scater::runUMAP()

      # only return the UMAP coordinates and metadata
      cbind(
        as.data.frame(reducedDim(spe, "UMAP")),
        as.data.frame(colData(spe))
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c120g")
    )
  ),
  # UMAP plot
  tar_target(
    umap_hd_ko_wt_8um,
    format = "file",
    {
      p <- umap_df_hd_ko_wt_8um |>
        mutate(
          cluster = fct_recode(cluster,
            "GC B Cells" = "Light zone",
            "GC B Cells" = "Dark zone",
            "GC B Cells" = "Germinal centre (unassigned)"
          )
        ) |>
        ggplot(aes(x = UMAP1, y = UMAP2, color = cluster)) +
        ggrastr::geom_point_rast(size = 0.01, alpha = 0.3) +
        scale_color_manual(values = ash_palette) +
        theme_classic() +
        labs(
          title = "UMAP of VisiumHD samples",
          x = "UMAP1",
          y = "UMAP2",
          color = "Cell Type"
        ) +
        guides(color = guide_legend(override.aes = list(size = 4)))
      file_path <- file.path("output", "umap_hd_ko_wt_8um.pdf")
      ggsave(file_path, p, width = 12, height = 9, dpi = 300)
      file_path 
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c40g")
    )
  ),

  # ---------------------------------------------------------------------------
  # Tbx21 spatial expression map
  # Raw count for the Tbx21 gene plotted across all four tissue sections to
  # visually verify that KO samples have substantially reduced expression.
  # ---------------------------------------------------------------------------
  # Tbet counts
  tar_target(
    spe_hd_ko_wt_tbx21_8um,
    format = "file",
    {
      tbx21 <- SingleCellExperiment::counts(
        spe_hd_ko_wt_cut_8um_labeled
      )["Tbx21", ]

      p <- plot_spatial_binned(
        spe = spe_hd_ko_wt_cut_8um_labeled,
        value = tbx21
      ) +
        scale_fill_gradientn(
          colours = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100)
        ) +
        guides(fill = guide_colorbar(title = "Raw counts")) +
        ggtitle("Tbx21")

      file_path <- file.path("output", "spe_hd_ko_wt_tbx21_8um.pdf")
      ggsave(file_path, p, width = 20, height = 20, dpi = 300, units = "cm")
      file_path
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c40g")
    )
  ),

  # Use dot plot instead
  # zone_hd_ko_wt_cut_8um does not exists anymore
  # tar_target(
  #   heatmap_hd_ko_wt_zones_8um,
  #   format = "file",
  #   {
  #     markers_df <- rbind(zone_markers_df, dz_lz_signatures) |>
  #       as.data.frame() |>
  #       mutate(
  #         zone = factor(case_match(zone,
  #           "LZ" ~ "Light zone",
  #           "DZ" ~ "Dark zone",
  #           .default = zone
  #         ))
  #       ) |>
  #       dplyr::filter(
  #         zone %in% zone_hd_ko_wt_cut_8um
  #       ) |>
  #       mutate(zone = forcats::fct_drop(zone))
  #     tb_wt <- logfc_tb_zone_rowsum(
  #       spe_hd_ko_wt_cut_8um[, spe_hd_ko_wt_cut_8um$group == "wt"],
  #       zone_hd_ko_wt_cut_8um[spe_hd_ko_wt_cut_8um$group == "wt"],
  #       markers_df
  #     ) |>
  #       mutate(group = "wt")
  #     tb_ko <- logfc_tb_zone_rowsum(
  #       spe_hd_ko_wt_cut_8um[, spe_hd_ko_wt_cut_8um$group == "ko"],
  #       zone_hd_ko_wt_cut_8um[spe_hd_ko_wt_cut_8um$group == "ko"],
  #       markers_df
  #     ) |>
  #       mutate(group = "ko")
  #     zones <- tb_wt |>
  #       dplyr::filter(cluster == zone) |>
  #       arrange(desc(logfc)) |>
  #       dplyr::pull(zone)
  #     tb <- rbind(tb_wt, tb_ko)
  #     p <- tb |>
  #       ggplot(aes(
  #         x = factor(cluster, levels = zones),
  #         y = factor(zone, levels = rev(zones)), fill = logfc
  #       )) +
  #       geom_tile(
  #         color = "white",
  #         lwd = 0.5,
  #         linetype = 1
  #       ) +
  #       theme_void() +
  #       theme(
  #         axis.title = element_blank(), plot.title = element_text(hjust = 0.5),
  #         axis.text = element_text(size = 12), legend.title = element_blank(),
  #         legend.text = element_text(size = 12),
  #         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  #       ) +
  #       geom_text(aes(label = format(round(logfc, digits = 1), nsmall = 1)), color = "black", size = 4) +
  #       # scale_fill_distiller(palette = "RdBu") +
  #       scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
  #       facet_wrap(~group, ncol = 2) +
  #       coord_fixed() +
  #       labs(fill = "logFC") +
  #       ggtitle("logFC of marker genes in annotated clusters") +
  #       xlab("Annotated cluster") +
  #       ylab("Marker genes")

  #     file_path <- file.path("output", "heatmap_hd_ko_wt_zones_8um.pdf")
  #     ggsave(file_path, p, width = 45, height = 15, dpi = 600, units = "cm")
  #     file_path
  #   },
  #   resources = tar_resources(
  #     crew = tar_resources_crew(controller = "slurm_1c80g")
  #   )
  # ),

  # ---------------------------------------------------------------------------
  # Cell type composition bar chart
  # Proportion of each RCTD-assigned cell type per sample, grouped by KO/WT,
  # as a grouped bar chart to summarise overall tissue composition.
  # ---------------------------------------------------------------------------
  tar_target(
    cell_type_composition_hd_ko_wt_ggplot_8um,
    format = "file",
    {
      tb <- SummarizedExperiment::colData(spe_hd_ko_wt_cut_8um_labeled) |>
        as.data.frame() |>
        as_tibble()

      p_bars <- tb |>
        dplyr::count(group, cluster) |>
        mutate(pct = n / sum(n), .by = group) |>
        ggplot(aes(x = cluster, y = pct, fill = group)) +
        geom_col(position = "dodge")
      file_path <- file.path("output", "cell_type_composition_hd_ko_wt_ggplot_8um.pdf")
      ggsave(file_path, p_bars, width = 15, height = 10, dpi = 300, units = "cm")
      file_path
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c120g")
    )
  ),

  # ---------------------------------------------------------------------------
  # GC marker gene heatmaps (Banksy DZ/LZ subclusters)
  # For each genotype (KO, WT), bins classified as GC B cells are extracted
  # and their Banksy subclusters are used as columns.  Rows are DZ/LZ marker
  # genes; expression is row-scaled z-scored mean counts per cluster.
  # ---------------------------------------------------------------------------
  tar_target(
    GC_markers_heatmap_hd_8um_ko,
    format = "file",
    {
      #spe <- spe_hd_ko_wt_8um[, as.numeric(colnames(spe_banksy_gc_subclustering_hd_8um_ko))]
      #spe$cluster <- spe_banksy_gc_subclustering_hd_8um_ko$clust_M0_lam0.2_k50_res0.7
      spe <- spe_hd_ko_wt_cut_8um[, as.numeric(colnames(rctd_hd_ko_wt_8um_ko_gc_subclusters))]
      spe$cluster <- rctd_hd_ko_wt_8um_ko_gc_subclusters$clust_M0_lam0.2_k50_res0.5

      mtx <- SingleCellExperiment::counts(spe[dz_lz_signatures$gene, ]) |>
        as("CsparseMatrix") |>
        S4Arrays::colsum(group = spe$cluster, reorder = TRUE)

      cell_numbers <- table(spe$cluster) |>
        as.data.frame() |>
        (\(x) setNames(x$Freq, x$Var1))()
      cell_numbers_elsewhere <- sum(cell_numbers) - cell_numbers

      expr <- sweep(mtx, 2, cell_numbers[colnames(mtx)], "/")
      expr_elsewhere <- rowSums(mtx) - mtx
      expr_elsewhere <- sweep(expr_elsewhere, 2, cell_numbers_elsewhere[colnames(expr_elsewhere)], "/")
      logfc <- log2(expr / expr_elsewhere)

      p <- ComplexHeatmap::Heatmap(
        t(scale(t(expr))),
        name = "expression",
        col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_split = dz_lz_signatures$zone,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_title = "",
        column_title = "Cluster"
      )

      file_path <- file.path("output", "GC_markers_heatmap_hd_8um_ko.pdf")
      pdf(file_path, width = 10, height = 7)
      print(p)
      dev.off()
      file_path
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c80g")
    )
  ),

  tar_target(
    GC_markers_heatmap_hd_8um_wt,
    format = "file",
    {
      #spe <- spe_hd_ko_wt_8um[, as.numeric(colnames(spe_banksy_gc_subclustering_hd_8um_wt))]
      #spe$cluster <- spe_banksy_gc_subclustering_hd_8um_wt$clust_M0_lam0.2_k50_res0.7
      spe <- spe_hd_ko_wt_cut_8um[, as.numeric(colnames(rctd_hd_ko_wt_8um_wt_gc_subclusters))]
      spe$cluster <- rctd_hd_ko_wt_8um_wt_gc_subclusters$clust_M0_lam0.2_k50_res0.5

      mtx <- SingleCellExperiment::counts(spe[dz_lz_signatures$gene, ]) |>
        as("CsparseMatrix") |>
        S4Arrays::colsum(group = spe$cluster, reorder = TRUE)

      cell_numbers <- table(spe$cluster) |>
        as.data.frame() |>
        (\(x) setNames(x$Freq, x$Var1))()
      cell_numbers_elsewhere <- sum(cell_numbers) - cell_numbers

      expr <- sweep(mtx, 2, cell_numbers[colnames(mtx)], "/")
      expr_elsewhere <- rowSums(mtx) - mtx
      expr_elsewhere <- sweep(expr_elsewhere, 2, cell_numbers_elsewhere[colnames(expr_elsewhere)], "/")
      logfc <- log2(expr / expr_elsewhere)

      p <- ComplexHeatmap::Heatmap(
        t(scale(t(expr))),
        name = "expression",
        col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_split = dz_lz_signatures$zone,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_title = "",
        column_title = "Cluster"
      )

      file_path <- file.path("output", "GC_markers_heatmap_hd_8um_wt.pdf")
      pdf(file_path, width = 10, height = 7)
      print(p)
      dev.off()
      file_path
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c80g")
    )
  ),

  tar_target(
    GC_markers_heatmap_hd_ko_wt_ggplot_8um,
    format = "file",
    {
      spe <- spe_hd_ko_wt_cut_8um
      spe <- spe[dz_lz_signatures$gene, str_detect(spe$cluster, "GC") | spe$cluster %in% c("11", "13", "15")]

      mtx <- SingleCellExperiment::counts(spe) |>
        as("CsparseMatrix") |>
        S4Arrays::colsum(group = spe$cluster, reorder = TRUE)

      cell_numbers <- table(spe$cluster) |>
        as.data.frame() |>
        (\(x) setNames(x$Freq, x$Var1))()
      cell_numbers_elsewhere <- sum(cell_numbers) - cell_numbers

      expr <- sweep(mtx, 2, cell_numbers[colnames(mtx)], "/")
      expr_elsewhere <- rowSums(mtx) - mtx
      expr_elsewhere <- sweep(expr_elsewhere, 2, cell_numbers_elsewhere[colnames(expr_elsewhere)], "/")
      logfc <- log2(expr / expr_elsewhere)

      p <- ComplexHeatmap::Heatmap(
        t(scale(t(expr))),
        name = "expression",
        col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_split = dz_lz_signatures$zone,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_title = "",
        column_title = "Cluster"
      )

      file_path <- file.path("output", "GC_markers_heatmap_hd_ko_wt_8um.pdf")
      pdf(file_path, width = 10, height = 7)
      print(p)
      dev.off()
      file_path
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c40g")
    )
  ),

  # ---------------------------------------------------------------------------
  # Volcano plots
  # Standard −log10(adjusted p-value) vs log2 fold-change plots.
  # HD_DZ_LZ_volcano_plot_8um   — DZ vs LZ in WT samples.
  # HD_KO_WT_volcano_plot_B_cell_8um — KO vs WT in naive B cells.
  # HD_KO_WT_volcano_plot_gc_8um     — KO vs WT in germinal centre B cells
  #                                    (Tbx21 highlighted in yellow).
  # ---------------------------------------------------------------------------
  # volcano plots
  ## DZ vs. LZ
  tar_target(
    HD_DZ_LZ_volcano_plot_8um,
    format = "file",
    (
      topTable(efits_dz_lz_cut_8um$wt, n = Inf) |>
        ggplot(aes(x = logFC, y = -log10(adj.P.Val))) +
        geom_point(
          aes(
            color = ifelse(
              adj.P.Val < 0.05 & abs(logFC) > 1,
              ifelse(
                logFC > 0,
                "Significant (Upregulated)",
                "Significant (Downregulated)"
              ),
              "Not Significant"
            )
          ), size = 2, alpha = 0.8
        ) +
        scale_color_manual(values = c("Significant (Downregulated)" = "#2166ac",
                                      "Significant (Upregulated)" = "#b2182b",
                                      "Not Significant" = "grey")) +
        # Add reference lines
        geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
        geom_vline(xintercept = -1, linetype = "dashed", color = "black") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        # Center the x-axis by setting limits from -4 to 4
        scale_x_continuous(limits = c(-5.1, 5.1)) +
        # Add plot title and axis labels
        labs(title = "Dark Zone vs. Light Zone in Wile Type samples",
             x = "Log2 Fold Change",
             y = "-Log10(Adjusted P-value)",
             color = "DE Status") +
        # Use a theme with a border
        theme_classic() +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
              plot.title = element_text(hjust = 0.5))
    ) |>
      ggsave(
        filename = file.path("output", "HD_DZ_LZ_volcano_plot_8um.pdf"),
        plot = _,
        width = 8, height = 6
      )
  ),

  # GC WT vs. KO volcano plot
  tar_target(
    HD_KO_WT_volcano_plot_B_cell_8um,
    format = "file",
    (
      topTable(efits_hd_ko_wt_cut_8um$`Naive B cells`, n = Inf) |>
        mutate(
          DE_Status = case_when(
          adj.P.Val < 0.05 & logFC > 1 ~ "Significant (Upregulated in KO)",
          adj.P.Val < 0.05 & logFC < -1 ~ "Significant (Downregulated in KO)",
          TRUE ~ "Not Significant"
        )
      ) |>
        ggplot(aes(x = logFC, y = -log10(adj.P.Val))) +
        geom_point(aes(color = DE_Status), size = 2, alpha = 0.8) +
        # highlight Tbx21
        geom_point(
          data = filter(topTable(efits_hd_ko_wt_cut_8um$`Germinal centre`, n = Inf), symbol == "Tbx21"),
          aes(x = logFC, y = -log10(adj.P.Val)),
          color = "yellow", size = 3
        ) +
        ggrepel::geom_text_repel(
          data = subset(topTable(efits_hd_ko_wt_cut_8um$`Germinal centre`, n = Inf), symbol == "Tbx21"),
          aes(label = "Tbx21"),
          box.padding = 0.35, point.padding = 0.5, size = 5, nudge_x = -2, nudge_y = 0.4
        ) +
        scale_color_manual(values = c("Significant (Upregulated in KO)" = "#b2182b",
                                      "Significant (Downregulated in KO)" = "#2166ac",
                                      "Not Significant" = "grey")) +
        # Add reference lines
        geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
        geom_vline(xintercept = -1, linetype = "dashed", color = "black") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        # Center the x-axis by setting limits from -4 to 4
        scale_x_continuous(limits = c(-9, 9)) +
        # Add plot title and axis labels
        labs(title = "B Cell KO vs. WT",
             x = "Log2 Fold Change",
             y = "-Log10(Adjusted P-value)",
             color = "DE Status") +
        # Use a theme with a border
        theme_classic() +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
              plot.title = element_text(hjust = 0.5))
    ) |>
      ggsave(
        filename = file.path("output", "HD_KO_WT_volcano_plot_B_cell_8um.pdf"),
        plot = _,
        width = 8, height = 6
      )
  ),

  # GC WT vs. KO volcano plot
  tar_target(
    HD_KO_WT_volcano_plot_gc_8um,
    format = "file",
    (
      topTable(efits_hd_ko_wt_cut_8um$`Germinal centre`, n = Inf) |>
        mutate(
          DE_Status = case_when(
          adj.P.Val < 0.05 & logFC > 1 ~ "Significant (Upregulated in KO)",
          adj.P.Val < 0.05 & logFC < -1 ~ "Significant (Downregulated in KO)",
          TRUE ~ "Not Significant"
        )
      ) |>
        ggplot(aes(x = logFC, y = -log10(adj.P.Val))) +
        geom_point(aes(color = DE_Status), size = 2, alpha = 0.8) +
        # highlight Tbx21
        geom_point(
          data = filter(topTable(efits_hd_ko_wt_cut_8um$`Germinal centre`, n = Inf), symbol == "Tbx21"),
          aes(x = logFC, y = -log10(adj.P.Val)),
          color = "yellow", size = 3
        ) +
        ggrepel::geom_text_repel(
          data = subset(topTable(efits_hd_ko_wt_cut_8um$`Germinal centre`, n = Inf), symbol == "Tbx21"),
          aes(label = "Tbx21"),
          box.padding = 0.35, point.padding = 0.5, size = 5, nudge_x = -2, nudge_y = 0.4
        ) +
        scale_color_manual(values = c("Significant (Upregulated in KO)" = "#b2182b",
                                      "Significant (Downregulated in KO)" = "#2166ac",
                                      "Not Significant" = "grey")) +
        # Add reference lines
        geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
        geom_vline(xintercept = -1, linetype = "dashed", color = "black") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        # Center the x-axis by setting limits from -4 to 4
        scale_x_continuous(limits = c(-9, 9)) +
        # Add plot title and axis labels
        labs(title = "Germinal Centre KO vs. WT",
             x = "Log2 Fold Change",
             y = "-Log10(Adjusted P-value)",
             color = "DE Status") +
        # Use a theme with a border
        theme_classic() +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
              plot.title = element_text(hjust = 0.5))
    ) |>
      ggsave(
        filename = file.path("output", "HD_KO_WT_volcano_plot_gc_8um.pdf"),
        plot = _,
        width = 8, height = 6
      )
  ),

  # Tbet (ko vs. wt) volcano plots - all other cell types
  tar_target(
    HD_KO_WT_volcano_plot_other,
    format = "file",
    {
      cell_types <- roast_table_tbet_hd_ko_wt_cut_8um |>
        filter(! `cell-type` %in% c("Germinal centre", "Naive B cells")) |>
        arrange(P.Value) |>
        pull(`cell-type`) |>
        as.character()
      efits <- efits_hd_ko_wt_cut_8um[cell_types]
      add_arial()
      p_all <- plot_volcano_facet(efits, hl_gene = "Tbx21", size = 3) +
        theme(
          text = element_text(family = "Arial")
        )
      out_file <- file.path("output", "HD_KO_WT_volcano_plot_other_cell_types_8um.pdf")
      ggsave(out_file, p_all, width = 15, height = 10, dpi = 300)
      out_file
    }
  ),

  # expression for marker gene dot plots
  tar_target(
    dotplot_tb_rctd_hd_ko_wt_8um,
    {
      markers_df <- rbind(zone_markers_df, dz_lz_signatures) |>
        as_tibble() |>
        select(-EnsembleID) |>
        mutate(
          zone = fct_recode(
            zone,
            "Light zone" = "LZ",
            "Dark zone" = "DZ",
            "Plasma B cells" = "Plasma cell",
            "Erythrocytes" = "Erythrocyte",
            "Macrophages" = "Macrophage",
            "T cells" = "T cell"
          )
        )

      spe <- spe_hd_ko_wt_cut_8um_labeled
      gc_spe <- spe[, spe$cluster %in% c("Dark zone", "Light zone", "Germinal centre (unassigned)")]
      gc_spe$cluster <- "Germinal centre"
      spe <- cbind(spe, gc_spe)

      tb <- counts(spe)[rownames(spe) %in% markers_df$gene, ] |>
        as.matrix() |>
        t() |>
        as.data.frame()
      tb <- sapply(
        split.default(tb, names(tb)),
        function(x) rowSums(as.matrix(x))
      ) |>
        as.data.frame() |>
        as_tibble() |>
        mutate(
          # genotype = spe$group,
          cell_type = spe$cluster
        ) |>
        pivot_longer(
          cols = all_of(markers_df$gene),
          names_to = "gene",
          values_to = "expression"
        ) |>
        left_join(markers_df,
          by = "gene", relationship = "many-to-many"
        )

      tb
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c120g")
    )
  ),

  # dot plot for marker genes
  tar_target(
    dotplot_rctd_hd_ko_wt_8um,
    format = "file",
    {
      tb_summary <- dotplot_tb_rctd_hd_ko_wt_8um |>
        group_by(cell_type, gene, zone) |>
        summarise(
          average_expression = mean(expression),
          pct_expressing = sum(expression > 0) / n() * 100,
          .groups = "drop"
        )

      tb <- tb_summary |>
        filter(
          ! cell_type %in% c("Germinal centre (unassigned)", "Dark zone", "Light zone"),
          ! zone %in% c(
            "Dark zone", "Light zone",
            "Red pulp", "Marginal zone", "White pulp"
          )
        ) |>
        # reorder levels in cell_type and zone
        mutate(
          cell_type = factor(cell_type, levels = c(
            "Plasma B cells", "Naive B cells", "Germinal centre",
            "T cells", "Erythrocytes", "Macrophages", "Granulocytes",
            "Monocytes", "DC", "NK cells", "Stem cells"
          )),
          zone = factor(zone, levels = c(
            "Plasma B cells", "B cell", "Germinal centre",
            "T cells", "Erythrocytes", "Macrophages", "Neutrophil"
            # "Red pulp", "Marginal zone", "White pulp"
          ))
        ) |>
        mutate(
          scaled_average_expression = as.vector(scale(average_expression)),
          .by = c(gene, zone)
        ) |>
        na.fail()

      facet_bg <- tibble(zone = factor(levels(tb$zone), levels = levels(tb$zone))) |>
        mutate(
          id = row_number(),
          fill_group = id %% 2 == 0
        ) |>
        filter(fill_group) # keep every second zone for background shading
      diag_bg <- expand_grid(
        zone = factor(levels(tb$zone), levels = levels(tb$zone)),
        cell_type = factor(levels(tb$cell_type), levels = levels(tb$cell_type))
      ) |>
        mutate(
          fill_group = as.character(zone) == as.character(cell_type),
          y_center = as.numeric(cell_type),
          ymin = y_center - 0.5,
          ymax = y_center + 0.5
        ) |>
        filter(fill_group)


      p <- ggplot() +
        geom_rect(
          data = facet_bg,
          aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
          fill = "grey95",
          inherit.aes = FALSE # Prevents it from looking for 'color' or 'size'
        ) +
        geom_point(data = tb, aes(
          x = cell_type, y = gene, size = pct_expressing, color = scaled_average_expression
        )) +
        facet_grid(zone ~ ., scales = "free_y", space = "free_y", switch = "y") +
        scale_color_viridis_c(name = "Scaled Avg Expression") +
        scale_size_continuous(
          trans = "sqrt", range = c(0.01, 15),
          breaks = c(1, 5, 25, 50, 75), name = "% Expressed"
        ) +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7),
          strip.text = element_text(size = 9, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = NA),
          strip.placement = "outside", # Moves zone label below the gene names
          panel.spacing = unit(0.1, "lines")
        ) +
        labs(
          x = "Annotated cell type",
          y = "Gene"
        )
      file_path <- file.path("output", "dotplot_rctd_hd_ko_wt_8um.pdf")
      ggsave(file_path, p, width = 10, height = 16, dpi = 300)
      file_path
    }
  ),

  # dot plot for DZ LZ marker genes
  tar_target(
    dotplot_dzlz_rctd_hd_ko_wt_8um,
    format = "file",
    {
      tb_summary <- dotplot_tb_rctd_hd_ko_wt_8um |>
        filter(
          cell_type %in% c("Germinal centre (unassigned)", "Dark zone", "Light zone"),
          zone %in% c("Dark zone", "Light zone")
        ) |>
        group_by(cell_type, gene, zone) |>
        summarise(
          average_expression = mean(expression),
          pct_expressing = sum(expression > 0) / n() * 100,
          .groups = "drop"
        )

      tb <- tb_summary |>
        mutate(
          cell_type = factor(
            cell_type,
            levels = rev(c("Light zone", "Dark zone", "Germinal centre (unassigned)"))
          ),
          zone = factor(zone, levels = c("Light zone", "Dark zone"))
        ) |>
        na.fail() |>
        mutate(
          scaled_average_expression = as.vector(scale(average_expression)),
          .by = c(gene, zone)
        )

      facet_bg <- tibble(zone = factor(levels(tb$zone), levels = levels(tb$zone))) |>
        mutate(
          id = row_number(),
          fill_group = id %% 2 == 0
        ) |>
        filter(fill_group) # keep every second zone for background shading

      p <- ggplot() +
        geom_rect(
          data = facet_bg,
          aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
          fill = "grey95", # A very light grey
          alpha = 0.5,
          inherit.aes = FALSE # Prevents it from looking for 'color' or 'size'
        ) +
        geom_point(data = tb, aes(
          x = gene, y = cell_type, size = pct_expressing, color = scaled_average_expression
        )) +
        facet_grid(. ~ zone, scales = "free_x", space = "free_x", switch = "x") +
        scale_color_viridis_c(name = "Scaled Avg Expression") +
        scale_size_continuous(
          # trans = "sqrt", range = c(0.01, 15),
          breaks = c(1, 5, 25, 50, 75), name = "% Expressed"
        ) +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          strip.text = element_text(size = 9, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = NA),
          strip.placement = "outside",
          panel.spacing = unit(0.1, "lines")
        ) +
        labs(
          x = "Gene",
          y = "Annotated cell type"
        )
      file_path <- file.path("output", "dotplot_dzlz_rctd_hd_ko_wt_8um.pdf")
      ggsave(file_path, p, width = 10, height = 4, dpi = 300)
      file_path
    }
  ),

  # dot plot of Tbx21 by genotype and cell type
  tar_target(
    dotplot_Tbx21_hd_ko_wt_8um,
    format = "file",
    {
      spe <- spe_hd_ko_wt_cut_8um_labeled

      tb <- tibble(
        Tbx21 = SingleCellExperiment::counts(spe)["Tbx21", ],
        genotype = spe$group,
        cell_type = spe$cluster
      )
      tb <- bind_rows(
        tb,
        filter(tb, cell_type %in% c("Dark zone", "Light zone", "Germinal centre (unassigned)")) |>
          mutate(cell_type = "Germinal centre")
      )

      p <- tb |>
        mutate(
          genotype = fct_recode(genotype, "KO" = "ko", "WT" = "wt"),
          cell_type = factor(
              cell_type,
              levels = arrange(roast_table_tbet_hd_ko_wt_cut_8um, P.Value)$`cell-type`)
        ) |>
        group_by(genotype, cell_type) |>
        summarise(
          average_expression = mean(Tbx21),
          pct_expressing = sum(Tbx21 > 0) / n() * 100,
          .groups = "drop"
        ) |>
        ggplot(aes(x = cell_type, y = genotype, color = average_expression, size = pct_expressing)) +
        geom_point() +
        # Use a perceptually uniform color scale (e.g., Viridis or Magma)
        scale_color_viridis_c(name = "Avg Expression") +
        # Control the range of dot sizes
        scale_size_continuous(range = c(1, 10), name = "% Expressing") +
        theme_classic() +
        labs(
          title = paste(""),
          x = "Cell Type",
          y = "Genotype"
        ) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "right",      # Move to the left side
          legend.box = "horizontal"      # Place the two legends side-by-side
        )

      file_path <- file.path("output", "dotplot_Tbx21_hd_ko_wt_8um.pdf")
      ggsave(p, filename = file_path, width = 8, height = 3)
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c80g")
    )
  ),

  # venn / euler diagram
  tar_target(
    venn_dz_lz,
    format = "file",
    {
      p <- strings_to_matrix(
        `Victora et al` = bulk_dz_lz_signatures$symbol,
        `Not in Visium probe set` = setdiff(
          bulk_dz_lz_signatures$symbol,
          toupper(rownames(pseudo_bulk_counts_hd_ko_wt_cut_8um$Erythrocytes))
        ),
        `VisiumHD` =  topTable(efits_dz_lz_cut_8um$wt, p = 0.05, n = Inf) |>
          dplyr::pull(symbol) |>
          toupper()
      ) |>
        eulerr::euler() |>
        plot(
          quantities = TRUE,
          edges = FALSE,
          fills = list(fill = c("red", "grey", "steelblue4"), alpha = 0.5)
        )
      file_path <- file.path("output", "venn_dz_lz_8um.pdf")
      pdf(file_path, width = 6, height = 6)
      print(p)
      dev.off()
      file_path
    }
  )
)
