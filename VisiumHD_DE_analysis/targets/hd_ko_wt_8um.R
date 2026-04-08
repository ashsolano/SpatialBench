# Analysis targets for the Tbx21 KO vs wild-type Visium HD experiment.
# Four mouse spleen sections (2 KO: 167, 168; 2 WT: 709, 713) sequenced at 8 µm resolution.
# Each tar_target() is a node in the dependency graph; targets is responsible for caching
# and re-running only what is outdated.  See README.md for a narrative overview.

list(
  # ---------------------------------------------------------------------------
  # Probe set annotation
  # Parse the SpaceRanger probe set CSV to obtain gene symbol ↔ Ensembl ID
  # mappings used for deduplication and annotation downstream.
  # The raw CSV uses "|"-separated fields, so these are converted to commas first.
  # ---------------------------------------------------------------------------
  tar_target(probe_set_v2_df,
    deployment = "main",
    read_file("/stornext/Projects/score/Analyses/G000218_spatial_benchmarking_study_batch33/extdata/SpaceRanger/H1-7HZCTJJ/G000218_batch33_168/outs/probe_set.csv") |>
      str_replace_all("\\|", ",") |>
      str_replace("probe_id", "id,symbol,id2") |>
      read_csv(comment = "#")
  ),

  # ---------------------------------------------------------------------------
  # Per-sample preprocessing (tar_map expands to one set of targets per sample)
  # For each of the four samples:
  #   seu  — load the pre-processed Seurat RDS, switch to the 8 µm assay, and
  #           apply relative-counts (RC) normalisation (scale factor 5 000).
  #   svg  — identify spatially variable genes with Seurat FindVariableFeatures.
  #   spe  — convert Seurat to SpatialExperiment for Bioconductor compatibility.
  #   spe_banksy_neighborhood — augment each bin's expression vector with a
  #           weighted average of its k_geom=6 spatial neighbours (Banksy step 1).
  # ---------------------------------------------------------------------------
  tar_map(
    # set up the Seurat objects paths and names
    values = tibble::tibble(
      seurat_files = c(
        "/vast/projects/SpatialBench/VisiumHD/normalized_obj_08016um/batch33_167.rds",
        "/vast/projects/SpatialBench/VisiumHD/normalized_obj_08016um/batch33_168.rds",
        "/vast/projects/SpatialBench/VisiumHD/normalized_obj_08016um/batch33_709.rds",
        "/vast/projects/SpatialBench/VisiumHD/normalized_obj_08016um/batch33_713.rds"
      ),
      target_names = c("hd_ko_8um_167", "hd_ko_8um_168", "hd_wt_8um_709", "hd_wt_8um_713")
    ),
    names = "target_names",

    # preprocess the Seurat objects
    tar_target(
      seu,
      readRDS(seurat_files) |>
        (\(x) {
          SeuratObject::DefaultAssay(x) <- "Spatial.008um"
          x
        })() |>
        Seurat::NormalizeData(scale.factor = 5000, normalization.method = "RC"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "slurm_1c40g")
      )
    ),
    tar_target(
      svg,
      SeuratObject::VariableFeatures(Seurat::FindVariableFeatures(seu)),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "slurm_1c20g")
      )
    ),

    # convert Seurat objects to SpatialExperiment objects
    tar_target(
      spe,
      seu_to_spe_visium_hd(seu, sample_id = target_names, assay = "Spatial.008um", image = "slice1.008um"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "slurm_1c40g")
      )
    ),

    # Banksy
    tar_target(
      spe_banksy_neighborhood,
      Banksy::computeBanksy(
        spe[svgs_hd_ko_wt_8um, ],
        assay_name = "normcounts", compute_agf = FALSE, k_geom = 6
      ),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "slurm_4c120g")
      )
    )
  ),
  # ---------------------------------------------------------------------------
  # Union of spatially variable genes across all four samples.
  # Using the union (rather than intersection) maximises the feature set fed
  # into Banksy PCA, retaining sample-specific spatial patterns.
  # ---------------------------------------------------------------------------
  tar_target(
    svgs_hd_ko_wt_8um,
    deployment = "main",
    Reduce(union, list(
      svg_hd_ko_8um_168,
      svg_hd_ko_8um_167,
      svg_hd_wt_8um_709,
      svg_hd_wt_8um_713
    ))
  ),
  # ---------------------------------------------------------------------------
  # Joint Banksy PCA and UMAP across all four samples.
  # All Banksy neighbourhood SPEs are merged, a KO/WT group label is added,
  # and group-aware PCA (λ=0.2) followed by UMAP is run to embed all bins in
  # a shared low-dimensional space while accounting for genotype batch effects.
  # ---------------------------------------------------------------------------
  tar_target(
    spe_umap_hd_ko_wt_8um,
    cbind(
      spe_banksy_neighborhood_hd_ko_8um_167,
      spe_banksy_neighborhood_hd_ko_8um_168,
      spe_banksy_neighborhood_hd_wt_8um_709,
      spe_banksy_neighborhood_hd_wt_8um_713
    ) |>
      (\(x) {
        x$group <- ifelse(grepl("ko", x$sample_id), "ko", "wt")
        x
      })() |>
      Banksy::runBanksyPCA(use_agf = FALSE, lambda = 0.2, group = "group", seed = 1000) |>
      Banksy::runBanksyUMAP(use_agf = FALSE, lambda = 0.2, seed = 1000),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c1000g")
    )
  ),

  # ---------------------------------------------------------------------------
  # Merge all four per-sample SPEs (without Banksy augmentation) into one
  # combined object and add a KO/WT group label.
  # cue = tar_cue(mode = "never") freezes this target so it is not re-computed
  # even if upstream inputs change (used to avoid re-running expensive steps
  # after minor modifications that do not affect the merged counts).
  # ---------------------------------------------------------------------------
  # merge all spe objects
  tar_target(
    spe_hd_ko_wt_8um,
    cbind(
      spe_hd_ko_8um_167,
      spe_hd_ko_8um_168,
      spe_hd_wt_8um_709,
      spe_hd_wt_8um_713
    ) |>
      (\(x) {
        x$group <- ifelse(grepl("ko", x$sample_id), "ko", "wt")
        # x$cluster <- spe_banksy_hd_ko_wt_8um$clust_M0_lam0.2_k50_res0.7
        x
      })(),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c120g")
    ),
    # commented out adding banksy stuff
    # no need to re-run to get the same result without $cluster label
    cue = tar_cue(mode = "never")
  ),

  # ---------------------------------------------------------------------------
  # Data cleaning: spatial trimming, gene deduplication, and annotation.
  # 1. Spatial trim — sample 709 has an artefact; bins below
  #    the line y = 2.13x − 16000 in pixel space are excluded.
  # 2. Gene deduplication — some probe IDs map to the same symbol; their raw
  #    counts and normalised counts are summed and the duplicates are collapsed.
  # 3. Ensembl IDs are matched from the probe set; two missing entries are
  #    patched manually (H2-M10, Tex19).
  # 4. Chromosome names are fetched from Ensembl via biomaRt, needed for the
  #    sex-linked gene set test downstream.
  # ---------------------------------------------------------------------------
  # remove low quality area
  # add chromosome names
  tar_target(spe_hd_ko_wt_cut_8um,
    (\(spe, probe_set) {
      tb <- spe |>
        spatialCoords() |>
        cbind(colData(spe))
      # spe$cluster_first <- spe$cluster
      # spe$cluster <- as.character(spe$cluster)
      # spe$cluster[spe$cluster %in% c(19)] <-
      #   paste0("GC-", spe_banksy_gc_hd_ko_wt_8um$clust_M0_lam0.2_k50_res0.6)
      spe <- spe[, tb$sample_id != "hd_wt_8um_709" | tb$y >= 2.13 * tb$x - 16000]

      # merge rows with the same symbol
      rowData(spe)$symbol <- rownames(spe) |>
        stringr::str_remove("\\.\\d+$")
      duplicated_rows <- rowData(spe)$symbol %in% rowData(spe)$symbol[duplicated(rowData(spe)$symbol)]
      dup_spe <- spe[duplicated_rows, ]
      rownames(dup_spe) <- rowData(dup_spe)$symbol
      dup_counts <- rowsum(SingleCellExperiment::counts(dup_spe), group = rowData(dup_spe)$symbol)
      dup_nromcounts <- rowsum(SingleCellExperiment::normcounts(dup_spe), group = rowData(dup_spe)$symbol)
      dup_spe <- SingleCellExperiment(
        assays = list(counts = dup_counts, normcounts = dup_nromcounts),
        rowData = rowData(dup_spe)[rownames(dup_counts), , drop = FALSE],
        colData = colData(dup_spe)
      )
      spe <- rbind(
        spe[!duplicated_rows, ],
        dup_spe
      )

      # add EnsembleID
      rowData(spe)$EnsembleID <- probe_set$gene_id[match(rownames(spe), probe_set$symbol)]
      rowData(spe[c("H2-M10", "Tex19"), ])$EnsembleID <- c("ENSMUSG00000023083", "ENSMUSG00000039337")

      # add chromosome names
      mart_table <-
        use_mart() |>
        biomaRt::select(
          keys = rowData(spe)$EnsembleID,
          keytype = "ensembl_gene_id",
          columns = c("ensembl_gene_id", "chromosome_name")
        )
      rowData(spe)$chromosome_name <- "NotFound"
      rowData(spe)[
        match(mart_table$ensembl_gene_id, rowData(spe)$EnsembleID),
        "chromosome_name"
      ] <- mart_table$chromosome_name

      spe
    })(spe_hd_ko_wt_8um, probe_set_v2_df),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c120g")
    ),
    # commented out adding banksy stuff
    # no need to re-run to get the same result without $cluster label
    cue = tar_cue(mode = "never")
  ),

  # ---------------------------------------------------------------------------
  # Single-cell reference preparation for RCTD deconvolution.
  # A Flex single-cell RNA-seq annoteted reference is loaded
  # and log-normalised.  SCDC additionally requires an ExpressionSet (eset).
  # ---------------------------------------------------------------------------
  # Deconvolution with RCTD using Flex single-cell reference
  tar_target(
    yining_flex_sce,
    {
      seu <- readRDS("/vast/projects/SpatialBench/scGEM_subcluster/scFlex_seu.rds")
      # Layers have different number of genes
      # need to be consistent to convert to SCE
      SeuratObject::LayerData(seu, "counts") <-
        SeuratObject::LayerData(seu, "counts")[rownames(SeuratObject::LayerData(seu, "data")), ]
      Seurat::as.SingleCellExperiment(seu) |>
        scuttle::computeLibraryFactors() |>
        scuttle::logNormCounts()
    }
  ),
  # SCDC wants eset object
  tar_target(
    yining_flex_eset,
    {
      mtx <- SingleCellExperiment::counts(yining_flex_sce)
      normed_mtx <- scuttle::normalizeCounts(
        mtx,
        size.factors = scuttle::librarySizeFactors(mtx),
        log = FALSE
      )
      SCDC::getESET(
        normed_mtx,
        fdata = rownames(yining_flex_sce),
        pdata = cbind(
          sample = colData(yining_flex_sce)$SampleID,
          # factors will give error for reasons beyond me
          cell_type = as.character(colData(yining_flex_sce)$cluster_label)
        )
      )
    }
  ),
  # ---------------------------------------------------------------------------
  # RCTD cell type deconvolution (spacexr, doublet mode).
  # Each 8 µm bin is assigned a primary cell type from the single-cell
  # reference.  run_RCTD() (R/clustering.R) splits by sample, runs RCTD
  # per-sample in parallel, then reassembles in the original bin order.
  # ---------------------------------------------------------------------------
  tar_target(
    rctd_results_hd_ko_wt_8um,
    run_RCTD(
      spe_hd_ko_wt_cut_8um,
      yining_flex_sce,
      cell_type_col = "cluster_label",
      max_cores = 32
    ),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_32c500g")
    )
  ),
  # ---------------------------------------------------------------------------
  # GC B cell subclustering to resolve dark zone (DZ) and light zone (LZ).
  # The single-cell reference does not distinguish DZ and LZ GC B cells, so
  # these are identified spatially:
  #   1. Bins classified as "GC B cells" by RCTD are extracted separately for
  #      KO and WT (keeping only the spatially trimmed region of sample 709).
  #   2. Banksy is re-run at resolution 0.5 on each subset to produce fine
  #      spatial subclusters.
  #   3. Subcluster-to-zone assignments (below) are based on expression of
  #      DZ/LZ marker genes.
  # ---------------------------------------------------------------------------
  # do not have DZ LZ in reference, use Banksy to cluster on GC B cells
  # to annotate DZ and LZ
  tar_target(
    rctd_hd_ko_wt_8um_ko_gc_subclusters,
    {
      spe <- spe_umap_hd_ko_wt_8um
      tb <- spe |>
        spatialCoords() |>
        cbind(colData(spe))
      spe <- spe[, tb$sample_id != "hd_wt_8um_709" | tb$y >= 2.13 * tb$x - 16000]

      colnames(spe) <- seq_len(ncol(spe))
      subset_cells <- rctd_results_hd_ko_wt_8um$colData$spot_class != "reject" &
        rctd_results_hd_ko_wt_8um$colData$first_type == "GC B cells" &
        spe$group == "ko"
      spe <- spe[, subset_cells]
      Banksy::clusterBanksy(
        spe,
        use_agf = FALSE, lambda = 0.2, resolution = 0.5, seed = 1000
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c500g")
    )
  ),
  tar_target(
    rctd_hd_ko_wt_8um_wt_gc_subclusters,
    {
      spe <- spe_umap_hd_ko_wt_8um
      tb <- spe |>
        spatialCoords() |>
        cbind(colData(spe))
      spe <- spe[, tb$sample_id != "hd_wt_8um_709" | tb$y >= 2.13 * tb$x - 16000]

      colnames(spe) <- seq_len(ncol(spe))
      subset_cells <- rctd_results_hd_ko_wt_8um$colData$spot_class != "reject" &
        rctd_results_hd_ko_wt_8um$colData$first_type == "GC B cells" &
        spe$group == "wt"
      spe <- spe[, subset_cells]
      Banksy::clusterBanksy(
        spe,
        use_agf = FALSE, lambda = 0.2, resolution = 0.5, seed = 1000
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c500g")
    )
  ),
  # ---------------------------------------------------------------------------
  # Merge GC subcluster DZ/LZ labels back into the full RCTD label vector.
  # GC B cells not assigned to a DZ/LZ subcluster (with mixed / unclear signatures)
  # are labelled "Germinal centre (unassigned)".
  # ---------------------------------------------------------------------------
  # combine subcluster labels with RCTD labels
  tar_target(
    rctd_subclustered_labels_hd_ko_wt_8um,
    {
      labels <- rctd_results_hd_ko_wt_8um$colData$first_type
      #labels[as.numeric(colnames(rctd_hd_ko_wt_8um_wt_gc_subclusters))] <-
      #  as.characte(rctd_hd_ko_wt_8um_wt_gc_subclusters$clust_M0_lam0.2_k50_res0.3)
      subclusters <- tibble(
        idx = as.numeric(colnames(rctd_hd_ko_wt_8um_wt_gc_subclusters)),
        cluster = as.character(rctd_hd_ko_wt_8um_wt_gc_subclusters$clust_M0_lam0.2_k50_res0.5),
        group = "wt"
      ) |>
        bind_rows(
          tibble(
            idx = as.numeric(colnames(rctd_hd_ko_wt_8um_ko_gc_subclusters)),
            cluster = as.character(rctd_hd_ko_wt_8um_ko_gc_subclusters$clust_M0_lam0.2_k50_res0.5),
            group = "ko"
          )
        ) |>
        mutate(
          annotation = case_when(
            group == "ko" & cluster %in% c("1", "4") ~ "Light zone",
            group == "ko" & cluster %in% c("2") ~ "Dark zone",
            group == "wt" & cluster %in% c("2") ~ "Light zone",
            group == "wt" & cluster %in% c("1", "3") ~ "Dark zone",
            .default = "Germinal centre (unassigned)"
          )
        )
      labels[subclusters$idx] <- subclusters$annotation
      # there were GC B cells not assigned (reject by RCTD)
      # rename GC B cells to Germinial centre (unassigned)
      labels <- forcats::fct_recode(
        labels,
        "Germinal centre (unassigned)" = "GC B cells"
      )
      labels
    }
  ),

  # ---------------------------------------------------------------------------
  # Assemble the final annotated SPE.
  # RCTD results (spot_class, deconvolution weights) and the final
  # cluster labels (DZ/LZ-resolved) are added to colData.  The Banksy UMAP
  # embedding is transferred from spe_umap_hd_ko_wt_8um.
  # ---------------------------------------------------------------------------
  # new spe object with cluster labels and QC
  tar_target(spe_hd_ko_wt_cut_8um_labeled,
    {
      spe <- spe_hd_ko_wt_cut_8um |>
        scuttle::addPerCellQC(
          subsets = list(
            mt = str_detect(rownames(spe_hd_ko_wt_cut_8um), "^mt-")
          )
        )

      colData(spe) <- cbind(
        colData(spe),
        rctd_results_hd_ko_wt_8um$colData
      )

      spe$cluster <- rctd_subclustered_labels_hd_ko_wt_8um
      spe$filter <- spe$sum >= 175 & spe$subsets_mt_percent < 2 & spe$spot_class != "reject"
      # spe$filter <- spe$sum > 0

      spe_umap <- spe_umap_hd_ko_wt_8um
      tb <- spe_umap |>
        spatialCoords() |>
        cbind(colData(spe_umap))
      spe_umap <- spe_umap[, tb$sample_id != "hd_wt_8um_709" | tb$y >= 2.13 * tb$x - 16000]

      stopifnot(all(colnames(spe) == colnames(spe_umap)))
      reducedDim(spe, "UMAP") <- reducedDim(spe_umap, "UMAP_M0_lam0.2")
      spe
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c120g")
    )
  ),

  # ---------------------------------------------------------------------------
  # Marker gene lists for heatmap annotation.
  # zone_markers_df  — curated markers for each splenic zone/cell type,
  #                    used to confirm cell type assignments.
  # dark_light_zone_markers — DZ and LZ marker genes
  # ---------------------------------------------------------------------------
  # marker heatmap
  tar_target(zone_markers_df, list(
    "Macrophage" = c("Cd274", "Marco", "Csf1r", "Adgre1", "Cd209b", "Cd206", "Cd80", "Mac1", "Cd68")[c(1, 3, 4, 5, 7, 9)],
    "B cell" = c("Cd19", "Cd22", "Ighd", "Cd5"),
    "Germinal centre" = c("Cxcr4", "Cd83", "Bcl6", "Rgs13", "Aicda"),
    "Neutrophil" = c("S100a9", "S100a8", "Ngp"),
    "Erythrocyte" = c("Car2", "Car1", "Klf1"),
    "Plasma cell" = c("Cd38", "Cd138", "Xbp1", "Irf4", "Prdm1", "Cd27", "Cd319", "Mum1")[c(1, 3, 4, 5, 6)],
    "T cell" = c("Trac", "Cd3d", "Cd4", "Cd3e", "Cd8a"),
    "Red pulp" = c("Ifitm3", "C1qc", "Hmox1", "Hba-a1", "Klf1"),
    "Marginal zone" = c("Marco", "Lyz2", "Ighd", "Igfbp7", "Igfbp3", "Ly6d"),
    "White pulp" = c("Ighd", "Cd19", "Trac", "Trbc2")
  ) |>
    stack() |>
    setNames(c("gene", "zone")) |>
    merge(genes_hd_ko_wt_cut_8um[, c("symbol", "EnsembleID")],
      by.x = "gene",
      by.y = "symbol"
    ) |>
    na.omit()
  ),

  tar_target(
    dark_light_zone_markers,
    tibble::tribble(
      ~gene, ~zone,
      "Rrm1", "dark_zone",
      "Cdca8", "dark_zone",
      "Lmnb1", "dark_zone",
      "Gpsm2", "dark_zone",
      "Nek2", "dark_zone",
      "Pfn2", "dark_zone",
      "Reln", "dark_zone",
      "Anxa2", "dark_zone",
      "Lmo4", "dark_zone",
      "Polh", "dark_zone",
      "Akap12", "dark_zone",
      "Ube2h", "dark_zone",
      "Otub2", "dark_zone",
      "Scn8a", "dark_zone",
      "Tifa", "dark_zone",
      # "Actb", "light_zone",
      "Cd38", "light_zone",
      "Stx11", "light_zone",
      "Nfkbie", "light_zone",
      "Samsn1", "light_zone",
      "Serpinb6b", "light_zone",
      "Cd83", "light_zone",
      "Cd86", "light_zone",
      "Cd40", "light_zone",
      "Gstt2", "light_zone",
      "B3gnt5", "light_zone",
      # "Gas5", "light_zone",
      "Fcer2a", "light_zone",
      "Ankrd33b", "light_zone",
      "Ncf1", "light_zone",
      "Nfkbid", "light_zone"
    ) |>
      mutate(EnsembleID = gene)
  ),

  # ---------------------------------------------------------------------------
  # Pseudo-bulk differential expression: KO vs WT
  # Step 1 — Aggregate raw counts per cell type × sample (pseudo_bulk_counts).
  #   GC "Germinal centre" counts are the sum of DZ + LZ + unassigned bins.
  # Step 2 — For each cell type (vooms): filter low-expressed genes, TMM-
  #   normalise, and apply voom to estimate precision weights.
  # Step 3 — Fit a means model (efits) and test the KO vs WT contrast with
  #   limma::eBayes and an empirical Bayes shrinkage of gene-wise variances.
  # ---------------------------------------------------------------------------
  # DE
  tar_target(
    pseudo_bulk_counts_hd_ko_wt_cut_8um,
    lapply(
      unique(rctd_subclustered_labels_hd_ko_wt_8um),
      function(cluster) {
        spe <- spe_hd_ko_wt_cut_8um_labeled[, spe_hd_ko_wt_cut_8um_labeled$cluster == cluster]
        # spe <- spe[, spe$filter]
        SingleCellExperiment::counts(spe) |>
        as("CsparseMatrix") |>
        S4Arrays::colsum(group = spe$sample_id, reorder = TRUE)
      }
    ) |>
      setNames(unique(rctd_subclustered_labels_hd_ko_wt_8um)) |>
      (\(x) {
        x[["Germinal centre"]] <- x[["Light zone"]] + x[["Dark zone"]] + x[["Germinal centre (unassigned)"]]
        x
      })(),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c80g")
    )
  ),

  tar_target(
    genes_hd_ko_wt_cut_8um,
    rowData(spe_hd_ko_wt_cut_8um),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c40g")
    )
  ),

  tar_target(
    vooms_hd_ko_wt_cut_8um,
    sapply(
      pseudo_bulk_counts_hd_ko_wt_cut_8um,
      function(mtx) {
        dgelist <- edgeR::DGEList(
          counts = mtx,
          samples = colnames(mtx),
          genes = genes_hd_ko_wt_cut_8um,
          group = factor(ifelse(
            grepl("ko", colnames(mtx)),
            "ko",
            "wt"
          ))
        )
        dgelist <- dgelist[
          filterByExpr(dgelist), , # bias?
          keep.lib.sizes = FALSE
        ]
        dgelist <- calcNormFactors(dgelist, method = "TMM")
        # means model
        design <- model.matrix(~ 0 + dgelist$samples$group, data = dgelist$samples)
        colnames(design) <- gsub(".*\\$", "", colnames(design)) |>
          gsub(" .*$", "", x = _)
        dge_v <- voom(dgelist, design, save.plot = T, plot = F, span = 0.2)
        return(dge_v)
      },
      simplify = FALSE
    ),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c40g")
    )
  ),
  tar_target(
    efits_hd_ko_wt_cut_8um,
    sapply(names(vooms_hd_ko_wt_cut_8um),
      function(x) {
        efit <- lmFit(vooms_hd_ko_wt_cut_8um[[x]], vooms_hd_ko_wt_cut_8um[[x]]$design) |>
          contrasts.fit(contrasts = makeContrasts(
            KovsWt = "groupko - groupwt",
            levels = vooms_hd_ko_wt_cut_8um[[x]]$design
          )) |>
          eBayes()
        efit$voom.line <- vooms_hd_ko_wt_cut_8um[[x]]$voom.line
        efit$voom.xy <- vooms_hd_ko_wt_cut_8um[[x]]$voom.xy
        efit$cluster <- x
        return(efit)
      },
      simplify = FALSE
    )
  ),
  # ---------------------------------------------------------------------------
  # Gene set tests — T-bet target genes (ROAST)
  # prev_de: significant DE genes from a prior bulk RNA-seq study of Tbx21 KO
  #   GC B cells (data/GCB-DE-D15.TbetKO_GCB.v.WT_GCB.csv), used as a
  #   reference gene set.  Gene weights are the t-statistics from that study,
  #   so ROAST tests whether the same directional signal appears in each
  #   cell type's KO vs WT comparison here.
  # roasts_tbet: rotation-based gene set test (10^5 rotations) applied per
  #   cell type pseudo-bulk.
  # ---------------------------------------------------------------------------
  tar_target(
    prev_de,
    read_csv("data/GCB-DE-D15.TbetKO_GCB.v.WT_GCB.csv") |>
      filter(adj.P.Val < 0.05)
  ),
  tar_target(
    roasts_tbet_hd_ko_wt_cut_8um,
    sapply(names(efits_hd_ko_wt_cut_8um),
      function(x) {
        roast(vooms_hd_ko_wt_cut_8um[[x]],
          design = efits_hd_ko_wt_cut_8um[[x]]$design,
          index = efits_hd_ko_wt_cut_8um[[x]]$gene %>%
            mutate(idx = seq.int(nrow(.))) %>%
            filter(symbol %in% prev_de$Symbol) %>%
            pull(idx),
          gene.weights = efits_hd_ko_wt_cut_8um[[x]]$gene %>%
            filter(symbol %in% prev_de$Symbol) %>%
            left_join(prev_de, by = c("symbol" = "Symbol")) %>%
            pull(t),
          contrast = efits_hd_ko_wt_cut_8um[[x]]$contrasts,
          nrot = 99999
        )
      },
      simplify = FALSE
    )
  ),
  # ---------------------------------------------------------------------------
  # Gene set test — sex-linked genes (ROAST)
  # xie_genes: X-inactivation-escape genes found in mouse spleen from Berletch et al. 2015
  #   (https://doi.org/10.1371/journal.pgen.1005079).
  # Weights: Y-chromosome genes = +1, XiE genes = −1
  # ---------------------------------------------------------------------------
  tar_target(
    xie_genes,
    c(
      "Cybb", "Ddx3x", "Kdm6a", "Cfp", "Utp14a", "Firre", "Bgn", "5430427O19Rik",
      "Eif2s3x", "Vsig4", "Xist", "Ftx", "5530601H04Rik", "Pbdc1", "5730416F02Rik",
      "Kdm5c", "Tmsb4x"
    )
  ),
  tar_target(
    roasts_sex_hd_ko_wt_cut_8um,
    sapply(names(efits_hd_ko_wt_cut_8um),
      function(x) {
        roast(vooms_hd_ko_wt_cut_8um[[x]],
          design = efits_hd_ko_wt_cut_8um[[x]]$design,
          index = efits_hd_ko_wt_cut_8um[[x]]$gene %>%
            mutate(idx = seq.int(nrow(.))) %>%
            filter(chromosome_name == "Y" | symbol %in% xie_genes) %>%
            pull(idx),
          gene.weights = efits_hd_ko_wt_cut_8um[[x]]$gene %>%
            mutate(weight = ifelse(chromosome_name == "Y", 1, -1)) %>%
            pull(weight),
          contrast = efits_hd_ko_wt_cut_8um[[x]]$contrasts,
          nrot = 99999
        )
      },
      simplify = FALSE
    )
  ),
  tar_target(
    roast_table_tbet_hd_ko_wt_cut_8um,
    sapply(roasts_tbet_hd_ko_wt_cut_8um, simplify = FALSE, function(x) {
      as_tibble(x$p.value["Up", ])
    }) |>
      bind_rows(.id = "cell-type") |>
      mutate(Number.of.genes = sapply(roasts_tbet_hd_ko_wt_cut_8um, function(x) {
        x$ngenes.in.set
      }))
  ),
  tar_target(
    roast_table_sex_hd_ko_wt_cut_8um,
    sapply(roasts_sex_hd_ko_wt_cut_8um, simplify = FALSE, function(x) {
      as_tibble(x$p.value["Up", ])
    }) |>
      bind_rows(.id = "cell-type") |>
      mutate(Number.of.genes = sapply(roasts_sex_hd_ko_wt_cut_8um, function(x) {
        x$ngenes.in.set
      }))
  ),

  # ---------------------------------------------------------------------------
  # Pseudo-bulk differential expression: Dark Zone vs Light Zone
  # Performed separately within KO and WT groups
  # The same voom-limma pipeline as for KO vs WT is applied; the contrast is
  # DZ vs LZ within each genotype.
  # Results are compared against a previous bulk DZ/LZ signature
  # (Victora et al., 10.1182/blood-2012-03-415380) via a weighted ROAST test.
  # ---------------------------------------------------------------------------
  # DZ vs. LZ
  tar_target(
    vooms_dz_lz_cut_8um,
    sapply(
      c("ko", "wt"),
      function(group) {
        mtx <- cbind(
          pseudo_bulk_counts_hd_ko_wt_cut_8um[["Dark zone"]],
          pseudo_bulk_counts_hd_ko_wt_cut_8um[["Light zone"]]
        )
        colnames(mtx) <- paste0(
          colnames(mtx), "_",
          c(rep("DZ", ncol(pseudo_bulk_counts_hd_ko_wt_cut_8um[["Dark zone"]])),
            rep("LZ", ncol(pseudo_bulk_counts_hd_ko_wt_cut_8um[["Light zone"]]))
          )
        )
        mtx <- mtx[, grepl(group, colnames(mtx))]
        dgelist <- edgeR::DGEList(
          counts = mtx,
          samples = colnames(mtx),
          genes = genes_hd_ko_wt_cut_8um,
          group = factor(ifelse(
            grepl("DZ", colnames(mtx)),
            "DZ",
            "LZ"
          ))
        )
        dgelist <- dgelist[
          filterByExpr(dgelist), , # bias?
          keep.lib.sizes = FALSE
        ]
        dgelist <- calcNormFactors(dgelist, method = "TMM")
        # means model
        design <- model.matrix(~ 0 + dgelist$samples$group, data = dgelist$samples)
        colnames(design) <- gsub(".*\\$", "", colnames(design)) |>
          gsub(" .*$", "", x = _)
        dge_v <- voom(dgelist, design, save.plot = T, plot = F, span = 0.2)
        return(dge_v)
      },
      simplify = FALSE
    ),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c40g")
    )
  ),
  tar_target(
    efits_dz_lz_cut_8um,
    sapply(names(vooms_dz_lz_cut_8um),
      function(x) {
        efit <- lmFit(vooms_dz_lz_cut_8um[[x]], vooms_dz_lz_cut_8um[[x]]$design) |>
          contrasts.fit(contrasts = makeContrasts(
            DZvsLZ = "groupDZ - groupLZ",
            levels = vooms_dz_lz_cut_8um[[x]]$design
          )) |>
          eBayes()
        efit$voom.line <- vooms_dz_lz_cut_8um[[x]]$voom.line
        efit$voom.xy <- vooms_dz_lz_cut_8um[[x]]$voom.xy
        efit$cluster <- x
        return(efit)
      },
      simplify = FALSE
    )
  ),
  tar_target(
    bulk_dz_lz_signatures,
    readxl::read_xlsx("data/Gabriel_etal_supp.xlsx", skip = 4) |>
      dplyr::filter(`Gene Symbol` != "---") |>
      dplyr::summarize(
        z = -mean(`Z-score`),
        .by = `Gene Symbol`
      ) |>
      dplyr::rename(
        symbol = `Gene Symbol`
      ) |>
      dplyr::mutate(
        symbol = sapply(str_split(symbol, " /// "), function(x) x[1])
      )
  ),
  tar_target(
    dz_lz_signatures,
    read_csv("data/DZLZ_signature.csv") |>
    # Cc56: not a valid gene symbol, meant Cd56?
    #   Cd56 not in the count matrix
    # Reln was meant to be removed
    # Gas5: not present in count matrix
      filter(!Gene %in% c("Reln", "Cd56", "Gas5")) |>
      dplyr::rename_with(tolower) |>
      dplyr::mutate(EnsembleID = NA)
  ),
  tar_target(
    roasts_dz_lz_hd_ko_wt_cut_8um,
    lapply(
      names(vooms_dz_lz_cut_8um),
      function(x) {
        roast(vooms_dz_lz_cut_8um[[x]],
          design = efits_dz_lz_cut_8um[[x]]$design,
          index = efits_dz_lz_cut_8um[[x]]$gene %>%
            mutate(idx = seq.int(nrow(.))) %>%
            filter(toupper(symbol) %in% bulk_dz_lz_signatures$symbol) %>%
            pull(idx),
          gene.weights = efits_dz_lz_cut_8um[[x]]$gene %>%
            filter(toupper(symbol) %in% bulk_dz_lz_signatures$symbol) %>%
            mutate(symbol = toupper(symbol)) %>%
            left_join(bulk_dz_lz_signatures, by = "symbol") %>%
            pull(z),
          contrast = efits_dz_lz_cut_8um[[x]]$contrasts,
          nrot = 99999
        )
      }) |>
      setNames(names(vooms_dz_lz_cut_8um))
  ),
  tar_target(
    roast_table_dz_lz_hd_ko_wt_cut_8um,
    sapply(roasts_dz_lz_hd_ko_wt_cut_8um, function(x) {
      x$p.value["Up", ]
    }) |>
      t() |>
      cbind(Number.of.genes = sapply(roasts_dz_lz_hd_ko_wt_cut_8um, function(x) {
        x$ngenes.in.set
      }))
  ),
  # ---------------------------------------------------------------------------
  # Export DE results to Excel.
  # dz_lz_cut_8um_xlsx  — DZ vs LZ top table (FDR < 5 %), annotated with
  #   the published DZ/LZ z-score and a column indicating directional agreement.
  # hd_ko_wt_cut_8um_xlsx — KO vs WT top table per cell type (FDR < 20 %),
  #   with the prior bulk RNA-seq logFC appended for comparison.
  # ---------------------------------------------------------------------------
  tar_target(
    dz_lz_cut_8um_xlsx,
    format = "file",
    {
      output_file <- "output/dz_lz_cut_8um.xlsx"
      lapply(
        efits_dz_lz_cut_8um,
        function(x) {
          limma::topTable(
            x,
            n = Inf,
            p.value = 0.05
          ) |>
            mutate(symbol_upper = toupper(symbol)) |> 
            left_join(bulk_dz_lz_signatures, by = c('symbol_upper' = 'symbol')) |>
            mutate(agree = (logFC > 0) == (z > 0))
        }
      ) |>
      writexl::write_xlsx(
        path = output_file,
        format_headers = TRUE,
        col_names = TRUE
      )
      output_file
    } 
  ),
  tar_target(
    hd_ko_wt_cut_8um_xlsx,
    format = "file",
    {
      output_file <- "output/hd_ko_wt_cut_8um.xlsx"
      lapply(
        efits_hd_ko_wt_cut_8um,
        function(x) {
          limma::topTable(
            x,
            n = Inf,
            p.value = 0.2
          ) |>
            left_join(
              dplyr::rename(
                dplyr::select(prev_de, Symbol, logFC, adj.P.Val),
                prev_logFC = logFC, prev_adj.P.Val = adj.P.Val
              ),
              by = c("symbol" = "Symbol")
            )
        }
      ) |>
      writexl::write_xlsx(
        path = output_file,
        format_headers = TRUE,
        col_names = TRUE
      )
      output_file
    }
  ),

  # ---------------------------------------------------------------------------
  # Sanity checks
  # Cell type proportion validation with SCDC.
  # SCDC::SCDC_prop deconvolves the pseudo-bulk count matrices using the same
  # single-cell reference, providing an independent estimate of cell type
  # composition to cross-check the spatial RCTD assignments.
  # Applying SCDC to pseudo-bulk matrices derived from RCTD labels is
  # somewhat circular; this is used as a sanity check only.
  # ---------------------------------------------------------------------------
  # run SCDC for each pseudo-bulk sample
  tar_target(
    scdc_props_hd_ko_wt_cut_8um,
    {
      combined_mat <- do.call(
        cbind,
        lapply(names(pseudo_bulk_counts_hd_ko_wt_cut_8um), function(nm) {
          mat <- pseudo_bulk_counts_hd_ko_wt_cut_8um[[nm]]
          colnames(mat) <- stringr::str_extract(colnames(mat), "\\d+$")
          colnames(mat) <- paste(nm, colnames(mat), sep = "_")
          mat
        })
      )
      normed_mat <- scuttle::normalizeCounts(
        combined_mat,
        size.factors = scuttle::librarySizeFactors(combined_mat),
        log = FALSE
      )
      bulk_eset <- SCDC::getESET(
        normed_mat,
        fdata = rownames(normed_mat),
        pdata = colnames(normed_mat)
      )
      SCDC::SCDC_prop(
        bulk.eset = bulk_eset,
        sc.eset = yining_flex_eset,
        sample = "sample",
        ct.varname =  "cell_type",
        ct.sub = unique(Biobase::pData(yining_flex_eset)$cell_type)
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_1c20g")
    )
  ),

  # heatmap of deconvolution results
  tar_target(
    scdc_props_heatmap,
    # convert matrix to long format tibble
    scdc_props_hd_ko_wt_cut_8um$prop.est.mvw |>
      as_tibble(rownames = "sample") |>
      pivot_longer(
        cols = -sample,
        names_to = "cell_type",
        values_to = "proportion"
      ) |>
      mutate(
        sample = factor(
          sample, levels = sort(unique(sample))
        ),
        cell_type = factor(
          cell_type,
          levels = sort(unique(cell_type))
        )
      ) |>
      ggplot(aes(x = cell_type, y = sample, fill = proportion)) +
      geom_tile() +
      # scale_fill_viridis_c() +
      scale_fill_gradient(low = "white", high = "red") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        title = "SCDC Deconvolution Proportions",
        x = "Cell Type",
        y = "Spatial Transcriptomics Pseudo-bulk Samples",
        fill = "Proportion"
      )
  )
)
