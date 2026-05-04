# Purpose:  Compute data outputs for Figure 3: segmentation quality comparison.
#           Tidies pre-computed summary metrics CSVs (cell count, area,
#           transcripts per cell, transcripts assigned) for MERSCOPE and Xenium.
#           Loads annotated Seurat objects for all six segmentation methods,
#           extracts UMAP embeddings, cell type compositions, per-sample MECR
#           scores, and negative marker purity scores (overall and per cell type).
# Inputs:   config/config.yaml
#             segmentation_comp.merscope_metrics_csv — pre-computed summary CSV
#             segmentation_comp.xenium_metrics_csv   — pre-computed summary CSV
#             segmentation_comp.scrna_path           — scGEM_seu.rds reference
#             segmentation_comp.scrna_celltype_col   — cell type column in scRNA obj
#             segmentation_comp.celltype_col         — cell type column in spatial objs
#             segmentation_comp.xenium_keep_samples  — sample IDs to retain from CSV
#           results/02_spatial_analysis/{method}/annotated_final.rds
#             xenium_batch34_default, xenium_batch34_cellpose, xenium_batch34_proseg
#             merscope_default,       merscope_cellpose,       merscope_proseg
# Outputs:  results/03_benchmarking/segmentation_quality/
#             metrics_long.rds      — long-format metrics (platform/method/sample/metric/value)
#             umap_coords.rds       — UMAP1/UMAP2 + cell type per platform/method/cell
#             cell_type_counts.rds  — per-sample cell type counts (n_cells)
#             mecr_table.rds        — per-sample MECR scores
#             purity_summary.rds    — overall purity + tx_assigned per platform/method
#             purity_ct.rds         — per-cell-type purity scores

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(yaml)
  library(optparse)
})

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--config"),  type = "character",
              default = "config/config.yaml",
              help    = "Path to config.yaml [default: %default]"),
  make_option(c("--out_dir"), type = "character",
              default = "results/03_benchmarking/segmentation_quality",
              help    = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg     <- yaml::read_yaml(opt$config)
seg_cfg <- cfg$segmentation_comp

if (is.null(seg_cfg)) {
  stop("config.yaml is missing the 'segmentation_comp' section. ",
       "Add it with keys: merscope_metrics_csv, xenium_metrics_csv, ",
       "scrna_path, scrna_celltype_col, celltype_col, xenium_keep_samples.")
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

celltype_col <- seg_cfg$celltype_col %||% "cell_type"
sample_col   <- "sample_id"

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Seurat v4/v5 compatible raw counts extraction
get_counts_mat <- function(seu, assay) {
  tryCatch(
    GetAssayData(seu, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(seu, assay = assay, slot = "counts")
  )
}

# Return the first available UMAP reduction from a ranked list of candidates
find_umap_reduction <- function(seu) {
  candidates <- c("umap.harmony_banksy", "umap_harmony", "umap")
  avail      <- Reductions(seu)
  match      <- intersect(candidates, avail)
  if (length(match) == 0) {
    stop("No recognised UMAP reduction found. Available: ", paste(avail, collapse = ", "))
  }
  match[1]
}

# Extract sparse count matrix as dgCMatrix
get_sparse <- function(seu, assay_use = NULL, slot_use = "counts") {
  if (is.null(assay_use)) assay_use <- DefaultAssay(seu)
  m <- get_counts_mat(seu, assay_use)
  if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
  m
}

# Resolve the sample_id column (pipeline uses sample_id; fall back to orig.ident)
get_sample_vector <- function(seu) {
  if (sample_col %in% colnames(seu@meta.data)) return(seu@meta.data[[sample_col]])
  if ("orig.ident" %in% colnames(seu@meta.data)) {
    message("  'sample_id' not found — using 'orig.ident'")
    return(seu@meta.data[["orig.ident"]])
  }
  stop("Neither 'sample_id' nor 'orig.ident' found in Seurat metadata.")
}

# ---------------------------------------------------------------------------
# MECR (Mixed-cell Expression Co-expression Rate) function
# Adapted from https://github.com/Center-for-Spatial-OMICs/SpatialQM
# Measures cross-cell-type co-expression of marker genes as a proxy for
# cell boundary bleed-through / merged cells. Higher = worse segmentation.
# ---------------------------------------------------------------------------
marker_df_spleen <- data.frame(
  gene = c(
    "Cd3d", "Cd3e", "Trac", "Cd4", "Cd8a", "Foxp3",
    "Cd19", "Cd22", "Cr2", "Ighd", "Fcer2a", "Cd72", "Bcr",
    "Bcl6", "Aicda", "Rgs13", "Cd83",
    "Jchain", "Cd38",
    "Csf1r", "Adgre1", "Cd68", "Ifi30", "Mmp12", "Cd209b", "Siglech",
    "Ngp", "S100a9", "Ffar2"
  ),
  cell_type = c(
    rep("T cells", 6),
    rep("B cells", 7),
    rep("GC cells", 4),
    rep("Plasma B cells", 2),
    rep("Macrophages", 7),
    rep("Granulocytes", 3)
  ),
  stringsAsFactors = FALSE
)
rownames(marker_df_spleen) <- marker_df_spleen$gene

calc_mecr <- function(seu, assay_use = NULL, slot_use = "counts",
                      marker_df = marker_df_spleen, max_genes = 25, seed = 1) {
  if (is.null(assay_use)) assay_use <- DefaultAssay(seu)
  exp   <- get_counts_mat(seu, assay_use)
  genes <- intersect(rownames(exp), marker_df$gene)

  if (length(genes) < 6) {
    return(list(MECR = NA_real_, n_genes_used = length(genes),
                note = "Too few marker genes."))
  }

  mtx <- as.matrix(exp[genes, , drop = FALSE])
  set.seed(seed)
  if (length(genes) > max_genes) genes <- sample(genes, max_genes)

  coexp_rates <- c()
  for (i in seq_along(genes)) {
    for (j in seq_len(i - 1)) {
      g1 <- genes[i]; g2 <- genes[j]
      if (marker_df[g1, "cell_type"] == marker_df[g2, "cell_type"]) next
      c1    <- mtx[g1, ]; c2 <- mtx[g2, ]
      denom <- sum((c1 > 0) | (c2 > 0))
      if (denom == 0) next
      coexp_rates <- c(coexp_rates, sum((c1 > 0) & (c2 > 0)) / denom)
    }
  }

  list(MECR = round(mean(coexp_rates), 3), n_genes_used = length(genes),
       n_pairs = length(coexp_rates))
}

calc_mecr_per_sample <- function(seu, assay_use) {
  samples <- unique(get_sample_vector(seu))
  objs    <- SplitObject(seu, split.by = sample_col)
  out <- lapply(names(objs), function(s) {
    res <- calc_mecr(objs[[s]], assay_use = assay_use)
    data.frame(sample = s, MECR = res$MECR, n_pairs = res$n_pairs,
               n_genes = res$n_genes_used)
  })
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# Negative marker purity (cell-level)
# Higher score = fewer cells expressing markers from the wrong cell type = better.
# ---------------------------------------------------------------------------
neg_marker_purity_cells <- function(seu_sp, seu_sc,
                                    key_sp = celltype_col,
                                    key_sc = seg_cfg$scrna_celltype_col,
                                    assay_sp = NULL, assay_sc = NULL,
                                    min_number_cells = 3,
                                    max_ratio_cells  = 0.005,
                                    pipeline_output  = TRUE) {
  Xsp <- get_sparse(seu_sp, assay_sp)
  Xsc <- get_sparse(seu_sc, assay_sc)

  common_genes <- intersect(rownames(Xsp), rownames(Xsc))
  Xsp <- Xsp[common_genes, , drop = FALSE]
  Xsc <- Xsc[common_genes, , drop = FALSE]

  sp_ct <- as.character(seu_sp[[key_sp]][, 1])
  sc_ct <- as.character(seu_sc[[key_sc]][, 1])

  shared  <- intersect(unique(sp_ct), unique(sc_ct))
  sp_tab  <- table(sp_ct)[shared]
  sc_tab  <- table(sc_ct)[shared]
  keep_ct <- shared[(sp_tab >= min_number_cells) & (sc_tab >= min_number_cells)]

  if (length(keep_ct) < 2) {
    if (pipeline_output) return(NaN)
    return(list(score = NaN, score_per_celltype = NULL))
  }

  Xsp <- Xsp[, sp_ct %in% keep_ct, drop = FALSE]
  Xsc <- Xsc[, sc_ct %in% keep_ct, drop = FALSE]
  sp_ct <- sp_ct[sp_ct %in% keep_ct]
  sc_ct <- sc_ct[sc_ct %in% keep_ct]

  # Fraction of positive cells per (gene × cell type)
  calc_ratio_pos <- function(X, ct_vec, cts) {
    out <- matrix(NA_real_, nrow = nrow(X), ncol = length(cts),
                  dimnames = list(rownames(X), cts))
    for (ct in cts) {
      idx <- which(ct_vec == ct)
      out[, ct] <- as.numeric(Matrix::rowSums(X[, idx, drop = FALSE] > 0)) / length(idx)
    }
    out
  }

  ratio_sc <- calc_ratio_pos(Xsc, sc_ct, keep_ct)
  ratio_sp <- calc_ratio_pos(Xsp, sp_ct, keep_ct)

  neg_mask <- ratio_sc < max_ratio_cells
  if (sum(neg_mask, na.rm = TRUE) < 1) {
    if (pipeline_output) return(NaN)
    return(list(score = NaN, score_per_celltype = NULL))
  }

  mean_sc <- mean(ratio_sc[neg_mask], na.rm = TRUE)
  mean_sp <- mean(ratio_sp[neg_mask], na.rm = TRUE)
  score   <- max(1 - (mean_sp - mean_sc), 0)

  if (pipeline_output) return(score)

  # Per cell-type score
  score_per_ct <- sapply(keep_ct, function(ct) {
    mask_ct <- neg_mask[, ct]
    if (sum(mask_ct, na.rm = TRUE) < 1) return(NA_real_)
    m_sc <- mean(ratio_sc[mask_ct, ct], na.rm = TRUE)
    m_sp <- mean(ratio_sp[mask_ct, ct], na.rm = TRUE)
    max(1 - (m_sp - m_sc), 0)
  })

  list(score = score, score_per_celltype = score_per_ct)
}

# ---------------------------------------------------------------------------
# 1. Tidy summary metrics CSVs
# ---------------------------------------------------------------------------
message("Loading and tidying summary metrics CSVs...")

# Shared metric renaming applied to both platforms
rename_metrics <- function(m) {
  dplyr::recode(m,
    "Cell count"                    = "Total cell count",
    "Cell Area - median (um^2)"     = "Cell Area (µm²)",
    "Transcripts within a cell (%)" = "Transcripts assigned (%)",
    "Transcripts per cell - median" = "Transcripts per cell"
  )
}

# Column pattern to pivot: "metric - method"
pivot_metrics_csv <- function(csv_path, vendor_label, default_label,
                              keep_samples = NULL) {
  df <- readr::read_csv(csv_path, na = c("", "NA", "N/A", "NULL"),
                        show_col_types = FALSE)

  if (!is.null(keep_samples)) {
    df <- dplyr::filter(df, Sample %in% keep_samples)
  }

  long <- df %>%
    dplyr::select(Sample, dplyr::matches(
      paste0("(Cell Area - median|Cell count|Transcripts per cell - median|",
             "Transcripts within a cell \\(%\\)) - (Cellpose|Proseg|", vendor_label, ")")
    )) %>%
    dplyr::filter(!is.na(Sample)) %>%
    tidyr::pivot_longer(
      cols         = -Sample,
      names_to     = c("metric", "method_raw"),
      names_pattern = paste0("^(.*) - (Cellpose|Proseg|", vendor_label, ")$"),
      values_to    = "value"
    ) %>%
    dplyr::mutate(
      method = dplyr::recode(method_raw, !!vendor_label := "Default"),
      method = factor(method, levels = c("Default", "Proseg", "Cellpose")),
      metric = rename_metrics(metric),
      metric = factor(metric, levels = c(
        "Total cell count", "Cell Area (µm²)",
        "Transcripts assigned (%)", "Transcripts per cell"
      )),
      # vendor CSVs store transcripts assigned as a ratio (0-1); convert to %
      value  = dplyr::if_else(metric == "Transcripts assigned (%)", value * 100, value)
    )

  long
}

mer_long <- pivot_metrics_csv(
  csv_path      = seg_cfg$merscope_metrics_csv,
  vendor_label  = "Vizgen",
  default_label = "Default"
) %>%
  dplyr::mutate(platform = "MERSCOPE")

xen_long <- pivot_metrics_csv(
  csv_path      = seg_cfg$xenium_metrics_csv,
  vendor_label  = "Xenium",
  default_label = "Default",
  keep_samples  = seg_cfg$xenium_keep_samples
) %>%
  dplyr::mutate(platform = "Xenium")

metrics_long <- dplyr::bind_rows(mer_long, xen_long) %>%
  dplyr::mutate(platform = factor(platform, levels = c("MERSCOPE", "Xenium")))

message("  Rows in metrics_long: ", nrow(metrics_long))
saveRDS(metrics_long, file.path(opt$out_dir, "metrics_long.rds"))
message("  Saved: metrics_long.rds")

# ---------------------------------------------------------------------------
# 2. Load annotated Seurat objects
# ---------------------------------------------------------------------------
message("Loading annotated Seurat objects...")

seg_methods <- list(
  xenium   = c("xenium_batch34_default", "xenium_batch34_cellpose", "xenium_batch34_proseg"),
  merscope = c("merscope_default",       "merscope_cellpose",       "merscope_proseg")
)

platform_assay <- list(xenium = "Xenium", merscope = "Vizgen")

# Build a flat list of all 6 objects with metadata
all_objs <- list()

for (platform in names(seg_methods)) {
  for (method in seg_methods[[platform]]) {
    rds_path <- file.path("results/02_spatial_analysis", method, "annotated_final.rds")
    if (!file.exists(rds_path)) {
      stop("annotated_final.rds not found: ", rds_path,
           "\n  Run the spatial_analysis pipeline for method '", method, "' first.")
    }
    message("  Loading: ", rds_path)
    seu <- readRDS(rds_path)

    # Standardise segmentation label: strip platform prefix from method name
    seg_label <- sub("^(xenium_batch34_|merscope_)", "", method)
    seg_label <- tools::toTitleCase(seg_label)
    seg_label <- dplyr::recode(seg_label, "Default" = "Default",
                               "Cellpose" = "Cellpose", "Proseg" = "Proseg")

    all_objs[[method]] <- list(
      seu       = seu,
      platform  = if (platform == "xenium") "Xenium" else "MERSCOPE",
      method    = method,
      seg_label = seg_label,
      assay     = platform_assay[[platform]]
    )
  }
}

# ---------------------------------------------------------------------------
# 3. Extract UMAP coordinates + cell type labels
# ---------------------------------------------------------------------------
message("Extracting UMAP coordinates...")

umap_list <- lapply(all_objs, function(obj) {
  seu    <- obj$seu
  red    <- find_umap_reduction(seu)
  coords <- as.data.frame(Embeddings(seu, reduction = red))
  colnames(coords) <- c("UMAP1", "UMAP2")

  ct_vec <- as.character(seu@meta.data[[celltype_col]])
  if (all(is.na(ct_vec))) {
    warning("'", celltype_col, "' is all NA in ", obj$method,
            ". Check the annotated_final.rds cell type column.")
  }

  coords %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::mutate(
      cell_type    = ct_vec,
      platform     = obj$platform,
      segmentation = obj$seg_label
    )
})

umap_coords <- dplyr::bind_rows(umap_list) %>%
  dplyr::mutate(
    platform     = factor(platform, levels = c("MERSCOPE", "Xenium")),
    segmentation = factor(segmentation, levels = c("Default", "Proseg", "Cellpose"))
  )

saveRDS(umap_coords, file.path(opt$out_dir, "umap_coords.rds"))
message("  Saved: umap_coords.rds")

# ---------------------------------------------------------------------------
# 4. Cell type composition counts
# ---------------------------------------------------------------------------
message("Computing cell type counts...")

ct_count_list <- lapply(all_objs, function(obj) {
  seu <- obj$seu
  tibble::tibble(
    platform     = obj$platform,
    segmentation = obj$seg_label,
    sample       = get_sample_vector(seu),
    cell_type    = as.character(seu@meta.data[[celltype_col]])
  )
})

cell_type_counts <- dplyr::bind_rows(ct_count_list) %>%
  dplyr::filter(!is.na(sample), !is.na(cell_type), cell_type != "") %>%
  dplyr::count(platform, segmentation, sample, cell_type, name = "n_cells") %>%
  dplyr::mutate(
    platform     = factor(platform, levels = c("MERSCOPE", "Xenium")),
    segmentation = factor(segmentation, levels = c("Default", "Proseg", "Cellpose"))
  )

saveRDS(cell_type_counts, file.path(opt$out_dir, "cell_type_counts.rds"))
message("  Saved: cell_type_counts.rds")

# ---------------------------------------------------------------------------
# 5. MECR per sample
# ---------------------------------------------------------------------------
message("Computing MECR per sample...")

mecr_list <- lapply(all_objs, function(obj) {
  message("  MECR: ", obj$method)
  df <- calc_mecr_per_sample(obj$seu, assay_use = obj$assay)
  dplyr::mutate(df,
    platform     = obj$platform,
    segmentation = obj$seg_label
  )
})

mecr_table <- dplyr::bind_rows(mecr_list) %>%
  dplyr::mutate(
    segmentation = factor(segmentation, levels = c("Default", "Proseg", "Cellpose")),
    platform     = factor(platform, levels = c("MERSCOPE", "Xenium"))
  )

saveRDS(mecr_table, file.path(opt$out_dir, "mecr_table.rds"))
message("  Saved: mecr_table.rds")

# ---------------------------------------------------------------------------
# 6. Negative marker purity
# ---------------------------------------------------------------------------
scrna_path <- seg_cfg$scrna_path
if (is.null(scrna_path) || !file.exists(scrna_path)) {
  message("scRNA reference not found at: ", scrna_path)
  message("Skipping purity computation. ",
          "Set segmentation_comp.scrna_path in config.yaml to enable this step.")
  # Write empty placeholders so downstream rules don't fail
  saveRDS(tibble::tibble(), file.path(opt$out_dir, "purity_summary.rds"))
  saveRDS(tibble::tibble(), file.path(opt$out_dir, "purity_ct.rds"))
  message("Done (purity skipped). Outputs written to: ", opt$out_dir)
  quit(save = "no", status = 0)
}

message("Loading scRNA reference: ", scrna_path)
sc_seu <- readRDS(scrna_path)

# Ensure the scRNA cell type column exists
scrna_ct_col <- seg_cfg$scrna_celltype_col %||% "Bcell_subcluster_combined"
if (!scrna_ct_col %in% colnames(sc_seu@meta.data)) {
  stop("Column '", scrna_ct_col, "' not found in scRNA metadata. ",
       "Update segmentation_comp.scrna_celltype_col in config.yaml.")
}

message("Computing negative marker purity (cells) per method...")

purity_summary_list <- list()
purity_ct_list      <- list()

for (method in names(all_objs)) {
  obj <- all_objs[[method]]
  message("  Purity: ", method)

  res <- tryCatch(
    neg_marker_purity_cells(
      seu_sp        = obj$seu,
      seu_sc        = sc_seu,
      key_sp        = celltype_col,
      key_sc        = scrna_ct_col,
      assay_sp      = obj$assay,
      assay_sc      = "RNA",
      pipeline_output = FALSE
    ),
    error = function(e) {
      message("    Error: ", conditionMessage(e))
      list(score = NaN, score_per_celltype = NULL)
    }
  )

  purity_summary_list[[method]] <- tibble::tibble(
    platform     = obj$platform,
    segmentation = obj$seg_label,
    purity       = res$score
  )

  if (!is.null(res$score_per_celltype)) {
    purity_ct_list[[method]] <- tibble::tibble(
      platform     = obj$platform,
      segmentation = obj$seg_label,
      cell_type    = names(res$score_per_celltype),
      purity       = as.numeric(res$score_per_celltype)
    )
  }
}

# ---------------------------------------------------------------------------
# 7. Join purity scores with median tx_assigned for dotplot
# ---------------------------------------------------------------------------
tx_assigned_summary <- metrics_long %>%
  dplyr::filter(metric == "Transcripts assigned (%)") %>%
  dplyr::group_by(platform, method) %>%
  dplyr::summarise(tx_assigned = median(value / 100, na.rm = TRUE), .groups = "drop") %>%
  dplyr::rename(segmentation = method)

purity_summary <- dplyr::bind_rows(purity_summary_list) %>%
  dplyr::mutate(
    segmentation = factor(segmentation, levels = c("Default", "Proseg", "Cellpose")),
    platform     = factor(platform, levels = c("MERSCOPE", "Xenium"))
  ) %>%
  dplyr::left_join(tx_assigned_summary, by = c("platform", "segmentation"))

purity_ct <- dplyr::bind_rows(purity_ct_list) %>%
  dplyr::mutate(
    segmentation = factor(segmentation, levels = c("Default", "Proseg", "Cellpose")),
    platform     = factor(platform, levels = c("MERSCOPE", "Xenium"))
  )

saveRDS(purity_summary, file.path(opt$out_dir, "purity_summary.rds"))
message("  Saved: purity_summary.rds")

saveRDS(purity_ct, file.path(opt$out_dir, "purity_ct.rds"))
message("  Saved: purity_ct.rds")

message("Done. All outputs written to: ", opt$out_dir)
