# Purpose:  Pseudobulk log10(CPM + 1) and Pearson r between single-cell
#           references (10X FLEX, 10X 3' GEX) and VisiumHD, MERSCOPE and Xenium
#           (8µm bins), matched WT animals (cfg$scrna$wt_ids) only. Common-gene
#           values are the pairwise values subset to the gene_lists.rds
#           intersection (not renormalised).
# Inputs:   config/config.yaml; cfg$scrna$path (FLEX), cfg$scrna$sc_path (3' GEX);
#           results/03_benchmarking/dataset_summary/gene_lists.rds;
#           results/01_preprocessing/{merscope,xenium}_8um_filtered/*.rds;
#           cfg$visiumhd$data_dir
# Outputs:  results/03_benchmarking/scrna_correlation/avg_expr.rds (FLEX, for fig1.R)
#           results/03_benchmarking/scrna_correlation/correlation_by_reference.rds
#             (data, summary, common_genes, samples; for fig1ext_scrna_correlation.R)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(purrr)
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
  make_option(c("--gene_lists"), type = "character",
              default = "results/03_benchmarking/dataset_summary/gene_lists.rds",
              help    = "Per-platform gene panels from dataset_summary.R [default: %default]"),
  make_option(c("--out_dir"), type = "character",
              default = "results/03_benchmarking/scrna_correlation",
              help    = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
if (!file.exists(opt$gene_lists)) stop("--gene_lists not found: ", opt$gene_lists)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
# Filtered (empty-bin + DBSCAN QC) binning object written by filter_binned.R
filtered_bin_path <- function(output_dir, platform, sample, res) {
  file.path(output_dir, "01_preprocessing",
            paste0(platform, "_", res, "um_filtered"),
            paste0(sample, "_", res, "um_filtered.rds"))
}

# Number-boundary regex for WT IDs so "709" cannot match "1709" or "7090"
wt_id_regex <- function(wt_ids) {
  paste0("(^|[^0-9])(", paste(wt_ids, collapse = "|"), ")([^0-9]|$)")
}

# WT animal ID contained in each sample name (NA if none matches)
extract_animal_id <- function(sample_names, wt_ids) {
  m <- regmatches(sample_names, regexec(wt_id_regex(wt_ids), sample_names))
  vapply(m, function(x) if (length(x) >= 3) x[3] else NA_character_, character(1))
}

# Stop unless the animals contributing to a dataset are exactly wt_ids,
# each represented once (i.e. one pseudobulk sample per biological animal)
check_matched_animals <- function(animal_ids, wt_ids, dataset) {
  animal_ids <- animal_ids[!is.na(animal_ids)]
  if (!setequal(animal_ids, wt_ids) || anyDuplicated(animal_ids) > 0) {
    stop(dataset, ": expected exactly one sample per animal in {",
         paste(wt_ids, collapse = ", "), "}, found {",
         paste(animal_ids, collapse = ", "), "}")
  }
  message(dataset, " animals: ", paste(sort(animal_ids), collapse = ", "))
  invisible(TRUE)
}

# Keep the matched WT samples (e.g. drops wt710, KO, ctrl) and log the rest
select_matched_samples <- function(sample_names, wt_ids, dataset) {
  animal_ids <- extract_animal_id(sample_names, wt_ids)
  keep       <- !is.na(animal_ids)
  check_matched_animals(animal_ids[keep], wt_ids, dataset)
  message(dataset, " excluded (not matched WT): ",
          if (any(!keep)) paste(sample_names[!keep], collapse = ", ") else "none")
  sample_names[keep]
}

# Retrieve the counts layer/slot, supporting both Seurat v4 and v5
get_counts <- function(seu, assay) {
  tryCatch(
    GetAssayData(seu, assay = assay, layer  = "counts"),
    error = function(e) GetAssayData(seu, assay = assay, slot = "counts")
  )
}

# Pseudobulk one ST sample: drop empty bins, then sum across bins
pseudobulk_one_sample <- function(seu, assay) {
  counts    <- get_counts(seu, assay)
  keep_bins <- Matrix::colSums(counts) > 0
  if (any(!keep_bins)) {
    message("    removing ", sum(!keep_bins), " empty bins")
    counts <- counts[, keep_bins, drop = FALSE]
  }
  Matrix::rowSums(counts)
}

# Named list of ST objects -> sparse gene × sample matrix
pseudobulk_samples <- function(seurat_list, assay, target_genes) {
  vec_list <- lapply(seurat_list, pseudobulk_one_sample, assay = assay)

  pb_mat <- Matrix::Matrix(
    0, nrow = length(target_genes), ncol = length(vec_list),
    dimnames = list(target_genes, names(vec_list)),
    sparse = TRUE
  )
  for (s in names(vec_list)) {
    g <- intersect(names(vec_list[[s]]), target_genes)
    pb_mat[g, s] <- vec_list[[s]][g]
  }
  pb_mat
}

detect_sample_col <- function(meta,
                              candidates = c("SampleID", "sample_id", "Sample",
                                             "sample", "donor", "donor_id",
                                             "orig.ident")) {
  for (col in candidates) if (col %in% colnames(meta)) return(col)
  NA_character_
}

# Pseudobulk scRNA-seq into a sparse gene × sample matrix, optionally keeping
# only WT sample IDs
pseudobulk_scRNA <- function(sc_seu, assay, target_genes,
                             wt_ids = NULL, sample_col = NULL,
                             sc_label = "scRNA") {
  counts <- get_counts(sc_seu, assay)
  meta   <- sc_seu@meta.data

  if (is.null(sample_col)) sample_col <- detect_sample_col(meta)

  if (is.na(sample_col)) {
    warning("No sample-ID column found in scRNA metadata. Pooling all cells.")
    sample_vec <- rep("scRNA_all", ncol(counts))
  } else {
    sample_vec <- as.character(meta[[sample_col]])
    message("scRNA sample column: '", sample_col,
            "'  unique values: ", paste(unique(sample_vec), collapse = ", "))
  }

  if (!is.null(wt_ids) && !is.na(sample_col)) {
    keep <- grepl(wt_id_regex(wt_ids), sample_vec)
    if (sum(keep) == 0) {
      stop("No WT cells matched in scRNA column '", sample_col, "'.\n",
           "  WT IDs:  ", paste(wt_ids,            collapse = ", "), "\n",
           "  Present: ", paste(unique(sample_vec), collapse = ", "))
    }
    # NA = undemultiplexed cells
    excluded <- table(sample_vec[!keep], useNA = "ifany")
    message("  ", sc_label, ": excluded ", sum(!keep), " cells (not matched WT): ",
            if (length(excluded) > 0)
              paste0(names(excluded), " (", excluded, ")", collapse = ", ")
            else "none")
    counts     <- counts[, keep, drop = FALSE]
    sample_vec <- sample_vec[keep]
    message("  ", sc_label, ": kept ", sum(keep), " WT cells: ",
            paste(unique(sample_vec), collapse = ", "))
    check_matched_animals(extract_animal_id(unique(sample_vec), wt_ids),
                          wt_ids, sc_label)
  }

  unique_samples <- unique(sample_vec)
  shared         <- intersect(rownames(counts), target_genes)

  pb_mat <- Matrix::Matrix(
    0, nrow = length(target_genes), ncol = length(unique_samples),
    dimnames = list(target_genes, unique_samples),
    sparse = TRUE
  )
  for (s in unique_samples) {
    idx  <- which(sample_vec == s)
    sums <- Matrix::rowSums(counts[shared, idx, drop = FALSE])
    pb_mat[shared, s] <- sums
  }
  pb_mat
}

# log10(CPM+1); empty samples are dropped before scaling. Scaling stays sparse;
# only the final (small) result is densified.
sparse_log10cpm <- function(pb_mat) {
  lib_sizes <- Matrix::colSums(pb_mat)
  keep <- lib_sizes > 0
  if (any(!keep)) {
    message("  dropping empty pseudobulk samples: ",
            paste(colnames(pb_mat)[!keep], collapse = ", "))
    pb_mat    <- pb_mat[, keep, drop = FALSE]
    lib_sizes <- lib_sizes[keep]
  }
  if (ncol(pb_mat) == 0) stop("Pseudobulk matrix has no non-empty samples.")

  cpm_mat <- Matrix::t(Matrix::t(pb_mat) * (1e6 / lib_sizes))

  log10(as.matrix(cpm_mat) + 1)
}

# Per-gene mean over animals of log10(CPM+1) for one platform; library size =
# the pairwise shared genes. Returns list(data, correlation (Pearson r), n_genes).
prepare_correlation_data <- function(sc_seu, sc_assay,
                                     st_list, st_assay,
                                     platform_name,
                                     wt_ids        = NULL,
                                     sc_sample_col = NULL,
                                     sc_label      = "scRNA") {

  # Restrict ST list to the matched WT animals before gene selection and
  # pseudobulk aggregation; stops unless exactly one sample per animal remains
  if (!is.null(wt_ids)) {
    matched <- select_matched_samples(names(st_list), wt_ids, platform_name)
    message(platform_name, " WT samples: ", paste(matched, collapse = ", "))
    st_list <- st_list[matched]
  }

  # Genes shared by the scRNA reference and every matched ST sample
  sc_counts    <- get_counts(sc_seu, sc_assay)
  st_gene_sets <- lapply(st_list, function(seu) rownames(get_counts(seu, st_assay)))
  shared_genes <- Reduce(intersect, c(list(rownames(sc_counts)), st_gene_sets))
  message(platform_name, " vs ", sc_label, ": ", length(shared_genes), " shared genes")

  st_pb <- pseudobulk_samples(st_list, st_assay, shared_genes)
  sc_pb <- pseudobulk_scRNA(sc_seu, sc_assay, shared_genes,
                            wt_ids = wt_ids, sample_col = sc_sample_col,
                            sc_label = sc_label)

  sc_norm <- sparse_log10cpm(sc_pb)
  st_norm <- sparse_log10cpm(st_pb)

  sc_avg <- rowMeans(sc_norm, na.rm = TRUE)
  st_avg <- rowMeans(st_norm, na.rm = TRUE)

  expr_data <- data.frame(
    Gene  = shared_genes,
    scRNA = as.numeric(sc_avg[shared_genes]),
    ST    = as.numeric(st_avg[shared_genes])
  ) |>
    dplyr::filter(is.finite(scRNA), is.finite(ST))

  cor_value <- cor(expr_data$scRNA, expr_data$ST,
                   method = "pearson", use = "complete.obs")

  list(data = expr_data, correlation = cor_value, n_genes = nrow(expr_data))
}

# Run before loading the ST objects so a wrong sample column fails fast
check_reference_samples <- function(sc_seu, sample_col, wt_ids, sc_label) {
  if (!sample_col %in% colnames(sc_seu@meta.data)) {
    stop(sc_label, ": sample column '", sample_col, "' not found in metadata")
  }
  present <- unique(as.character(sc_seu@meta.data[[sample_col]]))
  message(sc_label, " sample IDs present: ",
          paste(sort(present, na.last = TRUE), collapse = ", "))
  found <- extract_animal_id(present, wt_ids)
  check_matched_animals(unique(found[!is.na(found)]), wt_ids, sc_label)
}

# Stop if any common gene is missing, so every intersection panel has the same genes
subset_common_genes <- function(expr_data, common_genes, label) {
  missing <- setdiff(common_genes, expr_data$Gene)
  if (length(missing) > 0) {
    stop(label, ": ", length(missing), " common genes missing from pairwise data: ",
         paste(missing, collapse = ", "))
  }
  dplyr::filter(expr_data, Gene %in% common_genes)
}

wt_ids <- as.character(cfg$scrna$wt_ids)

# ---------------------------------------------------------------------------
# Load single-cell references (10X FLEX first: avg_expr.rds is FLEX only)
# ---------------------------------------------------------------------------
references <- list(
  FLEX = list(label = "10X FLEX",   path  = cfg$scrna$path,
              assay = cfg$scrna$assay,    sample_col = cfg$scrna$sample_col),
  GEX  = list(label = "10X 3' GEX", path  = cfg$scrna$sc_path,
              assay = cfg$scrna$sc_assay, sample_col = cfg$scrna$sc_sample_col)
)

sc_objs <- lapply(references, function(ref) {
  message("Loading ", ref$label, " reference...\n  ", ref$path)
  seu <- readRDS(ref$path)
  check_reference_samples(seu, ref$sample_col, wt_ids, ref$label)
  seu
})

# ---------------------------------------------------------------------------
# Common gene set (VisiumHD ∩ MERSCOPE ∩ Xenium), defined in dataset_summary.R
# ---------------------------------------------------------------------------
gene_lists   <- readRDS(opt$gene_lists)
common_genes <- Reduce(intersect, gene_lists[c("VisiumHD", "MERSCOPE", "Xenium")])
message("Common gene set (dataset_summary gene_lists.rds): ",
        length(common_genes), " genes")

# ---------------------------------------------------------------------------
# Load VisiumHD samples
# ---------------------------------------------------------------------------
message("Loading VisiumHD samples...")

visiumhd_dir     <- cfg$visiumhd$data_dir
visiumhd_samples <- cfg$visiumhd$samples
# Matched WT animals only (drops KO samples)
visiumhd_samples <- visiumhd_samples[
  select_matched_samples(names(visiumhd_samples), wt_ids, "VisiumHD")
]

visiumhd_objs <- lapply(names(visiumhd_samples), function(samp) {
  path <- file.path(visiumhd_dir, visiumhd_samples[[samp]])
  message("  ", samp, ": ", path)
  obj <- readRDS(path)
  DefaultAssay(obj) <- "Spatial.008um"
  obj
})
names(visiumhd_objs) <- names(visiumhd_samples)

# ---------------------------------------------------------------------------
# Load filtered MERSCOPE 8µm binning objects
# ---------------------------------------------------------------------------
message("Loading filtered MERSCOPE 8µm binning objects...")

# Matched WT animals only (drops wt710, ctrl and KO samples)
merscope_samples <- select_matched_samples(
  unlist(cfg$spatial_analysis$merscope_samples), wt_ids, "MERSCOPE"
)
bin_res          <- cfg$bin_resolutions[[1]]   # 8

merscope_objs <- lapply(merscope_samples, function(samp) {
  path <- filtered_bin_path(cfg$output_dir, "merscope", samp, bin_res)
  message("  ", samp, ": ", path)
  readRDS(path)
}) |> setNames(merscope_samples)

# ---------------------------------------------------------------------------
# Load filtered Xenium 8µm binning objects
# ---------------------------------------------------------------------------
message("Loading filtered Xenium 8µm binning objects...")

# Matched WT animals only (drops wt710, ctrl and KO samples)
xenium_samples <- select_matched_samples(
  unlist(cfg$spatial_analysis$xenium_default_samples), wt_ids, "Xenium"
)

xenium_objs <- lapply(xenium_samples, function(samp) {
  path <- filtered_bin_path(cfg$output_dir, "xenium", samp, bin_res)
  message("  ", samp, ": ", path)
  readRDS(path)
}) |> setNames(xenium_samples)

# ---------------------------------------------------------------------------
# Compute per-platform pseudobulk correlations for each reference
# ---------------------------------------------------------------------------
st_platforms <- list(
  VisiumHD = list(objs = visiumhd_objs, assay = "Spatial.008um"),
  MERSCOPE = list(objs = merscope_objs, assay = "Vizgen"),
  Xenium   = list(objs = xenium_objs,   assay = "Xenium")
)

results <- lapply(names(references), function(ref_name) {
  ref <- references[[ref_name]]
  lapply(names(st_platforms), function(plat) {
    message("Computing ", plat, " vs ", ref$label, " correlation...")
    prepare_correlation_data(
      sc_seu        = sc_objs[[ref_name]], sc_assay = ref$assay,
      st_list       = st_platforms[[plat]]$objs,
      st_assay      = st_platforms[[plat]]$assay,
      platform_name = plat,
      wt_ids        = wt_ids, sc_sample_col = ref$sample_col,
      sc_label      = ref$label
    )
  }) |> setNames(names(st_platforms))
}) |> setNames(names(references))

# ---------------------------------------------------------------------------
# Long gene-level table: pairwise and common gene sets
# ---------------------------------------------------------------------------
corr_data <- purrr::imap_dfr(results, function(ref_res, ref_name) {
  purrr::imap_dfr(ref_res, function(res, plat) {
    label <- paste0(plat, " vs ", references[[ref_name]]$label)
    dplyr::bind_rows(
      dplyr::mutate(res$data, gene_set = "pairwise"),
      dplyr::mutate(subset_common_genes(res$data, common_genes, label),
                    gene_set = "common")
    ) |>
      dplyr::mutate(reference = ref_name, platform = plat, .before = 1)
  })
})

corr_summary <- corr_data |>
  dplyr::group_by(reference, platform, gene_set) |>
  dplyr::summarise(
    r       = cor(scRNA, ST, method = "pearson", use = "complete.obs"),
    n_genes = dplyr::n(),
    .groups = "drop"
  )
message("Correlation summary:")
print(as.data.frame(corr_summary))

# The pairwise rows must reproduce the per-platform results exactly
pairwise_check <- corr_summary |>
  dplyr::filter(gene_set == "pairwise") |>
  dplyr::rowwise() |>
  dplyr::mutate(ok = isTRUE(all.equal(r, results[[reference]][[platform]]$correlation)) &&
                     n_genes == results[[reference]][[platform]]$n_genes)
stopifnot(all(pairwise_check$ok))

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------
# FLEX only, structure expected by fig1.R
avg_expr <- results$FLEX

out_file <- file.path(opt$out_dir, "avg_expr.rds")
saveRDS(avg_expr, out_file)
message("Saved: ", out_file)

# Both references (fig1ext_scrna_correlation.R)
corr_by_ref <- list(
  data         = corr_data,
  summary      = corr_summary,
  common_genes = common_genes,
  samples      = list(
    wt_ids   = wt_ids,
    VisiumHD = names(visiumhd_objs),
    MERSCOPE = merscope_samples,
    Xenium   = xenium_samples
  )
)
out_file <- file.path(opt$out_dir, "correlation_by_reference.rds")
saveRDS(corr_by_ref, out_file)
message("Saved: ", out_file)
