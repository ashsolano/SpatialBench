# Purpose:  Compute pseudobulk correlations between 10X FLEX scRNA-seq and
#           VisiumHD, MERSCOPE, and Xenium binning objects at 8µm resolution.
#           Produces per-platform averaged log10(CPM+1) expression data frames
#           and Pearson correlations suitable for fig1c density scatter plots.
#           Only matched WT animals (cfg$scrna$wt_ids, i.e. 709 and 713) are used
#           for every dataset; this is enforced before pseudobulk aggregation.
# Inputs:   config/config.yaml  (scrna, visiumhd, spatial_analysis, bin_resolutions)
#           cfg$scrna$path  (scFlex_seu.rds)
#           results/01_preprocessing/merscope_8um_filtered/{sample}_8um_filtered.rds
#           results/01_preprocessing/xenium_8um_filtered/{sample}_8um_filtered.rds
#           cfg$visiumhd$data_dir  (VisiumHD Seurat objects with Spatial.008um assay)
# Outputs:  results/03_benchmarking/scrna_correlation/avg_expr.rds


suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
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
              default = "results/03_benchmarking/scrna_correlation",
              help    = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Path to a filtered (post empty-bin + DBSCAN QC) binning object written by
# filter_binned.R, e.g. results/01_preprocessing/merscope_8um_filtered/
# wt709_batch13_8um_filtered.rds
filtered_bin_path <- function(output_dir, platform, sample, res) {
  file.path(output_dir, "01_preprocessing",
            paste0(platform, "_", res, "um_filtered"),
            paste0(sample, "_", res, "um_filtered.rds"))
}

# Number-boundary regex for WT IDs so "709" cannot match "1709" or "7090"
wt_id_regex <- function(wt_ids) {
  paste0("(^|[^0-9])(", paste(wt_ids, collapse = "|"), ")([^0-9]|$)")
}

# Extract the WT animal ID contained in each sample name (NA if none matches)
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

# Subset a vector of sample names to the matched WT animals and verify that
# exactly the expected animals are present (e.g. drops wt710 and KO/ctrl)
select_matched_samples <- function(sample_names, wt_ids, dataset) {
  animal_ids <- extract_animal_id(sample_names, wt_ids)
  keep       <- !is.na(animal_ids)
  check_matched_animals(animal_ids[keep], wt_ids, dataset)
  sample_names[keep]
}

# Retrieve the counts layer/slot, supporting both Seurat v4 and v5
get_counts <- function(seu, assay) {
  tryCatch(
    GetAssayData(seu, assay = assay, layer  = "counts"),
    error = function(e) GetAssayData(seu, assay = assay, slot = "counts")
  )
}

# Pseudobulk a single ST sample: drop empty bins first, then sum across bins.
# Returns a named numeric vector (genes).
pseudobulk_one_sample <- function(seu, assay) {
  counts    <- get_counts(seu, assay)
  keep_bins <- Matrix::colSums(counts) > 0
  if (any(!keep_bins)) {
    message("    removing ", sum(!keep_bins), " empty bins")
    counts <- counts[, keep_bins, drop = FALSE]
  }
  Matrix::rowSums(counts)
}

# Aggregate a named list of ST Seurat objects into a sparse gene × sample matrix.
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

# Detect the per-cell sample ID column in scRNA metadata
detect_sample_col <- function(meta,
                              candidates = c("SampleID", "sample_id", "Sample",
                                             "sample", "donor", "donor_id",
                                             "orig.ident")) {
  for (col in candidates) if (col %in% colnames(meta)) return(col)
  NA_character_
}

# Pseudobulk the scRNA-seq object into a sparse gene × sample matrix,
# optionally subsetting to WT sample IDs (matched by grepl against column values).
pseudobulk_scRNA <- function(sc_seu, assay, target_genes,
                             wt_ids = NULL, sample_col = NULL) {
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
    counts     <- counts[, keep, drop = FALSE]
    sample_vec <- sample_vec[keep]
    message("  scRNA: kept ", sum(keep), " WT cells: ",
            paste(unique(sample_vec), collapse = ", "))
    # Verify exactly the matched animals contribute to the scRNA pseudobulk
    check_matched_animals(extract_animal_id(unique(sample_vec), wt_ids),
                          wt_ids, "10x FLEX")
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

# Sparse-compatible log10(CPM+1) normalisation.
# Empty samples are dropped before scaling. Column scaling is done via sparse
# matrix arithmetic; the final log transform materialises the (small) result.
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

  # Scale each column to CPM without materialising a dense intermediate
  cpm_mat <- Matrix::t(Matrix::t(pb_mat) * (1e6 / lib_sizes))

  # log10(CPM+1) densifies here; pb_mat is genes × samples (small)
  log10(as.matrix(cpm_mat) + 1)
}

# Prepare per-gene average log10(CPM+1) for one platform.
# Returns a list: data (gene-level data frame), correlation (Pearson r), n_genes.
prepare_correlation_data <- function(sc_seu, sc_assay,
                                     st_list, st_assay,
                                     platform_name,
                                     wt_ids        = NULL,
                                     sc_sample_col = NULL) {

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
  message(platform_name, ": ", length(shared_genes), " shared genes")

  st_pb <- pseudobulk_samples(st_list, st_assay, shared_genes)
  sc_pb <- pseudobulk_scRNA(sc_seu, sc_assay, shared_genes,
                            wt_ids = wt_ids, sample_col = sc_sample_col)

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

# ---------------------------------------------------------------------------
# Load scRNA-seq FLEX reference
# ---------------------------------------------------------------------------
message("Loading scRNA-seq FLEX reference...")
sc_path <- cfg$scrna$path
message("  ", sc_path)
sc_rnaseq <- readRDS(sc_path)

sc_assay      <- cfg$scrna$assay
sc_sample_col <- cfg$scrna$sample_col
wt_ids        <- as.character(cfg$scrna$wt_ids)

# ---------------------------------------------------------------------------
# Load VisiumHD samples
# ---------------------------------------------------------------------------
message("Loading VisiumHD samples...")

visiumhd_dir     <- cfg$visiumhd$data_dir
visiumhd_samples <- cfg$visiumhd$samples
# Load only the matched WT animals (drops KO samples)
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

# Load only the matched WT animals (drops wt710, ctrl and KO samples)
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

# Load only the matched WT animals (drops wt710, ctrl and KO samples)
xenium_samples <- select_matched_samples(
  unlist(cfg$spatial_analysis$xenium_default_samples), wt_ids, "Xenium"
)

xenium_objs <- lapply(xenium_samples, function(samp) {
  path <- filtered_bin_path(cfg$output_dir, "xenium", samp, bin_res)
  message("  ", samp, ": ", path)
  readRDS(path)
}) |> setNames(xenium_samples)

# ---------------------------------------------------------------------------
# Compute per-platform pseudobulk correlations
# ---------------------------------------------------------------------------
message("Computing VisiumHD correlation...")
visiumhd_results <- prepare_correlation_data(
  sc_seu        = sc_rnaseq,     sc_assay = sc_assay,
  st_list       = visiumhd_objs, st_assay = "Spatial.008um",
  platform_name = "VisiumHD",
  wt_ids        = wt_ids,        sc_sample_col = sc_sample_col
)

message("Computing MERSCOPE correlation...")
merscope_results <- prepare_correlation_data(
  sc_seu        = sc_rnaseq,     sc_assay = sc_assay,
  st_list       = merscope_objs, st_assay = "Vizgen",
  platform_name = "MERSCOPE",
  wt_ids        = wt_ids,        sc_sample_col = sc_sample_col
)

message("Computing Xenium correlation...")
xenium_results <- prepare_correlation_data(
  sc_seu        = sc_rnaseq,    sc_assay = sc_assay,
  st_list       = xenium_objs,  st_assay = "Xenium",
  platform_name = "Xenium",
  wt_ids        = wt_ids,       sc_sample_col = sc_sample_col
)

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------
avg_expr <- list(
  VisiumHD = visiumhd_results,
  MERSCOPE = merscope_results,
  Xenium   = xenium_results
)

out_file <- file.path(opt$out_dir, "avg_expr.rds")
saveRDS(avg_expr, out_file)
message("Saved: ", out_file)
