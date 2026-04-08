# Helper functions for clustering and cell type deconvolution.

print_error <- function(expr) {
  tryCatch(
    expr,
    error = function(e) {
      message("Error: ", conditionMessage(e))
      message("Traceback:")
      traceback()   # prints last error’s stack
      stop(e)       # rethrow the original error
    }
  )
}

intersect_spe_rows <- function(samples_tb) {
  rows <- Reduce(intersect, lapply(samples_tb$spe, rownames))
  samples_tb$spe <- sapply(samples_tb$spe, function(x) {
    x[rows, ]
  })
  return(samples_tb)
}

combine_spe <- function(samples_tb) {
  mapply(function(spe, sample_id) {
    spe$sample_id <- sample_id
    rowData(spe) <- rowData(spe)[, c("symbol", "EnsembleID")]
    colnames(spe) <- paste0(gsub("-1$", "-", colnames(spe)), sample_id)
    return(spe)
  }, samples_tb$spe, samples_tb$sample) %>%
    do.call(cbind, .)
}

get_cluster_label <- function(
    combined_spe, k = 10, resolution_parameter = 0.06,
    use.dimred = "PCA", n_iterations = 100) {
  combined_spe |>
    buildSNNGraph(k = k, use.dimred = use.dimred) |>
    igraph::cluster_leiden(resolution_parameter = resolution_parameter, n_iterations = n_iterations) |>
    igraph::membership() |>
    factor()
}

# iSC.MEB requires Seurat 4.4.0 / SeuratObject 4.1.4 (not the latest Seurat).
# use_seurat_4_4() temporarily prepends a library path containing those exact
# versions so they take precedence over the default library.
# See https://github.com/XiaoZhangryy/iSC.MEB/issues/2
# requires Seurat 4.4.0 and SeuratObject 4.1.4
run_iSC_MEB <- function(
    seu_list, customGenelist, k, n_pca = 10, maxIter = 20) {
  iSC.MEB::CreateiSCMEBObject(seu_list, customGenelist = customGenelist) |>
    iSC.MEB::CreateNeighbors(platform = "Visium") |>
    iSC.MEB::runPCA(npcs = n_pca, pca.method = "APCA") |>
    iSC.MEB::SetModelParameters(maxIter = maxIter, coreNum = length(seu_list)) |>
    iSC.MEB::iSCMEB(K = k)
}

# install older version of Seurat
# remotes::install_version("SeuratObject",
#   version = "4.1.4",
#   lib = "/stornext/Home/data/allstaff/w/wang.ch/R/x86_64-pc-linux-gnu-library/Seurat_4.4"
# )
# remotes::install_version("Seurat",
#   version = "4.4.0",
#   lib = "/stornext/Home/data/allstaff/w/wang.ch/R/x86_64-pc-linux-gnu-library/Seurat_4.4"
# )
# install latest version of Seurat for default library
# install.packages("Seurat")
# install.packages("SeuratObject")
use_seurat_4_4 <- function(code) {
  withr::with_libpaths(
    new = "/stornext/Home/data/allstaff/w/wang.ch/R/x86_64-pc-linux-gnu-library/Seurat_4.4",
    action = "prefix",
    code = code
  )
}

# check if a vector is consecutive
# e.g. c("A", "A", "B", "B", "B") -> TRUE
#      c("A", "B", "A", "B") -> FALSE
is_consecutive <- function(vec) {
  anyDuplicated(rle(vec)$values) == 0
}

# run RCTD on a single sample
run_RCTD_sample <- function(spe, ref_sce, cell_type_col, rctd_mode = "doublet", max_cores = 1) {
  rctd_data <- spacexr::createRctd(
    spe, ref_sce, cell_type_col = cell_type_col,
    # stop RCTD from removing spots
    pixel_count_min = 0, UMI_min = 0
  )
  spacexr::runRctd(rctd_data, rctd_mode = rctd_mode, max_cores = max_cores)
}

# run_RCTD() splits a multi-sample SPE by sample_id, runs RCTD on each sample
# independently, then reassembles the
# weight matrix and colData in the original bin order.
# Requires spe$sample_id to be consecutive (all bins of the same sample
# grouped together) — validated with is_consecutive() above.
run_RCTD <- function(spe, ref_sce, cell_type_col, rctd_mode = "doublet", max_cores = 1) {

  if (!is_consecutive(spe$sample_id)) {
    stop("spe$sample_id must be consecutive")
  }
  if (!cell_type_col %in% colnames(SummarizedExperiment::colData(ref_sce))) {
    stop(paste0(cell_type_col, " not found in colData of ref_sce"))
  }

  # filter reference sce to only include cell types with > 25 cells
  # as per required by RCTD
  cell_types <- table(ref_sce[[cell_type_col]]) |>
    as.data.frame() |>
    dplyr::filter(Freq > 25) |>
    dplyr::pull(Var1)
  ref_sce <- ref_sce[, ref_sce[[cell_type_col]] %in% cell_types]

  spe_list <- split(seq_len(ncol(spe)), spe$sample_id)
  results_list <- lapply(spe_list, function(idx) {
    message("Running RCTD on sample: ", spe$sample_id[idx[1]])
    results_spe <- run_RCTD_sample(
      spe[, idx], ref_sce,
      cell_type_col, rctd_mode, max_cores
    )
    list(
      mat = SummarizedExperiment::assays(results_spe)$weights,
      cd = as.data.frame(SummarizedExperiment::colData(results_spe)),
      idx = idx
    )
  })

  mat <- do.call(cbind, lapply(results_list, `[[`, "mat"))
  cd <- dplyr::bind_rows(lapply(results_list, `[[`, "cd"))

  original_order <- order(unlist(lapply(results_list, `[[`, "idx")))
  mat <- mat[, original_order]
  cd <- cd[original_order, ]

  list(weights = mat, colData = cd)
}
