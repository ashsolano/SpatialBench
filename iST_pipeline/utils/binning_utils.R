# Purpose:  Utility functions for binning spatial transcriptomics data into
#           fixed-size spatial bins and loading as Seurat objects.
#           Supports Xenium (10x Genomics) and MERSCOPE (Vizgen) platforms.
# Inputs:   Raw data directories containing transcripts.parquet (Xenium) or
#           detected_transcripts.csv (Vizgen)
# Outputs:  Named list of sparse count matrices (Read*) or a Seurat object (Load*)

library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(arrow)
library(purrr)
library(Seurat)
library(dbscan)


ReadXenium_binned <- function(data.dir, resolution = 8, pixel_size = 1, mols.qv.threshold = 20) {

  transcripts <- arrow::read_parquet(file.path(data.dir, "transcripts.parquet"))
  transcripts <- transcripts %>%
    filter(qv >= mols.qv.threshold) %>%
    mutate(
      x_location_um = x_location * pixel_size,
      y_location_um = y_location * pixel_size,
      binx          = floor(x_location_um / resolution),
      biny          = floor(y_location_um / resolution),
      spatial_bin   = paste(binx, biny, sep = "_")
    )

  # Older Xenium data lacks codeword_category; infer it from feature name patterns
  if (!"codeword_category" %in% colnames(transcripts)) {
    message("'codeword_category' missing — inferring from feature_name patterns.")
    transcripts <- transcripts %>%
      mutate(
        codeword_category = case_when(
          grepl("^NegControlCodeword_", feature_name) ~ "negative_control_codeword",
          grepl("^NegControlProbe_",    feature_name) ~ "negative_control_probe",
          grepl("^UnassignedCodeword_", feature_name) ~ "unassigned_codeword",
          TRUE                                         ~ "custom_gene"
        )
      )
  }

  # Compute all_bins up front so every category matrix has the same column set
  all_bins <- sort(unique(transcripts$spatial_bin))

  binned_data <- transcripts %>%
    group_by(spatial_bin, feature_name, codeword_category) %>%
    summarise(count = n(), .groups = "drop") %>%
    split(.$codeword_category) %>%
    purrr::map(~ {
      feature_data <- .x %>%
        dplyr::select(spatial_bin, feature_name, count) %>%
        complete(spatial_bin = all_bins, feature_name, fill = list(count = 0))

      count_matrix <- feature_data %>%
        pivot_wider(names_from = spatial_bin, values_from = count, values_fill = 0) %>%
        column_to_rownames("feature_name") %>%
        as.matrix()

      as(count_matrix, "dgCMatrix")
    })

  # Rename to Seurat-compatible names; reorder list and names together to keep them in sync
  name_mapping <- c(
    "custom_gene"               = "Gene Expression",
    "negative_control_codeword" = "Negative Control Codeword",
    "negative_control_probe"    = "Negative Control Probe",
    "unassigned_codeword"       = "Unassigned Codeword"
  )

  new_names <- name_mapping[names(binned_data)]

  if (any(is.na(new_names))) {
    stop(
      "Unrecognised codeword categories: ",
      paste(names(binned_data)[is.na(new_names)], collapse = ", ")
    )
  }

  sorted_idx         <- order(new_names)
  binned_data        <- binned_data[sorted_idx]
  names(binned_data) <- new_names[sorted_idx]

  return(binned_data)
}


LoadXenium_binned <- function(data.dir, resolution = 8, fov = "fov", assay = "Xenium") {

  data <- ReadXenium_binned(data.dir = data.dir, resolution = resolution)

  xenium.obj <- CreateSeuratObject(counts = data[["Gene Expression"]], assay = assay)

  # Older Xenium output uses "Blank Codeword"; newer uses "Unassigned Codeword"
  if ("Blank Codeword" %in% names(data)) {
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data[["Blank Codeword"]])
  } else {
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data[["Unassigned Codeword"]])
  }
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]]    <- CreateAssayObject(counts = data[["Negative Control Probe"]])

  # Bin IDs are "binx_biny"; multiply indices back by resolution to get micron coordinates
  bins <- colnames(xenium.obj)
  bin_centroid_df <- data.frame(
    x = as.numeric(sub("_.*", "", bins)) * resolution,
    y = as.numeric(sub(".*_", "", bins)) * resolution,
    row.names = bins
  )

  segmentations.data <- list(
    "centroids" = CreateCentroids(coords = bin_centroid_df, nsides = 4)
  )

  coords <- CreateFOV(
    coords    = segmentations.data,
    type      = "centroids",
    molecules = NULL,
    assay     = assay
  )

  xenium.obj[[fov]] <- coords

  return(xenium.obj)
}


ReadVizgen_binned <- function(data.dir, resolution = 8, z = "all", filter = NA_character_) {

  mx <- data.table::fread(file.path(data.dir, "detected_transcripts.csv"), sep = ",", verbose = FALSE)

  # z = "all" keeps every z-plane; otherwise filter down to the requested plane
  if (!identical(z, "all")) {
    mx <- mx[mx$global_z == z, , drop = FALSE]
  }

  if (!is.na(filter)) {
    mx <- mx[!grepl(pattern = filter, x = mx$gene), , drop = FALSE]
  }

  # global_x / global_y are already in micron space for Vizgen data
  microns <- data.frame(
    x             = mx$global_x,
    y             = mx$global_y,
    gene          = mx$gene,
    transcript_id = mx$transcript_id,
    stringsAsFactors = FALSE
  )

  # Genes prefixed with "Blank-" are negative controls; separate them from real genes
  transcripts <- microns %>%
    mutate(
      binx          = floor(x / resolution),
      biny          = floor(y / resolution),
      spatial_bin   = paste(binx, biny, sep = "_"),
      gene_category = ifelse(grepl("^Blank-", gene), "Blanks", "Gene Expression")
    )

  all_bins <- sort(unique(transcripts$spatial_bin))

  binned_matrices <- transcripts %>%
    group_by(spatial_bin, gene, gene_category) %>%
    summarise(count = n(), .groups = "drop") %>%
    split(.$gene_category) %>%
    purrr::map(~ {
      feature_data <- .x %>%
        dplyr::select(spatial_bin, gene, count) %>%
        complete(spatial_bin = all_bins, gene, fill = list(count = 0))

      count_matrix <- feature_data %>%
        pivot_wider(names_from = spatial_bin, values_from = count, values_fill = 0) %>%
        column_to_rownames("gene") %>%
        as.matrix()

      as(count_matrix, "dgCMatrix")
    })

  return(list(binned_matrices = binned_matrices, microns = microns))
}


LoadVizgen_binned <- function(data.dir, resolution = 8, fov = "fov", assay = "Vizgen", z = "all") {

  # All z-planes are used by default to match the Xenium binning workflow,
  # which does not filter by z-plane; addresses reviewer concerns about
  # unmatched z-plane handling between platforms.
  data <- ReadVizgen_binned(data.dir, resolution = resolution, z = z)

  vizgen.obj <- CreateSeuratObject(counts = data$binned_matrices[["Gene Expression"]], assay = assay)

  # Add Blanks assay only if Blank- genes were present in the data
  if ("Blanks" %in% names(data$binned_matrices)) {
    vizgen.obj[["Blanks"]] <- CreateAssayObject(counts = data$binned_matrices[["Blanks"]])
  }

  # Bin IDs are "binx_biny"; multiply indices back by resolution to get micron coordinates
  bins <- colnames(vizgen.obj)
  bin_centroid_df <- data.frame(
    x = as.numeric(sub("_.*", "", bins)) * resolution,
    y = as.numeric(sub(".*_", "", bins)) * resolution,
    row.names = bins
  )

  segmentations.data <- list(
    "centroids" = CreateCentroids(coords = bin_centroid_df, nsides = 4)
  )

  coords <- CreateFOV(
    coords    = segmentations.data,
    type      = "centroids",
    molecules = data$microns,
    assay     = assay
  )

  vizgen.obj[[fov]] <- coords

  return(vizgen.obj)
}


# Detect the background/tissue nCount boundary from the local minimum between
# the background and signal peaks of the count distribution.
find_background_threshold <- function(counts, max_count = 500, min_thresh = 5) {
  counts <- counts[counts > 0 & counts <= max_count]
  if (length(counts) == 0) return(min_thresh)
  d <- density(counts, bw = 2)
  dy <- diff(d$y)
  minima_idx <- which(dy[-1] > 0 & dy[-length(dy)] < 0) + 1
  minima_x   <- d$x[minima_idx]
  thresh <- minima_x[minima_x > min_thresh][1]
  if (is.na(thresh)) thresh <- min_thresh
  thresh
}


# Post-processing QC: drop background/empty bins and remove spatially isolated
# bins (e.g. debris, edge artefacts) via DBSCAN. min_cluster_size lets more
# than one tissue cluster survive, for samples with multiple tissue pieces.
filter_binned_seurat <- function(obj, assay = "Vizgen",
                                  min_count = 1,
                                  use_adaptive_threshold = TRUE,
                                  eps = 15, minPts = 5,
                                  min_cluster_size = 1000) {

  fov_name <- names(obj@images)[1]
  fov      <- obj@images[[fov_name]]
  cents    <- fov@boundaries$centroids@coords
  rownames(cents) <- colnames(obj)

  df <- data.frame(
    x      = cents[, 1],
    y      = cents[, 2],
    nCount = obj@meta.data[[paste0("nCount_", assay)]],
    row.names = colnames(obj)
  )

  # Adaptive threshold finds the valley between the background and tissue
  # peaks of the count distribution; otherwise fall back to a fixed min_count
  if (use_adaptive_threshold) {
    min_count <- find_background_threshold(df$nCount)
    message("Adaptive background threshold: nCount >= ", round(min_count, 2))
  }

  df_nonempty  <- df[df$nCount >= min_count, , drop = FALSE]

  cl            <- dbscan::dbscan(as.matrix(df_nonempty[, c("x", "y")]), eps = eps, minPts = minPts)$cluster
  cluster_sizes <- table(cl[cl > 0])
  keep_clusters <- as.integer(names(cluster_sizes)[cluster_sizes >= min_cluster_size])

  # Fall back to the single largest cluster if none meet min_cluster_size
  if (length(keep_clusters) == 0) {
    keep_clusters <- as.integer(names(which.max(cluster_sizes)))
  }

  keep_cells <- rownames(df_nonempty)[cl %in% keep_clusters]

  centroids_keep <- CreateCentroids(data.frame(
    x    = df_nonempty[keep_cells, "x"],
    y    = df_nonempty[keep_cells, "y"],
    cell = keep_cells
  ))

  fov_keep <- CreateFOV(
    coords = centroids_keep,
    assay  = DefaultAssay(fov),
    key    = Key(fov),
    name   = names(fov@boundaries)[1]
  )

  obj_no_fov        <- obj
  obj_no_fov@images <- list()
  obj_filt          <- subset(obj_no_fov, cells = keep_cells)
  obj_filt[[fov_name]] <- fov_keep
  DefaultFOV(obj_filt) <- fov_name

  message("Filtered: ", ncol(obj), " -> ", ncol(obj_filt), " bins")
  obj_filt
}
