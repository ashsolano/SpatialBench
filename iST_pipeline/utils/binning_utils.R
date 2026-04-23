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


ReadVizgen_8um <- function(data.dir, resolution = 8, z = 3L, filter = NA_character_) {

  mx <- data.table::fread(file.path(data.dir, "detected_transcripts.csv"), sep = ",", verbose = FALSE)
  mx <- mx[mx$global_z == z, , drop = FALSE]

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


LoadVizgen_8um <- function(data.dir, resolution = 8, fov = "fov", assay = "Vizgen") {

  data <- ReadVizgen_8um(data.dir, resolution = resolution, z = 3L)

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
