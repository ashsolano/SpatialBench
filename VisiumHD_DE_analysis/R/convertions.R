# Conversion functions between Seurat and SpatialExperiment objects.
#
# seu_to_spe_visium_hd()        — primary converter used in the pipeline.
#   Accepts arbitrary assay/image slot names so it works with both 8 µm and
#   16 µm resolutions.  Transfers counts and normcounts if a "data" layer
#   exists in the Seurat assay.
# seu_to_spe_visium_hd_simple() — legacy version, hardcoded to 16 µm assay.

seu_to_spe_visium_hd_simple <- function(seu, sample_id = "sample01") {
  counts <- SeuratObject::GetAssayData(seu, assay = "Spatial.016um", layer = "counts")
  spatial_coords <- SeuratObject::GetTissueCoordinates(seu, image = "slice1.016um")[, c("x", "y")]

  if ("data" %in% SeuratObject::Layers(seu[["Spatial.016um"]])) {
    normcounts <- SeuratObject::GetAssayData(seu, assay = "Spatial.016um", layer = "data")
    return(
      SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts, normcounts = normcounts),
        spatialCoords = as.matrix(spatial_coords),
        sample_id = sample_id
      )
    )
  }

  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    spatialCoords = as.matrix(spatial_coords),
    sample_id = sample_id
  )
}

seu_to_spe_visium_hd <- function(
    seu, assay = "Spatial.016um", image = "slice1.016um",
    sample_id = "sample01") {

  counts <- SeuratObject::GetAssayData(seu, assay = assay, layer = "counts")
  spatial_coords <- SeuratObject::GetTissueCoordinates(seu, image = image)[, c("x", "y")]

  if ("data" %in% SeuratObject::Layers(seu[[assay]])) {
    normcounts <- SeuratObject::GetAssayData(seu, assay = assay, layer = "data")
    return(
      SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts, normcounts = normcounts),
        spatialCoords = as.matrix(spatial_coords),
        sample_id = sample_id
      )
    )
  }

  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    spatialCoords = as.matrix(spatial_coords),
    sample_id = sample_id
  )
}
