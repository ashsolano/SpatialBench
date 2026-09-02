# Purpose: Utility functions for converting aligned spatial bin centroids
# from pixel to physical micron coordinates, and extracting a fixed-size
# physical region of interest (ROI) from a Seurat object.

library(Seurat)

# Convert a Seurat object's bin centroid coordinates from pixels to microns.
# Works on VisiumHD reference objects or aligned MERSCOPE/Xenium objects —
# both share the same pixel coordinate space post-alignment (see
# align_binned.R), so the same microns_per_pixel value applies to either.
pixels_to_microns <- function(obj, mpp, image_name = NULL) {
  img_name <- if (is.null(image_name)) names(obj@images)[1] else image_name
  boundary <- obj@images[[img_name]]@boundaries$centroids
  cents <- boundary@coords
  cell_ids <- boundary@cells
  data.frame(
    cell_id = cell_ids,
    x_um = cents[,1] * mpp,
    y_um = cents[,2] * mpp
  )
}

# Compute an ROI center (microns) from a reference object's tissue extent —
# typically VisiumHD, since it's the un-warped target space every platform
# was registered into.
get_roi_center <- function(vis_obj, mpp, image_name = NULL) {
  coords_um <- pixels_to_microns(vis_obj, mpp, image_name)
  list(cx = mean(range(coords_um$x_um)), cy = mean(range(coords_um$y_um)))
}

# Subset a Seurat object to bins within a fixed physical ROI, returning a
# fully intact Seurat object (all metadata/counts preserved).

extract_roi <- function(obj, mpp, roi_center_um, roi_size_um, image_name = NULL, assay_name = NULL) {
  fov_name <- if (is.null(image_name)) names(obj@images)[1] else image_name
  fov <- obj@images[[fov_name]]

  boundary <- fov@boundaries$centroids
  cents <- boundary@coords
  cell_ids <- boundary@cells

  x_um <- cents[,1] * mpp
  y_um <- cents[,2] * mpp

  half <- roi_size_um / 2
  keep <- x_um >= roi_center_um$cx - half & x_um <= roi_center_um$cx + half &
          y_um >= roi_center_um$cy - half & y_um <= roi_center_um$cy + half
  keep_cells <- cell_ids[keep]

  if (!is.null(assay_name)) {
    DefaultAssay(obj) <- assay_name
  }

  obj_no_fov <- obj
  obj_no_fov@images <- list()
  obj_filt <- subset(obj_no_fov, cells = keep_cells)

  centroids_keep <- CreateCentroids(data.frame(
    x    = cents[keep, 1],
    y    = cents[keep, 2],
    cell = keep_cells
  ))

  fov_keep <- CreateFOV(
    coords = centroids_keep,
    assay  = DefaultAssay(obj_filt),
    key    = Key(fov),
    name   = names(fov@boundaries)[1]
  )

  obj_filt[[fov_name]] <- fov_keep
  DefaultFOV(obj_filt) <- fov_name

  obj_filt
}
