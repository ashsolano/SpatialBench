# Purpose:  Compute per-feature Moran's I (spatial autocorrelation) for
#           background vs target features on the filtered 8µm binned objects
#           of one platform (MERSCOPE or Xenium). Bin counts are rasterised
#           onto a hexagonal grid (SEraster, 100µm) and Moran's I is computed
#           with MERINGUE on the rasterised pixel means.
# Inputs:   config/config.yaml  (spatial_analysis, qc_backgrounds, output_dir)
#           results/01_preprocessing/{platform}_8um_filtered/{sample}_8um_filtered.rds
# Outputs:  --out_rds, e.g. results/03_benchmarking/moransi/{platform}_moransi.rds
#               (feature, observed, expected, sd, p.adj, type, sample)
#
# NOTE: requires R >= 4.5 (SEraster 0.99.5); run with the R/4.5.1 module plus
#       geos/proj/gdal (sf) and ImageMagick (magick, via SpatialExperiment).

suppressPackageStartupMessages({
  library(Seurat)
  library(SpatialExperiment)
  library(SEraster)
  library(MERINGUE)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(yaml)
  library(optparse)
})

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--config"), type = "character",
              default = "config/config.yaml",
              help    = "Path to config.yaml [default: %default]"),
  make_option(c("--platform"), type = "character",
              help    = "Platform: merscope or xenium"),
  make_option(c("--resolution"), type = "double", default = 100,
              help    = "Rasterisation resolution in µm [default: %default]"),
  make_option(c("--out_rds"), type = "character",
              help    = "Output RDS path")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
dir.create(dirname(opt$out_rds), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Platform settings
# ---------------------------------------------------------------------------
# MERSCOPE blank probes are in the "Blanks" assay; Xenium "Unassigned Codeword"
# features are loaded by Seurat into the "BlankCodeword" assay. These match
# the background definitions used for FDR in qc_backgrounds.R.
platform_settings <- list(
  merscope = list(
    samples   = cfg$spatial_analysis$merscope_samples,
    bg_assay  = "Blanks",
    tg_assay  = "Vizgen"
  ),
  xenium = list(
    samples   = cfg$spatial_analysis$xenium_default_samples,
    bg_assay  = "BlankCodeword",
    tg_assay  = "Xenium"
  )
)
if (!opt$platform %in% names(platform_settings)) stop("Unknown --platform: ", opt$platform)
settings <- platform_settings[[opt$platform]]

# Drop samples without background probe signal (config["qc_backgrounds"])
excluded <- unlist(cfg$qc_backgrounds$exclude_samples[[opt$platform]])
if (length(excluded) > 0) message("Excluding samples (no background signal): ", paste(excluded, collapse = ", "))
settings$samples <- setdiff(settings$samples, excluded)

# ---------------------------------------------------------------------------
# Moran's I for one assay
# ---------------------------------------------------------------------------
run_moran <- function(seurat_obj, assay_name, type_label, res = 100) {
  # Name coordinates by the centroids' own cell IDs so counts can be matched to them
  fov       <- names(seurat_obj@images)[1]
  centroids <- seurat_obj@images[[fov]]@boundaries[["centroids"]]
  pos_mat   <- matrix(
    centroids@coords[, 1:2], ncol = 2,
    dimnames = list(centroids@cells, c("x", "y"))
  )

  counts <- GetAssayData(seurat_obj, assay = assay_name, layer = "counts")
  counts <- counts[, rownames(pos_mat)]

  # Hexagonal grid, mean per pixel
  se <- SpatialExperiment(
    assays        = list(counts = counts),
    spatialCoords = pos_mat
  )
  rast <- rasterizeGeneExpression(
    se,
    assay_name = "counts",
    resolution = res,
    square     = FALSE,
    n_threads  = 1
  )

  # Neighbours = pixels within one pixel width
  pix  <- assay(rast, "pixelval")
  w    <- getSpatialNeighbors(spatialCoords(rast), filterDist = res)
  Ires <- getSpatialPatterns(pix, w, verbose = FALSE)

  tibble(
    feature  = rownames(Ires),
    observed = Ires[, "observed"],
    expected = Ires[, "expected"],
    sd       = Ires[, "sd"],
    p.adj    = Ires[, "p.adj"],
    type     = type_label
  )
}

# ---------------------------------------------------------------------------
# Run across samples (loaded one at a time to limit memory)
# ---------------------------------------------------------------------------
in_dir <- file.path(cfg$output_dir, "01_preprocessing", paste0(opt$platform, "_8um_filtered"))
paths  <- file.path(in_dir, paste0(settings$samples, "_8um_filtered.rds"))
missing <- paths[!file.exists(paths)]
if (length(missing) > 0) stop("Missing filtered 8um objects:\n  ", paste(missing, collapse = "\n  "))

moransi <- map2_dfr(settings$samples, paths, function(samp, path) {
  message("Processing ", samp, ": ", path)
  obj <- readRDS(path)
  if (!settings$bg_assay %in% SeuratObject::Assays(obj)) {
    stop(samp, " has no ", settings$bg_assay, " assay. ",
         "Add it to qc_backgrounds$exclude_samples in config.yaml if intended.")
  }
  bind_rows(
    run_moran(obj, settings$bg_assay, "Background", res = opt$resolution),
    run_moran(obj, settings$tg_assay, "Target",     res = opt$resolution)
  ) %>%
    mutate(sample = samp)
})

# Check that all samples were processed
print(count(moransi, sample, type), n = Inf)
message("Samples processed: ", n_distinct(moransi$sample))

saveRDS(moransi, opt$out_rds)
message("Saved: ", opt$out_rds)
