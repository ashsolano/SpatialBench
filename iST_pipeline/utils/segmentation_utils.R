# Purpose:  Utility functions for loading cell-segmented spatial transcriptomics data
#           into Seurat objects. Covers Xenium (myReadXenium, myLoadXenium),
#           MERSCOPE (myReadVizgen, myLoadVizgen), and Proseg (myLoadProseg).
#           Used for default vendor, Cellpose, and Proseg segmentation methods.
# Inputs:   Xenium: output directory containing cell_feature_matrix/, cells.csv.gz,
#                   cell_boundaries.csv.gz, and transcripts.parquet
#           MERSCOPE: Vizgen output directory (cell_by_gene.csv, cell_metadata.csv,
#                     detected_transcripts.csv, and either cell_boundaries/ HDF5 files
#                     or cell_boundaries.parquet with WKB-encoded polygon geometry)
#           Proseg: output directory (expected-counts.csv.gz, cell-metadata.csv.gz,
#                   transcript-metadata.csv.gz, cell-polygons.geojson.gz)
# Outputs:  Named list of data frames / matrices (myReadXenium, myReadVizgen) or
#           Seurat objects with spatial coordinates, segmentation boundaries, and
#           control assays (myLoadXenium, myLoadVizgen, myLoadProseg)

library(Seurat)
library(sf)       # WKB polygon decoding for cell_boundaries.parquet


myReadXenium <- function(data.dir, outs = c("matrix", "microns"),
                         type = "centroids", mols.qv.threshold = 20) {

  type <- match.arg(arg = type, choices = c("centroids", "segmentations"), several.ok = TRUE)
  outs <- match.arg(arg = outs, choices = c("matrix", "microns"), several.ok = TRUE)
  outs <- c(outs, type)

  has_dt <- requireNamespace("data.table", quietly = TRUE) &&
    requireNamespace("R.utils", quietly = TRUE)

  data <- sapply(outs, function(otype) {
    switch(EXPR = otype,
      matrix = {
        suppressWarnings(Read10X(data.dir = file.path(data.dir, "cell_feature_matrix/")))
      },
      centroids = {
        if (has_dt) {
          cell_info <- as.data.frame(data.table::fread(file.path(data.dir, "cells.csv.gz")))
        } else {
          cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
        }
        data.frame(
          x    = cell_info$x_centroid,
          y    = cell_info$y_centroid,
          cell = cell_info$cell_id,
          stringsAsFactors = FALSE
        )
      },
      segmentations = {
        if (has_dt) {
          cell_boundaries_df <- as.data.frame(
            data.table::fread(file.path(data.dir, "cell_boundaries.csv.gz"))
          )
        } else {
          cell_boundaries_df <- read.csv(
            file.path(data.dir, "cell_boundaries.csv.gz"),
            stringsAsFactors = FALSE
          )
        }
        names(cell_boundaries_df) <- c("cell", "x", "y")
        cell_boundaries_df
      },
      microns = {
        transcripts <- arrow::read_parquet(file.path(data.dir, "transcripts.parquet"))
        transcripts <- subset(transcripts, qv >= mols.qv.threshold)
        data.frame(
          x    = transcripts$x_location,
          y    = transcripts$y_location,
          gene = transcripts$feature_name,
          stringsAsFactors = FALSE
        )
      },
      stop("Unknown Xenium input type: ", otype)
    )
  }, USE.NAMES = TRUE)

  return(data)
}


myLoadXenium <- function(data.dir, fov = "fov", assay = "Xenium") {

  data <- myReadXenium(
    data.dir = data.dir,
    type     = c("centroids", "segmentations")
  )

  segmentations.data <- list(
    "centroids"    = CreateCentroids(data$centroids),
    "segmentation" = CreateSegmentation(data$segmentations)
  )

  coords <- CreateFOV(
    coords    = segmentations.data,
    type      = c("segmentation", "centroids"),
    molecules = data$microns,
    assay     = assay
  )

  xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)

  # Older Xenium output uses "Blank Codeword"; newer uses "Unassigned Codeword"
  if ("Blank Codeword" %in% names(data$matrix)) {
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Blank Codeword"]])
  } else {
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  }
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]]    <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])

  xenium.obj[[fov]] <- coords

  return(xenium.obj)
}


myReadVizgen <- function(data.dir, z = 3L,
                         type     = c("segmentations", "centroids"),
                         mol.type = "microns") {
  # Adapted from Seurat::ReadVizgen. Reads MERSCOPE/Vizgen output files into a named
  # list ready for CreateFOV(). The segmentations block is extended to fall back to
  # cell_boundaries.parquet (WKB format) when no HDF5 cell_boundaries/ directory exists.

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install 'data.table' for this function")
  }
  hdf5     <- requireNamespace("hdf5r", quietly = TRUE)
  type     <- match.arg(type,     c("segmentations", "centroids", "boxes"), several.ok = TRUE)
  mol.type <- match.arg(mol.type, c("microns", "pixels"),                   several.ok = TRUE)

  if (!z %in% seq.int(0L, 6L)) stop("The z-index must be in the range [0, 6]")
  if (!dir.exists(data.dir))   stop("Cannot find Vizgen directory ", data.dir)

  # Locate input files using the same filename patterns as Seurat::ReadVizgen
  find_file <- function(pattern) {
    hits <- list.files(data.dir, pattern = pattern, full.names = TRUE, recursive = FALSE)
    if (length(hits) == 0L) return(NA_character_)
    sort(hits, decreasing = TRUE)[1L]
  }
  f_tx   <- find_file("cell_by_gene[_a-zA-Z0-9]*.csv")
  f_sp   <- find_file("cell_metadata[_a-zA-Z0-9]*.csv")
  f_mols <- find_file("detected_transcripts[_a-zA-Z0-9]*.csv")
  h5dir  <- file.path(data.dir, "cell_boundaries")
  zidx   <- paste0("zIndex_", z)

  # Preload spatial metadata (shared by centroids, segmentations, and boxes)
  if (is.na(f_sp)) stop("Cannot find cell_metadata CSV in ", data.dir)
  message("Preloading cell spatial coordinates")
  sp <- data.table::fread(f_sp, sep = ",", data.table = FALSE, verbose = FALSE)
  rownames(sp) <- as.character(sp[[1]])
  sp <- sp[, -1, drop = FALSE]

  # Check which segmentation source is available: HDF5 directory or parquet fallback
  parquet_file <- file.path(data.dir, "cell_boundaries.parquet")
  use_hdf5     <- hdf5 && dir.exists(h5dir)
  use_parquet  <- !use_hdf5 && file.exists(parquet_file)
  if ("segmentations" %in% type && !use_hdf5 && !use_parquet) {
    warning(
      "No segmentation source found (no HDF5 directory and no cell_boundaries.parquet); ",
      "dropping segmentations from output",
      immediate. = TRUE
    )
    type <- setdiff(type, "segmentations")
  }

  # Preload molecule coordinates if needed
  if (length(mol.type) > 0L && !is.na(f_mols)) {
    message("Preloading molecule coordinates")
    mx <- data.table::fread(f_mols, sep = ",", data.table = FALSE, verbose = FALSE)
    mx <- mx[mx$global_z == z, , drop = FALSE]
  }

  outs <- list()

  # --- Counts matrix (genes x cells) ---
  if (!is.na(f_tx)) {
    message("Reading counts matrix")
    # Force the cell ID column to character at read time. 18-digit Vizgen cell IDs
    # exceed double precision (2^53); if fread reads them as double (no bit64) or
    # CreateSeuratObject later coerces integer64, IDs silently lose precision and
    # stop matching the IDs in cell_boundaries.parquet.
    tx <- data.table::fread(f_tx, sep = ",", data.table = FALSE, verbose = FALSE,
                            colClasses = list(character = 1))
    rownames(tx) <- tx[[1]]
    tx <- t(as.matrix(tx[, -1, drop = FALSE]))
    ratio <- getOption("Seurat.input.sparse_ratio", default = 0.4)
    if ((sum(tx == 0) / length(tx)) > ratio) {
      tx <- as.sparse(tx)
    }
    outs[["transcripts"]] <- tx
  }

  # --- Spatial coordinate types ---
  for (otype in type) {
    outs[[otype]] <- switch(otype,

      centroids = {
        message("Creating centroid coordinates")
        data.frame(x = sp$center_x, y = sp$center_y, cell = rownames(sp),
                   stringsAsFactors = FALSE)
      },

      segmentations = {
        if (use_hdf5) {
          # Read polygon vertices from per-FOV HDF5 files
          message("Creating polygon coordinates from HDF5")
          pg <- lapply(unique(sp$fov), function(f) {
            fname <- file.path(h5dir, paste0("feature_data_", f, ".hdf5"))
            if (!file.exists(fname)) {
              warning("Cannot find HDF5 file for FOV ", f, immediate. = TRUE)
              return(NULL)
            }
            hfile <- hdf5r::H5File$new(filename = fname, mode = "r")
            on.exit(hfile$close_all())
            cells <- rownames(sp)[sp$fov == f]
            df <- lapply(cells, function(x) {
              tryCatch({
                cc <- hfile[["featuredata"]][[x]][[zidx]][["p_0"]][["coordinates"]]$read()
                cc <- as.data.frame(t(cc))
                colnames(cc) <- c("x", "y")
                cc$cell <- x
                cc
              }, error = function(...) NULL)
            })
            do.call("rbind", df)
          })
          pg <- do.call("rbind", pg)
          npg <- length(unique(pg$cell))
          if (npg < nrow(sp)) {
            warning(nrow(sp) - npg, " cells missing polygon information", immediate. = TRUE)
          }
          pg
        } else {
          # Parquet fallback: polygon geometry stored as WKB blobs
          message("Creating polygon coordinates from cell_boundaries.parquet")
          # Extract EntityID as character directly from the arrow column BEFORE
          # as.data.frame(); arrow converts int64 to R double during that step,
          # which silently corrupts 18-digit IDs (beyond 2^53 precision)
          bnd_arrow       <- arrow::read_parquet(parquet_file)
          entity_ids_char <- as.character(bnd_arrow$EntityID)
          boundaries      <- as.data.frame(bnd_arrow)
          boundaries$EntityID <- entity_ids_char
          boundaries <- boundaries[boundaries$Type == "cell" & boundaries$ZIndex == z, ]
          geom       <- sf::st_as_sfc(boundaries$Geometry, crs = NA)
          coords_mat <- as.data.frame(sf::st_coordinates(geom))
          # L2 is the 1-based feature index: maps each vertex row back to its EntityID
          entity_ids <- boundaries$EntityID[coords_mat$L2]
          data.frame(x    = coords_mat$X,
                     y    = coords_mat$Y,
                     cell = entity_ids,
                     stringsAsFactors = FALSE)
        }
      },

      boxes = {
        # Bounding-box fallback (min_x, max_x, min_y, max_y from cell_metadata.csv)
        message("Creating bounding box coordinates")
        bx <- lapply(rownames(sp), function(cell) {
          row <- sp[cell, , drop = FALSE]
          df <- expand.grid(
            x    = c(row$min_x, row$max_x),
            y    = c(row$min_y, row$max_y),
            cell = cell,
            KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
          )
          df[c(1L, 3L, 4L, 2L), , drop = FALSE]
        })
        do.call("rbind", bx)
      }
    )
  }

  # --- Molecule coordinates ---
  if (length(mol.type) > 0L && !is.na(f_mols)) {
    for (mtype in mol.type) {
      outs[[mtype]] <- switch(mtype,
        microns = {
          message("Creating micron-level molecule coordinates")
          data.frame(x = mx$global_x, y = mx$global_y, gene = mx$gene,
                     stringsAsFactors = FALSE)
        },
        pixels = {
          message("Creating pixel-level molecule coordinates")
          data.frame(x = mx$x, y = mx$y, gene = mx$gene, stringsAsFactors = FALSE)
        }
      )
    }
  }

  return(outs)
}


myLoadVizgen <- function(data.dir, fov = "fov", assay = "Vizgen",
                         mol.type = "microns", z = 3L) {
  # Loads a MERSCOPE sample into a Seurat object with spatial coordinates and a
  # Blanks assay for Blank- negative control genes. Uses centroids only — sufficient
  # for BANKSY neighbourhood extraction and avoids ZIndex mismatches in
  # cell_boundaries.parquet across different segmentation methods.

  data <- myReadVizgen(
    data.dir = data.dir, z = z,
    type     = "centroids",
    mol.type = mol.type
  )

  # Create Seurat object from the full counts matrix (Blank- genes included at this stage)
  obj <- CreateSeuratObject(counts = data[["transcripts"]], assay = assay)

  # Build spatial FOV from centroids
  message("Building FOV from centroids")
  cents  <- CreateCentroids(data[["centroids"]])
  coords <- CreateFOV(
    coords    = list(centroids = cents),
    type      = "centroids",
    molecules = data[[mol.type]],
    assay     = assay
  )
  coords <- subset(coords,
                   cells = intersect(Cells(coords[["centroids"]]), Cells(obj)))

  # Separate Blank- negative control genes into a dedicated assay
  counts      <- GetAssayData(obj, assay = assay, slot = "counts")
  blank_genes <- grep("^Blank-", rownames(counts), value = TRUE)
  if (length(blank_genes) > 0) {
    message(length(blank_genes), " Blank- genes moved to 'Blanks' assay")
    obj[["Blanks"]] <- CreateAssayObject(
      counts = counts[blank_genes, , drop = FALSE]
    )
    obj[[assay]] <- CreateAssayObject(
      counts = counts[setdiff(rownames(counts), blank_genes), , drop = FALSE]
    )
  }

  # Seurat requires valid R names for FOV slots (no underscores or hyphens)
  fov      <- gsub("[_-]", ".", fov)
  obj[[fov]] <- coords

  return(obj)
}


myLoadProseg <- function(data.dir, fov = "fov", assay = c("Vizgen", "Xenium")) {
  # Loads a Proseg-segmented sample (Xenium or MERSCOPE) into a Seurat object.
  # Polygon boundaries are read from cell-polygons.geojson.gz; MULTIPOLYGON cells are
  # resolved by keeping the polygon with the most vertices. Blank genes are moved to a
  # dedicated Blanks assay for consistency with myLoadXenium and myLoadVizgen.
  # Adapted from ProsegToSeurat (proseg2seurat.r) by Ji Zhang, with parallelism replaced
  # by base-R lapply and sf-based polygon extraction.

  assay <- match.arg(assay)

  # Helper: locate csv.gz or parquet fallback for a given Proseg output basename
  find_proseg_file <- function(basename) {
    csv_path     <- file.path(data.dir, paste0(basename, ".csv.gz"))
    parquet_path <- file.path(data.dir, paste0(basename, ".parquet"))
    if (file.exists(csv_path))     return(list(path = csv_path,     fmt = "csv"))
    if (file.exists(parquet_path)) return(list(path = parquet_path, fmt = "parquet"))
    NULL
  }

  read_proseg_file <- function(basename, required = TRUE) {
    f <- find_proseg_file(basename)
    if (is.null(f)) {
      if (required) stop("Cannot find ", basename, " (.csv.gz or .parquet) in ", data.dir)
      return(NULL)
    }
    if (f$fmt == "csv") {
      return(as.data.frame(data.table::fread(f$path, data.table = FALSE)))
    }
    as.data.frame(arrow::read_parquet(f$path))
  }

  # 1. Expected counts: rows = cells, columns = genes (no explicit cell-ID column)
  message("Loading expected counts")
  counts_df <- read_proseg_file("expected-counts")

  # 2. Cell metadata: includes cell, centroid_x, centroid_y
  message("Loading cell metadata")
  cell_metadata <- read_proseg_file("cell-metadata")

  # 3. Transcript metadata (optional): provides molecule coordinates
  message("Loading transcript metadata")
  tx_meta   <- read_proseg_file("transcript-metadata", required = FALSE)
  molecules <- NULL
  if (!is.null(tx_meta)) {
    tx_meta   <- tx_meta[!grepl("^Blank", tx_meta$gene), , drop = FALSE]
    molecules <- data.frame(x = tx_meta$x, y = tx_meta$y, gene = tx_meta$gene,
                            stringsAsFactors = FALSE)
  } else {
    message("Transcript metadata not found; Seurat object will have no molecules slot")
  }

  # 4. Build counts matrix (genes x cells) and create Seurat object
  message("Creating Seurat object")
  counts_matrix <- Matrix::Matrix(t(as.matrix(counts_df)), sparse = TRUE)
  obj           <- CreateSeuratObject(counts = counts_matrix, meta.data = cell_metadata,
                                      assay = assay)
  # Rename cells from integer row indices to actual cell IDs from metadata
  colnames(obj) <- obj$cell

  # 5. Load cell polygon boundaries from compressed GeoJSON
  message("Loading cell polygons")
  polygons_path <- file.path(data.dir, "cell-polygons.geojson.gz")
  if (!file.exists(polygons_path)) {
    stop("Cannot find cell-polygons.geojson.gz in ", data.dir)
  }
  con      <- gzfile(polygons_path, "rt")
  geo_text <- paste(readLines(con, warn = FALSE), collapse = "\n")
  close(con)
  polygons_sf <- geojsonsf::geojson_sf(geo_text)

  # Cast to POLYGON (splits any MULTIPOLYGON into individual polygon rows)
  polygons_sf <- sf::st_cast(polygons_sf, "POLYGON", warn = FALSE)

  # For cells with multiple polygons (split from MULTIPOLYGON), keep the largest
  cell_ids <- as.character(polygons_sf$cell)
  n_pts    <- sapply(polygons_sf$geometry, function(g) nrow(g[[1]]))
  keep_idx <- tapply(seq_len(nrow(polygons_sf)), cell_ids,
                     function(idx) idx[which.max(n_pts[idx])])
  polygons_sf <- polygons_sf[unlist(keep_idx), ]

  # Extract x/y vertex coordinates: L1 = ring (1 = exterior), L2 = feature index
  coords_mat <- as.data.frame(sf::st_coordinates(polygons_sf$geometry))
  segs_df    <- data.frame(
    x    = coords_mat$X,
    y    = coords_mat$Y,
    cell = as.character(polygons_sf$cell[coords_mat$L2]),
    stringsAsFactors = FALSE
  )
  message("Cell polygons loaded")

  # 6. Build FOV with centroids and segmentation polygons
  centroids_df <- data.frame(
    x    = obj$centroid_x,
    y    = obj$centroid_y,
    cell = colnames(obj),
    stringsAsFactors = FALSE
  )
  cents <- CreateCentroids(centroids_df)
  segs  <- CreateSegmentation(segs_df)

  # Retain only cells present in both counts and segmentation
  shared_cells <- intersect(Cells(obj), Cells(segs))
  cents <- subset(cents, cells = shared_cells)
  segs  <- subset(segs,  cells = shared_cells)
  obj   <- subset(obj,   cells = shared_cells)

  coords <- CreateFOV(
    coords    = list(centroids = cents, segmentation = segs),
    type      = c("segmentation", "centroids"),
    molecules = molecules,
    assay     = assay
  )

  # 7. Separate Blank genes into a dedicated assay (consistent with other platforms)
  counts_full <- GetAssayData(obj, assay = assay, slot = "counts")
  blank_genes <- grep("^Blank", rownames(counts_full), value = TRUE)
  if (length(blank_genes) > 0) {
    message(length(blank_genes), " Blank genes moved to 'Blanks' assay")
    obj[["Blanks"]] <- CreateAssayObject(
      counts = counts_full[blank_genes, , drop = FALSE]
    )
    obj[[assay]] <- CreateAssayObject(
      counts = counts_full[setdiff(rownames(counts_full), blank_genes), , drop = FALSE]
    )
  }

  obj[[fov]] <- coords

  return(obj)
}
