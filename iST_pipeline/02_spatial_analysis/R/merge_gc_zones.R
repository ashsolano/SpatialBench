# Purpose: Merge GC B cell zone assignments (Dark Zone / Light Zone) back into
#          the full annotated Seurat object. Reads a gc_subcluster->zone mapping
#          from config["gc_zones"][method], assigns gc_zone to the GC subset,
#          injects it into the full object via AddMetaData (NA for non-GC cells),
#          and builds cell_type2 — identical to cell_type for all non-GC cells,
#          replaced with the zone label for GC B cells.
# Inputs:  annotated.rds   — full object from annotate_clusters.R
#                            (must contain cell_type column)
#          subclustered_gc.rds — GC subset from subcluster_gc.R
#                            (must contain gc_subcluster column)
#          --gc_zones: JSON string e.g. '{"0":"Dark Zone","2":"Light Zone"}'
#                      produced by json.dumps() in the Snakemake rule
# Outputs: annotated_final.rds — full object with gc_zone and cell_type2 added.
#          If gc_zones mapping is empty the script completes without error:
#          gc_zone is NA for all cells and cell_type2 == cell_type.


library(optparse)
library(Seurat)
library(jsonlite)
library(dplyr)

option_list <- list(
  make_option(c("--input_annotated"), type = "character", default = NULL,
              help = "Path to annotated.rds from annotate_clusters.R"),
  make_option(c("--input_gc"),        type = "character", default = NULL,
              help = "Path to subclustered_gc.rds from subcluster_gc.R"),
  make_option(c("--method"),          type = "character", default = NULL,
              help = "Method name (e.g. xenium_default); used for logging"),
  make_option(c("--gc_zones"),        type = "character", default = NULL,
              help = "JSON string mapping gc_subcluster IDs to zone labels"),
  make_option(c("--out_file"),        type = "character", default = NULL,
              help = "Full path for the output annotated_final.rds")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_annotated)) stop("--input_annotated is required")
if (is.null(opt$input_gc))        stop("--input_gc is required")
if (is.null(opt$method))          stop("--method is required")
if (is.null(opt$gc_zones))        stop("--gc_zones is required")
if (is.null(opt$out_file))        stop("--out_file is required")

message("Method:          ", opt$method)
message("Input annotated: ", opt$input_annotated)
message("Input GC subset: ", opt$input_gc)
message("Output file:     ", opt$out_file)

# --- Load full annotated object ---
message("\nLoading full annotated object...")
obj <- readRDS(opt$input_annotated)
message("  Loaded: ", ncol(obj), " cells")

if (!"cell_type" %in% colnames(obj@meta.data)) {
  stop("'cell_type' column not found. Ensure annotate_clusters.R ran successfully.")
}

# --- Load GC subset ---
message("Loading GC subset...")
gc_obj <- readRDS(opt$input_gc)
message("  Loaded: ", ncol(gc_obj), " GC B cells")

if (!"gc_subcluster" %in% colnames(gc_obj@meta.data)) {
  stop("'gc_subcluster' column not found. Ensure subcluster_gc.R ran successfully.")
}

# --- Parse gc_zones mapping ---
# fromJSON on '{"0":"Dark Zone"}' returns a named character vector;
# fromJSON on '{}' returns an empty named list — unlist() normalises both
gc_zones <- unlist(jsonlite::fromJSON(opt$gc_zones))
message("  GC zone entries provided: ", length(gc_zones))

# --- Assign gc_zone to the GC subset ---
subcluster_ids <- as.character(gc_obj$gc_subcluster)
n_subclusters  <- length(unique(subcluster_ids))
message("  GC subclusters present: ", n_subclusters,
        " (", paste(sort(unique(subcluster_ids)), collapse = ", "), ")")

if (length(gc_zones) == 0) {
  # No mapping yet — proceed with NA; cell_type2 will equal cell_type for GC cells
  message(
    "\n  Warning: No gc_zones mapping provided for method '", opt$method, "'.",
    "\n  gc_zone will be NA for all GC B cells and cell_type2 will equal cell_type.",
    "\n  Populate config['gc_zones']['", opt$method,
    "'] after inspecting subclustered_gc.rds and re-run."
  )
  gc_zone_vals <- rep(NA_character_, ncol(gc_obj))
} else {
  gc_zone_vals <- gc_zones[subcluster_ids]

  # Subclusters not in the mapping produce NA — flag them
  n_unmapped <- sum(is.na(gc_zone_vals))
  if (n_unmapped > 0) {
    unmapped_ids <- sort(unique(subcluster_ids[is.na(gc_zone_vals)]))
    message(
      "\n  Warning: ", n_unmapped, " GC cells belong to unmapped subclusters: ",
      paste(unmapped_ids, collapse = ", "),
      "\n  These cells will have gc_zone = NA and cell_type2 = 'GC B cells'."
    )
  }

  message("\nGC zone assignments:")
  zone_table <- sort(table(gc_zone_vals, useNA = "ifany"), decreasing = TRUE)
  for (z in names(zone_table)) {
    message("  ", z, ": ", zone_table[[z]])
  }
}

# --- Inject gc_zone into the full object ---
# Build a named vector aligned to GC cell barcodes; AddMetaData fills NA for
# all cells not in the vector (i.e. all non-GC cells)
gc_zone_vec        <- gc_zone_vals
names(gc_zone_vec) <- colnames(gc_obj)

obj <- AddMetaData(obj, metadata = gc_zone_vec, col.name = "gc_zone")
message(
  "\ngc_zone injected into full object:",
  "\n  Non-GC cells (NA):  ", sum(is.na(obj$gc_zone)),
  "\n  GC cells assigned:  ", sum(!is.na(obj$gc_zone))
)

# --- Build cell_type2 ---
# For non-GC cells: cell_type2 == cell_type
# For GC cells with a zone assignment: cell_type2 == zone label
# For GC cells with NA gc_zone (unmapped or empty config): cell_type2 == "GC B cells"
cell_type2 <- as.character(obj$cell_type)
is_gc      <- cell_type2 == "GC B cells"
has_zone   <- is_gc & !is.na(obj$gc_zone)
cell_type2[has_zone] <- obj$gc_zone[has_zone]

# Use sort(unique()) for levels; figure scripts will set display order explicitly
obj$cell_type2 <- factor(cell_type2, levels = sort(unique(cell_type2)))

message("\ncell_type2 distribution:")
ct2_table <- sort(table(obj$cell_type2), decreasing = TRUE)
for (ct in names(ct2_table)) {
  message("  ", ct, ": ", ct2_table[[ct]])
}

# --- Save ---
dir.create(dirname(opt$out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, file = opt$out_file)
message("\nSaved: ", opt$out_file)
message("Done.")
