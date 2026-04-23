# Purpose: Assign cell type annotations to clusters in an embedded Seurat
#          object. Maps cluster_harmony_sampleid labels to cell type names
#          using a user-supplied cluster->cell_type JSON mapping from
#          config["annotations"][method].
# Inputs:  Embedded Seurat object from embed_harmony.R (must contain
#          cluster_harmony_sampleid metadata column)
#          --annotations: JSON string e.g. '{"0":"GC B cells","1":"T cells"}'
#          produced by json.dumps() in the Snakemake rule
# Outputs: Annotated Seurat object with cell_type column added to metadata.
#          Unmapped clusters are labelled "Unannotated" with a warning.


library(optparse)
library(jsonlite)
library(dplyr)

option_list <- list(
  make_option(c("--input_rds"),   type = "character", default = NULL,
              help = "Path to embedded RDS from embed_harmony.R"),
  make_option(c("--method"),      type = "character", default = NULL,
              help = "Method name (e.g. xenium_default); used for logging"),
  make_option(c("--annotations"), type = "character", default = NULL,
              help = "JSON string mapping cluster IDs to cell type names"),
  make_option(c("--out_file"),    type = "character", default = NULL,
              help = "Full path for the annotated output RDS file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds))   stop("--input_rds is required")
if (is.null(opt$method))      stop("--method is required")
if (is.null(opt$annotations)) stop("--annotations is required")
if (is.null(opt$out_file))    stop("--out_file is required")

message("Method:      ", opt$method)
message("Input RDS:   ", opt$input_rds)
message("Output file: ", opt$out_file)

# --- Load embedded object ---
message("\nLoading embedded object...")
obj <- readRDS(opt$input_rds)
message("  Loaded: ", ncol(obj), " cells")

# --- Validate prerequisites ---
if (!"cluster_harmony_sampleid" %in% colnames(obj@meta.data)) {
  stop(
    "'cluster_harmony_sampleid' not found in metadata. ",
    "Ensure embed_harmony.R ran successfully."
  )
}

# --- Parse annotation mapping ---
# fromJSON on '{"0":"GC B cells"}' returns a named character vector;
# fromJSON on '{}' returns an empty named list — unlist() normalises both
annotations <- unlist(jsonlite::fromJSON(opt$annotations))
message("  Annotation entries provided: ", length(annotations))

cluster_ids <- as.character(obj$cluster_harmony_sampleid)
n_clusters  <- length(unique(cluster_ids))
message("  Clusters present in object: ", n_clusters,
        " (", paste(sort(unique(cluster_ids)), collapse = ", "), ")")

# --- Map clusters to cell types ---
if (length(annotations) == 0) {
  # No annotations supplied yet — flag and label all cells as Unannotated
  message(
    "\n  Warning: No annotations provided for method '", opt$method, "'.",
    "\n  All cells will be labelled 'Unannotated'.",
    "\n  Populate config['annotations']['", opt$method, "'] and re-run."
  )
  obj$cell_type <- "Unannotated"
} else {
  # Look up each cell's cluster in the annotations vector
  cell_types <- annotations[cluster_ids]

  # Clusters not present in the mapping produce NA -> flag and fill
  n_unmapped <- sum(is.na(cell_types))
  if (n_unmapped > 0) {
    unmapped_ids <- sort(unique(cluster_ids[is.na(cell_types)]))
    message(
      "\n  Warning: ", n_unmapped, " cells belong to unmapped clusters: ",
      paste(unmapped_ids, collapse = ", "),
      "\n  These cells will be labelled 'Unannotated'."
    )
    cell_types[is.na(cell_types)] <- "Unannotated"
  }

  obj$cell_type <- unname(cell_types)

  # Report final distribution
  message("\nCell type assignments:")
  ct_table <- sort(table(obj$cell_type), decreasing = TRUE)
  for (ct in names(ct_table)) {
    message("  ", ct, ": ", ct_table[[ct]])
  }
}

# --- Save ---
dir.create(dirname(opt$out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, file = opt$out_file)
message("\nSaved: ", opt$out_file)
message("Done.")
