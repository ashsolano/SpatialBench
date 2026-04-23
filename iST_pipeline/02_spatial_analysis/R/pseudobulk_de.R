# Purpose:  Pseudobulk differential expression (KO vs WT) per cell type using
#           a limma-voom pipeline. Aggregates raw counts per cell type per
#           biological replicate (animal_id), runs edgeR TMM normalisation,
#           and fits a voom-lmFit contrast model. Saves one MArrayLM fit object
#           per qualifying cell type (>= min_reps replicates in both groups).
# Inputs:   annotated.rds from 02_spatial_analysis; requires metadata columns:
#           cell_type (annotate_clusters.R), condition and animal_id (preprocess_sample.R)
# Outputs:  pseudobulk_de/de_results.rds — named list of MArrayLM objects, one
#           per cell type with sufficient replicates in both KO and WT groups


suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(edgeR)
  library(limma)
})

option_list <- list(
  make_option(c("--input_rds"),     type = "character", default = NULL,
              help = "Path to annotated.rds from 02_spatial_analysis"),
  make_option(c("--method"),        type = "character", default = NULL,
              help = "Method name (e.g. xenium_default); used for logging only"),
  make_option(c("--assay"),         type = "character", default = NULL,
              help = "Seurat assay containing raw counts (e.g. Xenium or Vizgen)"),
  make_option(c("--cell_type_col"), type = "character", default = "cell_type",
              help = "Metadata column with cell type labels [default: cell_type]"),
  make_option(c("--condition_col"), type = "character", default = "condition",
              help = "Metadata column with condition labels (wt/ko/ctrl) [default: condition]"),
  make_option(c("--animal_id_col"), type = "character", default = "animal_id",
              help = "Metadata column with biological replicate ID [default: animal_id]"),
  make_option(c("--min_reps"),      type = "integer",   default = 2L,
              help = "Minimum replicates per group required to run DE [default: 2]"),
  make_option(c("--contrast_name"), type = "character", default = "KovsWt",
              help = "Name for the KO-vs-WT contrast coefficient [default: KovsWt]"),
  make_option(c("--out_file"),      type = "character", default = NULL,
              help = "Full output path for de_results.rds")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds))  stop("--input_rds is required")
if (is.null(opt$method))     stop("--method is required")
if (is.null(opt$assay))      stop("--assay is required")
if (is.null(opt$out_file))   stop("--out_file is required")

# Conditions and contrast levels are fixed for this project
keep_conditions <- c("wt", "ko")
contrast_levels <- c("wt", "ko")

message("Method:        ", opt$method)
message("Assay:         ", opt$assay)
message("Cell type col: ", opt$cell_type_col)
message("Condition col: ", opt$condition_col)
message("Animal ID col: ", opt$animal_id_col)
message("Min reps:      ", opt$min_reps)
message("Contrast name: ", opt$contrast_name)
message("Input RDS:     ", opt$input_rds)
message("Output:        ", opt$out_file)

# --- Load annotated object ---
message("\nLoading annotated object...")
seu <- readRDS(opt$input_rds)
message("  Loaded: ", ncol(seu), " cells, ", nrow(seu), " features")

stopifnot(inherits(seu, "Seurat"))

if (!opt$assay %in% names(seu@assays)) {
  stop("Assay '", opt$assay, "' not found. Available: ",
       paste(names(seu@assays), collapse = ", "))
}

# --- Validate required metadata columns ---
meta   <- seu@meta.data
needed <- c(opt$cell_type_col, opt$condition_col, opt$animal_id_col)
miss   <- setdiff(needed, colnames(meta))
if (length(miss)) {
  stop(
    "Missing metadata columns: ", paste(miss, collapse = ", "),
    "\nDid preprocess_sample.R and annotate_clusters.R complete successfully?"
  )
}

# --- Build replicate map: one row per cell with condition + replicate ID ---
rep_map <- meta %>%
  transmute(
    cell      = rownames(meta),
    condition = tolower(as.character(.data[[opt$condition_col]])),
    rep_id    = as.character(.data[[opt$animal_id_col]]),
    ct        = as.character(.data[[opt$cell_type_col]])
  )

message("\nCondition counts (all cells):")
print(table(rep_map$condition, useNA = "ifany"))

# --- Filter to KO and WT cells (exclude ctrl) ---
keep_cells <- rep_map$cell[rep_map$condition %in% keep_conditions]
if (!length(keep_cells)) {
  stop("No KO or WT cells found. Check the condition column.")
}
message("KO/WT cells retained: ", length(keep_cells))

# --- Aggregate raw counts per cell type x animal_id ---
message("\nAggregating counts (this can be the slow step)...")
agg <- AggregateExpression(
  seu,
  assays        = opt$assay,
  slot          = "counts",
  group.by      = c(opt$cell_type_col, opt$animal_id_col),
  return.seurat = FALSE,
  cells         = keep_cells
)[[opt$assay]]
message("  Aggregated matrix: ", nrow(agg), " genes x ", ncol(agg), " pseudobulk samples")
message("  Example columns: ", paste(head(colnames(agg)), collapse = ", "))

# --- Build animal_id -> condition lookup (KO/WT only) ---
samp_cond <- rep_map %>%
  filter(condition %in% keep_conditions) %>%
  distinct(rep_id, condition)

# --- Parse AggregateExpression column names into cell_type + rep_id ---
# AggregateExpression formats columns as "<cell_type>_<animal_id>".
# Split on the last underscore; valid because animal_id values (wt709, ko166)
# contain no underscores.
cn      <- colnames(agg)
last_us <- regexpr("_[^_]+$", cn)
if (!all(last_us > 0)) {
  stop(
    "Failed to parse AggregateExpression column names via last underscore. ",
    "Check cell_type and animal_id values for unexpected underscores."
  )
}
celltype <- substr(cn, 1, last_us - 1)
rep_id   <- substr(cn, last_us + 1, nchar(cn))
col_map  <- tibble(col = cn, cell_type = celltype, rep_id = rep_id)

# --- Build pseudobulk list per cell type, retaining only KO/WT replicates ---
keep_reps <- samp_cond$rep_id
pb_list <- split(col_map, col_map$cell_type) |>
  purrr::imap(function(df, ct) {
    df  <- df[df$rep_id %in% keep_reps, , drop = FALSE]
    mat <- agg[, df$col, drop = FALSE]
    colnames(mat) <- df$rep_id
    mat
  })
pb_list <- pb_list[lengths(pb_list) > 0]

# --- Summarise replicate counts per cell type ---
rep_summary <- purrr::imap_dfr(pb_list, function(mat, ct) {
  reps <- colnames(mat)
  cond <- samp_cond$condition[match(reps, samp_cond$rep_id)]
  tibble(
    cell_type = ct,
    wt        = sum(cond == "wt", na.rm = TRUE),
    ko        = sum(cond == "ko", na.rm = TRUE),
    unlabeled = sum(is.na(cond))
  )
}) %>% arrange(desc(wt + ko))

message("\nReplicate summary:")
print(rep_summary)

# --- Select cell types with sufficient replicates in both groups ---
good_ct <- rep_summary %>%
  filter(wt >= opt$min_reps, ko >= opt$min_reps) %>%
  pull(cell_type)

if (!length(good_ct)) {
  stop(
    "No cell types pass min_reps = ", opt$min_reps,
    ". Check condition assignment and replicate summary printed above."
  )
}
message(
  "\nRunning limma-voom DE on ", length(good_ct), " cell type(s): ",
  paste(good_ct, collapse = ", ")
)

# --- limma-voom pseudobulk DE per qualifying cell type ---
efits <- list()

for (ct in good_ct) {
  counts <- pb_list[[ct]]
  reps   <- colnames(counts)
  cond   <- samp_cond$condition[match(reps, samp_cond$rep_id)]

  # Drop replicates with no condition assignment
  keep   <- !is.na(cond)
  counts <- counts[, keep, drop = FALSE]
  cond   <- factor(cond[keep], levels = contrast_levels)

  # edgeR: DGEList, filter low-count genes, TMM normalisation
  dge <- DGEList(counts = counts, group = cond)
  k   <- filterByExpr(dge, group = cond)
  dge <- dge[k, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")

  # Design matrix: one coefficient per condition, no intercept
  design <- model.matrix(~0 + cond)
  colnames(design) <- contrast_levels

  # Contrast matrix: KO - WT
  cm <- matrix(0, nrow = ncol(design), ncol = 1,
               dimnames = list(colnames(design), opt$contrast_name))
  cm["ko", opt$contrast_name] <-  1
  cm["wt", opt$contrast_name] <- -1

  # limma-voom pipeline
  v   <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- contrasts.fit(fit, cm)
  fit <- eBayes(fit, robust = FALSE)
  fit$cluster <- ct

  efits[[ct]] <- fit
  message("  ", ct, ": ", sum(k), " genes tested")
}

# --- Save list of MArrayLM fits ---
dir.create(dirname(opt$out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(efits, file = opt$out_file)
message("\nSaved ", length(efits), " cell type fit(s) to: ", opt$out_file)
message("Done.")
