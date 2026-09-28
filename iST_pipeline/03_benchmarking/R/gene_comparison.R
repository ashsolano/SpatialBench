# Purpose:  Pseudobulk 8µm-bin counts per animal and platform (VisiumHD, MERSCOPE,
#           Xenium) for the four animals on all three platforms (WT709, WT713,
#           KO167, KO168), restricted to genes present in all 12 datasets, and
#           export an edgeR DGEList for the gene-comparison panels.
# Inputs:   config/config.yaml  (gene_comparison$animals, visiumhd, output_dir)
#           cfg$visiumhd$data_dir / cfg$visiumhd$samples  (VisiumHD 8µm Seurat
#               objects, "Spatial.008um" assay)
#           results/01_preprocessing/merscope_8um_filtered/{sample}_8um_filtered.rds
#               ("Vizgen" assay)
#           results/01_preprocessing/xenium_8um_filtered/{sample}_8um_filtered.rds
#               ("Xenium" assay)
# Outputs:  results/03_benchmarking/gene_comparison/dge.rds
#           results/03_benchmarking/gene_comparison/counts_mat.rds
#               (genes x 12 samples; columns named "{Platform}_{SampleID}")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(yaml)
  library(optparse)
})

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--config"),  type = "character",
              default = "config/config.yaml",
              help    = "Path to config.yaml [default: %default]"),
  make_option(c("--out_dir"), type = "character",
              default = "results/03_benchmarking/gene_comparison",
              help    = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Sample table: one row per animal x platform
# ---------------------------------------------------------------------------
# config["gene_comparison"]["animals"] keys are {CONDITION}{ID} (e.g. "WT709");
# condition and animal ID come from the key, so metadata does not depend on
# parsing platform-specific sample names.
bin_res  <- cfg$bin_resolutions[[1]]   # 8
platform_info <- tribble(
  ~Platform,  ~cfg_key,    ~assay,
  "VisiumHD", "visiumhd",  "Spatial.008um",
  "MERSCOPE", "merscope",  "Vizgen",
  "Xenium",   "xenium",    "Xenium"
)

sample_path <- function(cfg_key, sample_id) {
  if (cfg_key == "visiumhd") {
    file.path(cfg$visiumhd$data_dir, cfg$visiumhd$samples[[sample_id]])
  } else {
    file.path(cfg$output_dir, "01_preprocessing",
              paste0(cfg_key, "_", bin_res, "um_filtered"),
              paste0(sample_id, "_", bin_res, "um_filtered.rds"))
  }
}

samples_tbl <- imap_dfr(cfg$gene_comparison$animals, function(plat_samples, animal) {
  platform_info %>%
    mutate(
      Animal   = animal,
      SampleID = map_chr(cfg_key, ~ plat_samples[[.x]] %||% NA_character_)
    )
}) %>%
  mutate(
    Type  = sub("^([A-Z]+)\\d+$", "\\1", Animal),
    IDnum = sub("^[A-Z]+(\\d+)$", "\\1", Animal),
    path  = map2_chr(cfg_key, SampleID, ~ if (is.na(.y)) NA_character_ else sample_path(.x, .y)),
    Platform = factor(Platform, levels = platform_info$Platform)
  ) %>%
  arrange(Platform, match(Animal, names(cfg$gene_comparison$animals)))

# ---------------------------------------------------------------------------
# Validate the sample table before loading anything
# ---------------------------------------------------------------------------
# Requires exactly 4 matched WT/KO animals per platform with existing inputs
validate_samples <- function(tbl) {
  if (anyNA(tbl$SampleID)) {
    stop("Missing sample name in gene_comparison$animals for: ",
         paste(tbl$Animal[is.na(tbl$SampleID)], tbl$Platform[is.na(tbl$SampleID)], collapse = ", "))
  }
  per_platform <- count(tbl, Platform)
  if (any(per_platform$n != 4)) stop("Expected exactly 4 animals per platform, got: ",
                                     paste(per_platform$Platform, per_platform$n, collapse = ", "))
  if (nrow(tbl) != 12) stop("Expected exactly 12 samples, got ", nrow(tbl))
  animal_sets <- split(tbl$Animal, tbl$Platform) %>% map(sort)
  if (!all(map_lgl(animal_sets, identical, animal_sets[[1]]))) {
    stop("Animal IDs differ between platforms")
  }
  if (!all(tbl$Type %in% c("WT", "KO"))) stop("Unexpected condition in animal keys: ",
                                              paste(unique(tbl$Animal[!tbl$Type %in% c("WT", "KO")]), collapse = ", "))
  col_names <- paste(tbl$Platform, tbl$SampleID, sep = "_")
  if (anyDuplicated(col_names)) stop("Duplicated sample IDs: ", paste(col_names[duplicated(col_names)], collapse = ", "))
  missing <- tbl$path[!file.exists(tbl$path)]
  if (length(missing) > 0) stop("Missing input objects:\n  ", paste(missing, collapse = "\n  "))
  invisible(TRUE)
}
validate_samples(samples_tbl)

message("Selected samples:")
walk(seq_len(nrow(samples_tbl)), function(i) {
  with(samples_tbl[i, ], message(sprintf("  %-8s %-6s %-16s %s", Platform, Animal, SampleID, path)))
})

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Retrieve the counts layer/slot, supporting both Seurat v4 and v5 APIs.
get_counts_compat <- function(sobj, assay) {
  tryCatch(
    GetAssayData(sobj, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(sobj, assay = assay, slot = "counts")
  )
}

# Raw per-gene counts summed across all bins (empty bins add zero, so no bin
# filtering is needed). Objects are loaded one at a time so only one
# (up to ~1.9 GB) is in memory at once.
pseudobulk_one <- function(path, assay, label) {
  message("  ", label, ": loading ", path)
  sobj <- readRDS(path)
  if (!assay %in% Assays(sobj)) stop(label, " has no ", assay, " assay")
  mat  <- get_counts_compat(sobj, assay)
  message("    ", ncol(mat), " bins, ", nrow(mat), " genes")
  cnts <- Matrix::rowSums(mat)
  rm(sobj, mat); gc(verbose = FALSE)
  cnts
}

# ---------------------------------------------------------------------------
# Build per-sample pseudobulk vectors
# ---------------------------------------------------------------------------
message("Building pseudobulk counts per sample...")

samples_tbl <- samples_tbl %>%
  mutate(ColName = paste(Platform, SampleID, sep = "_"))

pb_list <- pmap(samples_tbl, function(path, assay, ColName, ...) {
  pseudobulk_one(path, assay, ColName)
}) %>% setNames(samples_tbl$ColName)

# ---------------------------------------------------------------------------
# Identify genes present in all 12 selected datasets
# ---------------------------------------------------------------------------
message("Finding genes common to all ", length(pb_list), " datasets...")

genes_by_platform <- split(pb_list, samples_tbl$Platform) %>%
  map(~ Reduce(intersect, map(.x, names)))
iwalk(genes_by_platform, ~ message("  ", .y, " genes (all 4 animals): ", length(.x)))

common_genes <- Reduce(intersect, map(pb_list, names))
message("  Common genes: ", length(common_genes))
if (length(common_genes) == 0) stop("No genes common to all selected datasets")

counts_mat <- map(pb_list, ~ .x[common_genes]) %>%
  do.call(cbind, .)
rownames(counts_mat) <- common_genes

validate_counts <- function(mat, pb_list, genes) {
  if (ncol(mat) != 12) stop("Expected 12 pseudobulk samples, got ", ncol(mat))
  if (anyDuplicated(colnames(mat))) stop("Duplicated sample IDs in count matrix")
  gene_ok <- map_lgl(pb_list, ~ all(genes %in% names(.x)))
  if (!all(gene_ok)) stop("Common genes missing from: ", paste(names(pb_list)[!gene_ok], collapse = ", "))
  if (anyNA(mat)) stop("Count matrix contains missing values")
  if (any(mat < 0)) stop("Count matrix contains negative values")
  if (any(mat != round(mat))) stop("Count matrix contains non-integer values (expected raw counts)")
  invisible(TRUE)
}
validate_counts(counts_mat, pb_list, common_genes)

# ---------------------------------------------------------------------------
# Sample metadata
# ---------------------------------------------------------------------------
# Type keeps the WT/KO/CTRL levels used by the figure scripts; no CTRL animals
# are included.
col_info <- samples_tbl %>%
  transmute(
    ColName,
    SampleID,
    Platform = as.character(Platform),
    Type     = factor(Type, levels = c("WT", "KO", "CTRL")),
    IDnum
  ) %>%
  column_to_rownames("ColName")
stopifnot(identical(rownames(col_info), colnames(counts_mat)))

message("Sample metadata:")
print(col_info)

dge <- DGEList(counts = counts_mat, samples = col_info)
message("DGEList: ", nrow(dge), " genes x ", ncol(dge), " samples")

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
dge_path  <- file.path(opt$out_dir, "dge.rds")
mat_path  <- file.path(opt$out_dir, "counts_mat.rds")

saveRDS(dge,        dge_path)
saveRDS(counts_mat, mat_path)

message("Saved: ", dge_path)
message("Saved: ", mat_path)
message("Done.")
