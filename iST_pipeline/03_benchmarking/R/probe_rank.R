# Purpose:  Compute per-probe rank tables and background overlap gene lists for
#           MERSCOPE and Xenium 8µm binning objects. Identifies target probes
#           whose mean expression falls below the 95th percentile of blank/
#           background probes (the "overlap" set), and builds the label pool
#           used for S-curve annotations in figure generation.
#           Adapted from probe_rank_platform_v2.R.
# Inputs:   config/config.yaml  (spatial_analysis, output_dir)
#           results/01_preprocessing/merscope_8um/{sample}_8um.rds
#           results/01_preprocessing/xenium_8um/{sample}_8um.rds
# Outputs:  results/03_benchmarking/probe_rank/ranked_plat.rds
#               (per-probe table: platform, feature, Type, mean_count, rank,
#                rank_frac, bg95, y, is_overlap)
#           results/03_benchmarking/probe_rank/probe_overlap_targets.csv
#               (target probes falling within the bg95 threshold)
#           results/03_benchmarking/probe_rank/overlap_genes_merscope.txt
#           results/03_benchmarking/probe_rank/overlap_genes_xenium.txt
#           results/03_benchmarking/probe_rank/label_pool_genes.csv
#               (union of top-N and tail-N MERSCOPE overlap probes, used for
#                S-curve text annotations in fig2_background.R)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
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
  make_option(c("--out_dir"), type = "character",
              default = "results/03_benchmarking/probe_rank",
              help    = "Output directory [default: %default]"),
  make_option(c("--label_top_n"), type = "integer",
              default = 8L,
              help    = "Top N overlap probes by mean count for label pool [default: %default]"),
  make_option(c("--label_tail_n"), type = "integer",
              default = 6L,
              help    = "Right-tail N overlap probes by rank percentile for label pool [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load 8um binning objects
# ---------------------------------------------------------------------------
merscope_samples <- cfg$spatial_analysis$merscope_samples
xenium_samples   <- cfg$spatial_analysis$xenium_default_samples

message("Loading MERSCOPE 8um objects...")
mer_dir      <- file.path(cfg$output_dir, "01_preprocessing", "merscope_8um")
merscope_8um <- setNames(
  lapply(merscope_samples, function(samp) {
    path <- file.path(mer_dir, paste0(samp, "_8um.rds"))
    message("  ", samp, ": ", path)
    readRDS(path)
  }),
  merscope_samples
)

message("Loading Xenium 8um objects...")
xen_dir    <- file.path(cfg$output_dir, "01_preprocessing", "xenium_8um")
xenium_8um <- setNames(
  lapply(xenium_samples, function(samp) {
    path <- file.path(xen_dir, paste0(samp, "_8um.rds"))
    message("  ", samp, ": ", path)
    readRDS(path)
  }),
  xenium_samples
)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Seurat v4/v5 compatible raw counts extraction; returns NULL if assay absent
get_counts_safe <- function(sobj, assay, slot = "counts") {
  if (!assay %in% names(sobj@assays)) return(NULL)
  tryCatch(
    GetAssayData(sobj, assay = assay, layer = slot),
    error = function(e) tryCatch(
      GetAssayData(sobj, assay = assay, slot = slot),
      error = function(e2) NULL
    )
  )
}

# Row-sum dispatch for sparse and dense matrices
rsum <- function(m) {
  if (inherits(m, "dgCMatrix") || inherits(m, "dgRMatrix")) Matrix::rowSums(m)
  else base::rowSums(m)
}

# Element-wise addition of named vectors with union of keys (zero-fills missing)
sum_named_vectors <- function(x, y) {
  allg <- union(names(x), names(y))
  x2 <- setNames(numeric(length(allg)), allg); x2[names(x)] <- x
  y2 <- setNames(numeric(length(allg)), allg); y2[names(y)] <- y
  x2 + y2
}

# Detect the Vizgen/MERSCOPE target assay by name; falls back to largest assay
pick_vizgen_assay <- function(sobj) {
  nm   <- names(sobj@assays)
  cand <- nm[nm %in% c("Vizgen", "VIZGEN", "vizgen")]
  if (length(cand) > 0) return(cand[1])
  if (length(nm)   == 1) return(nm[1])
  sizes <- vapply(nm, function(a) {
    m <- get_counts_safe(sobj, a); if (is.null(m)) 0L else nrow(m)
  }, numeric(1))
  nm[which.max(sizes)]
}

# Per-sample target/blank probe tally for one Seurat object.
# Returns a tidy tibble: platform, feature, count, Sample, Type.
summarise_sample_probes <- function(sobj, sample_name, platform) {

  if (platform == "Xenium") {
    blank_assays <- c("Unassigned", "ControlCodeword", "ControlProbe")
    target_mat   <- get_counts_safe(sobj, "Xenium")
    target_df    <- NULL
    if (!is.null(target_mat)) {
      targ_counts <- rsum(target_mat)
      target_df <- tibble(
        feature  = names(targ_counts), count = as.numeric(targ_counts),
        platform = "Xenium", Sample = sample_name, Type = "Target"
      )
    }
    blank_mats <- purrr::compact(lapply(blank_assays, get_counts_safe, sobj = sobj))
    blank_df   <- NULL
    if (length(blank_mats) > 0) {
      bl_counts <- Reduce(
        function(u, v) sum_named_vectors(u, rsum(v)),
        blank_mats,
        init = setNames(numeric(0), character(0))
      )
      blank_df <- tibble(
        feature  = names(bl_counts), count = as.numeric(bl_counts),
        platform = "Xenium", Sample = sample_name, Type = "Blank"
      )
    }
    return(bind_rows(target_df, blank_df))
  }

  if (platform == "MERSCOPE") {
    viz_assay <- pick_vizgen_assay(sobj)
    viz_mat   <- get_counts_safe(sobj, viz_assay)
    blank_mat <- get_counts_safe(sobj, "BlankProbe")
    if (is.null(viz_mat)) return(NULL)

    feats_v      <- rownames(viz_mat)
    # Detect blank rows by name when no dedicated BlankProbe assay is present
    is_blank_row <- if (is.null(blank_mat))
      grepl("^(Blank|BlankProbe)", feats_v, ignore.case = TRUE)
    else
      rep(FALSE, length(feats_v))

    targ_counts <- rsum(viz_mat[!is_blank_row, , drop = FALSE])
    targ_df <- tibble(
      feature  = names(targ_counts), count = as.numeric(targ_counts),
      platform = "MERSCOPE", Sample = sample_name, Type = "Target"
    )

    blank_df_from_viz <- NULL
    if (any(is_blank_row)) {
      bl_viz <- rsum(viz_mat[is_blank_row, , drop = FALSE])
      blank_df_from_viz <- tibble(
        feature  = names(bl_viz), count = as.numeric(bl_viz),
        platform = "MERSCOPE", Sample = sample_name, Type = "Blank"
      )
    }
    blank_df_from_assay <- NULL
    if (!is.null(blank_mat)) {
      bl_assay <- rsum(blank_mat)
      blank_df_from_assay <- tibble(
        feature  = names(bl_assay), count = as.numeric(bl_assay),
        platform = "MERSCOPE", Sample = sample_name, Type = "Blank"
      )
    }
    return(bind_rows(targ_df, blank_df_from_viz, blank_df_from_assay))
  }

  NULL
}

# ---------------------------------------------------------------------------
# Build per-sample probe count tables
# ---------------------------------------------------------------------------
message("Summarising per-sample probe counts (MERSCOPE)...")
mer_tbl <- purrr::map_dfr(
  names(merscope_8um),
  ~ summarise_sample_probes(merscope_8um[[.x]], .x, "MERSCOPE")
)

message("Summarising per-sample probe counts (Xenium)...")
xen_tbl <- purrr::map_dfr(
  names(xenium_8um),
  ~ summarise_sample_probes(xenium_8um[[.x]], .x, "Xenium")
)

all_counts <- bind_rows(mer_tbl, xen_tbl)

# ---------------------------------------------------------------------------
# Compute platform-level ranked probe table
# ---------------------------------------------------------------------------
# Mean count per probe averaged across all samples within each platform
plat_counts <- all_counts %>%
  group_by(platform, feature, Type) %>%
  summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop")

# 95th percentile of blank/background probes per platform — the overlap threshold
thr_platform <- plat_counts %>%
  filter(Type == "Blank") %>%
  group_by(platform) %>%
  summarise(bg95 = quantile(mean_count, 0.95, na.rm = TRUE), .groups = "drop")

# Rank each probe within its platform (descending by mean count)
ranked_plat <- plat_counts %>%
  group_by(platform) %>%
  arrange(desc(mean_count), .by_group = TRUE) %>%
  mutate(
    rank      = row_number(),
    rank_frac = rank / max(rank, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  left_join(thr_platform, by = "platform") %>%
  mutate(
    y          = log10(pmax(mean_count, 0) + 1),
    is_overlap = (Type == "Target" & mean_count <= bg95)
  )

message("Probe rank table: ", nrow(ranked_plat), " rows")

# ---------------------------------------------------------------------------
# Build label pool for S-curve annotations
# ---------------------------------------------------------------------------
# Label pool is the union of:
#   - top N MERSCOPE overlap targets by mean count
#   - right-tail N MERSCOPE overlap targets by rank percentile
# Applied to both platforms for consistent annotation across facets.
mer_overlap <- ranked_plat %>%
  filter(platform == "MERSCOPE", Type == "Target", is_overlap)

top_set  <- if (nrow(mer_overlap) > 0)
  slice_max(mer_overlap, order_by = mean_count, n = opt$label_top_n,  with_ties = FALSE)
else
  mer_overlap

tail_set <- if (nrow(mer_overlap) > 0)
  slice_max(mer_overlap, order_by = rank_frac,  n = opt$label_tail_n, with_ties = FALSE)
else
  mer_overlap

label_pool  <- union(top_set$feature, tail_set$feature)
label_genes <- ranked_plat %>%
  filter(feature %in% label_pool) %>%
  distinct(feature) %>%
  arrange(feature)

# ---------------------------------------------------------------------------
# Build overlap target table
# ---------------------------------------------------------------------------
probe_overlap_targets <- ranked_plat %>%
  filter(Type == "Target", is_overlap) %>%
  arrange(platform, desc(mean_count)) %>%
  select(platform, feature, mean_count, rank, rank_frac, bg95)

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
saveRDS(ranked_plat, file.path(opt$out_dir, "ranked_plat.rds"))
message("Saved: ranked_plat.rds")

utils::write.csv(probe_overlap_targets,
                 file.path(opt$out_dir, "probe_overlap_targets.csv"),
                 row.names = FALSE)
message("Saved: probe_overlap_targets.csv")

utils::write.csv(label_genes,
                 file.path(opt$out_dir, "label_pool_genes.csv"),
                 row.names = FALSE)
message("Saved: label_pool_genes.csv")

# Per-platform gene lists: one gene per line, no header
probe_overlap_targets %>%
  group_by(platform) %>%
  summarise(feature = sort(unique(feature)), .groups = "drop") %>%
  split(.$platform) %>%
  purrr::iwalk(function(df, plat) {
    out_file <- file.path(opt$out_dir, paste0("overlap_genes_", tolower(plat), ".txt"))
    utils::write.table(df["feature"], file = out_file,
                       quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("Saved: overlap_genes_", tolower(plat), ".txt")
  })

message("Done. Outputs written to: ", opt$out_dir)
