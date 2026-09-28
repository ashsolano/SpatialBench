# Purpose:  Per-probe rank tables for MERSCOPE and Xenium 8µm bins, the target
#           probes at or below the 95th percentile of background probes (the
#           "overlap" set), and the S-curve label pool for fig2_background.R.
#           Background sets match qc_backgrounds.R / moransi.R: MERSCOPE "Blanks"
#           assay (89 blank probes); Xenium "BlankCodeword" assay (380
#           Unassigned codewords only; control probes/codewords excluded).
# Inputs:   config/config.yaml  (spatial_analysis, qc_backgrounds$exclude_samples,
#                                output_dir)
#           results/01_preprocessing/merscope_8um_filtered/{sample}_8um_filtered.rds
#           results/01_preprocessing/xenium_8um_filtered/{sample}_8um_filtered.rds
# Outputs:  results/03_benchmarking/probe_rank/ranked_plat.rds
#               (per-probe table: platform, feature, Type, mean_count, rank,
#                rank_frac, bg95, y, is_overlap; pooled: mean_count is the
#                mean over samples of per-sample totals)
#           results/03_benchmarking/probe_rank/ranked_sample.rds
#               (per-sample probe table: platform, feature, count, Sample,
#                Type, rank, rank_frac, bg95, y, is_overlap; count is the
#                per-sample total, bg95 the 95th percentile of that sample's
#                background probes, is_overlap = target with count <= bg95)
#           results/03_benchmarking/probe_rank/probe_overlap_targets.csv
#               (target probes falling within the bg95 threshold)
#           results/03_benchmarking/probe_rank/overlap_genes_merscope.txt
#           results/03_benchmarking/probe_rank/overlap_genes_xenium.txt
#           results/03_benchmarking/probe_rank/label_pool_genes.csv
#               (all MERSCOPE overlap probes, used for S-curve text
#                annotations in fig2_background.R)

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
              help    = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load 8um binning objects
# ---------------------------------------------------------------------------
# Samples without background probe signal are excluded via
# config["qc_backgrounds"]["exclude_samples"] (same exclusions as
# qc_backgrounds.R / moransi.R), so target and background means are computed
# over the same samples
bg_cfg           <- cfg$qc_backgrounds
merscope_samples <- setdiff(cfg$spatial_analysis$merscope_samples,      unlist(bg_cfg$exclude_samples$merscope))
xenium_samples   <- setdiff(cfg$spatial_analysis$xenium_default_samples, unlist(bg_cfg$exclude_samples$xenium))
walk(c("merscope", "xenium"), function(p) {
  excl <- unlist(bg_cfg$exclude_samples[[p]])
  if (length(excl) > 0) message("Excluding ", p, " samples (no background signal): ", paste(excl, collapse = ", "))
})

message("Loading MERSCOPE 8um objects...")
mer_dir      <- file.path(cfg$output_dir, "01_preprocessing", "merscope_8um_filtered")
merscope_8um <- setNames(
  lapply(merscope_samples, function(samp) {
    path <- file.path(mer_dir, paste0(samp, "_8um_filtered.rds"))
    message("  ", samp, ": ", path)
    readRDS(path)
  }),
  merscope_samples
)

message("Loading Xenium 8um objects...")
xen_dir    <- file.path(cfg$output_dir, "01_preprocessing", "xenium_8um_filtered")
xenium_8um <- setNames(
  lapply(xenium_samples, function(samp) {
    path <- file.path(xen_dir, paste0(samp, "_8um_filtered.rds"))
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

# Per-sample target/blank probe totals: platform, feature, count, Sample, Type
summarise_sample_probes <- function(sobj, sample_name, platform) {

  if (platform == "Xenium") {
    # Background = Unassigned codewords only ("BlankCodeword" assay), matching
    # qc_backgrounds.R (FDR) and moransi.R; ControlCodeword / ControlProbe are
    # deliberately excluded.
    blank_assays <- c("BlankCodeword")
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
    blank_mat <- get_counts_safe(sobj, "Blanks")
    if (is.null(viz_mat)) return(NULL)

    feats_v      <- rownames(viz_mat)
    # Detect blank rows by name when no dedicated Blanks assay is present
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
# Per-sample probe count tables
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

# Stop if any retained sample has no background probes (e.g. a sample with no
# background assay that has not been listed in exclude_samples)
no_bg <- all_counts %>%
  group_by(platform, Sample) %>%
  summarise(has_bg = any(Type == "Blank"), .groups = "drop") %>%
  filter(!has_bg)
if (nrow(no_bg) > 0) {
  stop("No background probes for: ", paste(no_bg$Sample, collapse = ", "),
       ". Add to qc_backgrounds$exclude_samples in config.yaml if intended.")
}

# ---------------------------------------------------------------------------
# Platform-level (pooled) ranked probe table
# ---------------------------------------------------------------------------
# mean_count = mean over samples of per-sample totals
plat_counts <- all_counts %>%
  group_by(platform, feature, Type) %>%
  summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop")

# 95th percentile of blank/background probes per platform — the overlap threshold
thr_platform <- plat_counts %>%
  filter(Type == "Blank") %>%
  group_by(platform) %>%
  summarise(bg95 = quantile(mean_count, 0.95, na.rm = TRUE), .groups = "drop")

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
# Per-sample ranked probe table
# ---------------------------------------------------------------------------
# Same definitions as the pooled table, applied within each sample (count =
# per-sample total, bg95 = that sample's threshold). The pooled table remains
# the source for main Fig 2.
thr_sample <- all_counts %>%
  filter(Type == "Blank") %>%
  group_by(platform, Sample) %>%
  summarise(bg95 = unname(quantile(count, 0.95, na.rm = TRUE)), .groups = "drop")

ranked_sample <- all_counts %>%
  group_by(platform, Sample) %>%
  arrange(desc(count), .by_group = TRUE) %>%
  mutate(
    rank      = row_number(),
    rank_frac = rank / max(rank, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  left_join(thr_sample, by = c("platform", "Sample")) %>%
  mutate(
    y          = log10(pmax(count, 0) + 1),
    is_overlap = (Type == "Target" & count <= bg95)
  )

message("Per-sample probe rank table: ", nrow(ranked_sample), " rows, ",
        n_distinct(ranked_sample$Sample), " samples")

# ---------------------------------------------------------------------------
# S-curve label pool
# ---------------------------------------------------------------------------
# Every MERSCOPE overlap target, applied to both platforms for consistent
# annotation across facets.
label_pool <- ranked_plat %>%
  filter(platform == "MERSCOPE", Type == "Target", is_overlap) %>%
  pull(feature)

label_genes <- ranked_plat %>%
  filter(feature %in% label_pool) %>%
  distinct(feature) %>%
  arrange(feature)

probe_overlap_targets <- ranked_plat %>%
  filter(Type == "Target", is_overlap) %>%
  arrange(platform, desc(mean_count)) %>%
  select(platform, feature, mean_count, rank, rank_frac, bg95)

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
saveRDS(ranked_plat, file.path(opt$out_dir, "ranked_plat.rds"))
message("Saved: ranked_plat.rds")

saveRDS(ranked_sample, file.path(opt$out_dir, "ranked_sample.rds"))
message("Saved: ranked_sample.rds")

utils::write.csv(probe_overlap_targets,
                 file.path(opt$out_dir, "probe_overlap_targets.csv"),
                 row.names = FALSE)
message("Saved: probe_overlap_targets.csv")

utils::write.csv(label_genes,
                 file.path(opt$out_dir, "label_pool_genes.csv"),
                 row.names = FALSE)
message("Saved: label_pool_genes.csv")

# One gene per line, no header. Written for every platform (empty if no
# overlaps) so no stale list from a previous run is left behind
walk(c("MERSCOPE", "Xenium"), function(plat) {
  genes    <- sort(unique(probe_overlap_targets$feature[probe_overlap_targets$platform == plat]))
  out_file <- file.path(opt$out_dir, paste0("overlap_genes_", tolower(plat), ".txt"))
  writeLines(genes, out_file)
  message("Saved: overlap_genes_", tolower(plat), ".txt (", length(genes), " genes)")
})

message("Done. Outputs written to: ", opt$out_dir)
