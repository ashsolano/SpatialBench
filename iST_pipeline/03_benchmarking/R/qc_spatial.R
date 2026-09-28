# Purpose:  Per-bin coordinates with nCount / nFeature (all genes) for the spatial
#           QC maps, on the same bins as qc_metrics.R (checked against
#           metadata_combined.rds). Matched animals (config stalign) get aligned
#           coordinates in the matched VisiumHD pixel space (µm); others their
#           native bin grid. VisiumHD bins outside the main tissue section are
#           flagged (DBSCAN), not removed.
# Supports: 04_manuscript/R/extended/fig2ext_qc_spatial.R
# Inputs:   config/config.yaml, config/qc_spatial.yaml
#           results/01_preprocessing/{merscope,xenium}_8um_filtered/*.rds
#           results/01_preprocessing/{merscope,xenium}_8um_aligned/*.rds
#           VisiumHD objects in config visiumhd$data_dir
#           results/03_benchmarking/qc_metrics/metadata_combined.rds
# Outputs:  --out_rds (platform, Sample, bin_id, x_um, y_um, aligned,
#           main_section, nCount, nFeature)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(optparse)
})

source("03_benchmarking/R/utils/qc_utils.R")   # get_counts_mat(), build_meta_platform()
source("01_preprocessing/R/roi_utils.R")       # pixels_to_microns()

option_list <- list(
  make_option("--config",      type = "character", default = "config/config.yaml"),
  make_option("--qc_spatial",  type = "character", default = "config/qc_spatial.yaml"),
  make_option("--metadata",    type = "character",
              default = "results/03_benchmarking/qc_metrics/metadata_combined.rds"),
  make_option("--bin_size",    type = "integer",   default = 8L),
  make_option("--out_rds",     type = "character",
              default = "results/03_benchmarking/qc_spatial/bins_8um.rds")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
db  <- yaml::read_yaml(opt$qc_spatial)$qc_spatial$dbscan
res <- opt$bin_size
dir.create(dirname(opt$out_rds), recursive = TRUE, showWarnings = FALSE)

binned_path <- function(platform, stage, sample) {
  file.path(cfg$output_dir, "01_preprocessing", sprintf("%s_%dum_%s", platform, res, stage),
            sprintf("%s_%dum_%s.rds", sample, res, stage))
}

# Sample -> matched VisiumHD sample (NA if the animal has no alignment)
matched_visium <- function(platform, samples) {
  if (platform == "VisiumHD") return(setNames(samples, samples))
  by_animal <- cfg$stalign$matched_samples[[tolower(platform)]]
  vis       <- unlist(cfg$stalign$visiumhd[names(by_animal)])
  setNames(vis[match(samples, unlist(by_animal))], samples)
}

um_per_px <- function(vis_sample) as.numeric(cfg$stalign$microns_per_pixel[[vis_sample]])

# DBSCAN main section, as filter_noise() in fig2_qc.R but keeping every
# cluster of at least min_cluster_size bins
main_section <- function(x, y) {
  cl    <- dbscan::dbscan(cbind(x, y), eps = db$eps, minPts = db$minPts)$cluster
  sizes <- table(cl[cl > 0])
  cl %in% as.integer(names(sizes)[sizes >= db$min_cluster_size])
}

sample_bins <- function(obj, platform, sample, assay, vis_sample) {
  counts <- build_meta_platform(setNames(list(obj), sample), sample, assay, platform, "All",
                                paste0(res, "um"))
  counts$bin_id <- rownames(counts)

  if (platform == "VisiumHD") {
    xy <- pixels_to_microns(obj, um_per_px(sample), sprintf("slice1.%03dum", res))
  } else if (!is.na(vis_sample)) {
    aligned_obj <- readRDS(binned_path(tolower(platform), "aligned", sample))
    xy <- pixels_to_microns(aligned_obj, um_per_px(vis_sample))
    if (!setequal(xy$cell_id, counts$bin_id)) {
      stop(platform, " ", sample, ": aligned object does not hold the filtered bins")
    }
  } else {
    # Native grid: centroid = bin index x bin size (lower corner), so shift to the centre
    xy <- pixels_to_microns(obj, 1) %>% mutate(x_um = x_um + res / 2, y_um = y_um + res / 2)
  }

  out <- counts %>%
    select(platform, Sample, bin_id, nCount, nFeature) %>%
    left_join(select(xy, bin_id = cell_id, x_um, y_um), by = "bin_id") %>%
    mutate(aligned = !is.na(vis_sample))
  if (anyNA(out$x_um)) stop(platform, " ", sample, ": bins without coordinates")
  out$main_section <- if (platform == "VisiumHD") main_section(out$x_um, out$y_um) else TRUE
  message(sprintf("  %s %s: %d bins, %d outside the main section", platform, sample,
                  nrow(out), sum(!out$main_section)))
  out
}

platforms <- list(
  VisiumHD = list(assay = sprintf("Spatial.%03dum", res),
                  paths = setNames(file.path(cfg$visiumhd$data_dir, unlist(cfg$visiumhd$samples)),
                                   names(cfg$visiumhd$samples))),
  MERSCOPE = list(assay = "Vizgen", samples = cfg$spatial_analysis$merscope_samples),
  Xenium   = list(assay = "Xenium", samples = cfg$spatial_analysis$xenium_default_samples)
)
for (plat in c("MERSCOPE", "Xenium")) {
  s <- platforms[[plat]]$samples
  platforms[[plat]]$paths <- setNames(binned_path(tolower(plat), "filtered", s), s)
}
if (anyNA(matched_visium("VisiumHD", names(platforms$VisiumHD$paths)))) stop("VisiumHD sample without stalign entry")

bins <- imap_dfr(platforms, function(p, plat) {
  message("Extracting ", plat, "...")
  vis <- matched_visium(plat, names(p$paths))
  imap_dfr(p$paths, function(path, sample) {
    obj <- readRDS(path)
    out <- sample_bins(obj, plat, sample, p$assay, vis[[sample]])
    rm(obj); gc(verbose = FALSE)
    out
  })
})

# Same bins and totals as the QC violins
per_sample <- function(df) {
  df %>%
    mutate(platform = as.character(platform)) %>%
    group_by(platform, Sample) %>%
    summarise(n = n(), nCount = sum(as.numeric(nCount)), nFeature = sum(as.numeric(nFeature)),
              .groups = "drop")
}
ref <- readRDS(opt$metadata) %>% filter(Subset == "All", bin_size == paste0(res, "um"))
check <- full_join(per_sample(bins), per_sample(ref), by = c("platform", "Sample"), suffix = c("", "_ref"))
bad   <- filter(check, is.na(n) | is.na(n_ref) | n != n_ref | nCount != nCount_ref | nFeature != nFeature_ref)
if (nrow(bad) > 0) {
  print(as.data.frame(bad))
  stop("Bins do not match metadata_combined.rds for ", nrow(bad), " samples")
}
message("Matches metadata_combined.rds: ", nrow(check), " samples, ", nrow(bins), " bins")

bins <- mutate(bins, platform = factor(platform, levels = c("VisiumHD", "MERSCOPE", "Xenium")))
saveRDS(bins, opt$out_rds)
message("Saved: ", opt$out_rds)
