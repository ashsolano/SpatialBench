# Purpose:  Post-processing QC for a binned Seurat object — drop empty bins
#           and remove spatially isolated bins (e.g. debris, edge artefacts)
#           via DBSCAN, keeping only the dominant tissue cluster.
# Inputs:   Binned Seurat object RDS (output of create_seurat_binned_xenium.R
#           or create_seurat_binned_merscope.R)
# Outputs:  Filtered Seurat object saved to --out_rds, plus a before/after
#           spatial QC plot (PNG) saved alongside it as
#           <sample>_<resolution>um_filter_qc.png
# Usage:    Rscript 01_preprocessing/R/filter_binned.R \
#             --input_rds <path> --out_rds <path> --assay <Vizgen|Xenium> \
#             [--eps <15>] [--minPts <5>] [--min_count <1>] \
#             [--min_cluster_size <1000>] \
#             [--use_adaptive_threshold | --no_adaptive_threshold]

library(optparse)
library(Seurat)
library(ggplot2)
library(patchwork)

source("utils/binning_utils.R")  # must be run from the project root

option_list <- list(
  make_option(c("--input_rds"), type = "character", default = NULL,
              help = "Path to input binned Seurat object RDS"),
  make_option(c("--out_rds"),   type = "character", default = NULL,
              help = "Output path for the filtered Seurat object RDS"),
  make_option(c("--assay"),     type = "character", default = "Vizgen",
              help = "Assay name to QC (nCount column) [default: %default]"),
  make_option(c("--eps"),       type = "double",    default = 15,
              help = "DBSCAN neighbourhood radius, in microns [default: %default]"),
  make_option(c("--minPts"),    type = "integer",   default = 5,
              help = "DBSCAN minimum points per cluster [default: %default]"),
  make_option(c("--min_count"), type = "integer",   default = 1,
              help = "Minimum transcript count for a bin to be kept; used only when adaptive thresholding is disabled [default: %default]"),
  make_option(c("--min_cluster_size"), type = "integer", default = 1000,
              help = "Minimum DBSCAN cluster size (bins) to be kept as tissue [default: %default]"),
  make_option(c("--use_adaptive_threshold"), action = "store_true", default = TRUE,
              dest = "use_adaptive_threshold",
              help = "Detect the background/tissue nCount threshold automatically [default]"),
  make_option(c("--no_adaptive_threshold"), action = "store_false", default = TRUE,
              dest = "use_adaptive_threshold",
              help = "Disable adaptive thresholding; use --min_count instead")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds)) stop("--input_rds is required")
if (is.null(opt$out_rds))   stop("--out_rds is required")

out_dir <- dirname(opt$out_rds)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Input RDS:              ", opt$input_rds)
message("Assay:                  ", opt$assay)
message("eps:                    ", opt$eps)
message("minPts:                 ", opt$minPts)
message("min_count:              ", opt$min_count)
message("min_cluster_size:       ", opt$min_cluster_size)
message("use_adaptive_threshold: ", opt$use_adaptive_threshold)

obj <- readRDS(opt$input_rds)

obj_filt <- filter_binned_seurat(
  obj,
  assay                  = opt$assay,
  min_count              = opt$min_count,
  use_adaptive_threshold = opt$use_adaptive_threshold,
  eps                    = opt$eps,
  minPts                 = opt$minPts,
  min_cluster_size       = opt$min_cluster_size
)

saveRDS(obj_filt, file = opt$out_rds)
message("Saved: ", opt$out_rds)

# --- Before/after spatial QC plot ---------------------------------------------

# Bin centroids before filtering, labelled by whether each bin survived
fov_name  <- names(obj@images)[1]
before_df <- as.data.frame(obj@images[[fov_name]]@boundaries$centroids@coords)
colnames(before_df) <- c("x", "y")
before_df$status <- ifelse(colnames(obj) %in% colnames(obj_filt), "kept", "removed")

fov_name_filt <- names(obj_filt@images)[1]
after_df      <- as.data.frame(obj_filt@images[[fov_name_filt]]@boundaries$centroids@coords)
colnames(after_df) <- c("x", "y")

sample_label <- tools::file_path_sans_ext(basename(opt$input_rds))

p_before <- ggplot(before_df, aes(x = x, y = y, color = status)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_manual(values = c(kept = "grey30", removed = "firebrick")) +
  coord_fixed() +
  labs(title = "Before filtering", subtitle = paste0(nrow(before_df), " bins")) +
  theme_minimal() +
  theme(legend.position = "bottom")

p_after <- ggplot(after_df, aes(x = x, y = y)) +
  geom_point(size = 0.3, alpha = 0.6, color = "grey30") +
  coord_fixed() +
  labs(title = "After filtering", subtitle = paste0(nrow(after_df), " bins")) +
  theme_minimal()

combined_plot <- (p_before + p_after) +
  plot_annotation(title = sample_label)

qc_plot_file <- file.path(out_dir, paste0(sample_label, "_filter_qc.png"))
ggsave(qc_plot_file, combined_plot, width = 10, height = 5, dpi = 150)
message("Saved: ", qc_plot_file)
