# Purpose:  Comparison of cell type composition across segmentation methods.
#           Assesses whether different segmentation approaches recover
#           consistent cell type proportions within and across samples,
#           using annotated_final.rds objects as input.
# Inputs:   Seurat objects (annotated_final.rds) per method from
#           results/02_spatial_analysis/{method}/
# Outputs:  results/03_benchmarking/cell_composition/
#           - composition_long.rds  — per-sample per-cell-type proportions
#           - composition_summary.rds — cross-method comparison statistics
# Author:   Ashleigh Solano
# Date:     2026-04-30

message("cell_composition.R — in development")
