# Purpose:  Benchmarking of germinal centre (GC) zone recovery across
#           segmentation methods. Evaluates whether light zone / dark zone
#           subclusters are consistently identified, and quantifies the
#           number of GC B cells and zone-assigned cells per method.
# Inputs:   Seurat objects (annotated_final.rds and subclustered_gc.rds)
#           per method from results/02_spatial_analysis/{method}/
# Outputs:  results/03_benchmarking/gc_zone_recovery/
#           - gc_zone_metrics.rds   — per-method GC zone recovery statistics
#           - gc_zone_summary.rds   — cross-method comparison table

message("gc_zone_recovery.R — in development")
