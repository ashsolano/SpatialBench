# Purpose:  Concordance analysis of differential expression results across
#           segmentation methods. Compares DE gene lists, log-fold changes,
#           and significance rankings between methods to assess whether
#           biological conclusions are robust to segmentation choice.
# Inputs:   DE results (de_results.rds) per method from
#           results/02_spatial_analysis/{method}/pseudobulk_de/
# Outputs:  results/03_benchmarking/de_concordance/
#           - de_concordance.rds    — pairwise gene-level comparison table
#           - de_overlap_summary.rds — Jaccard / correlation summary

message("de_concordance.R — in development")
