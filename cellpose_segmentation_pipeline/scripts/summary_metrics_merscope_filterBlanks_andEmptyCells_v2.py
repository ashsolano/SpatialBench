# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 10:43:17 2025

@author: zaman.i
"""

import os
import pandas as pd

parent_dir = "path/to/folder/with/segmentation_outputs"
segmentations = ["Cellpose", "Vizgen"]

samples = [f for f in os.listdir(parent_dir) if "summary" not in f]

# Prepare lists to collect results
cp_results_filtered, cp_results_unfiltered = [], []
viz_results_filtered, viz_results_unfiltered = [], []

for sample in samples:
    for segmentation in segmentations:
        folder = os.path.join(parent_dir, sample, segmentation)
        print(f"Processing {sample} - {segmentation}")

        # Load files
        if segmentation == "Vizgen":
            cell_by_gene = pd.read_csv(os.path.join(folder, "cell_by_gene.csv"))
            cell_metadata = pd.read_csv(os.path.join(folder, "cell_metadata.csv"))
            detected_transcripts = pd.read_csv(os.path.join(folder, "detected_transcripts.csv"))
        else:  # Cellpose
            cell_by_gene = pd.read_csv(os.path.join(folder, "cell_by_gene_repartitioned.csv"))
            cell_metadata = pd.read_csv(os.path.join(folder, "cell_metadata_nCounts.csv"))
            detected_transcripts = pd.read_csv(os.path.join(folder, "detected_transcripts.csv"))

        # --- Unfiltered metrics ---
        num_blank_genes = cell_by_gene.filter(like="Blank", axis=1).shape[1]
        num_blank_transcripts = detected_transcripts["gene"].str.lower().str.contains("blank").sum()
        total_cells = cell_by_gene.shape[0]
        total_transcripts = detected_transcripts.shape[0]

        transcripts_per_cell = cell_by_gene.drop(columns=["cell"]).sum(axis=1)
        transcripts_per_cell_mean = transcripts_per_cell.mean()
        transcripts_per_cell_median = transcripts_per_cell.median()

        unique_genes_per_cell = (cell_by_gene.drop(columns=["cell"]) > 0).sum(axis=1)
        unique_genes_per_cell_mean = unique_genes_per_cell.mean()
        unique_genes_per_cell_median = unique_genes_per_cell.median()

        cell_volume_mean = cell_metadata["volume"].mean()
        cell_volume_median = cell_metadata["volume"].median()

        unfiltered_dict = {
            "Sample": sample,
            f"Cell count - {segmentation}": total_cells,
            f"Transcripts per cell - mean - {segmentation}": transcripts_per_cell_mean,
            f"Transcripts per cell - median - {segmentation}": transcripts_per_cell_median,
            f"Transcripts within a cell (%) - {segmentation}": transcripts_per_cell.sum() / total_transcripts,
            f"Unique genes per cell - mean - {segmentation}": unique_genes_per_cell_mean,
            f"Unique genes per cell - median - {segmentation}": unique_genes_per_cell_median,
            f"Cell volume - mean (um^3) - {segmentation}": cell_volume_mean,
            f"Cell volume - median (um^3) - {segmentation}": cell_volume_median,
            f"Cell Area - mean (um^2) - {segmentation}": cell_volume_mean / (7*1.5),
            f"Cell Area - median (um^2) - {segmentation}": cell_volume_median / (7*1.5),
            f"Total Transcript Count - {segmentation}": total_transcripts,
            f"Number of Blank Genes - {segmentation}": num_blank_genes,
            f"Number of Blank Transcripts - {segmentation}": num_blank_transcripts
        }

        if segmentation == "Cellpose":
            cp_results_unfiltered.append(unfiltered_dict)
        else:
            viz_results_unfiltered.append(unfiltered_dict)

        # --- Filter out blank genes and empty cells ---
        cell_by_gene_filtered = cell_by_gene.loc[:, ~cell_by_gene.columns.str.contains("Blank", case=False)]
        detected_transcripts_filtered = detected_transcripts[~detected_transcripts["gene"].str.contains("blank", case=False)]

        empty_cells_idx = cell_by_gene_filtered.index[(cell_by_gene_filtered.drop(columns=["cell"]) == 0).all(axis=1)]
        print(f"Filtering {len(empty_cells_idx)} empty cells in {folder}")

        cell_by_gene_filtered = cell_by_gene_filtered.drop(index=empty_cells_idx)
        cell_metadata_filtered = cell_metadata.drop(index=empty_cells_idx)

        # ---- Filtered metrics ----
        total_cells_filtered = cell_by_gene_filtered.shape[0]
        total_transcripts_filtered = detected_transcripts_filtered.shape[0]

        transcripts_per_cell_filtered = cell_by_gene_filtered.drop(columns=["cell"]).sum(axis=1)
        transcripts_per_cell_mean_filtered = transcripts_per_cell_filtered.mean()
        transcripts_per_cell_median_filtered = transcripts_per_cell_filtered.median()

        unique_genes_per_cell_filtered = (cell_by_gene_filtered.drop(columns=["cell"]) > 0).sum(axis=1)
        unique_genes_per_cell_mean_filtered = unique_genes_per_cell_filtered.mean()
        unique_genes_per_cell_median_filtered = unique_genes_per_cell_filtered.median()

        cell_volume_mean_filtered = cell_metadata_filtered["volume"].mean()
        cell_volume_median_filtered = cell_metadata_filtered["volume"].median()

        filtered_dict = {
            "Sample": sample,
            f"Cell count - {segmentation}": total_cells_filtered,
            f"Transcripts per cell - mean - {segmentation}": transcripts_per_cell_mean_filtered,
            f"Transcripts per cell - median - {segmentation}": transcripts_per_cell_median_filtered,
            f"Transcripts within a cell (%) - {segmentation}": transcripts_per_cell_filtered.sum() / total_transcripts_filtered,
            f"Unique genes per cell - mean - {segmentation}": unique_genes_per_cell_mean_filtered,
            f"Unique genes per cell - median - {segmentation}": unique_genes_per_cell_median_filtered,
            f"Cell volume - mean (um^3) - {segmentation}": cell_volume_mean_filtered,
            f"Cell volume - median (um^3) - {segmentation}": cell_volume_median_filtered,
            f"Cell Area - mean (um^2) - {segmentation}": cell_volume_mean_filtered / (7*1.5),
            f"Cell Area - median (um^2) - {segmentation}": cell_volume_median_filtered / (7*1.5),
            f"Total Transcript Count - {segmentation}": total_transcripts_filtered,
            f"Number of Empty Cells Filtered Out - {segmentation}": len(empty_cells_idx)
        }

        if segmentation == "Cellpose":
            cp_results_filtered.append(filtered_dict)
        else:
            viz_results_filtered.append(filtered_dict)

# --- Save results ---
cp_summary_unfiltered = pd.DataFrame(cp_results_unfiltered)
cp_summary_filtered = pd.DataFrame(cp_results_filtered)
viz_summary_unfiltered = pd.DataFrame(viz_results_unfiltered)
viz_summary_filtered = pd.DataFrame(viz_results_filtered)

cp_summary_unfiltered.to_csv(os.path.join(parent_dir, "final_summary_metrics/merscope_summary_metrics_unfiltered_cp.csv"), index=False)
cp_summary_filtered.to_csv(os.path.join(parent_dir, "final_summary_metrics/merscope_summary_metrics_noBlanks_noEmpty_cp.csv"), index=False)

viz_summary_unfiltered.to_csv(os.path.join(parent_dir, "final_summary_metrics/merscope_summary_metrics_unfiltered_viz.csv"), index=False)
viz_summary_filtered.to_csv(os.path.join(parent_dir, "final_summary_metrics/merscope_summary_metrics_noBlanks_noEmpty_viz.csv"), index=False)

print("✅ Saved filtered and unfiltered summary metrics for both segmentations")
