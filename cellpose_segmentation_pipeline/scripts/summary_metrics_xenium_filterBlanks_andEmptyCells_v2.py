# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 11:53:29 2025

@author: zaman.i
"""


import os
import pandas as pd

cp_parent_dir = "path/to/resegmentation_cellpose_xenium"
xen_parent_dir = "path/to/xenium_default_segmentation"
results_dir = "path/to/results"
segmentations = ["Cellpose", "Xenium"]

cp_samples = [f for f in os.listdir(cp_parent_dir) if "Batch" not in f and "final" not in f] 
xen_samples = [f for f in os.listdir(xen_parent_dir) if "output" in f]
# Prepare lists to collect results
cp_results_filtered, cp_results_unfiltered = [], []
xen_results_filtered, xen_results_unfiltered = [], []

for segmentation in segmentations:
    if segmentation == "Cellpose":
        parent_folder = cp_parent_dir
        samples = cp_samples
    else: 
        parent_folder = xen_parent_dir
        samples = xen_samples
    
    for sample in samples: 
        folder = os.path.join(parent_folder, sample)
        sample_name = "sample_name"
        print(f"Processing {sample} - {segmentation}")
        
        # Load files
        if segmentation == "Xenium":
            cell_metadata = pd.read_parquet(os.path.join(folder, "cells.parquet"))
            detected_transcripts = pd.read_parquet(os.path.join(folder, "transcripts.parquet"))
        else:  # Cellpose
            cell_metadata = pd.read_parquet(os.path.join(folder, "outs/cells.parquet"))
            detected_transcripts = pd.read_parquet(os.path.join(folder, "outs/transcripts.parquet"))
    
        # --- Unfiltered metrics ---
        
        num_non_genes = (~detected_transcripts["codeword_category"].isin(["predesigned_gene", "custom_gene"])).sum()
        
        total_cells = cell_metadata.shape[0]
        total_transcripts = detected_transcripts.shape[0]

        transcripts_per_cell_mean = cell_metadata["total_counts"].mean()
        transcripts_per_cell_median = cell_metadata["total_counts"].median()

        cell_area_mean = cell_metadata["cell_area"].mean()
        cell_area_median = cell_metadata["cell_area"].median()
    
        unfiltered_dict = {
            "Sample": sample_name,
            f"Cell count - {segmentation}": total_cells,
            f"Transcripts per cell - mean - {segmentation}": transcripts_per_cell_mean,
            f"Transcripts per cell - median - {segmentation}": transcripts_per_cell_median,
            f"Transcripts within a cell (%) - {segmentation}": cell_metadata["total_counts"].sum() / total_transcripts,
            f"Cell Area - mean (um^2) - {segmentation}": cell_area_mean,
            f"Cell Area - median (um^2) - {segmentation}": cell_area_median,
            f"Total Transcript Count - {segmentation}": total_transcripts,
            f"Number of Non-Genes - {segmentation}": num_non_genes
            
        }
    
        if segmentation == "Cellpose":
            cp_results_unfiltered.append(unfiltered_dict)
        else:
            xen_results_unfiltered.append(unfiltered_dict)
    
        # --- Filter out blank genes and empty cells ---
        
        detected_transcripts_filtered = detected_transcripts[
            (detected_transcripts["is_gene"] == True)
            & (detected_transcripts["codeword_category"].isin(["predesigned_gene", "custom_gene"]))
            & (detected_transcripts["qv"] >= 20)
        ]
    
        empty_cells_idx = cell_metadata.index[cell_metadata["transcript_counts"] <1]
        print(f"Filtering {len(empty_cells_idx)} empty cells in {folder}")
    
        cell_metadata_filtered = cell_metadata.drop(index=empty_cells_idx)
    
        # ---- Filtered metrics ----
        total_cells_filtered = cell_metadata_filtered.shape[0]
        total_transcripts_filtered = detected_transcripts_filtered.shape[0]
    
        transcripts_per_cell_mean_filtered = cell_metadata_filtered["transcript_counts"].mean()
        transcripts_per_cell_median_filtered = cell_metadata_filtered["transcript_counts"].median()
    
        cell_area_mean_filtered = cell_metadata_filtered["cell_area"].mean()
        cell_area_median_filtered = cell_metadata_filtered["cell_area"].median()
    
        filtered_dict = {
            "Sample": sample,
            f"Cell count - {segmentation}": total_cells_filtered,
            f"Transcripts per cell - mean - {segmentation}": transcripts_per_cell_mean_filtered,
            f"Transcripts per cell - median - {segmentation}": transcripts_per_cell_median_filtered,
            f"Transcripts within a cell (%) - {segmentation}": cell_metadata_filtered["transcript_counts"].sum()/detected_transcripts_filtered.shape[0],
            f"Cell Area - mean (um^2) - {segmentation}": cell_area_mean_filtered,
            f"Cell Area - median (um^2) - {segmentation}": cell_area_median_filtered,
            f"Total Transcript Count - {segmentation}": total_transcripts_filtered,
            f"Number of Empty Cells Filtered Out - {segmentation}": len(empty_cells_idx)
        }
    
        if segmentation == "Cellpose":
            cp_results_filtered.append(filtered_dict)
        else:
            xen_results_filtered.append(filtered_dict)

# --- Save results ---
cp_summary_unfiltered = pd.DataFrame(cp_results_unfiltered)
cp_summary_filtered = pd.DataFrame(cp_results_filtered)
xen_summary_unfiltered = pd.DataFrame(xen_results_unfiltered)
xen_summary_filtered = pd.DataFrame(xen_results_filtered)

cp_summary_unfiltered.to_csv(os.path.join(results_dir, "xenium_summary_metrics_unfiltered_cp.csv"), index=False)
cp_summary_filtered.to_csv(os.path.join(results_dir, "xenium_summary_metrics_filtered_cp.csv"), index=False)

xen_summary_unfiltered.to_csv(os.path.join(results_dir, "xenium_summary_metrics_unfiltered_xen.csv"), index=False)
xen_summary_filtered.to_csv(os.path.join(results_dir, "xenium_summary_metrics_filtered_xen.csv"), index=False)

print("✅ Saved filtered and unfiltered summary metrics for both segmentations")