# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 13:22:00 2025

@author: zaman.i
"""

# STEP 1: Calculate areas of layer 2 and layer 3 (for MERSCOPE and Xenium)
# STEP 2: Pull out GeoJSON of layer 2 and 3 for visualisation (for MERSCOPE and Xenium)
# STEP 3: Use expected counts and detected transcripts file or transcripts.csv.gz to get % transcripts and transcripts p/cell 

import pandas as pd
import os
import gzip
import json
import geopandas as gpd
from shapely.geometry import shape


# %% MERSCOPE Transcripts Metrics + Layer 2 Areas

parent_dir = '/path/to/merscope_proseg_segmentation'
results_dir = "/path/to/results"
folders = [f for f in os.listdir(parent_dir) if "proseg_merscope" in f] 



results = []
proseg_results_filtered, proseg_results_unfiltered = [], []
for folder in folders:
    print(f"Processing: {folder}")
    
    sample_name = folder.split("merscope_")[-1]
    segmentation = "Proseg"
    
    detected_transcripts = pd.read_csv(os.path.join(parent_dir, folder, "output/detected_transcripts.csv"))
    expected_counts = pd.read_csv(os.path.join(parent_dir, folder, "output/expected-counts.csv.gz"), compression='gzip')
    cell_metadata = pd.read_csv(os.path.join(parent_dir, folder, "output/cell-metadata.csv.gz"), compression='gzip')
    cell_by_gene = expected_counts
    
    
    num_blank_genes = cell_by_gene.filter(like="Blank", axis=1).shape[1]
    num_blank_transcripts = detected_transcripts["gene"].str.lower().str.contains("blank").sum()
    total_cells = cell_by_gene.shape[0]
    total_transcripts = detected_transcripts.shape[0]

    transcripts_per_cell = cell_by_gene.sum(axis=1)
    transcripts_per_cell_mean = transcripts_per_cell.mean()
    transcripts_per_cell_median = transcripts_per_cell.median()

    unique_genes_per_cell = (cell_by_gene > 0).sum(axis=1)
    unique_genes_per_cell_mean = unique_genes_per_cell.mean()
    unique_genes_per_cell_median = unique_genes_per_cell.median()

    geojson_path = os.path.join(parent_dir, folder, "output", "cell-polygons-layers.geojson.gz")
    mean_area = None
    median_area = None

    with gzip.open(geojson_path, "rt", encoding="utf-8") as f:
        data = json.load(f)
    records = []
    for feat in data["features"]:
        props = feat["properties"].copy()
        props["geometry"] = shape(feat["geometry"])
        records.append(props)
    gdf = gpd.GeoDataFrame(records, geometry="geometry")

    # filter layer 2
    gdf_layer2 = gdf[gdf["layer"] == 2].copy()
    # compute area
    gdf_layer2["area_um2"] = gdf_layer2.geometry.area
    mean_area = gdf_layer2["area_um2"].mean()
    median_area= gdf_layer2["area_um2"].median()

    unfiltered_dict = {
        "Sample": sample_name,
        f"Cell count - {segmentation}": total_cells,
        f"Transcripts per cell - mean - {segmentation}": transcripts_per_cell_mean,
        f"Transcripts per cell - median - {segmentation}": transcripts_per_cell_median,
        f"Transcripts within a cell (%) - {segmentation}": transcripts_per_cell.sum() / total_transcripts,
        f"Unique genes per cell - mean - {segmentation}": unique_genes_per_cell_mean,
        f"Unique genes per cell - median - {segmentation}": unique_genes_per_cell_median,
        f"Cell Area - mean (um^2) - {segmentation}": mean_area,
        f"Cell Area - median (um^2) - {segmentation}": median_area,
        f"Total Transcript Count - {segmentation}": total_transcripts,
        f"Number of Blank Genes - {segmentation}": num_blank_genes,
        f"Number of Blank Transcripts - {segmentation}": num_blank_transcripts
    }
    
    proseg_results_unfiltered.append(unfiltered_dict)
    
    # --- Filter out blank genes and empty cells ---
    cell_by_gene_filtered = cell_by_gene.loc[:, ~cell_by_gene.columns.str.contains("Blank", case=False)]
    detected_transcripts_filtered = detected_transcripts[~detected_transcripts["gene"].str.contains("blank", case=False)]

    empty_cells_idx = cell_by_gene_filtered.index[(cell_by_gene_filtered == 0).all(axis=1)]
    print(f"Filtering {len(empty_cells_idx)} empty cells in {folder}")

    cell_by_gene_filtered = cell_by_gene_filtered.drop(index=empty_cells_idx)
    cell_metadata_filtered = cell_metadata.drop(index=empty_cells_idx)

    # ---- Filtered metrics ----
    total_cells_filtered = cell_by_gene_filtered.shape[0]
    total_transcripts_filtered = detected_transcripts_filtered.shape[0]

    transcripts_per_cell_filtered = cell_by_gene_filtered.sum(axis=1)
    transcripts_per_cell_mean_filtered = transcripts_per_cell_filtered.mean()
    transcripts_per_cell_median_filtered = transcripts_per_cell_filtered.median()

    unique_genes_per_cell_filtered = (cell_by_gene_filtered > 0).sum(axis=1)
    unique_genes_per_cell_mean_filtered = unique_genes_per_cell_filtered.mean()
    unique_genes_per_cell_median_filtered = unique_genes_per_cell_filtered.median()

    # ---- Layer 2 Area (µm²), excluding empty cells ----

    # drop empty cells by matching 'cell' to expected_counts index
    gdf_layer2_filtered = gdf_layer2[~gdf_layer2["cell"].isin(empty_cells_idx)]

    # compute area
    gdf_layer2_filtered["area_um2"] = gdf_layer2_filtered.geometry.area
    mean_area_filtered = gdf_layer2_filtered["area_um2"].mean()
    median_area_filtered = gdf_layer2_filtered["area_um2"].median()
    print(f"✅ Computed layer 2 area stats for {folder} (excluded {len(empty_cells_idx)} empty cells)")


    # Collect metrics
    filtered_dict = {
        "Sample": sample_name,
        f"Cell count - {segmentation}": total_cells_filtered,
        f"Transcripts per cell - mean - {segmentation}": transcripts_per_cell_mean_filtered,
        f"Transcripts per cell - median - {segmentation}": transcripts_per_cell_median_filtered,
        f"Transcripts within a cell (%) - {segmentation}": transcripts_per_cell_filtered.sum() / total_transcripts_filtered,
        f"Unique genes per cell - mean - {segmentation}": unique_genes_per_cell_mean_filtered,
        f"Unique genes per cell - median - {segmentation}": unique_genes_per_cell_median_filtered,
        f"Cell Area - mean (um^2) - {segmentation}": mean_area_filtered,
        f"Cell Area - median (um^2) - {segmentation}": median_area_filtered,
        f"Total Transcript Count - {segmentation}": total_transcripts_filtered,
        f"Number of Empty Cells Filtered Out - {segmentation}": len(empty_cells_idx)
    }


    proseg_results_filtered.append(filtered_dict)
    

# Save results summary
proseg_summary_unfiltered = pd.DataFrame(proseg_results_unfiltered)
proseg_summary_filtered = pd.DataFrame(proseg_results_filtered)
proseg_summary_unfiltered.to_csv(os.path.join(results_dir, "merscope_summary_metrics_unfiltered_proseg_updated_samples.csv"), index=False)
proseg_summary_filtered.to_csv(os.path.join(results_dir, "merscope_summary_metrics_filtered_proseg_updated_samples.csv"), index=False)

print(f"✅ Summary metrics saved")




# %% Xenium Transcripts Metrics + Layer 2 Areas

parent_dir = '/path/to/xenium_proseg_segmentation'
results_dir = '/path/to/results'
folders = [f for f in os.listdir(parent_dir) if "summary" not in f] 



results = []
proseg_results_filtered, proseg_results_unfiltered = [], []
for folder in folders:
    print(f"Processing: {folder}")
    
    sample_name = folder.split("xenium__")[-1]
    segmentation = "Proseg"
    
    detected_transcripts = pd.read_parquet(os.path.join(parent_dir, folder, "input/transcripts.parquet"))
    expected_counts = pd.read_csv(os.path.join(parent_dir, folder, "output/expected-counts.csv.gz"), compression='gzip')
    cell_metadata = pd.read_csv(os.path.join(parent_dir, folder, "output/cell-metadata.csv.gz"), compression='gzip')
    cell_by_gene = expected_counts
    
    
        
    num_non_genes = (~detected_transcripts["codeword_category"].isin(["predesigned_gene", "custom_gene"])).sum()
    
    total_cells = cell_metadata.shape[0]
    total_transcripts = detected_transcripts.shape[0]
    
    transcripts_per_cell = expected_counts.sum(axis=1)
    transcripts_per_cell_mean = transcripts_per_cell.mean()
    transcripts_per_cell_median = transcripts_per_cell.median()

    transcripts_total = transcripts_per_cell.sum()
    transcripts_within_cell_pct = transcripts_total / detected_transcripts.shape[0]

    
    geojson_path = os.path.join(parent_dir, folder, "output", "cell-polygons-layers.geojson.gz")
    mean_area = None
    median_area = None

    with gzip.open(geojson_path, "rt", encoding="utf-8") as f:
        data = json.load(f)
    records = []
    for feat in data["features"]:
        props = feat["properties"].copy()
        props["geometry"] = shape(feat["geometry"])
        records.append(props)
    gdf = gpd.GeoDataFrame(records, geometry="geometry")

    # filter layer 2
    gdf_layer2 = gdf[gdf["layer"] == 2].copy()
    # compute area
    gdf_layer2["area_um2"] = gdf_layer2.geometry.area
    mean_area = gdf_layer2["area_um2"].mean()
    median_area= gdf_layer2["area_um2"].median()
    
    unfiltered_dict = {
        "Sample": sample_name,
        f"Cell count - {segmentation}": total_cells,
        f"Transcripts per cell - mean - {segmentation}": transcripts_per_cell_mean,
        f"Transcripts per cell - median - {segmentation}": transcripts_per_cell_median,
        f"Transcripts within a cell (%) - {segmentation}":  transcripts_within_cell_pct,
        f"Cell Area - mean (um^2) - {segmentation}": mean_area,
        f"Cell Area - median (um^2) - {segmentation}": median_area,
        f"Total Transcript Count - {segmentation}": total_transcripts,
        f"Number of Non-Genes - {segmentation}": num_non_genes
        
    }
    

    proseg_results_unfiltered.append(unfiltered_dict)

    
    # --- Filter out blank genes and empty cells ---
    detected_transcripts_filtered = detected_transcripts[
        (detected_transcripts["is_gene"] == True)
        & (detected_transcripts["codeword_category"].isin(["predesigned_gene", "custom_gene"]))
        & (detected_transcripts["qv"] >= 20)
    ]
    
    empty_cells_idx = cell_by_gene.index[(expected_counts == 0).all(axis=1)]
    print(f"Filtering {len(empty_cells_idx)} empty cells in {folder}")
    
    #cell_by_gene_filtered = cell_metadata.drop(index=empty_cells_idx)
    expected_counts_filtered = expected_counts.drop(index=empty_cells_idx)
    cell_metadata_filtered = cell_metadata.drop(index=empty_cells_idx)
    # ---- Filtered metrics ----
    total_cells_filtered = cell_metadata_filtered.shape[0]
    total_transcripts_filtered = detected_transcripts_filtered.shape[0]
    
    transcripts_per_cell = expected_counts_filtered.sum(axis=1)
    transcripts_per_cell_mean_filtered = transcripts_per_cell.mean()
    transcripts_per_cell_median_filtered = transcripts_per_cell.median()

    transcripts_total = transcripts_per_cell.sum()
    transcripts_within_cell_pct = transcripts_total / detected_transcripts_filtered.shape[0]
    gdf_layer2_filtered = gdf_layer2[~gdf_layer2["cell"].isin(empty_cells_idx)]

    # compute area
    gdf_layer2_filtered["area_um2"] = gdf_layer2_filtered.geometry.area
    mean_area_filtered = gdf_layer2_filtered["area_um2"].mean()
    median_area_filtered = gdf_layer2_filtered["area_um2"].median()
    
    filtered_dict = {
        "Sample": sample_name,
        f"Cell count - {segmentation}": total_cells_filtered,
        f"Transcripts per cell - mean - {segmentation}": transcripts_per_cell_mean_filtered,
        f"Transcripts per cell - median - {segmentation}": transcripts_per_cell_median_filtered,
        f"Transcripts within a cell (%) - {segmentation}":  transcripts_within_cell_pct,
        f"Cell Area - mean (um^2) - {segmentation}": mean_area_filtered,
        f"Cell Area - median (um^2) - {segmentation}": median_area_filtered,
        f"Total Transcript Count - {segmentation}": total_transcripts_filtered,
        f"Number of Empty Cells Filtered Out - {segmentation}": len(empty_cells_idx)
    }

    proseg_results_filtered.append(filtered_dict)

# Save results summary
proseg_summary_unfiltered = pd.DataFrame(proseg_results_unfiltered)
proseg_summary_filtered = pd.DataFrame(proseg_results_filtered)
proseg_summary_unfiltered.to_csv(os.path.join(results_dir, "xenium_summary_metrics_unfiltered_proseg.csv"), index=False)
proseg_summary_filtered.to_csv(os.path.join(results_dir, "xenium_summary_metrics_filtered_proseg.csv"), index=False)

print(f"✅ Summary metrics saved")


