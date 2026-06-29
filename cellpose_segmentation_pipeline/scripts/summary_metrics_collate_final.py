# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 10:55:36 2025

@author: zaman.i
"""

import pandas as pd
import os

# %% MERSCOPE Samples 

parent_dir = "/path/to/summary_metrics_merscope"

filtered_samples = [f for f in os.listdir(parent_dir) if "_filtered" in f and "merged" not in f]
unfiltered_samples = [f for f in os.listdir(parent_dir) if "unfiltered" in f and "merged" not in f]

mapping_samples= {
    "G000218ZalcensteinBatch27_VMSC04501_region_1": "Batch27_region1_Ctrl172",
    "G000218BenchmarkingBatch26_VMSC04501_region_0": "Batch26_region0_Ctrl173",
    "G000218Batch22FFMuSpleen_VMSC04501_region_1": "Batch22_region1_Ctrl174",
    "G000218Batch22FFMuSpleen_VMSC04501_region_2": "Batch22_region2_KO166",
    "G000218ZalcensteinBatch27_VMSC04501_region_2": "Batch27_region2_KO167",
    "G000218BenchmarkingBatch26_VMSC04501_region_2": "Batch26_region2_KO168",
    "Sample-709": "Sample-709",
    "SpleenBenchmarkingslide1_VMSC04501_Sample-710" : "Sample-710",
    "SpleenBenchmarkingslide1_VMSC04501_Sample-713" : "Sample-713"
}


dfs = []

for f in filtered_samples:
    # Read CSV
    df = pd.read_csv(os.path.join(parent_dir,f))

    dfs.append(df)

# Merge all DataFrames on 'Sample'
from functools import reduce
merged_df = reduce(lambda left, right: pd.merge(left, right, on="Sample", how="outer"), dfs)

# Reorder columns alphabetically (Sample first)
cols = ["Sample"] + sorted([c for c in merged_df.columns if c != "Sample"])
merged_df = merged_df[cols]

print(merged_df.head())

mapping_samples_2 = {
    "Batch22_region1_Ctrl174": "Ctrl174_Batch22_region1",
    "Batch22_region2_KO166": "KO166_Batch22_region2",
    "Batch26_region0_Ctrl173": "Ctrl173_Batch26_region0",
    "Batch26_region2_KO168": "KO168_Batch26_region2",
    "Batch27_region1_Ctrl172": "Ctrl172_Batch27_region1",
    "Batch27_region2_KO167": "KO167_Batch27_region2",
    "Sample-709": "WT709",
    "Sample-710": "WT710",
    "Sample-713": "WT713",
    }

merged_df["Sample"] = merged_df["Sample"].map(mapping_samples_2)
merged_df.to_csv(os.path.join(parent_dir, "merscope_summary_metrics_merged_filtered_final.csv"))


# %% Xenium Samples 

parent_dir = "/path/to/summary_metrics_xenium"

filtered_samples = [f for f in os.listdir(parent_dir) if "_filtered" in f]
unfiltered_samples = [f for f in os.listdir(parent_dir) if "unfiltered" in f]


dfs_filtered = []

for f in filtered_samples:
    # Read CSV

    df = pd.read_csv(os.path.join(parent_dir,f))
    
    if "proseg" in f:
        df["Sample"] = df["Sample"].str.split("__", n=1).str[1]
    else:
        df["Sample"] = df["Sample"].str.split("__").apply(lambda x: "__".join(x[1:3]))
        
    dfs_filtered.append(df)


dfs_unfiltered = []

for f in unfiltered_samples:
    # Read CSV
    df = pd.read_csv(os.path.join(parent_dir,f))
    
    if "proseg" in f:
        df["Sample"] = df["Sample"].str.split("__", n=1).str[1]
    else:
        df["Sample"] = df["Sample"].str.split("__").apply(lambda x: "__".join(x[1:3]))
        
    dfs_unfiltered.append(df)


# Merge all DataFrames on 'Sample'
from functools import reduce
merged_df_filtered = reduce(lambda left, right: pd.merge(left, right, on="Sample", how="outer"), dfs_filtered)
# Reorder columns alphabetically (Sample first)
cols_filtered = ["Sample"] + sorted([c for c in merged_df_filtered.columns if c != "Sample"])
merged_df_filtered = merged_df_filtered[cols_filtered]

merged_df_filtered.to_csv(os.path.join(parent_dir, r"xenium_summary_metrics_merged_filtered.csv"))


merged_df_unfiltered = reduce(lambda left, right: pd.merge(left, right, on="Sample", how="outer"), dfs_unfiltered)
cols_unfiltered = ["Sample"] + sorted([c for c in merged_df_unfiltered.columns if c != "Sample"])
merged_df_unfiltered = merged_df_unfiltered[cols_unfiltered]
merged_df_unfiltered.to_csv(os.path.join(parent_dir, r"xenium_summary_metrics_merged_unfiltered.csv"))

