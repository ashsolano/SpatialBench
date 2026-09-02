import os


# Rule:   create_seurat_binned_xenium
# Purpose: Create a binned Seurat object for one Xenium sample at a given resolution.
#
# Config keys used:
#   config["xenium"]["samples"]   — dict of sample name -> full path to data directory
#   config["bin_resolutions"]     — list of resolutions, e.g. [8, 16]

rule create_seurat_binned_xenium:
    input:
        data_dir = lambda wc: config["xenium"]["samples"][wc.sample]
    output:
        rds = "results/01_preprocessing/xenium_{resolution}um/{sample}_{resolution}um.rds"
    log:
        "logs/01_preprocessing/xenium_{resolution}um/{sample}.log"
    benchmark:
        "benchmarks/01_preprocessing/xenium_{resolution}um/{sample}.txt"
    params:
        out_dir = "results/01_preprocessing/xenium_{resolution}um"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 60000,
        cpus_per_task = 4,
        runtime       = 600,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 01_preprocessing/R/create_seurat_binned_xenium.R \
            --data_dir    {input.data_dir} \
            --sample_name {wildcards.sample} \
            --resolution  {wildcards.resolution} \
            --out_dir     {params.out_dir} \
            > {log} 2>&1
        """


# Rule:   create_seurat_binned_merscope
# Purpose: Create a binned Seurat object for one MERSCOPE sample at a given resolution.
#
# Config keys used:
#   config["merscope"]["samples"]  — dict of sample name -> full path to data directory
#   config["bin_resolutions"]      — list of resolutions, e.g. [8, 16]

rule create_seurat_binned_merscope:
    input:
        data_dir = lambda wc: config["merscope"]["samples"][wc.sample]
    output:
        rds = "results/01_preprocessing/merscope_{resolution}um/{sample}_{resolution}um.rds"
    log:
        "logs/01_preprocessing/merscope_{resolution}um/{sample}.log"
    benchmark:
        "benchmarks/01_preprocessing/merscope_{resolution}um/{sample}.txt"
    params:
        out_dir = "results/01_preprocessing/merscope_{resolution}um"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 40000,
        cpus_per_task = 4,
        runtime       = 600,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 01_preprocessing/R/create_seurat_binned_merscope.R \
            --data_dir    {input.data_dir} \
            --sample_name {wildcards.sample} \
            --resolution  {wildcards.resolution} \
            --out_dir     {params.out_dir} \
            > {log} 2>&1
        """


# Rule:   filter_binned_xenium
# Purpose: Post-processing QC on a binned Xenium object — drop empty bins and
#          remove spatially isolated bins (DBSCAN) outside the dominant tissue
#          cluster. Also writes a before/after spatial QC plot alongside the
#          filtered RDS.
#
# Config keys used:
#   config["xenium"]["samples"]   — dict of sample name -> full path to data directory
#   config["bin_resolutions"]     — list of resolutions, e.g. [8, 16]

rule filter_binned_xenium:
    input:
        rds = "results/01_preprocessing/xenium_{resolution}um/{sample}_{resolution}um.rds"
    output:
        rds = "results/01_preprocessing/xenium_{resolution}um_filtered/{sample}_{resolution}um_filtered.rds"
    log:
        "logs/01_preprocessing/xenium_{resolution}um_filtered/{sample}.log"
    benchmark:
        "benchmarks/01_preprocessing/xenium_{resolution}um_filtered/{sample}.txt"
    params:
        assay = "Xenium",
        # DBSCAN neighbourhood radius scales with bin size: 15um for 8um bins,
        # 20um for 16um bins
        eps   = lambda wc: 15 if wc.resolution == "8" else 20,
        # Xenium keeps only substantial tissue clusters — 5000 drops small
        # fragments that were slipping through at the previous 1000 threshold
        min_cluster_size = 5000
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 20000,
        cpus_per_task = 2,
        runtime       = 120,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 01_preprocessing/R/filter_binned.R \
            --input_rds {input.rds} \
            --out_rds   {output.rds} \
            --assay     {params.assay} \
            --eps       {params.eps} \
            --min_cluster_size {params.min_cluster_size} \
            --use_adaptive_threshold \
            > {log} 2>&1
        """


# Rule:   filter_binned_merscope
# Purpose: Post-processing QC on a binned MERSCOPE object — drop empty bins
#          and remove spatially isolated bins (DBSCAN) outside the dominant
#          tissue cluster. Also writes a before/after spatial QC plot
#          alongside the filtered RDS.
#
# Config keys used:
#   config["merscope"]["samples"]  — dict of sample name -> full path to data directory
#   config["bin_resolutions"]      — list of resolutions, e.g. [8, 16]

rule filter_binned_merscope:
    input:
        rds = "results/01_preprocessing/merscope_{resolution}um/{sample}_{resolution}um.rds"
    output:
        rds = "results/01_preprocessing/merscope_{resolution}um_filtered/{sample}_{resolution}um_filtered.rds"
    log:
        "logs/01_preprocessing/merscope_{resolution}um_filtered/{sample}.log"
    benchmark:
        "benchmarks/01_preprocessing/merscope_{resolution}um_filtered/{sample}.txt"
    params:
        assay     = "Vizgen",
        # DBSCAN neighbourhood radius scales with bin size: 25um for 8um bins,
        # 30um for 16um bins
        eps       = lambda wc: 25 if wc.resolution == "8" else 30,
        # MERSCOPE has no comparable background signal to Xenium, so adaptive
        # thresholding is disabled in favour of a fixed, low min_count
        min_count = 1,
        # MERSCOPE has legitimate smaller tissue pieces, so keep the default
        # (lower) cluster size threshold rather than Xenium's 5000
        min_cluster_size = 1000
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 20000,
        cpus_per_task = 2,
        runtime       = 120,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 01_preprocessing/R/filter_binned.R \
            --input_rds {input.rds} \
            --out_rds   {output.rds} \
            --assay     {params.assay} \
            --eps       {params.eps} \
            --min_count {params.min_count} \
            --min_cluster_size {params.min_cluster_size} \
            --no_adaptive_threshold \
            > {log} 2>&1
        """


# Rule:   combine_filter_qc_xenium
# Purpose: Combine every sample's before/after filtering QC PNG (produced by
#          filter_binned_xenium) into a single multi-page PDF for one
#          resolution. Depends on all filter_binned_xenium outputs for that
#          resolution, so filtering finishes for every sample first.
#
# Config keys used:
#   config["xenium"]["samples"]   — dict of sample name -> full path to data directory

rule combine_filter_qc_xenium:
    input:
        rds = lambda wc: expand(
            "results/01_preprocessing/xenium_{resolution}um_filtered/{sample}_{resolution}um_filtered.rds",
            sample     = config["xenium"]["samples"],
            resolution = wc.resolution,
        )
    output:
        pdf = "results/01_preprocessing/xenium_{resolution}um_filtered/filter_qc_all.pdf"
    log:
        "logs/01_preprocessing/xenium_{resolution}um_filtered/combine_filter_qc.log"
    benchmark:
        "benchmarks/01_preprocessing/xenium_{resolution}um_filtered/combine_filter_qc.txt"
    envmodules:
        "ImageMagick/7.1.2-18"
    resources:
        mem_mb        = 4000,
        cpus_per_task = 1,
        runtime       = 30,
        slurm_partition     = "regular"
    shell:
        """
        convert results/01_preprocessing/xenium_{wildcards.resolution}um_filtered/*_filter_qc.png {output.pdf} \
            > {log} 2>&1
        """


# Rule:   combine_filter_qc_merscope
# Purpose: Combine every sample's before/after filtering QC PNG (produced by
#          filter_binned_merscope) into a single multi-page PDF for one
#          resolution. Depends on all filter_binned_merscope outputs for that
#          resolution, so filtering finishes for every sample first.
#
# Config keys used:
#   config["merscope"]["samples"]  — dict of sample name -> full path to data directory

rule combine_filter_qc_merscope:
    input:
        rds = lambda wc: expand(
            "results/01_preprocessing/merscope_{resolution}um_filtered/{sample}_{resolution}um_filtered.rds",
            sample     = config["merscope"]["samples"],
            resolution = wc.resolution,
        )
    output:
        pdf = "results/01_preprocessing/merscope_{resolution}um_filtered/filter_qc_all.pdf"
    log:
        "logs/01_preprocessing/merscope_{resolution}um_filtered/combine_filter_qc.log"
    benchmark:
        "benchmarks/01_preprocessing/merscope_{resolution}um_filtered/combine_filter_qc.txt"
    envmodules:
        "ImageMagick/7.1.2-18"
    resources:
        mem_mb        = 4000,
        cpus_per_task = 1,
        runtime       = 30,
        slurm_partition     = "regular"
    shell:
        """
        convert results/01_preprocessing/merscope_{wildcards.resolution}um_filtered/*_filter_qc.png {output.pdf} \
            > {log} 2>&1
        """


# Rule:   align_binned_xenium
# Purpose: Apply a previously-fitted STalign registration (affine + LDDMM) to
#          a filtered, binned Xenium object, mapping its bin centroids into
#          VisiumHD target coordinate space. Only runs for the 4 samples with
#          manually-picked landmarks (config["stalign"]["landmarks"]).
#
# Config keys used:
#   config["stalign"] — python_bin, matched_samples, visiumhd,
#                        visiumhd_scalefactors, landmarks
#   config["visiumhd"] — data_dir, samples

rule align_binned_xenium:
    wildcard_constraints:
        sample     = "ko167_batch24|ko168_batch27|wt709_batch27|wt713_batch24",
        resolution = "8|16"
    input:
        rds = "results/01_preprocessing/xenium_{resolution}um_filtered/{sample}_{resolution}um_filtered.rds"
    output:
        rds = "results/01_preprocessing/xenium_{resolution}um_aligned/{sample}_{resolution}um_aligned.rds"
    log:
        "logs/01_preprocessing/xenium_{resolution}um_aligned/{sample}.log"
    benchmark:
        "benchmarks/01_preprocessing/xenium_{resolution}um_aligned/{sample}.txt"
    params:
        animal_id = lambda wc: wc.sample.split("_batch")[0]
    resources:
        mem_mb        = 40000,
        runtime       = 30,
        slurm_partition     = "gpuq",
        gres          = "gpu:A30:1"
    shell:
        """
        module load R/4.4.1 && Rscript --vanilla 01_preprocessing/R/align_binned.R \
            --input_rds  {input.rds} \
            --platform   xenium \
            --sample     {params.animal_id} \
            --resolution {wildcards.resolution} \
            --config     config/config.yaml \
            --out_rds    {output.rds} \
            > {log} 2>&1
        """


# Rule:   align_binned_merscope
# Purpose: Apply a previously-fitted STalign registration (affine + LDDMM) to
#          a filtered, binned MERSCOPE object, mapping its bin centroids into
#          VisiumHD target coordinate space. Only runs for the 4 samples with
#          manually-picked landmarks (config["stalign"]["landmarks"]).
#
# Config keys used:
#   config["stalign"] — python_bin, matched_samples, visiumhd,
#                        visiumhd_scalefactors, landmarks
#   config["visiumhd"] — data_dir, samples

rule align_binned_merscope:
    wildcard_constraints:
        sample     = "ko167_batch10|ko168_batch9|wt709_batch13|wt713_batch13",
        resolution = "8|16"
    input:
        rds = "results/01_preprocessing/merscope_{resolution}um_filtered/{sample}_{resolution}um_filtered.rds"
    output:
        rds = "results/01_preprocessing/merscope_{resolution}um_aligned/{sample}_{resolution}um_aligned.rds"
    log:
        "logs/01_preprocessing/merscope_{resolution}um_aligned/{sample}.log"
    benchmark:
        "benchmarks/01_preprocessing/merscope_{resolution}um_aligned/{sample}.txt"
    params:
        animal_id = lambda wc: wc.sample.split("_batch")[0]
    resources:
        mem_mb        = 40000,
        runtime       = 30,
        slurm_partition     = "gpuq",
        gres          = "gpu:A30:1"
    shell:
        """
        module load R/4.4.1 && Rscript --vanilla 01_preprocessing/R/align_binned.R \
            --input_rds  {input.rds} \
            --platform   merscope \
            --sample     {params.animal_id} \
            --resolution {wildcards.resolution} \
            --config     config/config.yaml \
            --out_rds    {output.rds} \
            > {log} 2>&1
        """


# Rule:   create_seurat_segmented_xenium
# Purpose: Create a cell-segmented Seurat object for one Xenium sample using
#          a named segmentation method (default vendor or Cellpose).
#
# Config keys used:
#   config["xenium_default"]["samples"]   — dict of sample name -> data directory
#   config["xenium_cellpose"]["samples"]  — dict of sample name -> data directory

def _segmented_xenium_path(wc):
    # method=default reuses the shared xenium sample paths; other methods
    # (e.g. cellpose) have their own config section
    section = "xenium" if wc.method == "default" else f"xenium_{wc.method}"
    path = config[section]["samples"][wc.sample]
    if path is None:
        raise ValueError(
            f"Path not set for {section} sample '{wc.sample}'. "
            f"Please fill in the path in config/config.yaml."
        )
    return path


rule create_seurat_segmented_xenium:
    wildcard_constraints:
        method = "default|cellpose"
    input:
        data_dir = _segmented_xenium_path
    output:
        rds = "results/01_preprocessing/xenium_{method}/{sample}_{method}.rds"
    log:
        "logs/01_preprocessing/xenium_{method}/{sample}.log"
    benchmark:
        "benchmarks/01_preprocessing/xenium_{method}/{sample}.txt"
    params:
        out_dir = "results/01_preprocessing/xenium_{method}"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 50000,
        cpus_per_task = 4,
        runtime       = 600,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 01_preprocessing/R/create_seurat_segmented_xenium.R \
            --data_dir    {input.data_dir} \
            --sample_name {wildcards.sample} \
            --method      {wildcards.method} \
            --out_dir     {params.out_dir} \
            > {log} 2>&1
        """


# Rule:   create_seurat_segmented_merscope
# Purpose: Create a cell-segmented Seurat object for one MERSCOPE sample using
#          a named segmentation method (default vendor or Cellpose).
#
# Config keys used:
#   config["merscope"]["samples"]          — default: paths ending in /Vizgen
#   config["merscope_cellpose"]["samples"] — Cellpose: paths ending in /Cellpose

def _segmented_merscope_path(wc):
    # method=default reuses the shared merscope sample paths; other methods
    # (e.g. cellpose) have their own config section
    section = "merscope" if wc.method == "default" else f"merscope_{wc.method}"
    path = config[section]["samples"][wc.sample]
    if path is None:
        raise ValueError(
            f"Path not set for {section} sample '{wc.sample}'. "
            f"Please fill in the path in config/config.yaml."
        )
    return path


rule create_seurat_segmented_merscope:
    wildcard_constraints:
        method = "default|cellpose"
    input:
        data_dir = _segmented_merscope_path
    output:
        rds = "results/01_preprocessing/merscope_{method}/{sample}_{method}.rds"
    log:
        "logs/01_preprocessing/merscope_{method}/{sample}.log"
    benchmark:
        "benchmarks/01_preprocessing/merscope_{method}/{sample}.txt"
    params:
        out_dir = "results/01_preprocessing/merscope_{method}"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 50000,
        cpus_per_task = 4,
        runtime       = 600,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 01_preprocessing/R/create_seurat_segmented_merscope.R \
            --data_dir    {input.data_dir} \
            --sample_name {wildcards.sample} \
            --method      {wildcards.method} \
            --out_dir     {params.out_dir} \
            > {log} 2>&1
        """


# Rule:   create_seurat_segmented_proseg
# Purpose: Create a Proseg-segmented Seurat object for one sample on either platform.
#          A single rule covers both Xenium ({platform}=xenium, assay=Xenium) and
#          MERSCOPE ({platform}=merscope, assay=Vizgen).
#
# Config keys used:
#   config["xenium_proseg"]["samples"]   — dict of sample name -> Proseg output directory
#   config["merscope_proseg"]["samples"] — dict of sample name -> Proseg output directory

def _proseg_path(wc):
    section = f"{wc.platform}_proseg"
    path = config[section]["samples"][wc.sample]
    if path is None:
        raise ValueError(
            f"Path not set for {section} sample '{wc.sample}'. "
            f"Please fill in the path in config/config.yaml."
        )
    return path


rule create_seurat_segmented_proseg:
    input:
        data_dir = _proseg_path
    output:
        rds = "results/01_preprocessing/{platform}_proseg/{sample}_proseg.rds"
    log:
        "logs/01_preprocessing/{platform}_proseg/{sample}.log"
    benchmark:
        "benchmarks/01_preprocessing/{platform}_proseg/{sample}.txt"
    params:
        out_dir = "results/01_preprocessing/{platform}_proseg",
        assay   = lambda wc: "Xenium" if wc.platform == "xenium" else "Vizgen"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 60000,
        cpus_per_task = 4,
        runtime       = 600,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 01_preprocessing/R/create_seurat_segmented_proseg.R \
            --data_dir    {input.data_dir} \
            --sample_name {wildcards.sample} \
            --assay       {params.assay} \
            --out_dir     {params.out_dir} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# ROI extraction (01_preprocessing/R/compute_roi_box.R, extract_roi.R)
# ---------------------------------------------------------------------------
# Extracts a fixed physical region of interest (config["roi"]["size_um"]),
# centred on each animal's VisiumHD tissue extent, from VisiumHD and matched
# aligned MERSCOPE/Xenium objects. Only runs for the 4 matched animals
# (config["stalign"]["matched_samples"]).


def _roi_size_label(size_um):
    # e.g. 2000 -> "2mm"; 500 -> "500um" — keeps ROI size unambiguous from
    # the filename alone without a distracting "2000um" for round mm sizes.
    if size_um % 1000 == 0:
        return f"{size_um // 1000}mm"
    return f"{size_um}um"


ROI_SIZE_LABEL = _roi_size_label(config["roi"]["size_um"])

# VisiumHD sample ID (e.g. "batch33_709") -> animal ID (e.g. "wt709")
_VISIUMHD_TO_ANIMAL = {v: k for k, v in config["stalign"]["visiumhd"].items()}


# Rule:   compute_roi_box
# Purpose: Compute one animal's ROI box definition (center + bounds, in
#          microns) from its VisiumHD reference object. Run once per animal;
#          reused by extract_roi_visium/merscope/xenium below.
#
# Config keys used:
#   config["roi"]["size_um"]
#   config["stalign"] — visiumhd, microns_per_pixel
#   config["visiumhd"] — data_dir, samples

rule compute_roi_box:
    wildcard_constraints:
        sample = "ko167|ko168|wt709|wt713"
    input:
        rds = lambda wc: os.path.join(
            config["visiumhd"]["data_dir"],
            config["visiumhd"]["samples"][config["stalign"]["visiumhd"][wc.sample]],
        )
    output:
        rds = "results/01_preprocessing/roi_boxes/{sample}_roi_box_" + ROI_SIZE_LABEL + ".rds"
    log:
        "logs/01_preprocessing/roi_boxes/{sample}.log"
    benchmark:
        "benchmarks/01_preprocessing/roi_boxes/{sample}.txt"
    resources:
        mem_mb  = 16000,
        runtime = 15,
        slurm_partition = "regular"
    shell:
        """
        module load R/4.4.1 && Rscript --vanilla 01_preprocessing/R/compute_roi_box.R \
            --sample   {wildcards.sample} \
            --config   config/config.yaml \
            --out_rds  {output.rds} \
            > {log} 2>&1
        """


# Rule:   extract_roi_visium
# Purpose: Subset a VisiumHD reference object to one animal's ROI box.
#
# Config keys used:
#   config["visiumhd"] — data_dir, samples
#   config["stalign"]["visiumhd"] — animal -> VisiumHD sample ID

rule extract_roi_visium:
    wildcard_constraints:
        visium_sample = "batch33_167|batch33_168|batch33_709|batch33_713"
    input:
        rds         = lambda wc: os.path.join(
            config["visiumhd"]["data_dir"],
            config["visiumhd"]["samples"][wc.visium_sample],
        ),
        roi_box_rds = lambda wc: (
            "results/01_preprocessing/roi_boxes/"
            + _VISIUMHD_TO_ANIMAL[wc.visium_sample] + "_roi_box_" + ROI_SIZE_LABEL + ".rds"
        )
    output:
        rds = "results/01_preprocessing/visium_8um_roi/{visium_sample}_8um_roi_" + ROI_SIZE_LABEL + ".rds"
    log:
        "logs/01_preprocessing/visium_8um_roi/{visium_sample}.log"
    benchmark:
        "benchmarks/01_preprocessing/visium_8um_roi/{visium_sample}.txt"
    resources:
        mem_mb  = 40000,
        runtime = 30,
        slurm_partition = "regular"
    shell:
        """
        module load R/4.4.1 && Rscript --vanilla 01_preprocessing/R/extract_roi.R \
            --input_rds   {input.rds} \
            --platform    visium \
            --roi_box_rds {input.roi_box_rds} \
            --out_rds     {output.rds} \
            > {log} 2>&1
        """


# Rule:   extract_roi_merscope
# Purpose: Subset an aligned MERSCOPE object to its matched animal's ROI box.
#
# Config keys used:
#   config["stalign"]["matched_samples"]["merscope"]

rule extract_roi_merscope:
    wildcard_constraints:
        sample = "ko167_batch10|ko168_batch9|wt709_batch13|wt713_batch13"
    input:
        rds         = "results/01_preprocessing/merscope_8um_aligned/{sample}_8um_aligned.rds",
        roi_box_rds = lambda wc: (
            "results/01_preprocessing/roi_boxes/"
            + wc.sample.split("_batch")[0] + "_roi_box_" + ROI_SIZE_LABEL + ".rds"
        )
    output:
        rds = "results/01_preprocessing/merscope_8um_roi/{sample}_8um_roi_" + ROI_SIZE_LABEL + ".rds"
    log:
        "logs/01_preprocessing/merscope_8um_roi/{sample}.log"
    benchmark:
        "benchmarks/01_preprocessing/merscope_8um_roi/{sample}.txt"
    resources:
        mem_mb  = 40000,
        runtime = 30,
        slurm_partition = "regular"
    shell:
        """
        module load R/4.4.1 && Rscript --vanilla 01_preprocessing/R/extract_roi.R \
            --input_rds   {input.rds} \
            --platform    merscope \
            --roi_box_rds {input.roi_box_rds} \
            --out_rds     {output.rds} \
            > {log} 2>&1
        """


# Rule:   extract_roi_xenium
# Purpose: Subset an aligned Xenium object to its matched animal's ROI box.
#
# Config keys used:
#   config["stalign"]["matched_samples"]["xenium"]

rule extract_roi_xenium:
    wildcard_constraints:
        sample = "ko167_batch24|ko168_batch27|wt709_batch27|wt713_batch24"
    input:
        rds         = "results/01_preprocessing/xenium_8um_aligned/{sample}_8um_aligned.rds",
        roi_box_rds = lambda wc: (
            "results/01_preprocessing/roi_boxes/"
            + wc.sample.split("_batch")[0] + "_roi_box_" + ROI_SIZE_LABEL + ".rds"
        )
    output:
        rds = "results/01_preprocessing/xenium_8um_roi/{sample}_8um_roi_" + ROI_SIZE_LABEL + ".rds"
    log:
        "logs/01_preprocessing/xenium_8um_roi/{sample}.log"
    benchmark:
        "benchmarks/01_preprocessing/xenium_8um_roi/{sample}.txt"
    resources:
        mem_mb  = 40000,
        runtime = 30,
        slurm_partition = "regular"
    shell:
        """
        module load R/4.4.1 && Rscript --vanilla 01_preprocessing/R/extract_roi.R \
            --input_rds   {input.rds} \
            --platform    xenium \
            --roi_box_rds {input.roi_box_rds} \
            --out_rds     {output.rds} \
            > {log} 2>&1
        """


rule visium_binning_roi:
    input:
        expand(
            "results/01_preprocessing/visium_8um_roi/{visium_sample}_8um_roi_" + ROI_SIZE_LABEL + ".rds",
            visium_sample = config["stalign"]["visiumhd"].values(),
        )


rule merscope_binning_roi:
    input:
        expand(
            "results/01_preprocessing/merscope_8um_roi/{sample}_8um_roi_" + ROI_SIZE_LABEL + ".rds",
            sample = config["stalign"]["matched_samples"]["merscope"].values(),
        )


rule xenium_binning_roi:
    input:
        expand(
            "results/01_preprocessing/xenium_8um_roi/{sample}_8um_roi_" + ROI_SIZE_LABEL + ".rds",
            sample = config["stalign"]["matched_samples"]["xenium"].values(),
        )
