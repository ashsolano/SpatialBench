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
        partition     = "regular"
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
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 01_preprocessing/R/create_seurat_binned_merscope.R \
            --data_dir    {input.data_dir} \
            --sample_name {wildcards.sample} \
            --resolution  {wildcards.resolution} \
            --out_dir     {params.out_dir} \
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
        partition     = "regular"
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
        partition     = "regular"
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
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 01_preprocessing/R/create_seurat_segmented_proseg.R \
            --data_dir    {input.data_dir} \
            --sample_name {wildcards.sample} \
            --assay       {params.assay} \
            --out_dir     {params.out_dir} \
            > {log} 2>&1
        """
