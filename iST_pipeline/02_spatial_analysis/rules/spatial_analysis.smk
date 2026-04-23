import json

# Rules:   spatial_preprocess_sample, spatial_merge_samples, spatial_embed_harmony,
#          spatial_annotate, spatial_subcluster_gc, spatial_merge_gc_zones,
#          spatial_de_pseudobulk
# Purpose: Per-method spatial analysis pipeline — QC filtering, BANKSY spatial
#          feature extraction, cross-sample merging, Harmony-corrected
#          dimensionality reduction with Louvain clustering, cell type
#          annotation, GC B cell subclustering, and GC zone merging.
#
# Pipeline DAG per method:
#   spatial_preprocess_sample (one job per sample, run in parallel)
#     -> spatial_merge_samples  (one job per method)
#       -> spatial_embed_harmony (one job per method)
#
# Methods and sample sets:
#   xenium_default          — 9 non-batch34 Xenium samples (batch24 + batch27)
#   xenium_batch34_default  — 5 batch34 Xenium samples, default segmentation
#   xenium_batch34_cellpose — 5 batch34 Xenium samples, Cellpose segmentation
#   xenium_batch34_proseg   — 5 batch34 Xenium samples, Proseg segmentation
#   merscope_default        — 9 MERSCOPE samples, default segmentation
#   merscope_cellpose       — 9 MERSCOPE samples, Cellpose segmentation
#   merscope_proseg         — 9 MERSCOPE samples, Proseg segmentation
#
# Config keys used:
#   config["spatial_analysis"]["xenium_default_samples"] — 9 non-batch34 Xenium samples
#   config["spatial_analysis"]["xenium_batch34_samples"] — 5 batch34 Xenium samples
#   config["spatial_analysis"]["merscope_samples"]       — 9 MERSCOPE samples
#   config["spatial_analysis"]["platforms"]              — assay and banksy_lambda per platform
#   config["spatial_analysis"]["qc_min_counts"]          — minimum count filter threshold
#   config["spatial_analysis"]["scale_factor"]           — LogNormalize scale factor
#   config["spatial_analysis"]["banksy_k_geom"]          — BANKSY geometric neighbours
#   config["spatial_analysis"]["npcs"]                   — PCA components
#   config["spatial_analysis"]["dims"]                   — embedding dimensions (Harmony/UMAP)
#   config["spatial_analysis"]["resolution"]             — Leiden clustering resolution
#   config["spatial_analysis"]["umap_seed"]              — UMAP random seed
#   config["annotations"][method]                        — cluster->cell_type mapping per method
#   config["spatial_analysis"]["gc_dims"]                — dims for GC FindNeighbors and RunUMAP
#   config["spatial_analysis"]["gc_resolution"]          — Leiden resolution for GC subclustering
#   config["gc_zones"][method]                           — gc_subcluster->zone label mapping per method
#
# Named sub-targets (defined in Snakefile):
#   snakemake spatial_xenium_default
#   snakemake spatial_xenium_batch34_default
#   snakemake spatial_xenium_batch34_cellpose
#   snakemake spatial_xenium_batch34_proseg
#   snakemake spatial_merscope_default
#   snakemake spatial_merscope_cellpose
#   snakemake spatial_merscope_proseg
#   snakemake spatial_analysis              # all seven methods (preprocess->embed)
#   snakemake spatial_annotate_{method}        # annotate one method
#   snakemake spatial_annotate_all             # annotate all seven methods
#   snakemake spatial_subcluster_gc_{method}   # GC subclustering, one method
#   snakemake spatial_subcluster_gc_all        # GC subclustering, all seven methods
#   snakemake spatial_merge_gc_zones_{method}  # merge GC zones into full object, one method
#   snakemake spatial_merge_gc_zones_all       # merge GC zones, all seven methods
#   snakemake spatial_de_{method}              # pseudobulk DE for one method (4 DE methods only)
#   snakemake spatial_de_all                   # pseudobulk DE for all 4 DE methods


# ---------------------------------------------------------------------------
# Lookup: method -> (preprocessing_dir, seg_suffix)
# ---------------------------------------------------------------------------
# xenium_batch34_default reads from the xenium_default preprocessing outputs;
# there is no xenium_batch34_default directory under 01_preprocessing.

_METHOD_PREPROCESS = {
    "xenium_default":          ("xenium_default",    "default"),
    "xenium_batch34_default":  ("xenium_default",    "default"),
    "xenium_batch34_cellpose": ("xenium_cellpose",   "cellpose"),
    "xenium_batch34_proseg":   ("xenium_proseg",     "proseg"),
    "merscope_default":        ("merscope_default",  "default"),
    "merscope_cellpose":       ("merscope_cellpose", "cellpose"),
    "merscope_proseg":         ("merscope_proseg",   "proseg"),
}

# Constraint string shared across all three rules
_METHODS = (
    "xenium_default"
    "|xenium_batch34_default"
    "|xenium_batch34_cellpose"
    "|xenium_batch34_proseg"
    "|merscope_default"
    "|merscope_cellpose"
    "|merscope_proseg"
)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _platform_from_method(method):
    """Return platform from method string (e.g. 'xenium_batch34_default' -> 'xenium')."""
    return method.split("_")[0]


def _preprocess_dir(method):
    """Return the 01_preprocessing directory name for this method."""
    return _METHOD_PREPROCESS[method][0]


def _seg_from_method(method):
    """Return the seg suffix used in 01_preprocessing output filenames."""
    return _METHOD_PREPROCESS[method][1]


def _spatial_samples(method):
    """Return the configured sample list for a given method.

    xenium_default        -> xenium_default_samples (9 non-batch34 samples)
    xenium_batch34_*      -> xenium_batch34_samples (5 batch34 samples)
    merscope_*            -> merscope_samples (9 samples)
    """
    if method == "xenium_default":
        return config["spatial_analysis"]["xenium_default_samples"]
    elif method.startswith("xenium_batch34"):
        return config["spatial_analysis"]["xenium_batch34_samples"]
    return config["spatial_analysis"]["merscope_samples"]


def _platform_cfg(method, key):
    """Look up a platform-specific value (assay or banksy_lambda) from config."""
    platform = _platform_from_method(method)
    return config["spatial_analysis"]["platforms"][platform][key]


def _preprocess_inputs(wc):
    """Return the full list of 02_spatial_analysis/preprocessed RDS paths for one method's merge."""
    samples = _spatial_samples(wc.method)
    return expand(
        "results/02_spatial_analysis/{method}/preprocessed/{sample}.rds",
        method = wc.method,
        sample = samples,
    )


# ---------------------------------------------------------------------------
# Rule: spatial_preprocess_sample
# ---------------------------------------------------------------------------
# QC filter, LogNormalize, ScaleData, and RunBanksy for one sample.
# Input comes from 01_preprocessing; output is one RDS per sample per method.
#
# Config keys used:
#   config["spatial_analysis"]["platforms"][platform]["assay"]
#   config["spatial_analysis"]["platforms"][platform]["banksy_lambda"]
#   config["spatial_analysis"]["qc_min_counts"]
#   config["spatial_analysis"]["scale_factor"]
#   config["spatial_analysis"]["banksy_k_geom"]

rule spatial_preprocess_sample:
    wildcard_constraints:
        method = _METHODS
    input:
        rds = lambda wc: "results/01_preprocessing/{preprocess_dir}/{sample}_{seg}.rds".format(
            preprocess_dir = _preprocess_dir(wc.method),
            sample         = wc.sample,
            seg            = _seg_from_method(wc.method),
        )
    output:
        rds = "results/02_spatial_analysis/{method}/preprocessed/{sample}.rds"
    log:
        "logs/02_spatial_analysis/{method}/preprocess/{sample}.log"
    benchmark:
        "benchmarks/02_spatial_analysis/{method}/preprocess/{sample}.txt"
    params:
        platform      = lambda wc: _platform_from_method(wc.method),
        seg           = lambda wc: _seg_from_method(wc.method),
        assay         = lambda wc: _platform_cfg(wc.method, "assay"),
        banksy_lambda = lambda wc: _platform_cfg(wc.method, "banksy_lambda"),
        out_dir       = "results/02_spatial_analysis/{method}/preprocessed",
        qc_min_counts = config["spatial_analysis"]["qc_min_counts"],
        scale_factor  = config["spatial_analysis"]["scale_factor"],
        banksy_k_geom = config["spatial_analysis"]["banksy_k_geom"],
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0", 
        "ImageMagick/7.1.2-18"
        
    resources:
        mem_mb        = 200000,
        cpus_per_task = 16,
        runtime       = 1440,
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 02_spatial_analysis/R/preprocess_sample.R \
            --input_rds      {input.rds}           \
            --sample_name    {wildcards.sample}     \
            --platform       {params.platform}      \
            --seg            {params.seg}           \
            --assay          {params.assay}         \
            --qc_min_counts  {params.qc_min_counts} \
            --scale_factor   {params.scale_factor}  \
            --banksy_lambda  {params.banksy_lambda} \
            --banksy_k_geom  {params.banksy_k_geom} \
            --out_dir        {params.out_dir}       \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: spatial_merge_samples
# ---------------------------------------------------------------------------
# Merge all per-sample preprocessed objects for one method into a single
# combined Seurat object. Depends on all spatial_preprocess_sample outputs.
#
# Config keys used (indirectly, via _preprocess_inputs):
#   config["spatial_analysis"]["xenium_default_samples"]
#   config["spatial_analysis"]["xenium_batch34_samples"]
#   config["spatial_analysis"]["merscope_samples"]

rule spatial_merge_samples:
    wildcard_constraints:
        method = _METHODS
    input:
        preprocessed = _preprocess_inputs
    output:
        rds = "results/02_spatial_analysis/{method}/merged.rds"
    log:
        "logs/02_spatial_analysis/{method}/merge.log"
    benchmark:
        "benchmarks/02_spatial_analysis/{method}/merge.txt"
    params:
        input_dir = "results/02_spatial_analysis/{method}/preprocessed"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 200000,
        cpus_per_task = 16,
        runtime       = 600,
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 02_spatial_analysis/R/merge_samples.R \
            --input_dir {params.input_dir} \
            --method    {wildcards.method} \
            --out_file  {output.rds}       \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: spatial_embed_harmony
# ---------------------------------------------------------------------------
# Scale BANKSY, run PCA, Harmony batch correction (by sample_name), SNN graph,
# Louvain clustering, and UMAP. Produces the final analysis object per method.
#
# Config keys used:
#   config["spatial_analysis"]["npcs"]
#   config["spatial_analysis"]["dims"]
#   config["spatial_analysis"]["resolution"]
#   config["spatial_analysis"]["umap_seed"]

rule spatial_embed_harmony:
    wildcard_constraints:
        method = _METHODS
    input:
        rds = "results/02_spatial_analysis/{method}/merged.rds"
    output:
        rds = "results/02_spatial_analysis/{method}/embedded.rds"
    log:
        "logs/02_spatial_analysis/{method}/embed.log"
    benchmark:
        "benchmarks/02_spatial_analysis/{method}/embed.txt"
    params:
        npcs       = config["spatial_analysis"]["npcs"],
        dims       = config["spatial_analysis"]["dims"],
        resolution = config["spatial_analysis"]["resolution"],
        umap_seed  = config["spatial_analysis"]["umap_seed"],
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0",
        "ImageMagick/7.1.2-18"
    resources:
        mem_mb        = 300000,
        cpus_per_task = 16,
        runtime       = 900,
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 02_spatial_analysis/R/embed_harmony.R \
            --input_rds  {input.rds}          \
            --method     {wildcards.method}   \
            --npcs       {params.npcs}        \
            --dims       {params.dims}        \
            --resolution {params.resolution}  \
            --umap_seed  {params.umap_seed}   \
            --out_file   {output.rds}         \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: spatial_annotate
# ---------------------------------------------------------------------------
# Assigns cell type labels to clusters using a cluster->cell_type mapping
# from config["annotations"][method]. Produces a new cell_type metadata
# column; unmapped clusters are labelled "Unannotated".
#
# Populate config["annotations"][method] after inspecting embedded.rds UMAPs,
# then re-run this rule. Because the input (embedded.rds) is unchanged, only
# this rule and its downstream dependents will re-execute.
#
# Config keys used:
#   config["annotations"][method]  — dict of cluster_id (str) -> cell type name

rule spatial_annotate:
    wildcard_constraints:
        method = _METHODS
    input:
        rds = "results/02_spatial_analysis/{method}/embedded.rds"
    output:
        rds = "results/02_spatial_analysis/{method}/annotated.rds"
    log:
        "logs/02_spatial_analysis/{method}/annotate.log"
    benchmark:
        "benchmarks/02_spatial_analysis/{method}/annotate.txt"
    params:
        # Serialise the per-method annotation dict to a JSON string so it can
        # be passed as a single shell argument; .get() returns {} if the method
        # key is absent (e.g. before annotations have been filled in)
        annotations = lambda wc: json.dumps(
            config.get("annotations", {}).get(wc.method, {})
        )
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 100000,
        cpus_per_task = 4,
        runtime       = 120,
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 02_spatial_analysis/R/annotate_clusters.R \
            --input_rds   {input.rds}              \
            --method      {wildcards.method}       \
            --annotations '{params.annotations}'   \
            --out_file    {output.rds}             \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: spatial_subcluster_gc
# ---------------------------------------------------------------------------
# Subsets GC B cells (cell_type == "GC B cells") and runs an independent
# embedding: ScaleData -> PCA (pca.gc, npcs=30) -> Harmony (harmony.gc,
# grouped by sample_name) -> FindNeighbors (dims 1:gc_dims) -> Louvain
# clustering -> UMAP (umap.gc, dims 1:gc_dims).
# Cluster labels are stored in gc_subcluster.
#
# Config keys used:
#   config["spatial_analysis"]["gc_dims"]       — dims for FindNeighbors and RunUMAP
#   config["spatial_analysis"]["gc_resolution"] — Louvain resolution
#   config["spatial_analysis"]["umap_seed"]     — UMAP random seed

rule spatial_subcluster_gc:
    wildcard_constraints:
        method = _METHODS
    input:
        rds = "results/02_spatial_analysis/{method}/annotated.rds"
    output:
        rds = "results/02_spatial_analysis/{method}/subclustered_gc.rds"
    log:
        "logs/02_spatial_analysis/{method}/subcluster_gc.log"
    benchmark:
        "benchmarks/02_spatial_analysis/{method}/subcluster_gc.txt"
    params:
        dims       = config["spatial_analysis"]["gc_dims"],
        resolution = config["spatial_analysis"]["gc_resolution"],
        umap_seed  = config["spatial_analysis"]["umap_seed"],
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0",
        "ImageMagick/7.1.2-18"
    resources:
        mem_mb        = 200000,
        cpus_per_task = 16,
        runtime       = 600,
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 02_spatial_analysis/R/subcluster_gc.R \
            --input_rds  {input.rds}          \
            --method     {wildcards.method}   \
            --dims       {params.dims}        \
            --resolution {params.resolution}  \
            --umap_seed  {params.umap_seed}   \
            --out_file   {output.rds}         \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: spatial_merge_gc_zones
# ---------------------------------------------------------------------------
# Injects GC B cell zone labels (Dark Zone / Light Zone) from the GC subset
# back into the full annotated object. Adds two new metadata columns:
#   gc_zone    — zone label for GC cells, NA for all other cells
#   cell_type2 — cell_type with GC B cells replaced by their zone label;
#                GC cells with no zone mapping keep "GC B cells"
#
# Populate config["gc_zones"][method] after inspecting subclustered_gc.rds,
# then re-run this rule. Input files are unchanged so only this rule fires.
# If gc_zones is empty the rule completes without error and cell_type2 == cell_type.
#
# Pipeline position:
#   spatial_annotate -> spatial_subcluster_gc -> spatial_merge_gc_zones
#
# Config keys used:
#   config["gc_zones"][method]  — dict of gc_subcluster_id (str) -> zone label

rule spatial_merge_gc_zones:
    wildcard_constraints:
        method = _METHODS
    input:
        annotated = "results/02_spatial_analysis/{method}/annotated.rds",
        gc_subset = "results/02_spatial_analysis/{method}/subclustered_gc.rds"
    output:
        rds = "results/02_spatial_analysis/{method}/annotated_final.rds"
    log:
        "logs/02_spatial_analysis/{method}/merge_gc_zones.log"
    benchmark:
        "benchmarks/02_spatial_analysis/{method}/merge_gc_zones.txt"
    params:
        # Serialise the per-method gc_zones dict to a JSON string; .get() returns
        # {} if the method key is absent (before zones have been assigned)
        gc_zones = lambda wc: json.dumps(
            config.get("gc_zones", {}).get(wc.method, {})
        )
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 100000,
        cpus_per_task = 4,
        runtime       = 60,
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 02_spatial_analysis/R/merge_gc_zones.R \
            --input_annotated {input.annotated}      \
            --input_gc        {input.gc_subset}      \
            --method          {wildcards.method}     \
            --gc_zones        '{params.gc_zones}'    \
            --out_file        {output.rds}           \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: spatial_de_pseudobulk
# ---------------------------------------------------------------------------
# Pseudobulk KO vs WT differential expression per cell type using a
# limma-voom pipeline. Only runs on methods that have KO samples:
#   xenium_default, merscope_default, merscope_cellpose, merscope_proseg
# xenium_batch34_* methods are excluded — that cohort has no KO samples.
#
# Config keys used:
#   config["pseudobulk_de"]["assay"]          — platform -> assay name
#   config["pseudobulk_de"]["cell_type_col"]  — cell type metadata column
#   config["pseudobulk_de"]["condition_col"]  — condition metadata column
#   config["pseudobulk_de"]["animal_id_col"]  — replicate ID metadata column
#   config["pseudobulk_de"]["min_reps"]       — minimum replicates per group
#   config["pseudobulk_de"]["contrast_name"]  — contrast coefficient name

_DE_METHODS = (
    "xenium_default"
    "|merscope_default"
    "|merscope_cellpose"
    "|merscope_proseg"
)


rule spatial_de_pseudobulk:
    wildcard_constraints:
        method = _DE_METHODS
    input:
        rds = "results/02_spatial_analysis/{method}/annotated.rds"
    output:
        rds = "results/02_spatial_analysis/{method}/pseudobulk_de/de_results.rds"
    log:
        "logs/02_spatial_analysis/{method}/pseudobulk_de.log"
    benchmark:
        "benchmarks/02_spatial_analysis/{method}/pseudobulk_de.txt"
    params:
        assay         = lambda wc: config["pseudobulk_de"]["assay"][_platform_from_method(wc.method)],
        cell_type_col = config["pseudobulk_de"]["cell_type_col"],
        condition_col = config["pseudobulk_de"]["condition_col"],
        animal_id_col = config["pseudobulk_de"]["animal_id_col"],
        min_reps      = config["pseudobulk_de"]["min_reps"],
        contrast_name = config["pseudobulk_de"]["contrast_name"],
    envmodules:
        "R/4.4.1"
    resources:
        mem_mb        = 150000,
        cpus_per_task = 4,
        runtime       = 240,
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 02_spatial_analysis/R/pseudobulk_de.R \
            --input_rds      {input.rds}               \
            --method         {wildcards.method}        \
            --assay          {params.assay}            \
            --cell_type_col  {params.cell_type_col}    \
            --condition_col  {params.condition_col}    \
            --animal_id_col  {params.animal_id_col}    \
            --min_reps       {params.min_reps}         \
            --contrast_name  {params.contrast_name}    \
            --out_file       {output.rds}              \
            > {log} 2>&1
        """
