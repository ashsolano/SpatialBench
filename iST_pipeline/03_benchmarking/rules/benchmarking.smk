# Purpose: Cross-platform benchmarking metrics and quality checks.
# Rules:   dataset_summary, scrna_correlation, qc_metrics, qc_metrics_roi, qc_spatial, moransi,
#          qc_backgrounds, probe_rank, gene_comparison, variance_partition,
#          segmentation_quality
#
# DAG (from 01_preprocessing binned objects; single job unless noted):
#   dataset_summary      all platforms, all resolutions
#   scrna_correlation    8µm; FLEX + 3' GEX; needs gene_lists
#   qc_metrics           all resolutions; needs gene_lists
#   qc_spatial           8µm; filtered + aligned objects; needs qc_metrics
#   moransi              one job per platform (MERSCOPE, Xenium); filtered 8µm
#   qc_backgrounds       MERSCOPE + Xenium, filtered 8µm; needs moransi
#   probe_rank           MERSCOPE + Xenium, 8µm
#   gene_comparison      all platforms, 8µm -> variance_partition
#
# Config keys used:
#   visiumhd.data_dir, visiumhd.samples          VisiumHD objects (external) and mapping
#   spatial_analysis.merscope_samples            9 MERSCOPE samples
#   spatial_analysis.xenium_default_samples      9 Xenium samples
#   bin_resolutions                              [8, 16]
#   scrna.path, scrna.sc_path                    FLEX / 3' GEX references (external)
#   scrna.wt_ids                                 matched WT animals (709, 713)
#   output_dir                                   results root
#
# Named sub-targets (Snakefile): benchmarking_{dataset_summary, scrna_correlation,
#   qc_metrics, qc_backgrounds, probe_rank, gene_comparison, variance_partition,
#   segmentation_quality}, benchmarking_all

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _binning_inputs(platforms, resolutions, filtered = False):
    """Expand preprocessing binning paths for each platform/sample/resolution.

    filtered = True points at the post-QC (empty-bin + DBSCAN) objects written
    by filter_binned_{xenium,merscope}; False points at the raw binned objects.
    """
    suffix = "_filtered" if filtered else ""
    inputs = []
    for platform, samples in platforms.items():
        inputs += expand(
            "results/01_preprocessing/{platform}_{res}um{suffix}/{sample}_{res}um{suffix}.rds",
            platform = platform,
            sample   = samples,
            res      = resolutions,
            suffix   = suffix,
        )
    return inputs


def _qc_bg_samples(platform):
    """Samples used for qc_backgrounds / moransi: the platform's 8µm sample
    list minus config["qc_backgrounds"]["exclude_samples"][platform]."""
    samples = {
        "merscope": config["spatial_analysis"]["merscope_samples"],
        "xenium":   config["spatial_analysis"]["xenium_default_samples"],
    }[platform]
    excluded = config["qc_backgrounds"]["exclude_samples"].get(platform) or []
    return [s for s in samples if s not in excluded]


# ---------------------------------------------------------------------------
# Rule: dataset_summary
# ---------------------------------------------------------------------------
# Per-sample dataset metrics and gene-panel lists for fig1.R.
# VisiumHD objects are external (config["visiumhd"]["data_dir"]), not listed in input:.

rule dataset_summary:
    input:
        _binning_inputs(
            platforms = {
                "merscope": config["spatial_analysis"]["merscope_samples"],
                "xenium":   config["spatial_analysis"]["xenium_default_samples"],
            },
            resolutions = config["bin_resolutions"],
            filtered    = True,
        )
    output:
        metrics    = "results/03_benchmarking/dataset_summary/metrics.rds",
        gene_lists = "results/03_benchmarking/dataset_summary/gene_lists.rds"
    log:
        "logs/03_benchmarking/dataset_summary.log"
    benchmark:
        "benchmarks/03_benchmarking/dataset_summary.txt"
    params:
        out_dir = "results/03_benchmarking/dataset_summary"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 300000,
        cpus_per_task = 16,
        runtime       = 360,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/dataset_summary.R \
            --config  config/config.yaml \
            --out_dir {params.out_dir}   \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: scrna_correlation
# ---------------------------------------------------------------------------
# Pseudobulk log10(CPM + 1) of the matched WT animals (config["scrna"]["wt_ids"])
# vs the FLEX and 3' GEX references, with Pearson r:
#   avg_expr.rds                  FLEX only, for fig1c (fig1.R)
#   correlation_by_reference.rds  both references, pairwise and common gene
#                                 sets, for fig1ext_scrna_correlation
# scRNA references and VisiumHD objects are external, not listed in input:.

rule scrna_correlation:
    input:
        script     = "03_benchmarking/R/scrna_correlation.R",   # edits to the script trigger a rerun
        gene_lists = "results/03_benchmarking/dataset_summary/gene_lists.rds",
        binned     = _binning_inputs(
            platforms = {
                "merscope": config["spatial_analysis"]["merscope_samples"],
                "xenium":   config["spatial_analysis"]["xenium_default_samples"],
            },
            resolutions = [config["bin_resolutions"][0]],
            filtered    = True,
        )
    output:
        avg_expr    = "results/03_benchmarking/scrna_correlation/avg_expr.rds",
        corr_by_ref = "results/03_benchmarking/scrna_correlation/correlation_by_reference.rds"
    log:
        "logs/03_benchmarking/scrna_correlation.log"
    benchmark:
        "benchmarks/03_benchmarking/scrna_correlation.txt"
    params:
        out_dir = "results/03_benchmarking/scrna_correlation"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 400000,
        cpus_per_task = 16,
        runtime       = 480,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose {input.script} \
            --config     config/config.yaml  \
            --gene_lists {input.gene_lists}  \
            --out_dir    {params.out_dir}    \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: qc_metrics
# ---------------------------------------------------------------------------
# Per-bin nCount / nFeature (all genes and common genes) for fig2_qc.R, plus
# FLEX / scRNA-seq / VisiumHD intersect-gene metadata. VisiumHD objects are
# external (config["visiumhd"]["data_dir"]).

rule qc_metrics:
    input:
        gene_lists = "results/03_benchmarking/dataset_summary/gene_lists.rds",
        binning    = _binning_inputs(
            platforms = {
                "merscope": config["spatial_analysis"]["merscope_samples"],
                "xenium":   config["spatial_analysis"]["xenium_default_samples"],
            },
            resolutions = config["bin_resolutions"],
            filtered    = True,
        ),
        flex_rds   = config["scrna"]["path"],
        sc_rds     = config["scrna"]["sc_path"]
    output:
        metadata        = "results/03_benchmarking/qc_metrics/metadata_combined.rds",
        meta_flex       = "results/03_benchmarking/qc_metrics/metadata_flex_scrna.rds",
        genes_intersect = "results/03_benchmarking/qc_metrics/genes_intersect_flex.rds"
    log:
        "logs/03_benchmarking/qc_metrics.log"
    benchmark:
        "benchmarks/03_benchmarking/qc_metrics.txt"
    params:
        out_dir = "results/03_benchmarking/qc_metrics"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 200000,
        cpus_per_task = 16,
        runtime       = 360,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/qc_metrics.R \
            --config     config/config.yaml \
            --out_dir    {params.out_dir}   \
            --gene_lists {input.gene_lists} \
            --sc_rds     {input.sc_rds}     \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: qc_metrics_roi
# ---------------------------------------------------------------------------
# ROI-restricted qc_metrics on the fixed-size 8µm ROIs of the matched animals
# (config["stalign"]), for the Figure 2 ROI QC boxplots.

rule qc_metrics_roi:
    input:
        gene_lists = "results/03_benchmarking/dataset_summary/gene_lists.rds",
        visium     = expand(
            "results/01_preprocessing/visium_8um_roi/{sample}_8um_roi_" + ROI_SIZE_LABEL + ".rds",
            sample = config["stalign"]["visiumhd"].values(),
        ),
        merscope   = expand(
            "results/01_preprocessing/merscope_8um_roi/{sample}_8um_roi_" + ROI_SIZE_LABEL + ".rds",
            sample = config["stalign"]["matched_samples"]["merscope"].values(),
        ),
        xenium     = expand(
            "results/01_preprocessing/xenium_8um_roi/{sample}_8um_roi_" + ROI_SIZE_LABEL + ".rds",
            sample = config["stalign"]["matched_samples"]["xenium"].values(),
        )
    output:
        metadata = "results/03_benchmarking/qc_metrics/metadata_roi" + ROI_SIZE_LABEL + ".rds"
    log:
        "logs/03_benchmarking/qc_metrics_roi.log"
    benchmark:
        "benchmarks/03_benchmarking/qc_metrics_roi.txt"
    params:
        out_dir   = "results/03_benchmarking/qc_metrics",
        roi_label = ROI_SIZE_LABEL
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 32000,
        cpus_per_task = 2,
        runtime       = 60,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/qc_metrics_roi.R \
            --config     config/config.yaml \
            --out_dir    {params.out_dir}   \
            --gene_lists {input.gene_lists} \
            --roi_label  {params.roi_label} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: qc_spatial
# ---------------------------------------------------------------------------
# Per-bin coordinates and nCount / nFeature (all genes, 8µm) for the
# fig2ext_qc_spatial maps; same bins as qc_metrics (checked). Matched animals
# (config["stalign"]) use aligned coordinates. Map parameters are in
# config/qc_spatial.yaml.

rule qc_spatial:
    input:
        script     = "03_benchmarking/R/qc_spatial.R",
        helpers    = ["03_benchmarking/R/utils/qc_utils.R", "01_preprocessing/R/roi_utils.R"],
        config     = "config/config.yaml",
        qc_spatial = "config/qc_spatial.yaml",
        filtered   = _binning_inputs(
            platforms = {
                "merscope": config["spatial_analysis"]["merscope_samples"],
                "xenium":   config["spatial_analysis"]["xenium_default_samples"],
            },
            resolutions = [8],
            filtered    = True,
        ),
        aligned    = expand(
            "results/01_preprocessing/{platform}_8um_aligned/{sample}_8um_aligned.rds",
            zip,
            platform = ["merscope"] * len(config["stalign"]["matched_samples"]["merscope"])
                     + ["xenium"]   * len(config["stalign"]["matched_samples"]["xenium"]),
            sample   = list(config["stalign"]["matched_samples"]["merscope"].values())
                     + list(config["stalign"]["matched_samples"]["xenium"].values()),
        ),
        visium     = expand(
            config["visiumhd"]["data_dir"] + "/{file}",
            file = config["visiumhd"]["samples"].values(),
        ),
        metadata   = "results/03_benchmarking/qc_metrics/metadata_combined.rds"
    output:
        bins = "results/03_benchmarking/qc_spatial/bins_8um.rds"
    log:
        "logs/03_benchmarking/qc_spatial.log"
    benchmark:
        "benchmarks/03_benchmarking/qc_spatial.txt"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb          = 24000,
        cpus_per_task   = 2,
        runtime         = 60,
        slurm_partition = "regular"
    shell:
        """
        Rscript --vanilla --verbose {input.script} \
            --config     {input.config}     \
            --qc_spatial {input.qc_spatial} \
            --metadata   {input.metadata}   \
            --bin_size   8                  \
            --out_rds    {output.bins}      \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: moransi
# ---------------------------------------------------------------------------
# Moran's I of background (MERSCOPE Blanks / Xenium unassigned codewords) vs
# target features, filtered 8µm. Samples exclude
# config["qc_backgrounds"]["exclude_samples"] (MERSCOPE Batch22 samples have
# no blank-probe signal).
# R 4.5.1 because SEraster 0.99.5 requires R >= 4.5; sf and magick need the
# geos/proj/gdal and ImageMagick modules.

rule moransi:
    wildcard_constraints:
        platform = "merscope|xenium"
    input:
        lambda wc: _binning_inputs(
            platforms   = {wc.platform: _qc_bg_samples(wc.platform)},
            resolutions = [config["bin_resolutions"][0]],
            filtered    = True,
        )
    output:
        rds = "results/03_benchmarking/moransi/{platform}_moransi.rds"
    log:
        "logs/03_benchmarking/moransi/{platform}.log"
    benchmark:
        "benchmarks/03_benchmarking/moransi/{platform}.txt"
    params:
        resolution = 100
    envmodules:
        "R/4.5.1",
        "geos/3.12.1",
        "proj/9.4.0",
        "gdal/3.9.0",
        "ImageMagick/7.1.2-18"
    resources:
        mem_mb        = 64000,
        cpus_per_task = 2,
        runtime       = 240,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/moransi.R \
            --config     config/config.yaml \
            --platform   {wildcards.platform} \
            --resolution {params.resolution} \
            --out_rds    {output.rds} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: qc_backgrounds
# ---------------------------------------------------------------------------
# Background vs target count summaries and FDR, plus combined Moran's I, for
# fig2_background.R.

rule qc_backgrounds:
    input:
        rds = _binning_inputs(
            platforms = {
                "merscope": _qc_bg_samples("merscope"),
                "xenium":   _qc_bg_samples("xenium"),
            },
            resolutions = [config["bin_resolutions"][0]],
            filtered    = True,
        ),
        moransi_mer = "results/03_benchmarking/moransi/merscope_moransi.rds",
        moransi_xen = "results/03_benchmarking/moransi/xenium_moransi.rds"
    output:
        bg_per_sample = "results/03_benchmarking/qc_backgrounds/background_per_sample.rds",
        bg_summary    = "results/03_benchmarking/qc_backgrounds/background_summary.rds",
        fdr           = "results/03_benchmarking/qc_backgrounds/fdr_results.rds",
        moransi       = "results/03_benchmarking/qc_backgrounds/moransi_combined.rds"
    log:
        "logs/03_benchmarking/qc_backgrounds.log"
    benchmark:
        "benchmarks/03_benchmarking/qc_backgrounds.txt"
    params:
        out_dir = "results/03_benchmarking/qc_backgrounds"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 100000,
        cpus_per_task = 8,
        runtime       = 180,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/qc_backgrounds.R \
            --config      config/config.yaml   \
            --out_dir     {params.out_dir}     \
            --moransi_mer {input.moransi_mer}  \
            --moransi_xen {input.moransi_xen}  \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: probe_rank
# ---------------------------------------------------------------------------
# Per-probe rank tables, target probes within background levels (is_overlap),
# and the S-curve label pool for fig2_background.R. Same sample exclusions as
# qc_backgrounds / moransi.

rule probe_rank:
    input:
        _binning_inputs(
            platforms = {
                "merscope": _qc_bg_samples("merscope"),
                "xenium":   _qc_bg_samples("xenium"),
            },
            resolutions = [config["bin_resolutions"][0]],
            filtered    = True,
        )
    output:
        ranked        = "results/03_benchmarking/probe_rank/ranked_plat.rds",
        ranked_sample = "results/03_benchmarking/probe_rank/ranked_sample.rds",
        label_pool    = "results/03_benchmarking/probe_rank/label_pool_genes.csv"
    log:
        "logs/03_benchmarking/probe_rank.log"
    benchmark:
        "benchmarks/03_benchmarking/probe_rank.txt"
    params:
        out_dir      = "results/03_benchmarking/probe_rank"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 100000,
        cpus_per_task = 8,
        runtime       = 180,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/probe_rank.R \
            --config       config/config.yaml      \
            --out_dir      {params.out_dir}        \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: gene_comparison
# ---------------------------------------------------------------------------
# 8µm pseudobulk counts for the four matched animals
# (config["gene_comparison"]["animals"]) on genes present in all 12 datasets,
# for fig2_gene_comparison.R.
# External VisiumHD objects are listed as inputs so a missing file fails at
# DAG build time.

def _gene_comparison_samples(platform):
    """Sample names on one platform for the matched gene_comparison animals."""
    return [a[platform] for a in config["gene_comparison"]["animals"].values()]


def _gene_comparison_visiumhd_inputs():
    """Paths to the VisiumHD 8µm objects for the matched animals."""
    vhd = config["visiumhd"]
    return [f'{vhd["data_dir"].rstrip("/")}/{vhd["samples"][s]}'
            for s in _gene_comparison_samples("visiumhd")]


rule gene_comparison:
    input:
        binned = _binning_inputs(
            platforms = {
                "merscope": _gene_comparison_samples("merscope"),
                "xenium":   _gene_comparison_samples("xenium"),
            },
            resolutions = [config["bin_resolutions"][0]],
            filtered    = True,
        ),
        visiumhd = _gene_comparison_visiumhd_inputs()
    output:
        dge        = "results/03_benchmarking/gene_comparison/dge.rds",
        counts_mat = "results/03_benchmarking/gene_comparison/counts_mat.rds"
    log:
        "logs/03_benchmarking/gene_comparison.log"
    benchmark:
        "benchmarks/03_benchmarking/gene_comparison.txt"
    params:
        out_dir = "results/03_benchmarking/gene_comparison"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 200000,
        cpus_per_task = 16,
        runtime       = 240,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/gene_comparison.R \
            --config  config/config.yaml \
            --out_dir {params.out_dir}   \
            > {log} 2>&1
        """

# ---------------------------------------------------------------------------
# Rule: variance_partition
# ---------------------------------------------------------------------------
# PERMANOVA of the matched pseudobulk profiles (~ Platform + Group, sequential
# SS) for fig2_gene_comparison.R. Platform is permuted within Animal; Group
# (constant within Animal) is not tested.

rule variance_partition:
    input:
        dge = "results/03_benchmarking/gene_comparison/dge.rds"
    output:
        permanova          = "results/03_benchmarking/variance_partition/permanova_results.rds",
        variance_explained = "results/03_benchmarking/variance_partition/variance_explained.rds"
    log:
        "logs/03_benchmarking/variance_partition.log"
    benchmark:
        "benchmarks/03_benchmarking/variance_partition.txt"
    params:
        out_dir = "results/03_benchmarking/variance_partition"
    envmodules:
        "R/4.4.1"
    resources:
        mem_mb        = 8000,
        cpus_per_task = 1,
        runtime       = 30,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/variance_partition.R \
            --dge     {input.dge}      \
            --out_dir {params.out_dir} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: segmentation_quality
# ---------------------------------------------------------------------------
# Morphology, UMAP, composition, MECR and negative-marker purity for the six
# annotated segmentation objects (Xenium batch34 and MERSCOPE x 3 methods).
# Summary-metrics CSVs and scRNA reference are external; paths are set under
# config segmentation_comp (see segmentation_quality.R header).

rule segmentation_quality:
    input:
        expand(
            "results/02_spatial_analysis/{method}/annotated_final.rds",
            method = [
                "xenium_batch34_default", "xenium_batch34_cellpose", "xenium_batch34_proseg",
                "merscope_default",       "merscope_cellpose",       "merscope_proseg",
            ],
        )
    output:
        metrics_long     = "results/03_benchmarking/segmentation_quality/metrics_long.rds",
        umap_coords      = "results/03_benchmarking/segmentation_quality/umap_coords.rds",
        cell_type_counts = "results/03_benchmarking/segmentation_quality/cell_type_counts.rds",
        mecr_table       = "results/03_benchmarking/segmentation_quality/mecr_table.rds",
        purity_summary   = "results/03_benchmarking/segmentation_quality/purity_summary.rds",
        purity_ct        = "results/03_benchmarking/segmentation_quality/purity_ct.rds"
    log:
        "logs/03_benchmarking/segmentation_quality.log"
    benchmark:
        "benchmarks/03_benchmarking/segmentation_quality.txt"
    params:
        out_dir = "results/03_benchmarking/segmentation_quality"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 500000,
        cpus_per_task = 16,
        runtime       = 480,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/segmentation_quality.R \
            --config  config/config.yaml \
            --out_dir {params.out_dir}   \
            > {log} 2>&1
        """
