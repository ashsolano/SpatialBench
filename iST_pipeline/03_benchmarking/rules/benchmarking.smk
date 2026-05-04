# Rules:   dataset_summary, scrna_correlation, qc_metrics, qc_backgrounds,
#          probe_rank, gene_comparison, segmentation_quality
# Purpose: Cross-platform benchmarking metrics and quality checks.
#
# Pipeline DAG:
#   (01_preprocessing binning outputs)
#     -> dataset_summary     (single job — all platforms, all resolutions)
#     -> scrna_correlation   (single job — all platforms, 8µm only)
#     -> qc_metrics          (single job — all platforms, all resolutions; needs gene_lists)
#     -> qc_backgrounds      (single job — MERSCOPE + Xenium, 8µm only)
#     -> probe_rank          (single job — MERSCOPE + Xenium, 8µm only)
#     -> gene_comparison     (single job — all platforms, 8µm only)
#
# Config keys used:
#   config["visiumhd"]["data_dir"]                        — VisiumHD processed objects (external)
#   config["visiumhd"]["samples"]                         — VisiumHD sample mapping
#   config["spatial_analysis"]["merscope_samples"]        — 9 MERSCOPE samples
#   config["spatial_analysis"]["xenium_default_samples"]  — 9 Xenium samples
#   config["bin_resolutions"]                             — [8, 16]
#   config["scrna"]["path"]                               — scFlex_seu.rds (external)
#   config["output_dir"]                                  — results root
#
# Named sub-targets (defined in Snakefile):
#   snakemake benchmarking_dataset_summary
#   snakemake benchmarking_scrna_correlation
#   snakemake benchmarking_qc_metrics
#   snakemake benchmarking_qc_backgrounds
#   snakemake benchmarking_probe_rank
#   snakemake benchmarking_gene_comparison
#   snakemake benchmarking_segmentation_quality
#   snakemake benchmarking_all


# ---------------------------------------------------------------------------
# Helper: collect binning RDS inputs for a dict of {platform: [samples]}
# ---------------------------------------------------------------------------

def _binning_inputs(platforms, resolutions):
    """Expand preprocessing binning paths for each platform/sample/resolution."""
    inputs = []
    for platform, samples in platforms.items():
        inputs += expand(
            "results/01_preprocessing/{platform}_{res}um/{sample}_{res}um.rds",
            platform = platform,
            sample   = samples,
            res      = resolutions,
        )
    return inputs


# ---------------------------------------------------------------------------
# Rule: dataset_summary
# ---------------------------------------------------------------------------
# Loads all MERSCOPE and Xenium binning objects (8µm and 16µm) plus the
# VisiumHD processed objects. Computes per-sample QC metrics (bins,
# transcripts, sparsity, common-gene equivalents across the three platforms)
# and saves a long-format metrics table and gene-panel lists for fig1.R.
#
# VisiumHD objects are loaded from the external path config["visiumhd"]["data_dir"]
# and are not pipeline outputs, so they are not listed in input:.

rule dataset_summary:
    input:
        _binning_inputs(
            platforms = {
                "merscope": config["spatial_analysis"]["merscope_samples"],
                "xenium":   config["spatial_analysis"]["xenium_default_samples"],
            },
            resolutions = config["bin_resolutions"],
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
        partition     = "regular"
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
# Loads the 10X FLEX scRNA-seq reference and the 8µm binning objects for
# MERSCOPE and Xenium, plus VisiumHD processed objects. Pseudobulks WT
# samples, computes sparse log10(CPM+1), and saves per-platform averaged
# expression data frames with Pearson r values for fig1c.
#
# The scRNA-seq reference and VisiumHD objects are external (not pipeline
# outputs) so they are not listed in input:.

rule scrna_correlation:
    input:
        _binning_inputs(
            platforms = {
                "merscope": config["spatial_analysis"]["merscope_samples"],
                "xenium":   config["spatial_analysis"]["xenium_default_samples"],
            },
            resolutions = [config["bin_resolutions"][0]],
        )
    output:
        avg_expr = "results/03_benchmarking/scrna_correlation/avg_expr.rds"
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
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/scrna_correlation.R \
            --config  config/config.yaml \
            --out_dir {params.out_dir}   \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: qc_metrics
# ---------------------------------------------------------------------------
# Computes per-bin nCount / nFeature metadata across all platforms and bin
# resolutions, for both full and common-gene subsets. Saves
# metadata_combined.rds for use by fig2_qc.R.
#
# VisiumHD objects are loaded from config["visiumhd"]["data_dir"] (external).
# The --sc_rds flag for FLEX/scRNA comparison is omitted; re-run manually
# with that flag if metadata_flex_scrna.rds is needed.

rule qc_metrics:
    input:
        gene_lists = "results/03_benchmarking/dataset_summary/gene_lists.rds",
        binning    = _binning_inputs(
            platforms = {
                "merscope": config["spatial_analysis"]["merscope_samples"],
                "xenium":   config["spatial_analysis"]["xenium_default_samples"],
            },
            resolutions = config["bin_resolutions"],
        )
    output:
        metadata = "results/03_benchmarking/qc_metrics/metadata_combined.rds"
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
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/qc_metrics.R \
            --config     config/config.yaml \
            --out_dir    {params.out_dir}   \
            --gene_lists {input.gene_lists} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: qc_backgrounds
# ---------------------------------------------------------------------------
# Computes background vs target count summaries and FDR for MERSCOPE and
# Xenium 8µm binning objects. Saves background_per_sample.rds,
# background_summary.rds, and fdr_results.rds for fig2_background.R.
#
# The --moransi_mer / --moransi_xen flags for pre-computed Moran's I are
# omitted; re-run manually with those flags to generate moransi_combined.rds.

rule qc_backgrounds:
    input:
        _binning_inputs(
            platforms = {
                "merscope": config["spatial_analysis"]["merscope_samples"],
                "xenium":   config["spatial_analysis"]["xenium_default_samples"],
            },
            resolutions = [config["bin_resolutions"][0]],
        )
    output:
        bg_per_sample = "results/03_benchmarking/qc_backgrounds/background_per_sample.rds",
        bg_summary    = "results/03_benchmarking/qc_backgrounds/background_summary.rds",
        fdr           = "results/03_benchmarking/qc_backgrounds/fdr_results.rds"
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
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/qc_backgrounds.R \
            --config  config/config.yaml \
            --out_dir {params.out_dir}   \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: probe_rank
# ---------------------------------------------------------------------------
# Builds per-probe rank tables for MERSCOPE and Xenium 8µm objects, identifies
# target probes overlapping background probe levels (is_overlap), and generates
# the label pool CSV used for S-curve annotations in fig2_background.R.

rule probe_rank:
    input:
        _binning_inputs(
            platforms = {
                "merscope": config["spatial_analysis"]["merscope_samples"],
                "xenium":   config["spatial_analysis"]["xenium_default_samples"],
            },
            resolutions = [config["bin_resolutions"][0]],
        )
    output:
        ranked     = "results/03_benchmarking/probe_rank/ranked_plat.rds",
        label_pool = "results/03_benchmarking/probe_rank/label_pool_genes.csv"
    log:
        "logs/03_benchmarking/probe_rank.log"
    benchmark:
        "benchmarks/03_benchmarking/probe_rank.txt"
    params:
        out_dir = "results/03_benchmarking/probe_rank"
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
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/probe_rank.R \
            --config  config/config.yaml \
            --out_dir {params.out_dir}   \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: gene_comparison
# ---------------------------------------------------------------------------
# Builds pseudobulk count matrices across VisiumHD, MERSCOPE, and Xenium
# using 8µm bins restricted to the common gene set. Saves a DGEList (dge.rds)
# and the raw count matrix (counts_mat.rds) for fig2_gene_comparison.R.
#
# VisiumHD objects are loaded from config["visiumhd"]["data_dir"] (external).

rule gene_comparison:
    input:
        _binning_inputs(
            platforms = {
                "merscope": config["spatial_analysis"]["merscope_samples"],
                "xenium":   config["spatial_analysis"]["xenium_default_samples"],
            },
            resolutions = [config["bin_resolutions"][0]],
        )
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
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/gene_comparison.R \
            --config  config/config.yaml \
            --out_dir {params.out_dir}   \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: segmentation_quality
# ---------------------------------------------------------------------------
# Loads six annotated Seurat objects (xenium_batch34 and merscope × 3 methods)
# and pre-computed summary metrics CSVs (paths in config segmentation_comp).
# Computes:
#   - long-format cell morphology metrics (area, count, transcripts)
#   - UMAP coordinates + cell type labels extracted from each object
#   - cell type composition counts per sample
#   - per-sample MECR (Mixed-cell Expression Co-expression Rate)
#   - per-sample and per-cell-type negative marker purity (requires scRNA ref)
#
# The summary metrics CSVs and scRNA reference are external inputs; their paths
# must be set in config.yaml under the segmentation_comp section.
# See segmentation_quality.R header for required config keys.

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
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 03_benchmarking/R/segmentation_quality.R \
            --config  config/config.yaml \
            --out_dir {params.out_dir}   \
            > {log} 2>&1
        """
