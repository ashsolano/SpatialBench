# Rules:   fig1, fig2_qc, fig2_background, fig2_gene_comparison, fig3
# Purpose: Generate publication-ready manuscript figures from pre-computed
#          benchmarking outputs.
#
# Pipeline DAG:
#   results/03_benchmarking/dataset_summary/metrics.rds
#   results/03_benchmarking/dataset_summary/gene_lists.rds   -> fig1
#   results/03_benchmarking/scrna_correlation/avg_expr.rds   /
#
#   results/03_benchmarking/qc_metrics/metadata_combined.rds
#   results/01_preprocessing/merscope_8um/ + xenium_8um/     -> fig2_qc
#
#   results/03_benchmarking/qc_backgrounds/                  -> fig2_background
#   results/03_benchmarking/probe_rank/                      /
#
#   results/03_benchmarking/gene_comparison/dge.rds          -> fig2_gene_comparison
#
# Config keys used:
#   config["visiumhd"]["data_dir"]                       — VisiumHD objects (fig2_qc montage)
#   config["spatial_analysis"]["merscope_samples"]       — 8µm binning inputs (fig2_qc)
#   config["spatial_analysis"]["xenium_default_samples"] — 8µm binning inputs (fig2_qc)
#   config["bin_resolutions"]                            — first element = 8µm
#   config["output_dir"]                                 — results root
#
# Named sub-targets (defined in Snakefile):
#   snakemake manuscript_fig1
#   snakemake manuscript_fig2_qc
#   snakemake manuscript_fig2_background
#   snakemake manuscript_fig2_gene_comparison
#   snakemake manuscript_fig3
#   snakemake manuscript_all


# ---------------------------------------------------------------------------
# Rule: fig1
# ---------------------------------------------------------------------------
# Reads pre-computed benchmarking outputs and produces three panels:
#   Panel A — horizontal bar chart of bins/transcripts/sparsity (8µm vs 16µm)
#   Panel B — Venn diagram of gene panel overlap
#   Panel C — scRNA-seq vs ST pseudobulk correlation density plots

rule fig1:
    input:
        metrics    = "results/03_benchmarking/dataset_summary/metrics.rds",
        gene_lists = "results/03_benchmarking/dataset_summary/gene_lists.rds",
        avg_expr   = "results/03_benchmarking/scrna_correlation/avg_expr.rds"
    output:
        barplot = "figures/fig1/fig1_barplot.pdf",
        venn    = "figures/fig1/fig1_venn.pdf",
        scrna   = "figures/fig1/fig1_scrna_correlation.pdf"
    log:
        "logs/04_manuscript/fig1.log"
    benchmark:
        "benchmarks/04_manuscript/fig1.txt"
    params:
        out_dir = "figures/fig1"
    envmodules:
        "R/4.4.1"
    resources:
        mem_mb        = 32000,
        cpus_per_task = 4,
        runtime       = 60,
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 04_manuscript/R/fig1.R \
            --input_rds     {input.metrics}    \
            --gene_lists    {input.gene_lists} \
            --scrna_cor_rds {input.avg_expr}   \
            --out_dir       {params.out_dir}   \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: fig2_qc
# ---------------------------------------------------------------------------
# Generates platform QC panels for Figure 2:
#   - Spatial scatter montages of nCount and nFeature across 4 representative
#     samples per platform (loads 8µm Seurat objects for MERSCOPE and Xenium;
#     VisiumHD loaded from config["visiumhd"]["data_dir"], not listed in input)
#   - Per-platform boxplots of median counts/bin and genes/bin (90 common genes)
# Optional FLEX/scRNA comparison panels are skipped if metadata_flex_scrna.rds
# is absent; re-run qc_metrics.R with --sc_rds to enable them.

rule fig2_qc:
    input:
        metadata = "results/03_benchmarking/qc_metrics/metadata_combined.rds",
        binning  = _binning_inputs(
            platforms = {
                "merscope": config["spatial_analysis"]["merscope_samples"],
                "xenium":   config["spatial_analysis"]["xenium_default_samples"],
            },
            resolutions = [config["bin_resolutions"][0]],
        )
    output:
        spatial_ncount   = "figures/fig2/spatial_ncount.pdf",
        spatial_nfeature = "figures/fig2/spatial_nfeature.pdf",
        qc_counts        = "figures/fig2/qc_counts_spatial.pdf",
        qc_genes         = "figures/fig2/qc_genes_spatial.pdf"
    log:
        "logs/04_manuscript/fig2_qc.log"
    benchmark:
        "benchmarks/04_manuscript/fig2_qc.txt"
    params:
        input_dir = "results/03_benchmarking/qc_metrics",
        out_dir   = "figures/fig2"
    envmodules:
        "R/4.4.1",
        "geos/3.12.1",
        "hdf5/1.12.3",
        "proj/9.4.0",
        "gdal/3.9.0"
    resources:
        mem_mb        = 100000,
        cpus_per_task = 8,
        runtime       = 60,
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 04_manuscript/R/fig2_qc.R \
            --config    config/config.yaml \
            --input_dir {params.input_dir} \
            --out_dir   {params.out_dir}   \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: fig2_background
# ---------------------------------------------------------------------------
# Generates background QC panels for Figure 2:
#   Panel 1 (background_vs_target) — total counts by assay type, boxplot + jitter
#   Panel 2 (moransi)              — Moran's I by assay type; skipped if
#                                    moransi_combined.rds is absent
#   Panel 3 (fdr)                  — FDR bar chart with per-sample points
#   Panel 4 (probe_scurves)        — probe rank S-curves with overlap highlights

rule fig2_background:
    input:
        bg_per_sample = "results/03_benchmarking/qc_backgrounds/background_per_sample.rds",
        fdr           = "results/03_benchmarking/qc_backgrounds/fdr_results.rds",
        ranked_plat   = "results/03_benchmarking/probe_rank/ranked_plat.rds",
        label_pool    = "results/03_benchmarking/probe_rank/label_pool_genes.csv"
    output:
        background = "figures/fig2/background_vs_target.pdf",
        fdr        = "figures/fig2/fdr.pdf",
        scurves    = "figures/fig2/probe_scurves.pdf"
    log:
        "logs/04_manuscript/fig2_background.log"
    benchmark:
        "benchmarks/04_manuscript/fig2_background.txt"
    params:
        qc_backgrounds_dir = "results/03_benchmarking/qc_backgrounds",
        probe_rank_dir     = "results/03_benchmarking/probe_rank",
        out_dir            = "figures/fig2"
    envmodules:
        "R/4.4.1"
    resources:
        mem_mb        = 32000,
        cpus_per_task = 4,
        runtime       = 30,
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 04_manuscript/R/fig2_background.R \
            --qc_backgrounds_dir {params.qc_backgrounds_dir} \
            --probe_rank_dir     {params.probe_rank_dir}      \
            --out_dir            {params.out_dir}             \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: fig2_gene_comparison
# ---------------------------------------------------------------------------
# Generates gene-comparison panels for Figure 2:
#   - Pseudobulk MDS plot coloured by platform, shaped by sample type
#   - Average log10(CPM+1) scatter plots for each platform pair, with density
#     contours, Pearson R annotation, and top-5 divergent gene labels

rule fig2_gene_comparison:
    input:
        dge = "results/03_benchmarking/gene_comparison/dge.rds"
    output:
        mds             = "figures/fig2/pseudobulk_mds.pdf",
        scatter_vs_mer  = "figures/fig2/avgexpr_scatter_visiumhd_merscope.pdf",
        scatter_vs_xen  = "figures/fig2/avgexpr_scatter_visiumhd_xenium.pdf",
        scatter_mer_xen = "figures/fig2/avgexpr_scatter_merscope_xenium.pdf"
    log:
        "logs/04_manuscript/fig2_gene_comparison.log"
    benchmark:
        "benchmarks/04_manuscript/fig2_gene_comparison.txt"
    params:
        input_dir = "results/03_benchmarking/gene_comparison",
        out_dir   = "figures/fig2"
    envmodules:
        "R/4.4.1"
    resources:
        mem_mb        = 32000,
        cpus_per_task = 4,
        runtime       = 30,
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 04_manuscript/R/fig2_gene_comparison.R \
            --input_dir {params.input_dir} \
            --out_dir   {params.out_dir}   \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: fig3
# ---------------------------------------------------------------------------
# Reads all six outputs from segmentation_quality.R and produces eight panels:
#   metrics_merscope / metrics_xenium — boxplots of 4 QC metrics per platform
#   umap_merscope    / umap_xenium    — UMAP triptychs coloured by cell type
#   cell_counts                       — bar chart of mean cells per sample ± SEM
#   mecr                              — MECR boxplot faceted by platform
#   purity_dotplot                    — tx_assigned vs negative-marker purity
#   purity_heatmap                    — ComplexHeatmap of per-cell-type purity
#
# purity panels are skipped gracefully when purity_summary.rds or purity_ct.rds
# are empty (i.e. segmentation_quality.R ran without a valid scRNA reference).

rule fig3:
    input:
        metrics_long     = "results/03_benchmarking/segmentation_quality/metrics_long.rds",
        umap_coords      = "results/03_benchmarking/segmentation_quality/umap_coords.rds",
        cell_type_counts = "results/03_benchmarking/segmentation_quality/cell_type_counts.rds",
        mecr_table       = "results/03_benchmarking/segmentation_quality/mecr_table.rds",
        purity_summary   = "results/03_benchmarking/segmentation_quality/purity_summary.rds",
        purity_ct        = "results/03_benchmarking/segmentation_quality/purity_ct.rds"
    output:
        metrics_merscope = "figures/fig3/metrics_merscope.pdf",
        metrics_xenium   = "figures/fig3/metrics_xenium.pdf",
        umap_merscope    = "figures/fig3/umap_merscope.pdf",
        umap_xenium      = "figures/fig3/umap_xenium.pdf",
        cell_counts      = "figures/fig3/cell_counts.pdf",
        mecr             = "figures/fig3/mecr.pdf",
        purity_dotplot   = "figures/fig3/purity_dotplot.pdf",
        purity_heatmap   = "figures/fig3/purity_heatmap.pdf"
    log:
        "logs/04_manuscript/fig3.log"
    benchmark:
        "benchmarks/04_manuscript/fig3.txt"
    params:
        input_dir = "results/03_benchmarking/segmentation_quality",
        out_dir   = "figures/fig3"
    envmodules:
        "R/4.4.1"
    resources:
        mem_mb        = 32000,
        cpus_per_task = 4,
        runtime       = 30,
        partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose 04_manuscript/R/fig3.R \
            --input_dir {params.input_dir} \
            --out_dir   {params.out_dir}   \
            > {log} 2>&1
        """
