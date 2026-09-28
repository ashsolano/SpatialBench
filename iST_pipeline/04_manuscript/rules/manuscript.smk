# Purpose: Publication-ready manuscript figures (fig1, fig2_qc, fig2_background,
#          fig2_gene_comparison, fig3) from pre-computed benchmarking outputs.
# Inputs:  results/03_benchmarking/{dataset_summary,scrna_correlation}/  -> fig1
#          results/03_benchmarking/qc_metrics/ + 8um aligned objects    -> fig2_qc
#          results/03_benchmarking/{qc_backgrounds,probe_rank}/          -> fig2_background
#          results/03_benchmarking/{gene_comparison,variance_partition}/ -> fig2_gene_comparison
#          results/03_benchmarking/segmentation_quality/                 -> fig3
# Outputs: figures/fig{1,2,3}/*.pdf
# Config keys used:
#   config["stalign"]["matched_samples"] — aligned 8µm montage objects (fig2_qc)
#   config["visiumhd"]["data_dir"]       — VisiumHD objects, read by fig2_qc.R via config.yaml
#   ROI_SIZE_LABEL (Snakefile)           — ROI size suffix for fig2_qc files
# Targets (Snakefile): manuscript_{fig1,fig2_qc,fig2_background,fig2_gene_comparison,fig3,all}

# ---------------------------------------------------------------------------
# Rule: fig1
# ---------------------------------------------------------------------------
# Panels:
#   - Bins/transcripts/sparsity bar chart (8µm vs 16µm)
#   - Gene panel overlap Venn diagram
#   - scRNA-seq vs ST pseudobulk correlation density plots

rule fig1:
    input:
        script     = "04_manuscript/R/fig1.R",
        helpers    = "04_manuscript/R/utils/plot_helpers.R",   # generate_density_plot()
        palettes   = "04_manuscript/R/utils/palettes.R",       # pal_muted_light
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
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose {input.script} \
            --input_rds     {input.metrics}    \
            --gene_lists    {input.gene_lists} \
            --scrna_cor_rds {input.avg_expr}   \
            --out_dir       {params.out_dir}   \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: fig2_qc
# ---------------------------------------------------------------------------
# Panels:
#   - nCount / nFeature spatial montages (VisiumHD objects are read from
#     config["visiumhd"]["data_dir"] and are not listed as inputs)
#   - Median counts/bin and genes/bin per animal (90 common genes), full tissue and ROI
#   - FLEX / scRNA-seq / VisiumHD intersect-gene dot plots (qc_metrics run with --sc_rds)

rule fig2_qc:
    input:
        script   = "04_manuscript/R/fig2_qc.R",
        config   = "config/config.yaml",   # stalign mapping, visiumhd data_dir
        metadata = "results/03_benchmarking/qc_metrics/metadata_combined.rds",
        meta_roi = "results/03_benchmarking/qc_metrics/metadata_roi" + ROI_SIZE_LABEL + ".rds",
        meta_flex       = "results/03_benchmarking/qc_metrics/metadata_flex_scrna.rds",
        genes_intersect = "results/03_benchmarking/qc_metrics/genes_intersect_flex.rds",
        # STalign-aligned 8µm objects for the montage samples
        aligned  = expand(
            "results/01_preprocessing/{platform}_8um_aligned/{sample}_8um_aligned.rds",
            zip,
            platform = ["merscope"] * len(config["stalign"]["matched_samples"]["merscope"])
                     + ["xenium"]   * len(config["stalign"]["matched_samples"]["xenium"]),
            sample   = list(config["stalign"]["matched_samples"]["merscope"].values())
                     + list(config["stalign"]["matched_samples"]["xenium"].values()),
        )
    output:
        spatial_ncount   = "figures/fig2/spatial_ncount.pdf",
        spatial_nfeature = "figures/fig2/spatial_nfeature.pdf",
        qc_counts        = "figures/fig2/qc_counts_spatial.pdf",
        qc_genes         = "figures/fig2/qc_genes_spatial.pdf",
        qc_counts_roi    = "figures/fig2/qc_counts_spatial_roi" + ROI_SIZE_LABEL + ".pdf",
        qc_genes_roi     = "figures/fig2/qc_genes_spatial_roi" + ROI_SIZE_LABEL + ".pdf",
        qc_counts_flex_all      = "figures/fig2/qc_counts_flex_all.pdf",
        qc_genes_flex_all       = "figures/fig2/qc_genes_flex_all.pdf",
        qc_counts_flex_visiumhd = "figures/fig2/qc_counts_flex_visiumhd.pdf",
        qc_genes_flex_visiumhd  = "figures/fig2/qc_genes_flex_visiumhd.pdf"
    log:
        "logs/04_manuscript/fig2_qc.log"
    benchmark:
        "benchmarks/04_manuscript/fig2_qc.txt"
    params:
        input_dir = "results/03_benchmarking/qc_metrics",
        out_dir   = "figures/fig2",
        roi_label = ROI_SIZE_LABEL
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
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose {input.script} \
            --config    {input.config} \
            --input_dir {params.input_dir} \
            --out_dir   {params.out_dir}   \
            --roi_label {params.roi_label} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: fig2_background
# ---------------------------------------------------------------------------
# Panels:
#   - background_vs_target: total counts by assay type
#   - moransi: Moran's I by assay type (filtered 8µm bins)
#   - fdr: FDR bar chart with per-sample points
#   - probe_scurves: probe rank S-curves

rule fig2_background:
    input:
        script        = "04_manuscript/R/fig2_background.R",
        bg_per_sample = "results/03_benchmarking/qc_backgrounds/background_per_sample.rds",
        fdr           = "results/03_benchmarking/qc_backgrounds/fdr_results.rds",
        moransi       = "results/03_benchmarking/qc_backgrounds/moransi_combined.rds",
        ranked_plat   = "results/03_benchmarking/probe_rank/ranked_plat.rds",
        label_pool    = "results/03_benchmarking/probe_rank/label_pool_genes.csv"
    output:
        background = "figures/fig2/background_vs_target.pdf",
        moransi    = "figures/fig2/moransi.pdf",
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
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose {input.script} \
            --qc_backgrounds_dir {params.qc_backgrounds_dir} \
            --probe_rank_dir     {params.probe_rank_dir}      \
            --out_dir            {params.out_dir}             \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: fig2_gene_comparison
# ---------------------------------------------------------------------------
# Panels:
#   - Pseudobulk MDS
#   - Variance explained by Platform / Group / Residual (PERMANOVA R2, rule variance_partition)
#   - Average log10(CPM+1) scatter per platform pair

rule fig2_gene_comparison:
    input:
        script   = "04_manuscript/R/fig2_gene_comparison.R",
        dge      = "results/03_benchmarking/gene_comparison/dge.rds",
        var_part = "results/03_benchmarking/variance_partition/variance_explained.rds"
    output:
        mds             = "figures/fig2/pseudobulk_mds.pdf",
        var_explained   = "figures/fig2/variance_explained.pdf",
        scatter_vs_mer  = "figures/fig2/avgexpr_scatter_visiumhd_merscope.pdf",
        scatter_vs_xen  = "figures/fig2/avgexpr_scatter_visiumhd_xenium.pdf",
        scatter_mer_xen = "figures/fig2/avgexpr_scatter_merscope_xenium.pdf"
    log:
        "logs/04_manuscript/fig2_gene_comparison.log"
    benchmark:
        "benchmarks/04_manuscript/fig2_gene_comparison.txt"
    params:
        input_dir    = "results/03_benchmarking/gene_comparison",
        var_part_dir = "results/03_benchmarking/variance_partition",
        out_dir      = "figures/fig2"
    envmodules:
        "R/4.4.1"
    resources:
        mem_mb        = 32000,
        cpus_per_task = 4,
        runtime       = 30,
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose {input.script} \
            --input_dir    {params.input_dir}    \
            --var_part_dir {params.var_part_dir} \
            --out_dir      {params.out_dir}      \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: fig3
# ---------------------------------------------------------------------------
# Panels: metrics_{merscope,xenium}, umap_{merscope,xenium}, cell_counts, mecr,
#   purity_dotplot, purity_heatmap
# Purity panels are skipped when purity_summary.rds / purity_ct.rds are empty
# (segmentation_quality.R ran without a valid scRNA reference).

rule fig3:
    input:
        script           = "04_manuscript/R/fig3.R",
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
        slurm_partition     = "regular"
    shell:
        """
        Rscript --vanilla --verbose {input.script} \
            --input_dir {params.input_dir} \
            --out_dir   {params.out_dir}   \
            > {log} 2>&1
        """
