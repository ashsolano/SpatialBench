# Rules:   fig1ext_scrna_correlation, fig2ext_background_persample,
#          fig2ext_probe_rank_persample, fig2ext_qc_violins
# Purpose: Extended manuscript figures from pre-computed benchmarking outputs
#          (plotting only). Extended Data numbers are mapped in CLAUDE.md.
# Inputs:  results/03_benchmarking/scrna_correlation/correlation_by_reference.rds
#          results/03_benchmarking/qc_backgrounds/{background_per_sample,moransi_combined}.rds
#          results/03_benchmarking/probe_rank/{ranked_sample.rds,label_pool_genes.csv}
#          results/03_benchmarking/qc_metrics/metadata_combined.rds
# Config keys used (fig2ext rules; fig1ext WT animals are set upstream in
# scrna_correlation.R from config["scrna"]["wt_ids"]):
#   config["spatial_analysis"]["merscope_samples"]       — animal columns
#   config["spatial_analysis"]["xenium_default_samples"] — animal columns
#   config["qc_backgrounds"]["exclude_samples"]          — expected missing animals
#   config["gene_comparison"]["animals"]                 — VisiumHD animal IDs (qc_violins)
# Targets (Snakefile): extended_fig1ext_scrna_correlation,
#   extended_fig2ext_background, extended_fig2ext_probe_rank,
#   extended_fig2ext_qc_violins, extended_all

# ---------------------------------------------------------------------------
# Rule: fig1ext_scrna_correlation
# ---------------------------------------------------------------------------
# Supports: fig1.R — spatial vs single-cell agreement holds for a second
#           reference (3' GEX) and on the gene set common to all three platforms
# Panels (single PDF, 2 x 6 grid):
#   - rows: 10X FLEX, 10X 3' GEX; columns: Visium HD, MERSCOPE, Xenium on
#     pairwise gene sets, then on the 90-gene common set

rule fig1ext_scrna_correlation:
    input:
        script      = "04_manuscript/R/extended/fig1ext_scrna_correlation.R",   # edits to the script trigger a rerun
        helpers     = "04_manuscript/R/utils/plot_helpers.R",   # generate_density_plot()
        palettes    = "04_manuscript/R/utils/palettes.R",       # pal_density, pal_muted
        corr_by_ref = "results/03_benchmarking/scrna_correlation/correlation_by_reference.rds"
    output:
        grid = "extended_figures/fig1ext_scrna_correlation/scrna_correlation_grid.pdf"
    log:
        "logs/04_manuscript/extended/fig1ext_scrna_correlation.log"
    benchmark:
        "benchmarks/04_manuscript/extended/fig1ext_scrna_correlation.txt"
    params:
        out_dir = "extended_figures/fig1ext_scrna_correlation"
    envmodules:
        "R/4.4.1"
    resources:
        mem_mb          = 8000,
        cpus_per_task   = 1,
        runtime         = 15,
        slurm_partition = "regular"
    shell:
        """
        Rscript --vanilla --verbose {input.script} \
            --in_file {input.corr_by_ref} \
            --out_dir {params.out_dir}    \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: fig2ext_background_persample
# ---------------------------------------------------------------------------
# Supports: fig2_background.R — background signal is low and spatially
#           unstructured relative to genes, per sample (not only pooled)
# Panels (one row per platform, MERSCOPE top, one facet per animal):
#   - background_counts_persample: total counts per assay type per sample
#   - moransi_persample: Moran's I of gene vs background features per sample

rule fig2ext_background_persample:
    input:
        script        = "04_manuscript/R/extended/fig2ext_background_persample.R",
        helpers       = "04_manuscript/R/utils/extended_helpers.R",
        config        = "config/config.yaml",   # sample lists, exclude_samples
        bg_per_sample = "results/03_benchmarking/qc_backgrounds/background_per_sample.rds",
        moransi       = "results/03_benchmarking/qc_backgrounds/moransi_combined.rds"
    output:
        counts  = "extended_figures/fig2ext_background/background_counts_persample.pdf",
        moransi = "extended_figures/fig2ext_background/moransi_persample.pdf"
    log:
        "logs/04_manuscript/extended/fig2ext_background_persample.log"
    benchmark:
        "benchmarks/04_manuscript/extended/fig2ext_background_persample.txt"
    params:
        qc_backgrounds_dir = "results/03_benchmarking/qc_backgrounds",
        out_dir            = "extended_figures/fig2ext_background"
    envmodules:
        "R/4.4.1"
    resources:
        mem_mb          = 8000,
        cpus_per_task   = 1,
        runtime         = 15,
        slurm_partition = "regular"
    shell:
        """
        Rscript --vanilla --verbose {input.script} \
            --qc_backgrounds_dir {params.qc_backgrounds_dir} \
            --config             {input.config}               \
            --out_dir            {params.out_dir}             \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: fig2ext_probe_rank_persample
# ---------------------------------------------------------------------------
# Supports: fig2_background.R — low-rank MERSCOPE genes fall within the
#           background distribution in individual samples (not only pooled),
#           and their spatial structure differs between platforms
# Panels:
#   - probe_rank_persample: per-sample probe rank curves vs the sample's
#     background threshold (one row per platform, MERSCOPE top)
#   - moransi_low_rank_genes: Moran's I of the main Fig 2 low-rank genes,
#     MERSCOPE vs Xenium, one facet per gene

rule fig2ext_probe_rank_persample:
    input:
        script        = "04_manuscript/R/extended/fig2ext_probe_rank_persample.R",
        helpers       = "04_manuscript/R/utils/extended_helpers.R",
        config        = "config/config.yaml",   # sample lists, exclude_samples
        ranked_sample = "results/03_benchmarking/probe_rank/ranked_sample.rds",
        label_pool    = "results/03_benchmarking/probe_rank/label_pool_genes.csv",
        moransi       = "results/03_benchmarking/qc_backgrounds/moransi_combined.rds"
    output:
        probe_rank = "extended_figures/fig2ext_probe_rank/probe_rank_persample.pdf",
        moransi    = "extended_figures/fig2ext_probe_rank/moransi_low_rank_genes.pdf"
    log:
        "logs/04_manuscript/extended/fig2ext_probe_rank_persample.log"
    benchmark:
        "benchmarks/04_manuscript/extended/fig2ext_probe_rank_persample.txt"
    params:
        probe_rank_dir     = "results/03_benchmarking/probe_rank",
        qc_backgrounds_dir = "results/03_benchmarking/qc_backgrounds",
        out_dir            = "extended_figures/fig2ext_probe_rank"
    envmodules:
        "R/4.4.1"
    resources:
        mem_mb          = 8000,
        cpus_per_task   = 1,
        runtime         = 15,
        slurm_partition = "regular"
    shell:
        """
        Rscript --vanilla --verbose {input.script} \
            --probe_rank_dir     {params.probe_rank_dir}     \
            --qc_backgrounds_dir {params.qc_backgrounds_dir} \
            --config             {input.config}               \
            --out_dir            {params.out_dir}             \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule: fig2ext_qc_violins
# ---------------------------------------------------------------------------
# Supports: fig2_qc.R — the per-animal median counts/genes per bin summarise
#           per-bin distributions that are consistent across animals
# Panels (one PDF per bin size, 2 x 4 grid):
#   - rows: counts/bin, genes/bin; columns: Visium HD all genes, then Visium HD,
#     MERSCOPE and Xenium on the 90 common genes; one violin per animal

rule fig2ext_qc_violins:
    input:
        script   = "04_manuscript/R/extended/fig2ext_qc_violins.R",
        helpers  = "04_manuscript/R/utils/extended_helpers.R",
        config   = "config/config.yaml",   # sample lists, gene_comparison.animals
        metadata = "results/03_benchmarking/qc_metrics/metadata_combined.rds"
    output:
        violins_8um  = "extended_figures/fig2ext_qc_violins/qc_violins_8um.pdf",
        violins_16um = "extended_figures/fig2ext_qc_violins/qc_violins_16um.pdf"
    log:
        "logs/04_manuscript/extended/fig2ext_qc_violins.log"
    benchmark:
        "benchmarks/04_manuscript/extended/fig2ext_qc_violins.txt"
    params:
        out_dir = "extended_figures/fig2ext_qc_violins"
    envmodules:
        "R/4.4.1"
    resources:
        mem_mb          = 16000,
        cpus_per_task   = 2,
        runtime         = 30,
        slurm_partition = "regular"
    shell:
        """
        Rscript --vanilla --verbose {input.script} \
            --metadata {input.metadata} \
            --config   {input.config}   \
            --out_dir  {params.out_dir} \
            > {log} 2>&1
        """
