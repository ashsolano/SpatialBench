# SpatialBench - Imaging Spatial Transcriptomics Pipeline

## About the Snakemake workflow

This project uses [Snakemake](https://snakemake.readthedocs.io/) to manage a
multi-stage spatial transcriptomics analysis pipeline that tracks dependencies
between rules automatically, caches every result as an `.rds` binary under
`results/`, and dispatches heavy jobs to the HPC via SLURM.  Only rules whose
inputs have changed are re-run, so a partially complete run can always be
resumed.

All parameters (sample paths, QC thresholds, BANKSY settings, DE contrasts)
are controlled through `config/config.yaml`.

### Quickstart

```bash
# Dry-run to see what would be executed
snakemake -n

# Run (or resume) the full pipeline
snakemake --profile profiles/slurm

# Run a single preprocessing target
snakemake xenium_default --profile profiles/slurm
snakemake merscope_cellpose --profile profiles/slurm

# Run all seven spatial analysis methods
snakemake spatial_analysis --profile profiles/slurm

# Run annotation for one method (after filling in config["annotations"])
snakemake spatial_annotate_merscope_default --profile profiles/slurm

# Override bin resolution
snakemake xenium_binning --config bin_resolutions=[8] --profile profiles/slurm
```

---

## Platforms and segmentation methods

| Platform  | Segmentation methods                                                |
|-----------|---------------------------------------------------------------------|
| Xenium    | Default (XeniumRanger v2), Cellpose, Proseg, Binning (8 µm, 16 µm) |
| MERSCOPE  | Default (Vizgen), Cellpose, Proseg, Binning (8 µm, 16 µm)          |

The Xenium cohort is split into two batches. **Xenium default** (9 samples,
batch24 + batch27) is used for expression analysis and pseudobulk DE.
**Xenium batch34** (5 samples) is used for cross-segmentation benchmarking
only — it contains no KO samples, so DE is not run.

---

## Input data

Nine cryosectioned mouse spleen samples captured on both platforms:

| Sample ID | Genotype    |
|-----------|-------------|
| wt709     | Wild-type   |
| wt710     | Wild-type   |
| wt713     | Wild-type   |
| ctrl172   | Control     |
| ctrl173   | Control     |
| ctrl174   | Control     |
| ko166     | *Tbx21* KO  |
| ko167     | *Tbx21* KO  |
| ko168     | *Tbx21* KO  |

Raw data live under `/vast/projects/SpatialBench/data/` and are not tracked by
Git. Paths to each sample directory are defined per-method in
`config/config.yaml`.

---

## Analysis pipeline

### 01 — Preprocessing

One Seurat object is built per sample per method from raw platform output.
Named targets: `xenium_default`, `xenium_cellpose`, `xenium_proseg`,
`xenium_binning`, `merscope_default`, `merscope_cellpose`, `merscope_proseg`,
`merscope_binning`.

### 02 — Spatial analysis

Seven method-cohorts (e.g. `xenium_default`, `merscope_cellpose`) are
processed through the same analysis chain. Per-sample objects are QC-filtered
and LogNormalised, then merged and jointly embedded using
[BANKSY](https://github.com/prabhakarlab/Banksy) (spatial neighbourhood
augmentation) and Harmony batch correction. Louvain clusters are annotated
with cell type labels defined in `config["annotations"]`. GC B cells are
re-embedded at higher resolution to resolve Dark Zone and Light Zone
subpopulations, with zone assignments stored in `config["gc_zones"]`. Finally,
limma-voom pseudobulk DE (KO vs WT) is run per cell type on the four methods
that include KO samples.

Named targets follow the pattern:
`spatial_analysis`, `spatial_annotate_all`, `spatial_subcluster_gc_all`,
`spatial_merge_gc_zones_all`, `spatial_de_all` — or per-method variants such
as `spatial_de_merscope_default`.

### 03 — Benchmarking

Evaluates and compares segmentation methods across metrics (cell count, size
distribution, transcript capture, spatial coherence). In development.

### 04 — Manuscript figures

Generates publication-ready figures and tables for the methods manuscript.
In development.

---

## Repository layout

```
Snakefile                        Pipeline entry point — rule all and named sub-targets
config/
  config.yaml                    Sample paths, QC thresholds, BANKSY/Harmony/DE parameters
profiles/
  slurm/                         Snakemake SLURM profile
01_preprocessing/
  R/                             Per-platform and per-method Seurat object creation scripts
  rules/
    preprocessing.smk            Snakemake rules for all preprocessing targets
02_spatial_analysis/
  R/
    preprocess_sample.R          QC, normalisation, SVG identification
    merge_samples.R              Combine per-sample objects into one cohort object
    embed_harmony.R              BANKSY embedding + Harmony batch correction
    annotate_clusters.R          Assign cell type labels from config["annotations"]
    subcluster_gc.R              GC B cell re-embedding and subclustering
    merge_gc_zones.R             Inject DZ/LZ zone labels into the full object
    pseudobulk_de.R              limma-voom pseudobulk KO vs WT DE per cell type
  rules/
    spatial_analysis.smk         Snakemake rules for all spatial analysis targets
03_benchmarking/
  R/                             (in development)
  rules/                         (in development)
04_manuscript/
  R/                             (in development)
  rules/                         (in development)
utils/
  binning_utils.R                Shared helpers for binning workflows
  segmentation_utils.R           Shared helpers for segmentation workflows
data/
  raw/                           Links to /vast/projects/SpatialBench/data (not in Git)
results/                         Pipeline outputs, one subdirectory per stage (not in Git)
benchmarks/                      Snakemake benchmark TSVs (wall time, peak memory)
logs/                            Rule-level and SLURM worker logs
figures/                         Publication-ready figures (not in Git)
```
