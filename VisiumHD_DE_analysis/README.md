# Spatial Transcriptomics Analysis of VisiumHD samples

## About the `targets` workflow

This project uses the [`targets`](https://books.ropensci.org/targets/) R package, a pipeline
toolkit for reproducible analysis that:

- Tracks dependencies between analysis steps automatically — only re-runs targets whose
  inputs have changed.
- Caches every result as a `.qs` binary file in `_targets/`, so a partially complete run
  can be resumed without redoing finished steps.
- Runs independent targets in parallel using the `crew` package, dispatching heavy jobs
  to the HPC via SLURM.

### Quickstart

```r
# Run (or resume) the full pipeline
targets::tar_make()

# Inspect which targets are out of date
targets::tar_outdated()

# Load a cached result into the R session
targets::tar_load(target_name)   # assign to its original name
targets::tar_read(target_name)   # return the value

# Visualise the dependency graph
targets::tar_visnetwork()
```

For more details, refer to the [`targets`](https://books.ropensci.org/targets/) documentation.

## Input data

Four cryosectioned mouse spleen samples were captured at 8 µm bin resolution:

| Sample ID       | Genotype       |
|-----------------|----------------|
| batch33_167     | *Tbx21* KO     |
| batch33_168     | *Tbx21* KO     |
| batch33_709     | Wild-type      |
| batch33_713     | Wild-type      |

Pre-processed Seurat objects (SpaceRanger output → Seurat normalisation) are read
directly by the pipeline (`/vast/projects/SpatialBench/VisiumHD/normalized_obj_08016um/`); they are not included in this repository.

## Analysis pipeline (`targets/hd_ko_wt_8um.R`)

### Per-sample preprocessing

Each Seurat object is loaded, the 8 µm assay (`Spatial.008um`) is set as the default, and
relative-counts (RC) normalisation is applied (scale factor 5 000). Spatially variable genes
(SVGs) are identified with Seurat's `FindVariableFeatures`, and the object is converted to a
`SpatialExperiment` for downstream Bioconductor-compatible steps.

### Spatial neighbourhood embedding — Banksy

[Banksy](https://github.com/prabhakarlab/Banksy) augments each bin's expression with a
summary of its spatial neighbours (k_geom = 6 nearest neighbours, λ = 0.2 neighbourhood
weight). PCA and UMAP are computed jointly across all four samples, using group-aware
normalisation to reduce the KO/WT batch effect. The union of all per-sample SVGs is used as
the feature set.

### Data cleaning

- A low-quality tissue edge in sample 709 is removed by a linear spatial cut
  (`y ≥ 2.13x − 16 000` in pixel coordinates).
- Probe entries mapping to the same gene symbol are merged by summing raw counts and
  normalised counts.
- Ensembl gene IDs and chromosome names are appended via `biomaRt` to support downstream
  gene set analyses.

### Cell type deconvolution — RCTD

[RCTD](https://github.com/dmcable/spacexr) is run in doublet mode using a Flex single-cell
RNA-seq annotated reference. It assigns each 8 µm bin a primary cell type.
Because the reference does not distinguish dark zone (DZ) and light zone (LZ) GC B cells,
those are resolved separately:

- Bins classified as GC B cells by RCTD are extracted for KO and WT separately.
- Banksy is re-run on each subset (resolution 0.5) to produce spatial subclusters.
- Subclusters are manually mapped to DZ or LZ based on the expression of
  zone-specific marker genes (`data/DZLZ_signature.csv`).

### Pseudo-bulk differential expression — KO vs WT

Raw counts are aggregated (summed) per cell type × sample to form pseudo-bulk libraries.
For each cell type independently:

1. Low-expressed genes are removed with `edgeR::filterByExpr`.
2. Libraries are TMM-normalised.
3. Precision weights are estimated with `limma::voom`.
4. A means model is fitted and the KO vs WT contrast is tested with `limma::eBayes`.

### Gene set tests — T-bet targets and sex-linked genes

Two `limma::roast` rotation gene set tests are applied to each cell type's voom object:

- **T-bet target genes** — Genes found differentially expressed in *Tbx21* KO GC B cells
  by a previous bulk RNA-seq study (`data/GCB-DE-D15.TbetKO_GCB.v.WT_GCB.csv`), weighted
  by their prior fold-change direction.
- **Sex-linked genes** — Y-chromosome genes (up in males) and X-inactivation-escape genes
  found in mouse spleen tissues reported by
  [Berletch et al. (2015)](https://doi.org/10.1371/journal.pgen.1005079),
  used to confirm expected male-vs-female signal between samples.

### Dark zone vs light zone differential expression

The same voom-limma pipeline is applied within each genotype group (KO, WT), comparing DZ
and LZ pseudo-bulk counts.  Results are validated against a published DZ/LZ bulk RNA-seq
signature ([Gabriel et al.](data/Gabriel_etal_supp.xlsx)) using a ROAST gene set test.

## Visualisation (`targets/hd_ko_wt_plots.R`)

| Plot target | Description |
|---|---|
| QC spatial maps | UMI counts, mitochondrial %, and filter status across all four sections |
| Cluster / zone maps | RCTD cell type labels plotted on tissue; DZ/LZ-resolved GC regions |
| UMAP plots | Banksy joint UMAP and standard expression UMAP, coloured by cell type |
| Tbx21 spatial map | Raw *Tbx21* counts across tissue to verify knockout efficiency |
| MA and barcode plots | T-bet target gene and sex-linked gene enrichment in each cell type |
| Volcano plots | KO vs WT and DZ vs LZ, per cell type |
| GC marker heatmaps | DZ/LZ marker gene expression in Banksy subclusters |
| Cell type composition | Bar charts of RCTD cell type proportions per sample and group |

All plots are saved as PDFs in `output/`.

## Repository layout

```
_targets.R          Pipeline entry point — crew/SLURM controller setup and targets/ discovery
targets/            Target list files (each returns a list() of tar_target() calls)
  hd_ko_wt_8um.R   Main analysis: preprocessing, deconvolution, DE, gene set tests
  hd_ko_wt_plots.R Visualisation targets
R/                  Helper functions (auto-sourced by tar_source())
  pre-processing.R  Sample loading, QC, normalisation, biomaRt utilities
  clustering.R      SNN/Leiden clustering, iSC.MEB, RCTD wrappers
  visualization.R   Spatial plots, heatmaps, volcano plots
  DE.R              Pseudo-bulk aggregation, GO/KEGG enrichment
  DA.R              sPLS-DA discriminant analysis (exploratory)
  convertions.R     Seurat ↔ SpatialExperiment conversion
data/               Sample metadata, probe sets, published gene signatures
output/             Generated plots (PDF) and result tables (XLSX/CSV)
logs/               SLURM worker logs
```
