<p align="center">
  <img src="img/logo.png" width="450" alt="SpatialBench logo">
</p>

# SpatialBench : Comparative cross-platform benchmarking of high-resolution spatial transcriptomics using matched mouse lymphoid tissue

Spatial transcriptomics (ST) has rapidly expanded with the introduction of multiple high-resolution platforms, yet cross-platform benchmarking remains limited and largely focused on technical performance. SpatialBench is a matched multi-platform resource comprising Visium HD, Xenium and MERSCOPE data together with single-cell and single-nucleus references from a malaria-challenged wild-type and B cell-specific *Tbx21* knockout mouse spleen model. We systematically evaluate ST platform performance using technical and biological readouts, providing a biologically defined reference dataset for evaluation of ST technologies, method development, and computational benchmarking.

---

## Repository structure

| Folder | Description |
|--------|-------------|
| [`VisiumHD/`](VisiumHD/) | `targets`-based R pipeline for Visium HD analysis |
| [`iST_pipeline/`](iST_pipeline/) | Snakemake pipeline for Xenium and MERSCOPE analysis |

---

## Data availability

Raw data are deposited at BioStudies
([S-BSST2361](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/)).

### Imaging-based spatial data 


| Batch | Technology | Condition | Sample ID |
|-------|------------|-----------|-----------|
| batch 1 | 10X Xenium | Wild Type | [WT713](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/Xenium/WT713) |
| batch 2 | 10X Xenium | Wild Type | [WT709-rep1](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/Xenium/WT709-rep1) |
| batch 2 | 10X Xenium | Wild Type | [WT710-rep1](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/Xenium/WT710-rep1) |
| batch 3 | 10X Xenium | Wild Type | [WT709-rep2](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/Xenium/WT709-rep2) |
| batch 3 | 10X Xenium | Wild Type | [WT710-rep2](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/Xenium/WT710-rep2) |
| batch 1 | 10X Xenium | Control | [Ctrl173-rep1](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/Xenium/Ctrl173-rep1) |
| batch 2 | 10X Xenium | Control | [Ctrl172-rep1](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/Xenium/Ctrl172-rep1) |
| batch 2 | 10X Xenium | Control | [Ctrl174-rep1](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/Xenium/Ctrl174-rep1) |
| batch 3 | 10X Xenium | Control | [Ctrl172-rep2](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/Xenium/Ctrl172-rep2) |
| batch 3 | 10X Xenium | Control | [Ctrl173-rep2](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/Xenium/Ctrl173-rep2) |
| batch 3 | 10X Xenium | Control | [Ctrl174-rep2](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/Xenium/Ctrl174-rep2) |
| batch 1 | 10X Xenium | Knock Out | [KO166](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/Xenium/KO166) |
| batch 1 | 10X Xenium | Knock Out | [KO167](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/Xenium/KO167) |
| batch 2 | 10X Xenium | Knock Out | [KO168](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/Xenium/KO168) |


| Batch | Technology | Condition | Sample ID |
|-------|------------|-----------|-----------|
| batch 13 | MERSCOPE | Wild Type | [WT709](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/MERSCOPE/WT709/WT709) |
| batch 13 | MERSCOPE | Wild Type | [WT710](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/MERSCOPE/WT710) |
| batch 13 | MERSCOPE | Wild Type | [WT713](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/MERSCOPE/WT713) |
| batch 6 | MERSCOPE | Control | [Ctrl174a](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/MERSCOPE/Ctrl174a) |
| batch 9 | MERSCOPE | Control | [Ctrl173](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/MERSCOPE/Ctrl173) |
| batch 10 | MERSCOPE | Control | [Ctrl172](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/MERSCOPE/Ctrl172) |
| batch 6 | MERSCOPE | Knock Out | [KO166b](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/MERSCOPE/KO166b) |
| batch 9 | MERSCOPE | Knock Out | [KO168](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/MERSCOPE/KO168) |
| batch 10 | MERSCOPE | Knock Out | [KO167](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/MERSCOPE/KO167) |

### Sequencing data

| Batch | Technology | Condition | Sample ID |
|-------|------------|-----------|-----------|
| batch 1 | 10X Visium HD | Wild Type | [WT709](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/VisiumHD/WT709.tar.gz) |
| batch 1 | 10X Visium HD | Wild Type | [WT713](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/VisiumHD/WT713.tar.gz) |
| batch 1 | 10X Visium HD | Knock Out | [KO167](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/VisiumHD/KO167.tar.gz) |
| batch 1 | 10X Visium HD | Knock Out | [KO168](https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BSST/361/S-BSST2361/Files/VisiumHD/KO168.tar.gz) |
| batch 29 | 10X Chromium | Wild Type | [WT709](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8047887) |
| batch 29 | 10X Chromium | Wild Type | [WT713](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8047887) |

---

## Acknowledgements

We thank the WEHI Advanced Genomics Facility, Center for Dynamic Imaging and 
Advanced Histotechnology facility for supporting the generation and analysis of 
single-cell and spatial transcriptomics data for this project.

### Authors

Ashleigh Solano, Raymond K. H. Yip, Changqing Wang, Daniela Amann-Zalcenstein,
Pradeep Rajasekhar, Ishrat Zaman, Allan Motyer, Marek Cmero, Yang Xu, Yining Pan,
Casey J. A. Anttila, Stephanie I. Studniberg, Peter F. Hickey, Layla Wang,
Callum J. Sargeant, Ling Ling, Yunshun Chen, Ruvimbo D. Mishi, Lisa J. Ioannidis,
Kim L. Good-Jacobson, Hamish W. King, Kelly L. Rogers, Diana S. Hansen,
Rory Bowden and Matthew E. Ritchie