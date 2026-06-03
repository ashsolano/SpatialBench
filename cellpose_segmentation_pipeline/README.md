# Cellpose Segmentation Pipeline for SpatialBench Datasets

## Overview

This repository contains the training data, trained Cellpose models, configuration files, and analysis scripts used to develop and evaluate custom Cellpose segmentation models for MERSCOPE and Xenium spatial transcriptomics datasets.

The workflow includes:

- Training data preparation
- Data augmentation
- Cellpose model training and evaluation
- QuPath-based inference
- Transcript repartitioning
- Segmentation benchmarking against vendor default and Proseg outputs

This repository accompanies the Cellpose benchmarking component of the SpatialBench project and focuses specifically on segmentation model development, inference, and benchmarking.

---

## Input Data

Training and test datasets are provided for both technologies.

| Platform | Contents |
|-----------|-----------|
| MERSCOPE | Training and test ROIs with Cellpose-compatible masks |
| Xenium | Training and test ROIs with Cellpose-compatible masks |

Data are organised as:

```text
data/
├── MERSCOPE/
│   ├── train/
│   └── test/
│
└── Xenium/
    ├── train/
    └── test/
```

---

## Configuration Files

Example configuration files used for running the SpatialVPT workflow and Cellpose segmentation are provided in:

```text
configs/
```

These files are included as examples and were used during the development and benchmarking of the segmentation workflows.

The full NextFlow spatialVPT workflow is documented [here](https://github.com/WEHI-SODA-Hub/spatialvpt).

---

## Model Development Workflow

### 1. Data Augmentation

Training images are augmented using:

- Gaussian blur
- Gaussian noise
- Colour jitter

using:

```text
scripts/augment_dataset.py
```

This workflow was used to increase training set diversity and improve model robustness.

---

### 2. Cellpose Training and Evaluation

Models are trained and evaluated using:

```text
scripts/cellpose_training_and_evaluation.py
```

This workflow:

- Loads train and test datasets
- Trains Cellpose models
- Evaluates performance on held-out test data
- Saves trained model weights

---

### 3. Training for Xenium

To create a model for Xenium data, the MERSCOPE training images were resized to match the resolution of the Xenium platforms. 

```text
scripts/resize_and_train.py
```

This workflow was used to generate the final Xenium Cellpose model included in this repository.

---

## Trained Models

Final trained Cellpose models are provided in:

```text
models/
```

### cyto2_GaussianBlur_ColorJitter_GaussNoise

Custom Cellpose model trained using MERSCOPE ROIs with augmentation consisting of:

- Gaussian blur
- Gaussian noise
- Colour jitter

### cyto2_SpleenXenium_XeniumROIs_v2

Custom Cellpose model trained for Xenium segmentation using manually annotated Xenium ROIs and additional downsampled MERSCOPE training data.

---

## Segmentation Workflows

### MERSCOPE

#### Training

1. Generate augmented training data:

```bash
python scripts/augment_dataset.py
```

2. Train and evaluate the Cellpose model:

```bash
python scripts/cellpose_training_and_evaluation.py
```

3. Save trained model weights to:

```text
models/
```

#### Inference

1. Run Cellpose segmentation using the trained model.
2. Repartition transcripts to segmented cells using the SpatialVPT workflow.
3. Generate segmentation summary metrics.

```bash
python scripts/summary_metrics_merscope_filterBlanks.py
```

---

### Xenium

#### Training

1. Prepare Xenium training ROIs.
2. Incorporate downsampled MERSCOPE ROIs:

```bash
python scripts/resize_and_train.py
```

3. Train and evaluate the Cellpose model:

```bash
python scripts/cellpose_training_and_evaluation.py
```

#### Inference

1. Load Xenium images into QuPath.
2. Run as a batch script:

```text
scripts/run_cellpose_spleen.groovy
```

3. Repartition transcripts to the exported objects (resegmented cells) using [XeniumRanger - Import Segmentation](https://www.10xgenomics.com/support/software/xenium-ranger/latest/analysis/running-pipelines/XR-import-segmentation)
4. Generate segmentation summary metrics:

```bash
python scripts/summary_metrics_xenium_filterBlanks.py
```

---

## Benchmarking

Segmentation methods compared:

| Platform | Methods |
|-----------|-----------|
| MERSCOPE | Vizgen Default, Cellpose, Proseg |
| Xenium | Xenium Default, Cellpose, Proseg |


### Scripts

```text
summary_metrics_merscope_filterBlanks.py
summary_metrics_xenium_filterBlanks.py
summary_metrics_proseg_v2.py
summary_metrics_collate_final.py
```

### Summary Metric Scripts

#### summary_metrics_merscope_filterBlanks.py

Calculates segmentation quality metrics for MERSCOPE datasets.

#### summary_metrics_xenium_filterBlanks.py

Calculates segmentation quality metrics for Xenium datasets.

#### summary_metrics_proseg_v2.py

Calculates equivalent metrics for Proseg outputs to enable comparison with vendor default and Cellpose segmentation.

#### summary_metrics_collate_final.py

Collates all generated metrics into final benchmarking tables used for downstream analysis and figure generation.

---

## Repository Layout

```text
cellpose_segmentation_pipeline/
│
├── configs/
│   ├── Batch26_region2_KO168_NextFlow...
│   └── cellpose2_custom.json
│
├── data/
│   ├── MERSCOPE/
│   │   ├── train/
│   │   └── test/
│   │
│   └── Xenium/
│       ├── train/
│       └── test/
│
├── models/
│   ├── cyto2_GaussianBlur_ColorJitter_GaussNoise
│   └── cyto2_SpleenXenium_XeniumROIs_v2
│
├── scripts/
│   ├── augment_dataset.py
│   ├── cellpose_training_and_evaluation.py
│   ├── resize_and_train.py
│   ├── run_cellpose_spleen.groovy
│   ├── summary_metrics_collate_final.py
│   ├── summary_metrics_merscope_filterBlanks.py
│   ├── summary_metrics_proseg_v2.py
│   └── summary_metrics_xenium_filterBlanks.py
│
└── README.md
```
