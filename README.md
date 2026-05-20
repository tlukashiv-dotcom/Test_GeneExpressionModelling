# Memory-Aware Semi-Markov Framework for Longitudinal Single-Cell Transcriptomics

## Overview 

This repository contains the R implementation of a memory-aware Semi-Markov framework for modelling, validation, and forecasting of longitudinal single-cell RNA sequencing (scRNA-seq) data.

The framework combines:

* Gaussian Mixture Models (GMMs);
* Wasserstein-based transition matrices;
* Semi-Markov stochastic dynamics;
* Brownian Bridge trajectory simulation;
* genome-scale validation utilities;
* longitudinal prediction diagnostics.

The methodology was developed for modelling dopaminergic neuron differentiation dynamics in the human ventral midbrain using the GSE76381 dataset.

---

## Repository Structure

```text
.
├── requirements.R
├── process_cef_data.R
├── simulation_functions.R
├── parameter_estimation.R
├── prediction_functions.R
├── validation_functions.R
├── extended_validation_functions.R
├── plotting_functions.R
├── trajectory_analysis.R
├── main_forward_validation.R
└── README.md
```

---

## Main Components

### Core Modelling

* `parameter_estimation.R`

  * Gaussian mixture parameter estimation
  * Wasserstein transition matrices

* `simulation_functions.R`

  * Semi-Markov simulation
  * Brownian Bridge interpolation

* `prediction_functions.R`

  * Forward prediction pipeline
  * Gene-level stochastic forecasting

---

### Validation

* `validation_functions.R`

  * validation metrics
  * accuracy summaries
  * marker integrity analysis

* `extended_validation_functions.R`

  * genome-scale diagnostics
  * stratified deviation analysis
  * differential expression agreement

---

### Visualization

* `plotting_functions.R`

  * publication-ready plotting utilities
  * volcano plots
  * trajectory plots
  * density comparisons
  * validation diagnostics

---

## Data

The code was developed using:

* GEO accession: GSE76381
* Human ventral midbrain developmental scRNA-seq data

Input data should contain:

* a `TP` column with numeric developmental timepoints;
* gene expression columns.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/YOUR_REPOSITORY.git
cd YOUR_REPOSITORY
```

Install and load dependencies:

```r
source("requirements.R")
```

---

## Running Forward Validation

Run:

```r
source("main_forward_validation.R")
```

Outputs will be saved automatically into:

```text
Results_Forward_Validation_New/
```

---

## Main Validation Metrics

The framework evaluates:

* Student's t-test agreement;
* Wilcoxon agreement;
* Wasserstein distances;
* log2 fold-change deviations;
* prediction calibration;
* marker preservation.

---

## Authors

* Taras Lukashiv
* Igor Malyk
* Mathias Galati
* Ahmed Hemedan
* Venkata Satagopam

---

## Status

Current repository status:

* Forward validation: stable
* Forecasting: stable
* Backward validation: experimental / under active development
