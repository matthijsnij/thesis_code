# thesis_code

This repository contains the R code and thesis document for my MSc thesis:

**S-MPBART: Soft Multinomial Probit Bayesian Additive Regression Trees**

The full thesis is available here: [MSc_Thesis_MatthijsNijeboer.pdf](MSc_Thesis_MatthijsNijeboer.pdf)

---

## Repository Structure

```
thesis_code/
├── MSc_Thesis_MatthijsNijeboer.pdf  # Full thesis document
├── scripts/               # All R source files
├── data/                  # Preprocessed real-world datasets (CSV)
├── output/                # Excel files with prediction results per method/dataset
└── calibration_plots/     # Saved RDS objects used for calibration plot generation
```

---

## Scripts

| File | Description |
|------|-------------|
| `soft_mpbart.R` | Core implementation of S-MPBART: Gibbs sampler with SoftBART trees, truncated MVN sampling for latent variables, and covariate rank-normalization |
| `mpbart.R` | Wrapper around `GcompBART::model_bart` for running MPBART (hard splits) as in Xu et al. (2024) |
| `random_forest.R` | Random forest classifier with 5-fold cross-validation for tuning `mtry` |
| `simulation_functions.R` | Data-generating functions for DGP 1–3, and `run_method()` helper for simulation experiments |
| `dgp1.R` | Simulation study for **DGP 1** — piecewise sinusoids with a discontinuity (3-class, 2 predictors) |
| `dgp2.R` | Simulation study for **DGP 2** — Friedman-style latent functions (3-class, 10 predictors; also with 50 extra noise predictors) |
| `dgp3.R` | Simulation study for **DGP 3** — Gaussian process-inspired, choice-specific predictors (3-class, 16 predictors) |
| `data_prep.R` | Loads raw datasets, preprocesses and normalizes them, and writes `*_preprocessed.csv` files to `data/` |
| `data_application.R` | Runs a chosen method on a chosen real-world dataset using cross-validation folds |
| `data_application_functions.R` | Helper functions: `read_data()` and `run_method()` for the real-data experiments |
| `calibration_plots.R` | Generates calibration plots from saved RDS files in `calibration_plots/` |

---

## Datasets

Nine real-world multiclass datasets are included (preprocessed):

| Dataset | Classes |
|---------|---------|
| `glass` | 6 |
| `vertebral` | 3 |
| `iris` | 3 |
| `wine` | 3 |
| `vehicle` | 4 |
| `vowel` | 11 |
| `waveform` | 3 |
| `travel` | 3 |
| `fishing` | 4 |

---

## Output

The `output/` folder contains Excel files named `{method}_{dataset/dgp}_output.xlsx` with per-fold or per-replication predictions and accuracy metrics for all three methods across all datasets and simulation settings.

---

## Dependencies

Install the required R packages before running any scripts:

```r
install.packages(c(
  "GcompBART", "SoftBart", "randomForest", "caret",
  "MASS", "MCMCpack", "tmvmixnorm", "coda",
  "Rfast", "LaplacesDemon", "fields",
  "glue", "openxlsx", "writexl", "ggplot2",
  "AER", "mlogit"
))
```

> **Note:** `GcompBART` and `SoftBart` may need to be installed from GitHub.

---

## Usage

### Simulation Studies

Set the `method` variable (`"smpbart"`, `"mpbart"`, or `"rf"`) and run the corresponding DGP script:

```r
# e.g., DGP 1 with S-MPBART
source("scripts/dgp1.R")
```

### Real-Data Experiments

Set `dataset` and `method` in `data_application.R`, then run it:

```r
dataset <- "wine"   # one of the nine datasets above
method  <- "smpbart"
source("scripts/data_application.R")
```

Results are written to `output/{method}_{dataset}_output.xlsx`.

