## CHIMERA generates synthetic health data with balanced fidelity, utility and privacy

This repository contains the code and workflow for the paper "CHIMERA generates synthetic health data with balanced
fidelity, utility and privacy" by Amadou Tidiane Niang (University of Lille).

##  Getting Started

The repository is designed to run the analysis workflow in a highly automated manner. To ensure proper execution after cloning:

1.  **Install all required packages**.\
    Each script attempts to install missing packages automatically.

2.  **Download and place all datasets** in the `data/` folder. the real REIN registry dataset are not included in the repository due to privacy.

3.  **Run the scripts sequentially** from `0_functions.R` to `7_Figures.R`.\
    Each script depends on outputs from the previous one, so maintaining this order is essential.

---

## Overview

CHIMERA is an imputation-based framework for generating fully synthetic tabular health data. It extends the Multiple Imputation by Chained Equations (MICE) framework to address a key limitation of standard conditional imputation: the tendency to generate individually plausible values that are less coherent in the joint multivariate space.

The method combines four sequential components:

1. **Iterative MCAR masking** — progressive, controlled introduction of artificial missing values.
2. **Conditional MICE imputation** — variable-type-specific imputation (predictive mean matching, logistic regression, multinomial logistic regression).
3. **Mahalanobis-distance realignment** — post-imputation one-to-one matching to anchor synthetic profiles in realistic multivariate neighbourhoods of the real data.
4. **Survival data support** — Nelson–Aalen cumulative hazard integration and spline-based reconstruction of survival times.

The objective is to generate synthetic datasets that simultaneously preserve statistical fidelity, analytical and predictive utility, and privacy protection — three dimensions evaluated jointly using nine complementary metrics.

---
##  Essential Data for Reproducing the Analysis

To successfully run the analysis and reproduce the results, all input datasets must be available locally. Create a folder named `data/` in the same directory as the `scripts/` folder, and place all required datasets inside it.

| Dataset | N | Variables | Outcome | Access |
|---|---|---|---|---|
| [PIMA](https://www.kaggle.com/datasets/uciml/pima-indians-diabetes-database) | 768 | 8 continuous | Binary (diabetes) | Public |
| AIDS (Hammer et al., 1996) | 2,139 | 26 mixed | Survival | Public |
| REIN (French national registry) | 42,176 | 16 mixed | Survival (3-month mortality) | Restricted |

---

## Evaluation Framework

CHIMERA is evaluated using nine complementary metrics covering three dimensions. All metrics are implemented in `R/assessment.R`.

#### Fidelity

| Metric | Symbol | Description | Optimal value |
|---|---|---|---|
| Jensen–Shannon divergence | T1 — JSD | Mean marginal distributional similarity (Freedman–Diaconis binning) | → 0 |
| Correlation Frobenius Distance | T2 — CFD | Bivariate structure preservation (Pearson + Cramer's V) | → 0 |
| Discriminative performance | T3 — DiscPred | AUC of a real-vs-synthetic logistic classifier (5-fold CV) | → 0.5 |

#### Utility

| Metric | Symbol | Description | Optimal value |
|---|---|---|---|
| Standardized effect difference | T4 — SDiff | Discrepancy in regression coefficients (OR or HR) | → 0 |
| TSTR AUC gap | T5 | Train on synthetic, test on real vs. train on real | → 0 |
| TSRTR AUC gap | T6 | Train on synthetic + real, test on real vs. train on real | → 0 |

#### Privacy

| Metric | Symbol | Description | Optimal value |
|---|---|---|---|
| Record-matching risk | T7 — MIR-RM | Mean Gower distance: real to nearest synthetic record | → 1 |
| Holdout membership inference | T8 — MIR-Holdout | P(synthetic closer to training than to holdout) | → 0.5 |
| Attribute inference risk | T9 — AIR (F1) | F1 score of a 1-NN attribute reconstruction attack | → 0 |

---

## REIN Prognostic Workflow

Beyond distributional metrics, CHIMERA is evaluated through the full replication of a published clinical prognostic score for 3-month mortality in dialysis patients aged 75 years or older (Couchoud et al., 2015).

The workflow implemented in `6_Score.R` includes:

- 50 multiple imputations
- 100 bootstrap resamples per imputation → 5,000 logistic regression models
- Variable selection: retained if statistically significant in ≥ 70% of models (Wald test, *p* < 0.05)
- Fisher's method for p-value aggregation across dummy variables of multi-level categorical predictors
- Coefficient pooling using Rubin's rules
- Discrimination (AUC) and calibration (restricted cubic splines) on an independent test set (30%)

Agreement with the real-data reference model is quantified using sensitivity, specificity, and Cohen's kappa.

---

## R Scripts Overview

#### `0_functions.R`

This script contains all custom functions used throughout the CHIMERA pipeline, organized into four functional blocks.

**Data preparation.** `prepare_survival_data_time()` converts raw date columns into numeric follow-up times in days, applies administrative censoring, and optionally computes a censoring indicator. `fd_bins()` implements the Freedman–Diaconis rule for histogram binning of continuous variables. `JSD()` computes the Jensen–Shannon divergence between two discrete probability distributions.

**Synthetic data generation.** `CHIMERA_generate()` is the core generation function. It implements the full iterative MCAR masking, MICE imputation, and Mahalanobis-distance realignment pipeline, with optional support for censored survival outcomes via Nelson–Aalen estimation and spline-based time reconstruction. `Ablation_Analysis_generate_synthetic_data()` implements two simplified MICE-based synthesis strategies — direct synthesis and iterative masking without realignment — used in the complementary ablation analysis to isolate the contribution of the matching step. `predictive_utility_surv_curve()` computes time-dependent AUC with bootstrap confidence intervals for Cox models.

**Fidelity, utility, and privacy evaluation.** `Assessment_function_unified()` computes the full set of nine benchmarking metrics (T1–T9) for a given real–synthetic dataset pair, supporting both logistic and Cox model settings and three synthesizer types (CHIMERA, SYNTHPOP, CTGAN). `Assessment_function_unified_Ablation_Analysis()` provides the same evaluation pipeline adapted to the ablation synthesis strategies. `privacy_AIR()` implements the attribute inference risk metric using a 1-nearest-neighbour attack based on Gower distance.

**REIN prognostic score workflow.** `preprocessing()` prepares the REIN survival dataset for multiple imputation, including variable selection, renaming, and construction of the 90-day mortality endpoint. `Imputed_data_function()` performs multiple imputation using MICE with Nelson–Aalen cumulative hazard as an auxiliary predictor, then reattaches survival outcomes to each completed dataset. `bootstrap_model()` fits a logistic regression model on a bootstrap resample and returns Wald-test p-values for all coefficients. `Score()` applies this bootstrap procedure across all imputed datasets, aggregates p-values using Fisher's method for multi-level categorical predictors, and computes variable selection frequencies across the 5,000 resulting models. `compute_selection_metrics()` quantifies agreement between real- and synthetic-data variable selection using sensitivity, specificity, and Cohen's kappa. `clean_variable_names()` maps coefficient names from model output back to original column names in the dataset. `compute_calibration()` and compute_calibration_synth() evaluate model calibration on an external test set using restricted cubic spline calibration curves, Brier score, calibration intercept, and calibration slope, for real imputed datasets and synthetic datasets respectively.

---

## Reproducibility

All analyses were performed using:

| Software | Version |
|---|---|
| R | 4.4.1 |
| Python | 3.10.12 |
| SDV (CTGAN) | 1.20.0 |
| synthpop | 1.8-0 |

All evaluation metrics, the best-of-*M* selection procedure, and the REIN prognostic workflow are fully described in the paper (Supplementary Materials, Sections S2–S5) and implemented in this repository.

---

## Contact

**Amadou Tidiane Niang** — `amadou-tidiane.niang@univ-lille.fr`  
Univ Lille · CHU Lille · ULR 2694 – METRICS · F-59000 Lille, France
