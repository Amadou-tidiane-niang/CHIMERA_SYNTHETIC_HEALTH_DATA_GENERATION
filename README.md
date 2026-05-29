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

### R Scripts Overview

This section provides a brief description of each R script in the repository and its role in the analysis pipeline.

#### 3.1 `00_functions.R`

This script contains all the custom functions used throughout the analysis. These functions are primarily focused on post-processing outputs from the Bayesian hierarchical spatiotemporal models fitted using INLA.

#### 3.2 `01_Models_selection.R`

This script fits a series of candidate spatio-temporal Poisson BYM2 models to the ESKD incidence data and performs model selection based on fit criteria such as `DIC`, `WAIC`, and `log-CPO`.

#### 3.3 `02_Disease_mapping.R`

This script fits the core spatio-temporal BYM2 model (with a Type I space-time interaction and an unstructured temporal effect) to the ESKD data. It extracts area-specific and temporal relative risks, including posterior estimates, credible intervals, and exceedance probabilities, which are saved for subsequent mapping and analysis.

#### 3.4 `03_NonLinear_Analysis.R`

This script fits non-linear ecological regressions using spatio-temporal BYM2 models with RW2 (second-order random walk) splines for key covariates, including `deprivation (EDI)`, `diabetes and hypertension` prevalence, and `PM2.5 exposure (lagged)`. The output consists of posterior relative risks and credible intervals for each covariate, which are visualized in combined plots to illustrate potential non-linear associations.

#### 3.5 `04_EcoReg_Multivariable.R`

This script fits multivariable spatio-temporal BYM2 models including key covariates such as `deprivation (EDI)`, `diabetes prevalence`, `healthcare access`, and `PM2.5 exposure`. It extracts adjusted temporal and spatial relative risks, as well as posterior probabilities for space–time interactions. Additionally, the script evaluates the proportion of spatial variance explained by each covariate through ablation models and produces forest plots for visualizing adjusted effects.

#### 3.6 `05_Interactions_Analysis.R`

This script evaluates space-time ecological interactions in the multivariable BYM2 framework. It fits models testing interactions between deprivation (EDI) and key covariates, including PM2.5 exposure (quartiles), healthcare access (binary), and diabetes prevalence (continuous and quartiles). Posterior estimates and relative risks for each interaction term are extracted and tabulated, allowing assessment of effect modification across subgroups.

#### 3.7 `06_Mediation_Analysis.R`

This script performs a Bayesian mediation analysis using the spatio-temporal BYM2 model. It evaluates the indirect effect of social deprivation (EDI) on ESKD incidence mediated through diabetes prevalence, adjusting for PM2.5 exposure and healthcare access. Posterior samples are drawn from the mediator and outcome models to compute total, direct, and indirect effects, as well as the proportion mediated. Results are tabulated and visualized on an annotated DAG for intuitive interpretation.

#### 3.8 `07_Maps.R`

This script implements the full mapping workflow for end-stage kidney disease (ESKD) and ecological covariates. It produces spatial relative risk maps, exceedance probability maps, and space-time interaction visualizations. Additionally, it generates maps and density plots for ecological covariates (diabetes, PM2.5, deprivation index) and healthcare access. All maps include embedded bar or density insets for improved interpretability, and outputs are saved as high-resolution images.

#### 3.9 `08_PAF.R`

This script computes population attributable fractions (PAFs) and attributable numbers for end-stage kidney disease (ESKD) associated with multiple exposures, including PM2.5, diabetes prevalence, deprivation index (EDI), and healthcare access. Using the multivariable INLA spatio-temporal model, it calculates PAFs relative to various reference levels, accounting for uncertainty bounds. Results are compiled into a unified table with confidence intervals and exported in LaTeX format.

#### 3.10 `09_Sensitivity_Analyses.R`

This script performs multiple sensitivity analyses for the multivariable spatio-temporal INLA models of end-stage kidney disease. Analyses include:\
- Using alternative priors for BYM and IID components.\
- Multivariable ecological regressions with different covariates (diabetes, PM2.5, deprivation, healthcare access).\
- Estimation of temporal and spatial random slopes for deprivation (French-EDI) to assess variation in effects across time and space.\
- Visualization of model results via forest plots, spatial maps, and density plots with exceedance probabilities.

Outputs include forest plots for adjusted relative risks, temporal random slope plots, and spatial random slope maps, all saved in the `results/` folder.

To ensure proper execution and reproducibility, all scripts should be run sequentially, beginning with `00_Functions.R` and proceeding in order through `09_Sensitivity_Analyses.R`. Each script builds upon the outputs of the previous scripts, so maintaining this sequence is essential for correct results
