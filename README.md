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

This script contains all custom functions used throughout the CHIMERA pipeline, organized into four blocks.

**Data preparation.** `prepare_survival_data_time()` converts raw dates into numeric follow-up times and applies administrative censoring. `fd_bins()` computes Freedman–Diaconis histogram breaks. `JSD()` computes Jensen–Shannon divergence between two distributions.

**Synthetic data generation.** `CHIMERA_generate()` implements the full iterative MCAR masking, MICE imputation, and Mahalanobis realignment pipeline, with optional survival data support via Nelson–Aalen estimation and spline-based time reconstruction. `Ablation_Analysis_generate_synthetic_data()` provides two simplified MICE-based alternatives — direct synthesis and iterative masking without realignment — used in the ablation analysis. `predictive_utility_surv_curve()` computes time-dependent AUC with bootstrap confidence intervals for Cox models.

**Fidelity, utility, and privacy evaluation.** `Assessment_function_unified()` computes the nine benchmarking metrics (T1–T9) for a real–synthetic dataset pair under logistic or Cox settings. `Assessment_function_unified_Ablation_Analysis()` provides the same evaluation for the ablation strategies. `privacy_AIR()` implements the attribute inference risk metric via a 1-nearest-neighbour Gower distance attack.

**REIN prognostic score workflow.** `preprocessing()` prepares the REIN dataset for imputation and constructs the 90-day mortality endpoint. `Imputed_data_function()` performs MICE-based multiple imputation using the Nelson–Aalen hazard as an auxiliary predictor. `bootstrap_model()` fits logistic regression on bootstrap resamples and returns Wald-test p-values. `Score()` aggregates these p-values across 5,000 models — combining multi-level dummy variables via Fisher's method — and returns variable selection frequencies. `compute_selection_metrics()` computes sensitivity, specificity, and Cohen's kappa between real- and synthetic-data variable selections. `compute_calibration()` and `compute_calibration_synth()` evaluate model calibration on an external test set using restricted cubic spline curves, Brier score, and calibration slope, for real imputed and synthetic datasets respectively.

#### `1_preprocessing.R`

This script handles data loading, preprocessing, and splitting for the three benchmark datasets used in the study.

**PIMA.** The raw CSV is imported, columns are renamed, and the binary outcome is converted to a labelled factor. The dataset is split 50/50 into training and holdout sets for the membership inference privacy evaluation.

**AIDS.** The dataset is imported, identifiers and unused variables are dropped, columns are renamed for readability, and binary variables are converted to labelled factors. The same 50/50 train–holdout split is applied.

**REIN.** Two Excel files are merged and filtered. Raw date columns are converted to numeric follow-up times in days using `prepare_survival_data_time()`, with administrative censoring set to 31 December 2022. Variables with more than 50% missing values are excluded. The remaining variables are renamed, recoded into labelled factors (including ordered categories for heart failure, peripheral artery disease, and walking autonomy), and the cohort is restricted to patients aged 75 years or older. The resulting dataset is saved both with and without missing values. Two separate splits are then produced: a 50/50 train–holdout split for the privacy evaluation, and a 70/30 train–test split for the prognostic score replication workflow. Missing values in the training sets are handled by a single MICE imputation using predictive mean matching for continuous variables.

#### `2_BestOf_M_DatasetGenerator.R`

This script generates and evaluates the 50 synthetic datasets produced by each method (CHIMERA, SYNTHPOP, CTGAN) across the three benchmark datasets (PIMA, AIDS, REIN).

**Synthetic data generation.** For CHIMERA and SYNTHPOP, 50 independent datasets are generated using distinct random seeds, with execution times recorded for each run. PIMA is treated as a standard tabular dataset (logistic setting, masking rate 25%), while AIDS and REIN use the survival extension (Cox setting, masking rates 15% and 30% respectively). CTGAN datasets are pre-generated externally in Python and imported from CSV files; two sets are loaded for each dataset — one generated on the full dataset, one on the training set only, for use in the holdout-based privacy metric.

**Evaluation.** For each synthetic dataset, `Assessment_function_unified()` computes the nine fidelity–utility–privacy metrics (T1–T9). The attribute inference risk scenarios are constructed by incrementally adding the top five predictors (ranked by |z-value| from a logistic regression on the real data) to the attacker's information set.

**Best-of-M selection.** For each method–dataset combination, the composite objective function is computed by normalizing each metric relative to its ideal target vector, grouping metrics into three families (fidelity, utility, privacy), and combining the family-level losses with equal weights (1/3 each). The resulting ranked metric tables are saved for downstream selection of the representative synthetic dataset reported in the paper.

#### `3_BestOf_M_Utility.R`

This script selects the representative synthetic dataset for each method–dataset combination and computes the analytical utility outputs used in the paper.

**Best-of-M selection.** For each dataset (PIMA, AIDS, REIN) and each method (CHIMERA, SYNTHPOP, CTGAN), the composite loss score computed in `2_BestOf_M_DatasetGenerator.R` is used to identify the synthetic dataset minimizing `L_total`. The corresponding dataset is extracted from the 50-run list and saved for downstream analyses.

**Analytical utility — PIMA.** Continuous variables are rescaled for interpretability before fitting logistic regression models on the real dataset and each selected synthetic dataset. Odds ratios, 95% confidence intervals, and standard errors are extracted for all predictors. Standardized differences (SDiff) between log odds ratios are computed for each synthetic method relative to the real-data estimates. ROC curves are computed and saved for each dataset.

**Analytical utility — AIDS and REIN.** The same workflow is applied using Cox proportional hazards models, producing hazard ratios, 95% confidence intervals, standard errors, and SDiff values. Kaplan–Meier survival curves are fitted across all four datasets (real, CHIMERA, SYNTHPOP, CTGAN) and saved. Time-dependent AUC curves with 200-bootstrap confidence intervals are computed using `predictive_utility_surv_curve()` over a predefined grid of time points, for both within-dataset discrimination (70/30 split) and cross-dataset transfer scenarios.

#### `4_BestOf_M_Dataset_Generator_Ablation_Analysis.R`

This script implements the complementary ablation analysis designed to isolate the contribution of the Mahalanobis-distance matching step in CHIMERA.

**Synthetic data generation.** For each of the three benchmark datasets, 50 independent synthetic datasets are generated under two simplified MICE-based strategies using `Ablation_Analysis_generate_synthetic_data() `: direct MICE synthesis (iterative_masking = FALSE), which applies a single full-dataset imputation without masking or matching; and CHIMERA without matching (iterative_masking = TRUE), which applies the iterative MCAR masking and progressive reconstruction but omits the post-imputation Mahalanobis realignment step. The same seeds, masking rates, and survival settings as in the main benchmark are used throughout.

**Evaluation.** For each strategy and dataset, `Assessment_function_unified_Ablation_Analysis()` computes the nine fidelity–utility–privacy metrics (T1–T9). The attribute inference risk scenarios are constructed identically to those in `2_BestOf_M_DatasetGenerator.R`. The composite loss score is then computed using the same target vector, family grouping, and equal weighting scheme, producing ranked metric tables for each strategy–dataset combination. All results are saved for downstream best-of-M selection and comparison with full CHIMERA.

#### `5_BestOf_M_Dataset_Generator_Score_Analysis.R`

This script replicates the synthetic data generation and best-of-M selection pipeline specifically for the REIN prognostic score workflow, using the 70% training split of the REIN registry rather than the full dataset.

**Synthetic data generation.** Fifty independent synthetic datasets are generated from the REIN training set `df_rein_train_with_na_score.csv` using CHIMERA (masking rate 30%, survival mode) and SYNTHPOP (default settings). CTGAN datasets are pre-generated externally and imported from CSV, with a separate set of training-only CTGAN datasets also loaded for the holdout-based privacy metric.

**Evaluation and best-of-M selection.** For each method, `Assessment_function_unified()` computes the nine fidelity–utility–privacy metrics (T1–T9) on the training split, using the same Cox model formula, attribute inference scenarios, and composite loss function as in the main benchmark. The resulting ranked metric tables are saved and used downstream to select the representative synthetic dataset for the prognostic score replication workflow.

#### `6_Score.R`

This script implements the full REIN prognostic score development and validation pipeline on the real and selected synthetic datasets.

**Preprocessing and multiple imputation.** The original REIN training and test sets are loaded and the positions of real missing values are recorded. These positions are reintroduced into the CHIMERA and CTGAN synthetic datasets to reproduce the original missingness structure before analysis. The `preprocessing()` function is then applied to the real training set and to each of the 50 synthetic datasets from all three methods, generating 50 multiply imputed datasets per source.

**Bootstrap-based variable selection.** `Score()` is applied to each set of imputed datasets, yielding bootstrap selection frequencies for all candidate predictors across the 5,000 logistic regression models (50 imputations × 100 bootstrap resamples). Results are saved for the real data and for every synthetic dataset from each method.

**Best-of-M selection.** For each method, the composite fidelity–utility–privacy loss computed in `5_BestOf_M_Dataset_Generator_Score_Analysis.R` is attached to the variable-selection agreement metrics (sensitivity, specificity, Cohen's kappa). The synthetic dataset minimizing `L_total` is selected as the representative dataset for downstream reporting.

**Discrimination.** Pooled ROC curves and AUC values are computed for the real data and each selected synthetic dataset by averaging predicted probabilities across the 50 imputed datasets and evaluating discrimination on the independent test set. AUC results are also computed across all 50 synthetic datasets per method to characterize run-to-run variability.

**Calibration.** `compute_calibration()` produces restricted cubic spline calibration curves, Brier scores, calibration intercepts, and calibration slopes for the real imputed datasets on the test set. `compute_calibration_synth()` replicates this evaluation for each of the 50 synthetic datasets from CHIMERA, SYNTHPOP, and CTGAN, enabling method-level comparison of calibration variability.

#### `7_Figures.R`

This script produces all publication-ready figures reported in the main manuscript and supplementary materials.

**Metric distribution boxplots (Supplementary Figs. S2–S4 and S5–S7).** For the main benchmark, T1–T9 metrics are reshaped into long format and displayed as grouped boxplots across the 50 runs for each method–dataset combination, with a gold diamond marking the best-of-M run and a dashed line indicating the ideal target value. The same layout is reproduced for the ablation analysis, replacing SYNTHPOP and CTGAN with CHIMERA without matching and direct MICE synthesis.

**Univariate distribution comparison (Fig. 1).** For each dataset, continuous variables are displayed as kernel density estimates and categorical variables as grouped bar charts of empirical proportions.

**PIMA — ROC curve and forest plot (Fig. 2).** A ROC curve panel shows AUC values with 95% confidence intervals for the real dataset and each synthetic method. A multi-row forest plot displays odds ratios with 95% confidence intervals from logistic regression, along with standardized differences (SDiff) between synthetic and real-data estimates. The two panels are combined vertically.

**AIDS and REIN — Kaplan–Meier, time-dependent AUC, and forest plot (Figs. 3–4).** Kaplan–Meier survival curves and time-dependent AUC curves (with integrated AUC annotations) are placed side by side in the upper panel. A multi-row Cox model forest plot displaying hazard ratios and SDiff values occupies the lower panel.

**REIN score — selection performance boxplots (Supplementary Fig. S8).** Sensitivity, specificity, Cohen's kappa, and AUC across the 50 synthetic datasets are displayed as boxplots for each method, with the best-of-M run highlighted.

**REIN score — ROC and calibration curves (Figs. 5 and S9).** The final model ROC curve is plotted for all four datasets. Restricted cubic spline calibration curves are shown for the best synthetic run of each method alongside the real-data reference. The supplementary version additionally overlays all 50 synthetic calibration curves per method in light grey to characterize between-run variability.

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
