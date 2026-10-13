## CHIMERA generates synthetic health data with balanced fidelity, utility and privacy

This repository contains the code and workflow for the paper "CHIMERA generates synthetic health data with balanced
fidelity, utility and privacy" by Amadou Tidiane Niang (University of Lille).

##  Getting Started

The repository is designed to run the analysis workflow in a highly automated manner. To ensure proper execution after cloning:

1.  **Install all required packages**.\
    Each script attempts to install missing packages automatically.

2.  **Download and place all datasets** in the `data/` folder. the real REIN registry dataset are not included in the repository due to privacy.

3.  **Run the scripts sequentially** from `0_Methods_and_metrics_functions.R` to `6_REIN_Mortality_Score_Validation_Across_Synthesis_Methods.R`.\
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
##  Essential Data And Benchmark Overview for Reproducing the Analysis

To successfully run the analysis and reproduce the results, all input datasets must be available locally. Create a folder named `data/` in the same directory as the `scripts/` folder, and place all required datasets inside it.

| Dataset | N | Variables | Outcome | Access |
|---|---|---|---|---|
| [PIMA](https://www.kaggle.com/datasets/uciml/pima-indians-diabetes-database) | 768 | 8 continuous | Binary (diabetes) | Public |
| AIDS (Hammer et al., 1996) | 2,139 | 26 mixed | Survival | Public |
| REIN (French national registry) | 42,176 | 16 mixed | Survival (3-month mortality) | Restricted |

---

CHIMERA is compared with five synthetic data generators on three datasets:

| | |
|---|---|
| **Methods** | CHIMERA, SYNTHPOP, CTGAN, Gaussian Copula, TVAE, TabDiff |
| **Datasets** | PIMA (logistic setting), AIDS (Cox setting), REIN (Cox setting) |
| **Replicates** | 50 independent synthetic datasets per method and dataset |
| **Metrics** | fidelity / utility / privacy metrics (T1–T13) |


The best-of-*M* selection uses T1–T9 only (T10–T13 are reported but not used for selection).

---

## Repository Structure

```
.
├── data/        # input datasets, train/test/holdout splits, imported CTGAN/TVAE/GC/TabDiff CSVs
├── results/     # .Rdata  outputs of every script
└── scripts/
    ├── 0_Methods_and_metrics_functions.R
    ├── 1_Preprocessing.R
    ├── 2_Generate_Synthetic_Data.R
    ├── 3_Synthetic_Data_Generation_Evaluation.R
    ├── 4_Comprehensive_Synthetic_Data_Evaluation.R
    ├── 5_CHIMERA_Ablation_Study_And_Sensitivity_Analysis.R
    └── 6_REIN_Mortality_Score_Validation_Across_Synthesis_Methods.R
```

Scripts are meant to be run in numerical order from the project root (paths are handled with `here`). Each script sources `scripts/0_Methods_and_metrics_functions.R`, reads from `data/` and writes to `results/`.

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
| Distance to the Closest Record | T10 - (DCR)| The proximity of synthetic records to real records | → 1 |
| Nearest Neighbor Distance Ratio | T11 - (NNDR)| Measures the relative proximity of an observation to its first and second nearest neighbors | → 1 |
| Linkability Attack Strength | T12 — \(s\) | Difference between the reference linkability success rate and the naive random-matching baseline | → 0 |
| Signed Specific Linkability Risk | T13 — \(R\) | Normalized excess linkability success rate in the reference condition relative to the control condition | → 0 |
---

## REIN Prognostic Workflow

Beyond distributional metrics, CHIMERA is evaluated through the full replication of a published clinical prognostic score for 3-month mortality in dialysis patients aged 75 years or older (Couchoud et al., 2015).

The workflow implemented in `6_REIN_Mortality_Score_Validation_Across_Synthesis_Methods.R` includes:

- 16 candidate predictors (3 continuous, 13 categorical); real training set n = 29,523, independent test set n = 12,653 (70/30 split)
- 50 multiple imputations (MICE, Nelson–Aalen hazard as auxiliary predictor)
- 100 bootstrap resamples per imputation → 5,000 logistic regression models
- Variable selection: retained if statistically significant in > 70% of models (Wald test, *p* < 0.05)
- Fisher's method for p-value aggregation across dummy variables of multi-level categorical predictors (heart failure, peripheral artery disease, walking autonomy)
- Final logistic models fitted in each imputed dataset with the retained predictors; predicted probabilities averaged patient by patient across the 50 imputations
- Discrimination (AUC, DeLong 95% CI) and calibration (restricted cubic splines, Brier score, intercept, slope) on the independent real test set (30%), itself multiply imputed

The workflow is replicated on each of the 50 synthetic datasets generated by each of the six methods. Agreement with the real-data reference model is quantified using sensitivity, specificity, and Cohen's kappa over the 16 candidate predictors.

---

## R Scripts Overview

#### `0_Methods_and_metrics_functions.R`

This script contains all custom functions used throughout the CHIMERA pipeline, organized into six blocks.

**Data preparation and import.** `prepare_survival_data_time()` converts raw dates into numeric follow-up times and applies administrative censoring. `import_all_synthetic_datasets()` and its helpers (`default_file_paths()`, `tabdiff_file_paths()`, `read_synthetic_datasets()`, `save_synthetic_datasets()`) read the externally generated CTGAN, TVAE, Gaussian Copula and TabDiff CSV files and store them as `.Rdata` lists. `load_first()` loads the first object of an `.Rdata` file, and `impute_mice()` performs single MICE imputation of a data frame or of a list of data frames (optionally in parallel with `future`).

**Synthetic data generation.** `CHIMERA_generate()` implements the full iterative MCAR masking, MICE imputation, Mahalanobis realignment (via `MatchIt`) and selective transfer pipeline, with optional survival support via Nelson–Aalen estimation and spline-based time reconstruction. `Ablation_Analysis_generate_synthetic_data()` provides the two simplified MICE-based alternatives used in the ablation analysis: direct MICE synthesis (`iterative_masking = FALSE`) and CHIMERA without matching (`iterative_masking = TRUE`).

**Fidelity, utility and privacy evaluation.** `Assessment_function_unified()` computes the 13 benchmarking metrics (T1–T13) for a real–synthetic dataset pair under logistic or Cox settings. Supporting functions include `fd_bins()` (Freedman–Diaconis breaks), `JSD()`, `privacy_AIR()` (1-nearest-neighbour Gower attribute inference attack), `compute_dcr()` and `compute_rrd_thresholds()` (DCR / NNDR and real-to-real reference thresholds), and `linkability_attack()`. `run_method_assessment()`, `run_dataset_evaluation()` and `run_all_evaluations()` drive the evaluation of every dataset–method combination from a single configuration list.

**Best-of-*M* selection and statistical comparison.** `compute_metrics_loss()` normalizes T1–T9 against their ideal targets, groups them into fidelity, utility and privacy families, and combines them with equal weights (1/3) into `L_total`. `analyse_M()` summarizes metric distributions and the best run for increasing numbers of candidates (*M* = 1, 5, 10, 25, 50). `bootstrap_compare_methods()` compares CHIMERA with a comparator using an independent two-sample bootstrap (10,000 resamples, percentile 95% CI) on medians, or on median absolute deviation from the target for target-based metrics.

**Analytical and predictive utility.** `fit_and_summarise_classic()` fits logistic regression and returns odds ratios, 95% CIs, standard errors, SDiff and ROC curves (PIMA). `fit_and_summarise_survival()` fits Cox models and returns hazard ratios, CIs, SDiff and time-dependent AUC curves (AIDS, REIN), using `predictive_utility_surv_curve()` for the bootstrap confidence intervals.

**REIN prognostic score workflow.** `preprocessing()` prepares the REIN dataset for imputation and constructs the 90-day mortality endpoint. `Imputed_data_function()` performs MICE-based multiple imputation using the Nelson–Aalen hazard as auxiliary predictor. `bootstrap_model()` fits logistic regression on bootstrap resamples and returns Wald-test p-values. `Score()` aggregates these p-values across 5,000 models, combining multi-level dummy variables via Fisher's method, and returns variable selection frequencies. `compute_selection_metrics()` computes sensitivity, specificity and Cohen's kappa between real- and synthetic-data selections. `clean_variable_names()`, `compute_auc_ci()`, `compute_pooled_roc()` and `compute_method_rocs()` handle AUC and pooled ROC computation, while `compute_calibration()` and `compute_calibration_synth()` evaluate calibration (restricted cubic splines, Brier score, intercept and slope) on the external test set for the real and synthetic datasets respectively.

#### `1_Preprocessing.R`

This script handles data loading, preprocessing, and splitting for the three benchmark datasets.

**PIMA.** The raw CSV is imported, columns are renamed, and the binary outcome is converted to a labelled factor (`No` / `Yes`). The full dataset is saved as `PIMA.csv`.

**AIDS.** The semicolon-delimited file is imported, identifiers and unused variables are dropped, columns are renamed for readability (e.g. `Karnofsky_score`, `CD4_baseline`, `times`, `Censored`), and binary variables are converted to labelled factors. The full dataset is saved as `AIDS.csv`.

**REIN.** The preprocessed registry file `REIN.csv` is loaded (character columns converted to factors).

**Data splits.** For each dataset, two splits are produced with `set.seed(123)`:

- a **70/30 train–test split** (`*_train.csv`, `*_test.csv`), used for generator fitting and predictive utility (and for the REIN score workflow);
- a **50/50 train–holdout split** (`*_trn.csv`, `*_hol.csv`), used for the holdout-based membership inference and linkability assessments.

#### `2_Generate_Synthetic_Data.R`

This script generates the synthetic datasets for CHIMERA, SYNTHPOP and the ablation variants.

**CHIMERA.** For each dataset, 50 independent datasets are generated with distinct seeds (`set.seed(123)`, `sample(1:1000, 50)`), with execution times recorded. PIMA is treated as a standard tabular dataset (logistic setting), while AIDS and REIN use the survival extension (`timevar = "times"`, `statusvar = "Censored"`). The masking rate is drawn once per cohort from U(0.10, 0.30) (values used in the paper: 0.25 / 0.15 / 0.30 for PIMA / AIDS / REIN), with `nb_imput = 1` and `nb_max_it = 5`. A second set of datasets is generated from the 50/50 training split (`*_trn`) for the holdout-based privacy metrics.

**SYNTHPOP.** 50 datasets per cohort are generated with the default `syn()` settings (CART), for both the 70% training set and the 50% holdout training set. Execution times are saved in `execution_time_results.csv`.

**Ablation variants.** For each dataset, 50 datasets are generated with `Ablation_Analysis_generate_synthetic_data()` under two strategies: direct MICE synthesis, and CHIMERA without Mahalanobis matching (same seeds and masking rate as CHIMERA).

**Sensitivity analysis.** CHIMERA is run 50 times for each masking rate *p* ∈ {0.05, 0.10, 0.15, 0.20, 0.25} on PIMA and AIDS (training and holdout versions).

**REIN preprocessing.** The final block regenerates the complete REIN registry through CHIMERA, then reintroduces the original missingness mask.

CTGAN, TVAE, Gaussian Copula and TabDiff datasets are generated externally in Python (SDV 1.20.0 and the TabDiff implementation) and imported in script 3.

#### `3_Synthetic_Data_Generation_Evaluation.R`

This script evaluates the six synthesis methods on the three datasets.

**Import.** The externally generated CTGAN, TVAE, Gaussian Copula and TabDiff CSV files (main and `_trn` versions) are imported with `import_all_synthetic_datasets()`, and the CHIMERA and SYNTHPOP lists generated in script 2 are loaded. REIN datasets (real and synthetic) are imputed with `impute_mice()` since they contain missing values.

**Evaluation configuration.** The attribute inference risk (AIR) scenarios are built by incrementally adding the top predictors of the target variable to the attacker's information set (`Diabetes_diagnosis` for PIMA, `Homosexuality` for AIDS, `Diabetes` for REIN). DCR and NNDR reference thresholds are computed from the real train–test distances (`compute_rrd_thresholds()`).

**Evaluation.** `run_all_evaluations()` computes the 13 metrics (T1–T13) for each of the 50 synthetic datasets of each method and dataset, including the holdout membership score and the linkability attack. The metrics and composite loss for each combination are saved as `{DATASET}_FUP_Metrics_{method}.Rdata`.

#### `4_Comprehensive_Synthetic_Data_Evaluation.R`

This script aggregates the evaluation results, selects the representative datasets and computes the analytical utility outputs.

**Aggregation and best-of-*M* selection.** `compute_metrics_loss()` produces the composite loss for all 18 method–dataset combinations. `analyse_M()` summarizes metric distributions for *M* = 1, 5, 10, 25, 50 and identifies the best run (minimum `L_total`) for each method and dataset; the results at *M* = 50 are saved in `best_all_by_M50.Rdata`.

**Bootstrap comparison.** `bootstrap_compare_methods()` compares CHIMERA with each of the five comparators on T1–T13 (10,000 resamples), separately for each dataset.

**Ablation.** The same loss, *M*-analysis and bootstrap comparison are applied to CHIMERA, CHIMERA without matching, and direct MICE synthesis.

**Analytical utility — PIMA.** Continuous variables are rescaled for interpretability, and logistic regression is fitted on the real training set and on the best synthetic dataset of each method (and on all 50 replicates). Odds ratios, 95% CIs, standard errors, SDiff, and ROC curves are saved.

**Analytical utility — AIDS and REIN.** The same workflow is applied with Cox models, producing hazard ratios, 95% CIs, SDiff values, Kaplan–Meier curves, and time-dependent AUC curves with 200-bootstrap confidence intervals on the independent test set. Results are computed for the best dataset of each method and across the 50 replicates (including Kaplan–Meier curves for all replicates).

#### `5_CHIMERA_Ablation_Study_And_Sensitivity_Analysis.R`

This script evaluates the complementary analyses designed to isolate the contribution of the Mahalanobis matching step and to study the sensitivity to the masking rate.

**Ablation evaluation.** For each dataset, `run_all_evaluations()` computes the 13 metrics for CHIMERA without matching and for direct MICE synthesis, using the same target variables, AIR scenarios, RRD thresholds and seeds as in script 3. Results are combined with those of full CHIMERA in script 4.

**Sensitivity to the masking rate.** On PIMA and AIDS, the 50 datasets generated for each masking rate *p* ∈ {0.05, 0.10, 0.15, 0.20, 0.25} are evaluated with `Assessment_function_unified()`. Results are saved as `PIMA_assessment_results_CHIMERA_sensitivity_by_p.rds` and `AIDS_assessment_results_CHIMERA_sensitivity_by_p.rds`.

#### `6_REIN_Mortality_Score_Validation_Across_Synthesis_Methods.R`

This script implements the full REIN prognostic score development and validation pipeline on the real data and on the 50 synthetic datasets of each of the six methods.

**Multiple imputation.** `preprocessing()` is applied to the real training set and to each synthetic dataset, generating 50 imputed datasets per source (parallelised via `parlmice`-style internal parallelism).

**Bootstrap-based variable selection.** `Score()` yields the selection frequencies of all candidate predictors across the 5,000 logistic regression models (50 imputations × 100 bootstrap resamples) for the real data and every synthetic dataset.

**Selection agreement.** `compute_selection_metrics()` computes sensitivity, specificity and Cohen's kappa against the real-data selection (threshold: 70%). The best-of-*M* datasets are those minimizing `L_total` in `best_all_by_M50.Rdata` (REIN rows).

**Discrimination.** AUC with DeLong 95% CI is computed for each real and synthetic dataset on the independent test set (multiply imputed, outcome excluded from the imputation model); pooled ROC curves are obtained by averaging predicted probabilities across the 50 imputations. AUC distributions across the 50 synthetic datasets per method characterize run-to-run variability.

**Calibration.** `compute_calibration()` produces restricted cubic spline calibration curves, Brier scores, intercepts and slopes for the real imputed datasets; `compute_calibration_synth()` replicates this for each of the 50 synthetic datasets of every method.

---

## Reproducibility

All analyses were performed using:

| Software | Version |
|---|---|
| R | 4.4.1 |
| Python | 3.10.12 |
| SDV (CTGAN, TVAE, Gaussian Copula) | 1.20.0 |
| synthpop | 1.8-0 |
| mice | 3.19.0 |
| MatchIt | 4.7.2 |

Run-specific seeds: for the R-based generators, a master seed (`set.seed(123)`) generates 50 seeds with `sample(1:1000, 50)`, shared by all R-based methods.

All evaluation metrics, the best-of-*M* selection procedure, and the REIN prognostic workflow are fully described in the paper (Materials and Methods and Supplementary Materials, Sections S1.2–S1.5) and implemented in this repository.

**Data availability.** PIMA is publicly available (Kaggle); AIDS is available through the original ACTG trial publications. The REIN registry data cannot be shared due to French privacy regulations.

---

## Contact

**Amadou Tidiane Niang** — `amadou-tidiane.niang@univ-lille.fr`  
Univ Lille · CHU Lille · ULR 2694 – METRICS · F-59000 Lille, France

