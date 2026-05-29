## \# CHIMERA generates synthetic health data with balanced fidelity, utility and privacy

This repository contains the code and workflow for the paper "CHIMERA generates synthetic health data with balanced
fidelity, utility and privacy" by Amadou Tidiane Niang (University of Lille).

##  Getting Started

The repository is designed to run the analysis workflow in a highly automated manner. To ensure proper execution after cloning:

1.  **Install all required packages**.\
    Each script attempts to install missing packages automatically.

2.  **Download and place all datasets** in the `data/` folder. the real REIN registry dataset are not included in the repository due to privacy.

3.  **Run the scripts sequentially** from `0_functions.R` to `7_Figures.R`.\
    Each script depends on outputs from the previous one, so maintaining this order is essential.

## 1. 📦 Installing Required Packages

Each `.R` script automatically installs or loads missing packages.\
👉However, one special package,**R-INLA**,is not on CRAN and must be installed manually.

#### Installation de R-INLA

``` r
if (!require('INLA')) {
  install.packages(
    'INLA',
    repos = c(getOption('repos'), INLA = 'https://inla.r-inla-download.org/R/testing'),
    dependencies = TRUE
  )
}
```

If you encounter installation issues, please refer to the [R-INLA Download & Install guide](https://www.r-inla.org/download-install).

## 2. 📥 Essential Data for Reproducing the Analysis

To successfully run the analysis and reproduce the results, all input datasets must be available locally. Create a folder named `data/` in the same directory as the `scripts/` folder, and place all required datasets inside it.\
These datasets include `municipal boundaries`, `ESKD incidence`, `social deprivation indices`, `clinical covariates`, `healthcare access`, and `air pollution measurements`.

#### 2.1 🗺️ `Shapefile` – Administrative Boundaries of French Municipalities

-   **Description:** This dataset provides the administrative boundaries of all French municipalities in 2021, covering a total of `34,969 communes`.For the purpose of this project, we restrict the analysis to metropolitan France. After excluding the `overseas territories (DOMs/TOMs)` and performing data management procedures, the final dataset includes `34,830 municipalities`, which are used in the analyses.

-   **Source:** [Administrative boundaries of French municipalities from OpenStreetMap (data.gouv.fr)](https://www.data.gouv.fr/datasets/decoupage-administratif-communal-francais-issu-d-openstreetmap/)

-   **Format:** Shapefile (.shp)

-   **Naming Convention:** `shape_commune_2021ssDOMTOM.RData`

-   **Coordinate Reference System (CRS):** Projection **WGS 84** (latitude/longitude in degrees)

#### 2.2 🏚️ `Covariate Data` – Social Deprivation (EDI)

-   **Description:** Social deprivation was measured using the `French European Deprivation Index (EDI)`, an ecological index reflecting material and social disadvantage at the `municipal level`. The EDI is a weighted composite of census-based indicators such as housing conditions, household structure, education, employment, and nationality. `Higher EDI values indicate greater deprivation.` Municipal-level EDI data were available for `2011, 2015, and 2017`. For intermediate years (2012–2014 and 2016), values were interpolated using interannual averages, while values from `2018–2021` were assumed constant based on 2017 estimates.

-   **Source:** [MapInMed](https://mapinmed.unicaen.fr/)

-   **Format:** Municipality-level index (period: 2011–2021, interpolated where necessary)

-   **Naming Convention:** `EDI_scale`

#### 2.3 🩺 `Covariate Data` – Diabetes and Hypertension Prevalence

-   **Description:** Annual municipal-level prevalences of `diabetes` and `hypertension` were estimated using the `pathology and expenditure mapping system` of the French `National Health Data System (SNDS)`, which covers \~99% of the French population.

    -   `Diabetes prevalence`: based on individuals receiving diabetes care (any type).\
    -   `Hypertension prevalence`: based on individuals with at least three antihypertensive drug deliveries within the same year.

-   **Source:** [SNDS – French National Health Data System](https://snds.gouv.fr)

-   **Format:** Annual prevalence rates at the municipality level (latest version: G10).

-   **Naming Convention:** `prevalence_Diab_scale`

#### 2.4 🏥 `Covariate Data` – Dialysis Center Accessibility

-   **Description:** Healthcare access was assessed for each municipality and year (2012–2021) by measuring the presence of at least one `dialysis center` within a `30-minute drive` from the municipality center. This measure included all types of dialysis units, such as `in-center dialysis` and `self-care facilities`.

-   **Source:** National REIN Registry & geospatial accessibility analysis.

-   **Format:** Annual binary indicator (accessible / not accessible) at the municipality level.

-   **Naming Convention:** `accSoins_binaire`

#### 2.5 🌍 `Covariate Data` – PM2.5 Exposure

-   **Description:** Long-term air pollution data were obtained from the `National Institute for the Industrial Environment and Risks (INERIS)`. Using the `CHIMERE Chemistry Transport Model` combined with `kriging techniques` and meteorological parameters (temperature, humidity, precipitation), mean annual `PM2.5 concentrations (μg/m³)` were reconstructed over France.

-   **Time Coverage:** 21 years of modeled data.

-   **Metric:** Municipality-level `3-year lagged average` of PM2.5, weighted by population density and adjusted for municipal area.

-   **Spatial Resolution:** \~4 km grid.

-   **Source:** [INERIS](https://www.ineris.fr/fr)

-   **Naming Convention:** `PM25_lag3`

#### 2.6 💉 `Outcome Data` – Incidence of End-Stage Kidney Disease (ESKD)

-   **Description:** Data on the incidence of `End-Stage Kidney Disease (ESKD)` were obtained from the `national REIN registry`, which is the official database monitoring all patients under `Kidney Replacement Therapy (KRT)` in France. The registry records detailed patient-level information at the initiation of dialysis or transplantation and comprehensively covers the entire French population. For this study, we included `all individuals who initiated maintenance KRT between January 1, 2012, and December 31, 2021`. Based on the patient’s `place of residence` at the time of KRT initiation, the `annual number of incident ESKD cases` was aggregated at the `municipality level` for the period `2012–2021`.

    In the final dataset, we include:

    -   All previously presented **covariates** (Social Deprivation (EDI),Diabetes and Hypertension Prevalence, Healthcare Access,Air Pollution, etc.).\
    -   The **number of observed incident ESKD cases**.\
    -   The **number of expected cases** (based on standardization procedures).\
    -   Spatial, temporal, and spatio-temporal indices for modeling:
        -   `ID.space` (municipality identifier)\
        -   `ID.time`, `ID.time2` (temporal identifiers)\
        -   `ID.space.time` (spatio-temporal identifier)

-   **Source:** REIN Registry (Réseau Epidémiologie et Information en Néphrologie) – [Official Website](https://www.agence-biomedecine.fr/fr/observatoire-de-la-maladie-renale-chronique/le-registre-rein)

-   **Format:** Aggregated data (municipality-level counts, annual time series 2012–2021)

-   **Naming Convention:** `REIN_DATA_ST_2012_2021.RData`

Although raw ESKD incidence data cannot be shared for confidentiality reasons, we provide spatial relative risks (RRs) for every municipality in metropolitan France. Municipalities are identified by `SP_ID` , corresponding to the official INSEE geographic codes (2021 administrative boundaries).

Spatial RRs are calculated as the posterior median of the exponentiated spatially structured random effect estimated using a BYM2 model within a spatio-temporal disease-mapping framework without covariates (see the publication’s supplementary materials for methodological details).

In addition to the continuous estimates (`resRR_unadjust`), we include:

-   A categorized version of the spatial RRs (`resRRcat_unadjust`), based on the thresholds defined in the publication.

-   The 95% credible interval bounds for each RR (`resRR_unadjust_low` and `resRR_unadjust_high`).

-   The exceedance probability (`PP_unadjust`: probability that the estimated random effect exceeds zero), provided both as raw values and as a categorized variable (`PPcat_unadjust`) using thresholds of 0.2 and 0.8 (see supplementary materials for details).

-   All data are stored in the file `ESKD_spatialRR.RData` (an `sf` object), located in the `data/` directory.

### 3. 📂 R Scripts Overview

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
