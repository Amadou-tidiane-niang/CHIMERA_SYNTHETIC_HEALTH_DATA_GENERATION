# ===============================================================================================
# INITIAL SETUP
# ===============================================================================================

# Clear environment
rm(list = ls())

# -----------------------------------------------------------------------------------------------
# LOAD LIBRARIES
# -----------------------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)         # File path management
  library(tidyverse)    # Data manipulation
  library(mice)         # Multiple imputation (CHIMERA)
  library(missMethods)  # Missing data simulation
  library(MatchIt)      # Matching (not directly used)
  library(readxl)       # Excel import (not directly used)
  library(synthpop)     # Synthetic data generation
  library(pROC)         # ROC / AUC
  library(fitdistrplus) # Distribution fitting
  library(caret)        # ML utilities
  library(survivalROC)  # Survival AUC
  library(survAUC)      # Time-dependent AUC
  library(gower)        # Distance metrics
  library(FNN)          # Nearest neighbors
  library(vcd)          # Cramer's V
  library(future.apply)
  library(future)
  
  
})


# Load custom functions
source(here("scripts", "0_Methods_and_metrics_functions.R"))

# parameter
n_synth <- 50
set.seed(123)
seeds <- sample(1:1000, n_synth)

# ===============================================================================================
# Loading data
# ===============================================================================================
args <- c("data", "results")

data_dir <- normalizePath(args[1], mustWork = TRUE)

if (!dir.exists(args[2])) {
  dir.create(args[2], recursive = TRUE)
}

results_dir <- normalizePath(args[2], mustWork = TRUE)


PIMA <- read.csv(here(data_dir, "PIMA.csv")) |> mutate_if(is.character, as.factor)
AIDS <- read.csv(here(data_dir, "AIDS.csv")) |> mutate_if(is.character, as.factor)
REIN <- read.csv(here(data_dir, "REIN.csv")) |> mutate_if(is.character, as.factor)

PIMA_train <- read.csv(here(data_dir, "PIMA_train.csv")) |> mutate_if(is.character, as.factor)
AIDS_train <- read.csv(here(data_dir, "AIDS_train.csv")) |> mutate_if(is.character, as.factor)
REIN_train <- read.csv(here(data_dir, "REIN_train.csv")) |> mutate_if(is.character, as.factor)

PIMA_trn <- read.csv(here(data_dir, "PIMA_trn.csv")) |> mutate_if(is.character, as.factor)
AIDS_trn <- read.csv(here(data_dir, "AIDS_trn.csv")) |> mutate_if(is.character, as.factor)
REIN_trn <- read.csv(here(data_dir, "REIN_trn.csv")) |> mutate_if(is.character, as.factor)

PIMA_test <- read.csv(here(data_dir, "PIMA_test.csv")) |> mutate_if(is.character, as.factor)
PIMA_val  <- read.csv(here(data_dir, "PIMA_hol.csv"))  |> mutate_if(is.character, as.factor)

AIDS_test <- read.csv(here(data_dir, "AIDS_test.csv")) |> mutate_if(is.character, as.factor)
AIDS_val  <- read.csv(here(data_dir, "AIDS_hol.csv"))  |> mutate_if(is.character, as.factor)

REIN_test <- read.csv(here(data_dir, "REIN_test.csv")) |> mutate_if(is.character, as.factor)
REIN_val  <- read.csv(here(data_dir, "REIN_hol.csv"))  |> mutate_if(is.character, as.factor)



## loading synthetic data for Gaussian Copula, ctgan ,tvae, tabdiff
datasets     <- c("PIMA", "AIDS", "REIN")
synthesizers <- c("CTGAN","TVAE","GaussianCopula","TABDIFF")

# Named vector: label used in object/file names -> folder/file suffix
suffixes <- c(main = "", holdout = "_trn")

# Registry: synthesizer -> file-path builder. Add new synthesizers here.
file_path_builders <- list(
  GaussianCopula = default_file_paths,
  CTGAN          = default_file_paths,
  TVAE           = default_file_paths,
  TABDIFF        = tabdiff_file_paths
)

## Run
import_all_synthetic_datasets(
  datasets     = datasets,
  synthesizers = synthesizers,
  suffixes     = suffixes,
  n_synth      = n_synth
)

# ===============================================================================================
# impute data
# ===============================================================================================
REIN <- impute_mice(data = REIN,seed = 123)
REIN_train <- impute_mice(data = REIN_train,seed = 123)
REIN_trn <- impute_mice(data = REIN_trn,seed = 123)
REIN_test <- impute_mice(data = REIN_test,seed = 123)
REIN_val <- impute_mice(data = REIN_val,seed = 123)


#--------------------------------------------------------------------
# Loading synthetic list data
#----------------------------------------------------------------------
# PIMA synthetic data
load(here(results_dir, "pima_synthetic_list_data_chimera.Rdata"))
load(here(results_dir, "pima_synthetic_list_data_synthpop.Rdata"))
load(here(results_dir, "pima_synthetic_list_data_gaussiancopula.Rdata"))
load(here(results_dir, "pima_synthetic_list_data_ctgan.Rdata"))
load(here(results_dir, "pima_synthetic_list_data_tvae.Rdata"))
load(here(results_dir, "pima_synthetic_list_data_tabdiff.Rdata"))

load(here(results_dir, "pima_synthetic_list_data_chimera_trn.Rdata"))
load(here(results_dir, "pima_synthetic_list_data_synthpop_trn.Rdata"))
load(here(results_dir, "pima_synthetic_list_data_gaussiancopula_trn.Rdata"))
load(here(results_dir, "pima_synthetic_list_data_ctgan_trn.Rdata"))
load(here(results_dir, "pima_synthetic_list_data_tvae_trn.Rdata"))
load(here(results_dir, "pima_synthetic_list_data_tabdiff_trn.Rdata"))


# AIDS synthetic data
load(here(results_dir, "aids_synthetic_list_data_chimera.Rdata"))
load(here(results_dir, "aids_synthetic_list_data_synthpop.Rdata"))
load(here(results_dir, "aids_synthetic_list_data_gaussiancopula.Rdata"))
load(here(results_dir, "aids_synthetic_list_data_ctgan.Rdata"))
load(here(results_dir, "aids_synthetic_list_data_tvae.Rdata"))
load(here(results_dir, "aids_synthetic_list_data_tabdiff.Rdata"))

load(here(results_dir, "aids_synthetic_list_data_chimera_trn.Rdata"))
load(here(results_dir, "aids_synthetic_list_data_synthpop_trn.Rdata"))
load(here(results_dir, "aids_synthetic_list_data_gaussiancopula_trn.Rdata"))
load(here(results_dir, "aids_synthetic_list_data_ctgan_trn.Rdata"))
load(here(results_dir, "aids_synthetic_list_data_tvae_trn.Rdata"))
load(here(results_dir, "aids_synthetic_list_data_tabdiff_trn.Rdata"))

# REIN synthetic data
load(here(results_dir, "rein_synthetic_list_data_chimera.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_synthpop.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_gaussiancopula.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_ctgan.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_tvae.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_tabdiff.Rdata"))

load(here(results_dir, "rein_synthetic_list_data_chimera_trn.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_synthpop_trn.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_gaussiancopula_trn.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_ctgan_trn.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_tvae_trn.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_tabdiff_trn.Rdata"))


N_CORES <- max(1, parallel::detectCores() - 1)

# impute synthetic data
rein_synthetic_list_data_chimera <- impute_mice(rein_synthetic_list_data_chimera,seed = 123,parallel = TRUE,workers = N_CORES)
rein_synthetic_list_data_chimera_trn <- impute_mice(rein_synthetic_list_data_chimera_trn,seed = 123,parallel = TRUE,workers = N_CORES)

rein_synthetic_list_data_synthpop <- impute_mice(rein_synthetic_list_data_synthpop,seed = 123,parallel = TRUE,workers = N_CORES)
rein_synthetic_list_data_synthpop_trn <- impute_mice(rein_synthetic_list_data_synthpop_trn,seed = 123,parallel = TRUE,workers = N_CORES)

rein_synthetic_list_data_gaussiancopula <- impute_mice(rein_synthetic_list_data_gaussiancopula,seed = 123,parallel = TRUE,workers = N_CORES)
rein_synthetic_list_data_gaussiancopula_trn <- impute_mice(rein_synthetic_list_data_gaussiancopula_trn,seed = 123,parallel = TRUE,workers = N_CORES)

rein_synthetic_list_data_ctgan <- impute_mice(rein_synthetic_list_data_ctgan,seed = 123,parallel = TRUE,workers = N_CORES)
rein_synthetic_list_data_ctgan_trn <- impute_mice(rein_synthetic_list_data_ctgan_trn,seed = 123,parallel = TRUE,workers = N_CORES)

rein_synthetic_list_data_tvae <- impute_mice(rein_synthetic_list_data_tvae,seed = 123,parallel = TRUE,workers = N_CORES)
rein_synthetic_list_data_tvae_trn <- impute_mice(rein_synthetic_list_data_tvae_trn,seed = 123,parallel = TRUE,workers = N_CORES)

rein_synthetic_list_data_tabdiff <- impute_mice(rein_synthetic_list_data_tabdiff,seed = 123,parallel = TRUE,workers = N_CORES)
rein_synthetic_list_data_tabdiff_trn <- impute_mice(rein_synthetic_list_data_tabdiff_trn,seed = 123,parallel = TRUE,workers = N_CORES)


#================================================================================================
## Compute FUP evaluation metrics for each synthetic dataset and generation method
#================================================================================================
# Fit logistic regression
logit_model <- glm(Diabetes_diagnosis ~ ., data = PIMA_train, family = binomial)
coef_table <- summary(logit_model)$coefficients[-1, ] |> as.data.frame() |> rownames_to_column("Variable")
top_features <- coef_table |> arrange(desc(abs(`z value`))) |> slice(1:5) |> pull(Variable)

feature_associate_pima <- lapply(seq_along(top_features), function(i) top_features[1:i])
feature_associate_pima <- c(list(character(0)), feature_associate_pima)

# Fit logistic regression
logit_model <- glm(Homosexuality ~ ., data = AIDS_train|>dplyr::select(-c(times,Censored)), family = binomial)
coef_table <- summary(logit_model)$coefficients[-1, ] |>as.data.frame() |> rownames_to_column("Variable")
top_features <- coef_table |> arrange(desc(abs(`z value`))) |> slice(1:5) |> pull(Variable)

top_features <- c("Sex", "Hemophilia", "Intravenous_drug_use", "Race", "Symptomatic_indicator")
feature_associate_aids <- lapply(seq_along(top_features), function(i) top_features[1:i])
feature_associate_aids <- c(list(character(0)), feature_associate_aids)


# Fit logistic regression
logit_model <- glm(Diabetes ~ ., data = REIN_train|>dplyr::select(-c(times,Censored)), family = binomial)
coef_table <- summary(logit_model)$coefficients[-1, ] |>as.data.frame() |> rownames_to_column("Variable")
top_features <- coef_table |> arrange(desc(abs(`z value`))) |> slice(1:6) |> pull(Variable)

top_features <- c("Body_mass_index", "Age", "Peripheral_artery_disease", "Malignancy","Coronary_artery_disease")
feature_associate_rein <- lapply(seq_along(top_features), function(i) top_features[1:i])
feature_associate_rein <- c(list(character(0)), feature_associate_rein)

# for DCR and NNDR 
pima_rrd_thresholds <- compute_rrd_thresholds(PIMA_train,PIMA_test)
pima_rrd_thresholds_dcr <- pima_rrd_thresholds$rrd_alpha_dcr
pima_rrd_thresholds_nndr <- pima_rrd_thresholds$rrd_alpha_nndr
pima_n_rrd <- pima_rrd_thresholds$n_rrd

aids_rrd_thresholds <- compute_rrd_thresholds(AIDS_train,AIDS_test)
aids_rrd_thresholds_dcr <- aids_rrd_thresholds$rrd_alpha_dcr
aids_rrd_thresholds_nndr <- aids_rrd_thresholds$rrd_alpha_nndr
aids_n_rrd <- aids_rrd_thresholds$n_rrd

rein_rrd_thresholds <- compute_rrd_thresholds(REIN_train,REIN_test)
rein_rrd_thresholds_dcr <- rein_rrd_thresholds$rrd_alpha_dcr
rein_rrd_thresholds_nndr <- rein_rrd_thresholds$rrd_alpha_nndr
rein_n_rrd <- rein_rrd_thresholds$n_rrd

## Evaluate synthetic-data generation methods across datasets (PIMA, AIDS, REIN)

method_display_names <- list(
  chimera        = "CHIMERA",
  synthpop       = "SYNTHPOP",
  gaussiancopula = "GAUSSIANCOPULA",
  ctgan          = "CTGAN",
  tvae           = "TVAE",
  tabdiff        = "TABDIFF"
)

## Configuration

dataset_configs <- list(
  
  PIMA = list(
    X_real            = PIMA,
    X_train           = PIMA_train,
    X_test            = PIMA_test,
    trn               = PIMA_trn,
    val               = PIMA_val,
    feature_associate = feature_associate_pima,
    seeds             = seeds,
    formula           = Diabetes_diagnosis ~ .,
    target_col        = "Diabetes_diagnosis",
    modeling          = "logistic",
    is_survival       = FALSE,
    
    alpha             = pima_rrd_thresholds_dcr ,
    beta              = pima_rrd_thresholds_nndr,
    n_rrd             = pima_n_rrd,
    
    methods           = c("chimera","synthpop","ctgan","gaussiancopula","tvae","tabdiff")
  ),
  AIDS = list(
    X_real            = AIDS,
    X_train           = AIDS_train,
    X_test            = AIDS_test,
    trn               = AIDS_trn,
    val               = AIDS_val,
    feature_associate = feature_associate_aids,
    seeds             = seeds,
    formula           = Surv(times, Censored) ~ .,
    target_col        = "Homosexuality",
    modeling          = "cox",
    timevar           = "times",
    statusvar         = "Censored",
    is_survival       = TRUE,
    
    alpha             = aids_rrd_thresholds_dcr ,
    beta              = aids_rrd_thresholds_nndr,
    n_rrd             = aids_n_rrd,
    
    methods           = c("chimera","synthpop","ctgan","gaussiancopula","tvae","tabdiff")
  ),
  
  REIN = list(
    X_real            = REIN,
    X_train           = REIN_train,
    X_test            = REIN_test,
    trn               = REIN_trn,
    val               = REIN_val,
    feature_associate = feature_associate_rein,
    seeds             = seeds,
    formula           = Surv(times, Censored) ~ .,
    target_col        = "Diabetes",
    modeling          = "cox",
    timevar           = "times",
    statusvar         = "Censored",
    is_survival       = TRUE,
    
    alpha             = rein_rrd_thresholds_dcr ,
    beta              = rein_rrd_thresholds_nndr,
    n_rrd             = rein_n_rrd,
    
    methods           =c("chimera","synthpop","ctgan","gaussiancopula","tvae","tabdiff")
  )
)



## ------------------------------------------------------------------------
## Run everything with a single call
## ------------------------------------------------------------------------

#--------------------------------------------------------------------
# Loading synthetic list data
#----------------------------------------------------------------------
all_evaluation_results <- run_all_evaluations(dataset_configs)


