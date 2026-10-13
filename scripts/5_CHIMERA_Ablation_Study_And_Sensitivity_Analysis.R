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
  
  
})


# Load custom functions
source(here("scripts", "0_Methods_and_metrics_functions.R"))

args <- c("data", "results")

data_dir <- normalizePath(args[1], mustWork = TRUE)

if (!dir.exists(args[2])) {
  dir.create(args[2], recursive = TRUE)
}

results_dir <- normalizePath(args[2], mustWork = TRUE)


# ===============================================================================================
# Loading data
# ===============================================================================================
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
load(here(results_dir, "pima_synthetic_list_data_chimera_no_matching.Rdata"))
load(here(results_dir, "pima_synthetic_list_data_mice.Rdata"))

load(here(results_dir, "pima_synthetic_list_data_chimera_no_matching_trn.Rdata"))
load(here(results_dir, "pima_synthetic_list_data_mice_trn.Rdata"))


# AIDS synthetic data
load(here(results_dir, "aids_synthetic_list_data_chimera_no_matching.Rdata"))
load(here(results_dir, "aids_synthetic_list_data_mice.Rdata"))

load(here(results_dir, "aids_synthetic_list_data_chimera_no_matching_trn.Rdata"))
load(here(results_dir, "aids_synthetic_list_data_mice_trn.Rdata"))


# REIN synthetic data
load(here(results_dir, "rein_synthetic_list_data_chimera_no_matching.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_mice.Rdata"))

load(here(results_dir, "rein_synthetic_list_data_chimera_no_matching_trn.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_mice_trn.Rdata"))

N_CORES <- max(1, parallel::detectCores() - 1)

# impute synthetic data
rein_synthetic_list_data_chimera_no_matching <- impute_mice(rein_synthetic_list_data_chimera_no_matching,seed = 123,parallel = TRUE,workers = N_CORES)
rein_synthetic_list_data_chimera_no_matching_trn <- impute_mice(rein_synthetic_list_data_chimera_no_matching_trn,seed = 123,parallel = TRUE,workers = N_CORES)

rein_synthetic_list_data_mice <- impute_mice(rein_synthetic_list_data_mice,seed = 123,parallel = TRUE,workers = N_CORES)
rein_synthetic_list_data_mice_trn <- impute_mice(rein_synthetic_list_data_mice_trn,seed = 123,parallel = TRUE,workers = N_CORES)

#================================================================================================
## Compute FUP evaluation metrics for each synthetic dataset and generation method
#================================================================================================
n_synth <- 50
set.seed(123)
seeds <- sample(1:1000, n_synth)

logit_model <- glm(Diabetes_diagnosis ~ ., data = PIMA_train, family = binomial)
coef_table <- summary(logit_model)$coefficients[-1, ] |> as.data.frame() |> rownames_to_column("Variable")
top_features <- coef_table |> arrange(desc(abs(`z value`))) |> slice(1:5) |> pull(Variable)

feature_associate_pima <- lapply(seq_along(top_features), function(i) top_features[1:i])
feature_associate_pima <- c(list(character(0)), feature_associate_pima)

# Fit logistic regression
logit_model <- glm(Homosexuality ~ ., data = AIDS_train|>dplyr::select(-c(times,Censored)), family = binomial)
coef_table <- summary(logit_model)$coefficients[-1, ] |>as.data.frame() |> rownames_to_column("Variable")
top_features <- coef_table |> arrange(desc(abs(`z value`))) |> slice(1:5) |> pull(Variable)

top_features <- c("Sex", "Hemophilia", "Intravenous_drug_use", "Race", "Weight")
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
  chimera_no_matching = "CHIMERA_NO_MATCHING",
  mice                = "MICE"
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

    methods           = c("chimera_no_matching","mice")
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

    methods           = c("chimera_no_matching","mice")
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

    methods           = c("chimera_no_matching","mice")
  )
)



## ------------------------------------------------------------------------
## Run everything with a single call
## ------------------------------------------------------------------------

all_evaluation_results <- run_all_evaluations(dataset_configs)





# ============================================================
# CHIMERA sensitivity analysis across missingness rates
# Datasets: PIMA and AIDS
# ============================================================


# ------------------------------------------------------------
# Define missingness-rate settings
# ------------------------------------------------------------
p_values <- c("p_005", "p_010", "p_015", "p_020", "p_025")


# ============================================================
# PIMA: Evaluate synthetic datasets
# ============================================================

# Load synthetic datasets generated at each missingness rate
load(here(results_dir,"pima_synthetic_chimera_sensitivity_by_p.Rdata"))
load(here(results_dir,"pima_synthetic_chimera_sensitivity_by_p_trn.Rdata"))

# Initialize a list to store evaluation results for each rate
pima_results_by_p <- setNames(
  vector("list", length(p_values)),
  p_values
)


# ------------------------------------------------------------
# Evaluate all synthetic datasets for each missingness rate
# ------------------------------------------------------------

for (p in p_values) {
  
  message("\n------------------------------------------------------------")
  message("Dataset: PIMA | Missingness-rate setting: ", p)
  message("------------------------------------------------------------")
  
  # Retrieve the corresponding synthetic datasets
  X_syn_list <- pima_synthetic_chimera_sensitivity_by_p[[p]]
  X_syn_list_trn <- pima_synthetic_chimera_sensitivity_by_p_trn[[p]]
  
  # Check that both lists contain the same number of datasets
  if (length(X_syn_list) != length(X_syn_list_trn)) {
    stop(
      "The numbers of synthetic datasets and corresponding training ",
      "datasets do not match for PIMA at setting ", p, "."
    )
  }
  
  # Check that the expected number of datasets is available
  if (length(X_syn_list) == 0L) {
    stop("No synthetic datasets found for PIMA at setting ", p, ".")
  }
  
  # Evaluate each synthetic dataset
  results_list <- lapply(seq_along(X_syn_list), function(i) {
    
    message(
      "Evaluating PIMA | Setting: ", p,
      " | Synthetic dataset: ", i, "/", length(X_syn_list)
    )
    
    Assessment_function_unified(
      X_real = PIMA,
      X_train = PIMA_train,
      X_syn = X_syn_list[[i]],
      X_test = PIMA_test,
      S = 5,
      formula = Diabetes_diagnosis ~ .,
      modeling = "logistic",
      timevar = NULL,
      statusvar = NULL,
      target_col = "Diabetes_diagnosis",
      feature_associate = feature_associate_pima,
      is_survival = FALSE,
      trn = PIMA_trn,
      val = PIMA_val,
      syn_trn = X_syn_list_trn[[i]],
      seed = seeds[i],
      alpha = pima_rrd_thresholds_dcr,
      beta = pima_rrd_thresholds_nndr,
      n_rrd = pima_n_rrd,
      k = 1,
      n_workers = 8,
      run_linkability = TRUE
    )
  })
  
  # Combine evaluation results into a single data frame
  results_df <- as.data.frame(do.call(rbind, results_list))
  
  # Record the synthetic dataset identifier and missingness-rate setting
  results_df$generation <- seq_along(X_syn_list)
  results_df$missing_rate_setting <- p
  
  # Store results under the corresponding setting
  pima_results_by_p[[p]] <- results_df
  
  message("Completed PIMA evaluation for setting: ", p)
}


# ------------------------------------------------------------
# Save PIMA evaluation results
# ------------------------------------------------------------

saveRDS(
  pima_results_by_p,
  file = here(
    results_dir,
    "PIMA_assessment_results_CHIMERA_sensitivity_by_p.rds"
  )
)


# ============================================================
# AIDS: Evaluate synthetic datasets
# ============================================================

# Load synthetic datasets generated at each missingness rate
load(here(results_dir,"aids_synthetic_chimera_sensitivity_by_p.Rdata"))

load(here(results_dir,"aids_synthetic_chimera_sensitivity_by_p_trn.Rdata"))

# Initialize a list to store evaluation results for each rate
aids_results_by_p <- setNames(vector("list", length(p_values)),p_values)


# ------------------------------------------------------------
# Evaluate all synthetic datasets for each missingness rate
# ------------------------------------------------------------

for (p in p_values) {
  
  message("\n------------------------------------------------------------")
  message("Dataset: AIDS | Missingness-rate setting: ", p)
  message("------------------------------------------------------------")
  
  # Retrieve the corresponding synthetic datasets
  X_syn_list <- aids_synthetic_chimera_sensitivity_by_p[[p]]
  X_syn_list_trn <- aids_synthetic_chimera_sensitivity_by_p_trn[[p]]
  
  # Check that both lists contain the same number of datasets
  if (length(X_syn_list) != length(X_syn_list_trn)) {
    stop(
      "The numbers of synthetic datasets and corresponding training ",
      "datasets do not match for AIDS at setting ", p, "."
    )
  }
  
  # Check that the expected number of datasets is available
  if (length(X_syn_list) == 0L) {
    stop("No synthetic datasets found for AIDS at setting ", p, ".")
  }
  
  # Evaluate each synthetic dataset
  results_list <- lapply(seq_along(X_syn_list), function(i) {
    
    message(
      "Evaluating AIDS | Setting: ", p,
      " | Synthetic dataset: ", i, "/", length(X_syn_list)
    )
    
    Assessment_function_unified(
      X_real = AIDS,
      X_train = AIDS_train,
      X_syn = X_syn_list[[i]],
      X_test = AIDS_test,
      S = 5,
      formula = Surv(times, Censored) ~ .,
      modeling = "cox",
      timevar = "times",
      statusvar = "Censored",
      target_col = "Homosexuality",
      feature_associate = feature_associate_aids,
      is_survival = TRUE,
      trn = AIDS_trn,
      val = AIDS_val,
      syn_trn = X_syn_list_trn[[i]],
      seed = seeds[i],
      alpha = aids_rrd_thresholds_dcr,
      beta = aids_rrd_thresholds_nndr,
      n_rrd = aids_n_rrd,
      k = 1,
      n_workers = 8,
      run_linkability = TRUE
    )
  })
  
  # Combine evaluation results into a single data frame
  results_df <- as.data.frame(do.call(rbind, results_list))
  
  # Record the synthetic dataset identifier and missingness-rate setting
  results_df$generation <- seq_along(X_syn_list)
  results_df$missing_rate_setting <- p
  
  # Store results under the corresponding setting
  aids_results_by_p[[p]] <- results_df
  
  message("Completed AIDS evaluation for setting: ", p)
}


# ------------------------------------------------------------
# Save AIDS evaluation results
# ------------------------------------------------------------

saveRDS(
  aids_results_by_p,
  file = here(
    results_dir,
    "AIDS_assessment_results_CHIMERA_sensitivity_by_p.rds"
  )
)
