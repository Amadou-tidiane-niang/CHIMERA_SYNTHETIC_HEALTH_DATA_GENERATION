# ===============================================================================================
# INITIAL SETUP
# ===============================================================================================

# Clear workspace
rm(list = ls())

# -----------------------------------------------------------------------------------------------
# LOAD REQUIRED LIBRARIES
# -----------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(mice)
  library(gtsummary)
  library(pROC)
  library(survival)
  library(survminer)
  library(vcd)
  library(cluster)
  library(summarytools)
  library(readxl)
  library(dplyr)
  library(forestplot)
  library(grid)
  library(gridExtra)
  library(FNN)
  library(caret)
  library(fastDummies)
  library(stringr)
  library(boot)
  library(knitr)
  library(CalibrationCurves)
  
})



# Load custom utility functions
source(here("scripts", "0_functions.R"))

###############################################################################################
## REIN DATASET
###############################################################################################

# -----------------------------------------------------------------------------------------------
# LOAD ORIGINAL TRAINING AND TEST DATASETS
# -----------------------------------------------------------------------------------------------
REIN_train <- read.csv(here("data", "df_rein_train_with_na_score.csv"))

# Store the original missing-value positions
# These indices will later be reapplied to synthetic datasets
# in order to reproduce the original missing-data structure.
na_indices_initial <- which(is.na(REIN_train), arr.ind = TRUE)

# -----------------------------------------------------------------------------------------------
# LOAD SYNTHETIC DATASETS
# -----------------------------------------------------------------------------------------------

# CHIMERA synthetic datasets
load(here("results","rein_synthetic_list_data_chimera_score.Rdata"))

# SYNTHPOP synthetic datasets
load(here("results","rein_synthetic_list_data_synthpop_score.Rdata"))

# CTGAN synthetic datasets
load(here("results","rein_synthetic_list_data_ctgan_score.Rdata"))


# -----------------------------------------------------------------------------------------------
# REINTRODUCE ORIGINAL MISSING-DATA PATTERN
# -----------------------------------------------------------------------------------------------

# Reapply the original missing-value structure to CHIMERA datasets
for (i in seq_along(rein_synthetic_list_data_chimera_score)) {
  rein_synthetic_list_data_chimera_score[[i]][na_indices_initial] <- NA
}

# Reapply the original missing-value structure to CTGAN datasets
for (i in seq_along(rein_synthetic_list_data_ctgan_score)) {
  rein_synthetic_list_data_ctgan_score[[i]][na_indices_initial] <- NA
}

###############################################################################################
## MULTIPLE IMPUTATION PREPROCESSING
###############################################################################################

# -----------------------------------------------------------------------------------------------
# ORIGINAL DATASET
# -----------------------------------------------------------------------------------------------

nb_imputations <- 50

# Generate multiple imputed datasets from the original REIN dataset
results.original.imp <- preprocessing(
  REIN_train,
  nb_imputations = nb_imputations
)

# Save imputed datasets
save(
  results.original.imp,
  file = here(
    "results",
    "results.original.imp.Rdata"
  )
)

# -----------------------------------------------------------------------------------------------
# CHIMERA SYNTHETIC DATASETS
# -----------------------------------------------------------------------------------------------

# Apply preprocessing and multiple imputation to each CHIMERA dataset
results.chimera.imp <- lapply(
  rein_synthetic_list_data_chimera_score,
  preprocessing,
  nb_imputations = nb_imputations
)

# Save results
save(
  results.chimera.imp,
  file = here(
    "results",
    "results.chimera.imp.Rdata"
  )
)

# -----------------------------------------------------------------------------------------------
# SYNTHPOP SYNTHETIC DATASETS
# -----------------------------------------------------------------------------------------------

# Apply preprocessing and multiple imputation to each SYNTHPOP dataset
results.synthpop.imp <- lapply(
  rein_synthetic_list_data_synthpop_score,
  preprocessing,
  nb_imputations = nb_imputations
)

# Save results
save(
  results.synthpop.imp,
  file = here(
    "results",
    "results.synthpop.imp.Rdata"
  )
)

# -----------------------------------------------------------------------------------------------
# CTGAN SYNTHETIC DATASETS
# -----------------------------------------------------------------------------------------------

# Apply preprocessing and multiple imputation to each CTGAN dataset
results.ctgan.imp <- lapply(
  rein_synthetic_list_data_ctgan_score,
  preprocessing,
  nb_imputations = nb_imputations
)

# Save results
save(
  results.ctgan.imp,
  file = here(
    "results",
    "results.ctgan.imp.Rdata"
  )
)

###############################################################################################
## RISK SCORE ANALYSIS
## DEATH WITHIN 3 MONTHS AFTER INITIATION OF RENAL REPLACEMENT THERAPY
###############################################################################################

# -----------------------------------------------------------------------------------------------
# ORIGINAL DATASET
# -----------------------------------------------------------------------------------------------

# Compute bootstrap-based variable selection stability
# using the original imputed datasets
results.original <- Score(results.original.imp)

# Save results
save(
  results.original,
  file = here(
    "results",
    "results.original.Rdata"
  )
)

# -----------------------------------------------------------------------------------------------
# CHIMERA SYNTHETIC DATASETS
# -----------------------------------------------------------------------------------------------

# Load imputed CHIMERA datasets
load(
  here(
    "results",
    "results.chimera.imp.Rdata"
  )
)

# Initialize result container
results.chimera <- list()

# Compute variable-selection stability for each synthetic dataset
for (i in seq_along(results.chimera.imp)) {
  
  results.chimera[[i]] <- Score(results.chimera.imp[[i]])
  
  message("CHIMERA iteration ", i, " completed")
}

# Save results
save(
  results.chimera,
  file = here(
    "results",
    "results.chimera.Rdata"
  )
)

# -----------------------------------------------------------------------------------------------
# SYNTHPOP SYNTHETIC DATASETS
# -----------------------------------------------------------------------------------------------

# Initialize result container
results.synthpop <- list()

# Compute variable-selection stability for each synthetic dataset
for (i in seq_along(results.synthpop.imp)) {
  
  results.synthpop[[i]] <- Score(results.synthpop.imp[[i]])
  
  message("SYNTHPOP iteration ", i, " completed")
}

# Save results
save(
  results.synthpop,
  file = here(
    "results",
    "results.synthpop.Rdata"
  )
)

# -----------------------------------------------------------------------------------------------
# CTGAN SYNTHETIC DATASETS
# -----------------------------------------------------------------------------------------------

# Initialize result container
results.ctgan <- list()

# Compute variable-selection stability for each synthetic dataset
for (i in seq_along(results.ctgan.imp)) {
  
  results.ctgan[[i]] <- Score(results.ctgan.imp[[i]])
  
  message("CTGAN iteration ", i, " completed")
}

# Save results
save(
  results.ctgan,
  file = here(
    "results",
    "results.ctgan.Rdata"
  )
)


###############################################################################################
## VARIABLE SELECTION PERFORMANCE METRICS
## Sensitivity, Specificity, and Cohen's Kappa
###############################################################################################

# -----------------------------------------------------------------------------------------------
# CHIMERA
# -----------------------------------------------------------------------------------------------

# Reference results from the original dataset
tbl_df_real <- results.original
threshold <- 70

# Compute selection metrics for each CHIMERA synthetic dataset
chimera_results <- lapply(results.chimera, function(df_synth) {
  
  compute_selection_metrics(
    tbl_df_real,
    df_synth,
    real_col  = "nb_significative_coef.x",
    synth_col = "nb_significative_coef.y",
    threshold = threshold
  )
})

# Convert list output into a dataframe
chimera_selection_metrics_df <- do.call(
  rbind,
  lapply(chimera_results, as.data.frame)
)

# Save results
save(
  chimera_selection_metrics_df,
  file = here(
    "results",
    "chimera_selection_metrics_df.Rdata"
  )
)

# -----------------------------------------------------------------------------------------------
# SYNTHPOP
# -----------------------------------------------------------------------------------------------

# Compute selection metrics for each SYNTHPOP synthetic dataset
synthpop_results <- lapply(results.synthpop, function(df_synth) {
  
  compute_selection_metrics(
    tbl_df_real,
    df_synth,
    real_col  = "nb_significative_coef.x",
    synth_col = "nb_significative_coef.y",
    threshold = threshold
  )
})

# Convert list output into a dataframe
synthpop_selection_metrics_df <- do.call(
  rbind,
  lapply(synthpop_results, as.data.frame)
)

# Save results
save(
  synthpop_selection_metrics_df,
  file = here(
    "results",
    "synthpop_selection_metrics_df.Rdata"
  )
)

# -----------------------------------------------------------------------------------------------
# CTGAN
# -----------------------------------------------------------------------------------------------

# Compute selection metrics for each CTGAN synthetic dataset
ctgan_results <- lapply(results.ctgan, function(df_synth) {
  
  compute_selection_metrics(
    tbl_df_real,
    df_synth,
    real_col  = "nb_significative_coef.x",
    synth_col = "nb_significative_coef.y",
    threshold = threshold
  )
})

# Convert list output into a dataframe
ctgan_selection_metrics_df <- do.call(
  rbind,
  lapply(ctgan_results, as.data.frame)
)

# Save results
save(
  ctgan_selection_metrics_df,
  file = here(
    "data",
    "ctgan_selection_metrics_df.Rdata"
  )
)

###############################################################################################
## ADD GLOBAL PERFORMANCE SCORE
###############################################################################################

# -----------------------------------------------------------------------------------------------
# LOAD GLOBAL FIDELITY SCORES
# -----------------------------------------------------------------------------------------------
load(
  here(
    "results",
    "REIN_Metrics_T1_T9_chimera_score.Rdata"
  )
)

load(
  here(
    "results",
    "REIN_Metrics_T1_T9_synthpop_score.Rdata"
  )
)

load(
  here(
    "results",
    "REIN_Metrics_T1_T9_ctgan_score.Rdata"
  )
)


# -----------------------------------------------------------------------------------------------
# ATTACH GLOBAL DISTANCE METRIC
# -----------------------------------------------------------------------------------------------

# The distance corresponds to the global ranking criterion
# based on the aggregated fidelity / utility / privacy score.

chimera_selection_metrics_df$distance <-
  REIN_Metrics_T1_T9_chimera$L_total

save(
  chimera_selection_metrics_df,
  file = here(
    "results",
    "chimera_selection_metrics_df.Rdata"
  )
)

synthpop_selection_metrics_df$distance <-
  REIN_Metrics_T1_T9_synthpop$L_total

save(
  synthpop_selection_metrics_df,
  file = here(
    "results",
    "synthpop_selection_metrics_df.Rdata"
  )
)


ctgan_selection_metrics_df$distance <-
  REIN_Metrics_T1_T9_ctgan$L_total

save(
  ctgan_selection_metrics_df,
  file = here(
    "results",
    "ctgan_selection_metrics_df.Rdata"
  )
)


###############################################################################################
## BEST SYNTHETIC DATASET SELECTION
###############################################################################################

# Select the dataset minimizing the global distance score
# (lower distance = better overall performance)

best_row_chimera  <- which.min(chimera_selection_metrics_df$distance)
best_row_synthpop <- which.min(synthpop_selection_metrics_df$distance)
best_row_ctgan    <- which.min(ctgan_selection_metrics_df$distance)

# -----------------------------------------------------------------------------------------------
# EXTRACT BEST METRICS
# -----------------------------------------------------------------------------------------------

chimera_best_metrics  <- chimera_selection_metrics_df[best_row_chimera, ]
synthpop_best_metrics <- synthpop_selection_metrics_df[best_row_synthpop, ]
ctgan_best_metrics    <- ctgan_selection_metrics_df[best_row_ctgan, ]

# -----------------------------------------------------------------------------------------------
# EXTRACT BEST SYNTHETIC DATASETS
# -----------------------------------------------------------------------------------------------

best_chimera  <- results.chimera[[best_row_chimera]]
best_synthpop <- results.synthpop[[best_row_synthpop]]
best_ctgan    <- results.ctgan[[best_row_ctgan]]

###############################################################################################
## LOAD BEST IMPUTED DATASETS
###############################################################################################

# Original imputed datasets
df.imp.original <- results.original.imp

# Best CHIMERA synthetic dataset
df.imp.chimera <- results.chimera.imp[[best_row_chimera]]

# Best SYNTHPOP synthetic dataset
df.imp.synthpop <- results.synthpop.imp[[best_row_synthpop]]

# Best SYNTHPOP synthetic dataset
df.imp.ctgan <- results.synthpop.imp[[best_row_ctgan]]

###############################################################################################
## VARIABLE SELECTION
###############################################################################################

# -----------------------------------------------------------------------------------------------
# ORIGINAL DATASET
# -----------------------------------------------------------------------------------------------
df_test <- read.csv(here("data", "df_rein_test_with_na_score.csv"))
df_test <- df_test|>mutate_if(is.character,as.factor)

# --- Initial MICE setup
ini <- mice(df_test, maxit = 0)
meth <- ini$method
pred <- ini$predictorMatrix
meth[names(df_test)[sapply(df_test, is.numeric)]] <- "pmm"

# --- Impute missing data
df_test <- complete(mice(df_test, m = 1, maxit = 5, seed = 123,method = meth, predictorMatrix = pred, printFlag = TRUE))
df_test <- df_test|>
  rename(
    Censorship    = Censored,
    Follow_up_time = times
  )

df_test$death <- as.factor(ifelse(df_test$Follow_up_time < 90 & df_test$Censorship == 1, "Yes", "No"))

# Select variables identified in more than 70% of bootstrap samples
df.nbSignCoef.original <-
  results.original |>
  filter(nb_significative_coef > threshold)

# Clean variable names for model fitting
name_variable_original <-
  clean_variable_names(
    df.nbSignCoef.original$facteur_risque,
    df.imp.original[[1]]
  )

# Compute pooled AUC and confidence intervals
auc_original <- compute_auc_ci(
  df.imp.original,
  df_test,
  outcome     = "death",
  method_name = "REAL",
  name_variable = name_variable_original
)

# -----------------------------------------------------------------------------------------------
# CHIMERA
# -----------------------------------------------------------------------------------------------

df.nbSignCoef.chimera <-
  results.chimera[[best_row_chimera]] |>
  filter(nb_significative_coef > threshold)

name_variable_chimera <-
  clean_variable_names(
    df.nbSignCoef.chimera$facteur_risque,
    df.imp.chimera[[1]]
  )

auc_chimera <- compute_auc_ci(
  df.imp.chimera,
  df_test,
  outcome       = "death",
  method_name   = "CHIMERA",
  name_variable = name_variable_chimera
)

# -----------------------------------------------------------------------------------------------
# SYNTHPOP
# -----------------------------------------------------------------------------------------------

df.nbSignCoef.synthpop <-
  results.synthpop[[best_row_synthpop]] |>
  filter(nb_significative_coef > threshold)

name_variable_synthpop <-
  clean_variable_names(
    df.nbSignCoef.synthpop$facteur_risque,
    df.imp.synthpop[[1]]
  )

auc_synthpop <- compute_auc_ci(
  df.imp.synthpop,
  df_test,
  outcome       = "death",
  method_name   = "SYNTHPOP",
  name_variable = name_variable_synthpop
)

# -----------------------------------------------------------------------------------------------
# CTGAN
# -----------------------------------------------------------------------------------------------

df.nbSignCoef.ctgan <-
  results.ctgan[[best_row_ctgan]] |>
  filter(nb_significative_coef > threshold)

name_variable_ctgan <-
  clean_variable_names(
    df.nbSignCoef.ctgan$facteur_risque,
    df.imp.ctgan[[1]]
  )

auc_ctgan <- compute_auc_ci(
  df.imp.ctgan,
  df_test,
  outcome       = "death",
  method_name   = "CTGAN",
  name_variable = name_variable_ctgan
)

###############################################################################################
## COMBINE AUC RESULTS
###############################################################################################

# Merge all AUC results
auc_all <- bind_rows(
  auc_original,
  auc_chimera,
  auc_synthpop,
  auc_ctgan
)

# Define plotting order
auc_all$source <- factor(
  toupper(auc_all$source),
  levels = c("REAL", "CHIMERA", "SYNTHPOP", "CTGAN")
)

# Final AUC results
auc_final_model_rein <- auc_all

# Save results
save(
  auc_final_model_rein,
  file = here(
    "results",
    "auc_final_model_rein.Rdata"
  )
)

###############################################################################################
## ROC CURVE COMPUTATION
###############################################################################################

# -----------------------------------------------------------------------------------------------
# ORIGINAL DATASET
# -----------------------------------------------------------------------------------------------

# Initialize pooled predictions
predictions_pooled <- rep(0, nrow(df_test))

for (i in 1:nb_imputations) {
  
  # Logistic regression model
  formula <- as.formula(
    paste(
      "death ~",
      paste(name_variable_original, collapse = " + ")
    )
  )
  
  reg.log <- glm(
    formula,
    data   = df.imp.original[[i]],
    family = "binomial"
  )
  
  # Accumulate predictions
  predictions_pooled <-
    predictions_pooled +
    predict(reg.log, newdata = df_test, type = "response")
}

# Average predictions across imputations
predictions_pooled <- predictions_pooled / nb_imputations

# Compute pooled ROC curve
rein.roc.original <- roc(df_test$death, predictions_pooled)

# Save ROC object
save(
  rein.roc.original,
  file = here(
    "results",
    "rein.roc.original.Rdata"
  )
)

# -----------------------------------------------------------------------------------------------
# CHIMERA
# -----------------------------------------------------------------------------------------------

predictions_pooled <- rep(0, nrow(df_test))

for (i in 1:nb_imputations) {
  
  formula <- as.formula(
    paste(
      "death ~",
      paste(name_variable_chimera, collapse = " + ")
    )
  )
  
  reg.log <- glm(
    formula,
    data   = df.imp.chimera[[i]],
    family = "binomial"
  )
  
  predictions_pooled <-
    predictions_pooled +
    predict(reg.log, newdata = df_test, type = "response")
}

predictions_pooled <- predictions_pooled / nb_imputations

rein.roc.chimera <- roc(df_test$death, predictions_pooled)

save(
  rein.roc.chimera,
  file = here(
    "results",
    "rein.roc.chimera.Rdata"
  )
)

# -----------------------------------------------------------------------------------------------
# SYNTHPOP
# -----------------------------------------------------------------------------------------------

predictions_pooled <- rep(0, nrow(df_test))

for (i in 1:nb_imputations) {
  
  formula <- as.formula(
    paste(
      "death ~",
      paste(name_variable_synthpop, collapse = " + ")
    )
  )
  
  reg.log <- glm(
    formula,
    data   = df.imp.synthpop[[i]],
    family = "binomial"
  )
  
  predictions_pooled <-
    predictions_pooled +
    predict(reg.log, newdata = df_test, type = "response")
}

predictions_pooled <- predictions_pooled / nb_imputations

rein.roc.synthpop <- roc(df_test$death, predictions_pooled)

save(
  rein.roc.synthpop,
  file = here(
    "results",
    "rein.roc.synthpop.Rdata"
  )
)

# -----------------------------------------------------------------------------------------------
# CTGAN
# -----------------------------------------------------------------------------------------------

predictions_pooled <- rep(0, nrow(df_test))

for (i in 1:nb_imputations) {
  
  formula <- as.formula(
    paste(
      "death ~",
      paste(name_variable_ctgan, collapse = " + ")
    )
  )
  
  reg.log <- glm(
    formula,
    data   = df.imp.ctgan[[i]],
    family = "binomial"
  )
  
  predictions_pooled <-
    predictions_pooled +
    predict(reg.log, newdata = df_test, type = "response")
}

predictions_pooled <- predictions_pooled / nb_imputations

rein.roc.ctgan <- roc(df_test$death, predictions_pooled)

save(
  rein.roc.ctgan,
  file = here(
    "results",
    "rein.roc.ctgan.Rdata"
  )
)

###############################################################################################
## AUC EVALUATION ACROSS ALL SYNTHETIC DATASETS
###############################################################################################

# -----------------------------------------------------------------------------------------------
# CHIMERA
# -----------------------------------------------------------------------------------------------

synthetic_data_list <- lapply(
  rein_synthetic_list_data_chimera_score,
  preprocessing
)

df.nbSignCoef.chimera <-
  results.chimera[[best_row_chimera]] |>
  filter(nb_significative_coef > threshold)

name_variable_chimera <-
  clean_variable_names(
    df.nbSignCoef.chimera$facteur_risque,
    df.imp.chimera[[1]]
  )

auc_chimera <- compute_auc_ci(
  synthetic_data_list,
  df_test,
  outcome       = "death",
  method_name   = "CHIMERA",
  name_variable = name_variable_chimera
)

# -----------------------------------------------------------------------------------------------
# SYNTHPOP
# -----------------------------------------------------------------------------------------------

synthetic_data_list <- lapply(
  rein_synthetic_list_data_synthpop_score,
  preprocessing
)

df.nbSignCoef.synthpop <-
  results.synthpop[[best_row_synthpop]] |>
  filter(nb_significative_coef > threshold)

name_variable_synthpop <-
  clean_variable_names(
    df.nbSignCoef.synthpop$facteur_risque,
    df.imp.synthpop[[1]]
  )

auc_synthpop <- compute_auc_ci(
  synthetic_data_list,
  df_test,
  outcome       = "death",
  method_name   = "SYNTHPOP",
  name_variable = name_variable_synthpop
)

# -----------------------------------------------------------------------------------------------
# CTGAN
# -----------------------------------------------------------------------------------------------

synthetic_data_list <- lapply(
  rein_synthetic_list_data_ctgan_score,
  preprocessing
)

df.nbSignCoef.ctgan <-
  results.ctgan[[best_row_ctgan]] |>
  filter(nb_significative_coef > threshold)

name_variable_ctgan <-
  clean_variable_names(
    df.nbSignCoef.ctgan$facteur_risque,
    df.imp.ctgan[[1]]
  )

auc_ctgan <- compute_auc_ci(
  synthetic_data_list,
  df_test,
  outcome       = "death",
  method_name   = "CTGAN",
  name_variable = name_variable_ctgan
)

###############################################################################################
## FINAL SYNTHETIC AUC RESULTS
###############################################################################################

# Merge synthetic AUC results
auc_all <- bind_rows(
  auc_chimera,
  auc_synthpop,
  auc_ctgan
)

# Define plotting order
auc_all$source <- factor(
  toupper(auc_all$source),
  levels = c("CHIMERA", "SYNTHPOP", "CTGAN")
)

# Final synthetic AUC summary
auc_final_model_rein_syn <- auc_all

# Save results
save(
  auc_final_model_rein_syn,
  file = here(
    "results",
    "auc_final_model_rein_syn.Rdata"
  )
)


synthetic_data_list <- rein_synthetic_list_data_chimera_score
res.chimera.cal <- compute_calibration_synth(results.chimera, df_test, "death", synthetic_data_list)

# Save results
save(
  res.chimera.cal,
  file = here(
    "results",
    "res.chimera.cal.Rdata"
  )
)

# synthpop
synthetic_data_list <- rein_synthetic_list_data_synthpop_score
res.synthpop.cal <- compute_calibration_synth(results.synthpop, df_test, "death", synthetic_data_list)
# Save results
save(
  res.synthpop.cal,
  file = here(
    "results",
    "res.synthpop.cal.Rdata"
  )
)
# ctgan
synthetic_data_list <- rein_synthetic_list_data_ctgan_score
res.ctgan.cal <- compute_calibration_synth(results.ctgan, df_test, "death", synthetic_data_list)
# Save results
save(
  res.ctgan.cal,
  file = here(
    "results",
    "res.ctgan.cal.Rdata"
  )
)


# Original
res <- compute_calibration(results.original, df_test, "death",df.imp.original)
res.original.cal <- res
# Save results
save(
  res.original.cal,
  file = here(
    "results",
    "res.original.cal.Rdata"
  )
)



