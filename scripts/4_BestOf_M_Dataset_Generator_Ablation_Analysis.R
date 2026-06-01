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
  
})


# Load custom functions
source(here("scripts", "0_functions.R"))


# ===============================================================================================
# DATA LOADING : PIMA
# ===============================================================================================
PIMA <- read.csv(here("data", "PIMA.csv"))
PIMA <- PIMA|> mutate_if(is.character,as.factor)

# MICE
n_synth <- 50
set.seed(123)
seeds <- sample(1:1000, n_synth)

pima_synthetic_list_data_mice <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = PIMA,
      missing_rate = NULL,
      iterative_masking = FALSE,
      survival_data = FALSE,
      time_variable = NULL,
      event_variable = NULL,
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      continuous_method = "pmm",
      verbose = TRUE
    )
    
  })["elapsed"]
  
  pima_synthetic_list_data_mice[[i]] <- synthetic_data[["synthetic_data"]]
}

# -----------------------------------------------------------------------------------------------
# EXECUTION TIME SUMMARY
# -----------------------------------------------------------------------------------------------
execution_time_mean  <- mean(execution_times) # 2.988 s
execution_time_upper <- execution_time_mean + sd(execution_times) # 2.499 s
execution_time_lower <- execution_time_mean - sd(execution_times) # 3.478

# Save synthetic datasets
save(pima_synthetic_list_data_mice,
     file = here("results", "pima_synthetic_list_data_mice.Rdata")
)


# CHIMERA NO MATCHING
pima_synthetic_list_data_chimera_no_matching <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = PIMA,
      missing_rate = 0.25,
      iterative_masking = TRUE,
      survival_data = FALSE,
      time_variable = NULL,
      event_variable = NULL,
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      continuous_method = "pmm",
      verbose = TRUE
    )
    
  })["elapsed"]
  
  pima_synthetic_list_data_chimera_no_matching[[i]] <- synthetic_data[["synthetic_data"]]
}

# -----------------------------------------------------------------------------------------------
# EXECUTION TIME SUMMARY
# -----------------------------------------------------------------------------------------------
execution_time_mean  <- mean(execution_times) # 2.988 s
execution_time_upper <- execution_time_mean + sd(execution_times) # 2.499 s
execution_time_lower <- execution_time_mean - sd(execution_times) # 3.478

# Save synthetic datasets
save(pima_synthetic_list_data_chimera_no_matching,
     file = here("results", "pima_synthetic_list_data_chimera_no_matching.Rdata")
)


# ===============================================================================================
# DATA LOADING : AIDS
# ===============================================================================================
AIDS <- read.csv(here("data", "AIDS.csv"))
AIDS <- AIDS|> mutate_if(is.character,as.factor)

aids_synthetic_list_data_mice <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = AIDS,
      missing_rate = NULL,
      iterative_masking = FALSE,
      survival_data = TRUE,
      time_variable = "times",
      event_variable = "Censored",
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      continuous_method = "pmm",
      verbose = TRUE
    )
    
  })["elapsed"]
  
  aids_synthetic_list_data_mice[[i]] <- synthetic_data[["synthetic_data"]]
}

# -----------------------------------------------------------------------------------------------
# EXECUTION TIME SUMMARY
# -----------------------------------------------------------------------------------------------
execution_time_mean  <- mean(execution_times) # 2.988 s
execution_time_upper <- execution_time_mean + sd(execution_times) # 2.499 s
execution_time_lower <- execution_time_mean - sd(execution_times) # 3.478

# Save synthetic datasets
save(aids_synthetic_list_data_mice,
     file = here("results", "aids_synthetic_list_data_mice.Rdata")
)


# CHIMERA NO MATCHING
aids_synthetic_list_data_chimera_no_matching <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = AIDS,
      missing_rate = 0.15,
      iterative_masking = TRUE,
      survival_data = TRUE,
      time_variable = "times",
      event_variable = "Censored",
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      continuous_method = "pmm",
      verbose = TRUE
    )
    
  })["elapsed"]
  
  aids_synthetic_list_data_chimera_no_matching[[i]] <- synthetic_data[["synthetic_data"]]
}

# -----------------------------------------------------------------------------------------------
# EXECUTION TIME SUMMARY
# -----------------------------------------------------------------------------------------------
execution_time_mean  <- mean(execution_times) # 2.988 s
execution_time_upper <- execution_time_mean + sd(execution_times) # 2.499 s
execution_time_lower <- execution_time_mean - sd(execution_times) # 3.478

# Save synthetic datasets
save(aids_synthetic_list_data_chimera_no_matching,
     file = here("results", "aids_synthetic_list_data_chimera_no_matching.Rdata")
)


# ===============================================================================================
# DATA LOADING : REIN
# ===============================================================================================
df_rein <- read.csv(here("data", "df_rein_with_na.csv"))
df_rein <- df_rein|> mutate_if(is.character,as.factor)

rein_synthetic_list_data_mice <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = df_rein,
      missing_rate = NULL,
      iterative_masking = FALSE,
      survival_data = TRUE,
      time_variable = "times",
      event_variable = "Censored",
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      continuous_method = "pmm",
      verbose = TRUE
    )
    
  })["elapsed"]
  
  rein_synthetic_list_data_mice[[i]] <- synthetic_data[["synthetic_data"]]
}

# -----------------------------------------------------------------------------------------------
# EXECUTION TIME SUMMARY
# -----------------------------------------------------------------------------------------------
execution_time_mean  <- mean(execution_times) # 2.988 s
execution_time_upper <- execution_time_mean + sd(execution_times) # 2.499 s
execution_time_lower <- execution_time_mean - sd(execution_times) # 3.478

# Save synthetic datasets
save(rein_synthetic_list_data_mice,
     file = here("results", "rein_synthetic_list_data_mice.Rdata")
)


# CHIMERA NO MATCHING
rein_synthetic_list_data_chimera_no_matching <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = df_rein,
      missing_rate = 0.30,
      iterative_masking = TRUE,
      survival_data = TRUE,
      time_variable = "times",
      event_variable = "Censored",
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      continuous_method = "pmm",
      verbose = TRUE
    )
    
  })["elapsed"]
  
  rein_synthetic_list_data_chimera_no_matching[[i]] <- synthetic_data[["synthetic_data"]]
}

# -----------------------------------------------------------------------------------------------
# EXECUTION TIME SUMMARY
# -----------------------------------------------------------------------------------------------
execution_time_mean  <- mean(execution_times) # 2.988 s
execution_time_upper <- execution_time_mean + sd(execution_times) # 2.499 s
execution_time_lower <- execution_time_mean - sd(execution_times) # 3.478

# Save synthetic datasets
save(rein_synthetic_list_data_chimera_no_matching,
     file = here("results", "rein_synthetic_list_data_chimera_no_matching.Rdata")
)



#================================================================================================
# Select the best dataset use for downstream analyses
#================================================================================================
X_real <- read.csv(here("data", "PIMA.csv"))
X_real <- X_real|> mutate_if(is.character,as.factor)

trn <- read.csv(here("data", "PIMA_trn.csv"))
trn <- trn|> mutate_if(is.character,as.factor)
val <- read.csv(here("data", "PIMA_hol.csv"))
val <- val|> mutate_if(is.character,as.factor)

# Fit logistic regression
logit_model <- glm(Diabetes_diagnosis ~ ., data = X_real, family = binomial)

coef_table <- summary(logit_model)$coefficients[-1, ] %>%
  as.data.frame() %>%
  rownames_to_column("Variable")

# Select top variables based on |z-value|
top_features <- coef_table %>%
  arrange(desc(abs(`z value`))) %>%
  slice(1:5) %>%
  pull(Variable)


# GENERATE FEATURE COMBINATIONS
feature_associate <- lapply(seq_along(top_features), function(i) {
  top_features[1:i]
})

# Add empty set explicitly (useful for AIR baseline)
feature_associate <- c(list(character(0)), feature_associate)


# mice
load(here("results", "pima_synthetic_list_data_mice.Rdata"))
results_all <- lapply(seq_along(pima_synthetic_list_data_mice), function(i) {
  Assessment_function_unified_Ablation_Analysis(
    X_real = X_real,
    X_syn  = pima_synthetic_list_data_mice[[i]],
    S = 5,
    seed = seeds[i],
    formula = Diabetes_diagnosis ~ .,
    modeling = "logistic",
    timevar = NULL,
    statusvar = NULL,
    missing_rate = NULL,
    iterative_masking=FALSE,
    target_col = "Diabetes_diagnosis",
    feature_associate = feature_associate,
    is_survival = FALSE,
    trn = trn,
    val = val
  )
})

pima_results_all_mice <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(pima_results_all_mice,
     file = here("results", "pima_results_all_mice_min_max.Rdata")
)


# metrics
metrics <- paste0("T", 1:9)

#target vector
T_star <- c(
  T1 = 0,
  T2 = 0,
  T3 = 0.5,
  T4 = 0,
  T5 = 0,
  T6 = 0,
  T7 = 1,
  T8 = 0.5,
  T9 = 0
)

PIMA_Metrics_T1_T9_mice <- pima_results_all_mice|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(PIMA_Metrics_T1_T9_mice[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(PIMA_Metrics_T1_T9_mice[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

PIMA_Metrics_T1_T9_mice$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
PIMA_Metrics_T1_T9_mice$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
PIMA_Metrics_T1_T9_mice$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)

# Poids (modifiable)
w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

PIMA_Metrics_T1_T9_mice$L_total <-
  w_fid  * PIMA_Metrics_T1_T9_mice$L_fid +
  w_util * PIMA_Metrics_T1_T9_mice$L_util +
  w_conf * PIMA_Metrics_T1_T9_mice$L_conf

PIMA_Metrics_T1_T9_mice$dataset <- paste0("syn_", seq_len(n_synth))


# Save synthetic datasets
save(PIMA_Metrics_T1_T9_mice,
     file = here("results", "PIMA_Metrics_T1_T9_mice.Rdata")
)


# mice
load(here("results", "pima_synthetic_list_data_chimera_no_matching.Rdata"))
results_all <- lapply(seq_along(pima_synthetic_list_data_chimera_no_matching), function(i) {
  Assessment_function_unified_Ablation_Analysis(
    X_real = X_real,
    X_syn  = pima_synthetic_list_data_chimera_no_matching[[i]],
    S = 5,
    seed = seeds[i],
    formula = Diabetes_diagnosis ~ .,
    modeling = "logistic",
    timevar = NULL,
    statusvar = NULL,
    missing_rate = 0.25,
    iterative_masking=TRUE,
    target_col = "Diabetes_diagnosis",
    feature_associate = feature_associate,
    is_survival = FALSE,
    trn = trn,
    val = val
  )
})

pima_results_all_chimera_no_matching <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(pima_results_all_chimera_no_matching,
     file = here("results", "pima_results_all_chimera_no_matching_min_max.Rdata")
)


# metrics
metrics <- paste0("T", 1:9)

#target vector
T_star <- c(
  T1 = 0,
  T2 = 0,
  T3 = 0.5,
  T4 = 0,
  T5 = 0,
  T6 = 0,
  T7 = 1,
  T8 = 0.5,
  T9 = 0
)

PIMA_Metrics_T1_T9_chimera_no_matching <- pima_results_all_chimera_no_matching|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(PIMA_Metrics_T1_T9_chimera_no_matching[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(PIMA_Metrics_T1_T9_chimera_no_matching[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

PIMA_Metrics_T1_T9_chimera_no_matching$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
PIMA_Metrics_T1_T9_chimera_no_matching$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
PIMA_Metrics_T1_T9_chimera_no_matching$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)

# Poids (modifiable)
w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

PIMA_Metrics_T1_T9_chimera_no_matching$L_total <-
  w_fid  * PIMA_Metrics_T1_T9_chimera_no_matching$L_fid +
  w_util * PIMA_Metrics_T1_T9_chimera_no_matching$L_util +
  w_conf * PIMA_Metrics_T1_T9_chimera_no_matching$L_conf

PIMA_Metrics_T1_T9_chimera_no_matching$dataset <- paste0("syn_", seq_len(n_synth))

# Save synthetic datasets
save(PIMA_Metrics_T1_T9_chimera_no_matching,
     file = here("results", "PIMA_Metrics_T1_T9_chimera_no_matching.Rdata")
)


#================================================================================================
# Select the best dataset use for downstream analyses
#================================================================================================
X_real <- read.csv(here("data", "AIDS.csv"))
X_real <- X_real|> mutate_if(is.character,as.factor)

trn <- read.csv(here("data", "AIDS_trn.csv"))
trn <- trn|> mutate_if(is.character,as.factor)
val <- read.csv(here("data", "AIDS_hol.csv"))
val <- val|> mutate_if(is.character,as.factor)

# Fit logistic regression
logit_model <- glm(Homosexuality ~ ., data = X_real|>dplyr::select(-c(times,Censored)), family = binomial)

coef_table <- summary(logit_model)$coefficients[-1, ] %>%
  as.data.frame() %>%
  rownames_to_column("Variable")

# Select top variables based on |z-value|
top_features <- coef_table %>%
  arrange(desc(abs(`z value`))) %>%
  slice(1:5) %>%
  pull(Variable)

top_features <- c("Sex","Hemophilia","Intravenous_drug_use","Race","Weight")
# GENERATE FEATURE COMBINATIONS
feature_associate <- lapply(seq_along(top_features), function(i) {
  top_features[1:i]
})

# Add empty set explicitly (useful for AIR baseline)
feature_associate <- c(list(character(0)), feature_associate)


# mice
load(here("results", "aids_synthetic_list_data_mice.Rdata"))
results_all <- lapply(seq_along(aids_synthetic_list_data_mice), function(i) {
  Assessment_function_unified_Ablation_Analysis(
    X_real = X_real,
    X_syn  = aids_synthetic_list_data_mice[[i]],
    S = 5,
    seed = seeds[i],
    formula = Surv(times, Censored) ~ .,
    modeling = "cox",
    timevar = "times",
    statusvar = "Censored",
    missing_rate = NULL,
    iterative_masking=FALSE,
    target_col = "Homosexuality",
    feature_associate = feature_associate,
    is_survival = TRUE,
    trn = trn,
    val = val
  )
})

aids_results_all_mice <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(aids_results_all_mice,
     file = here("results", "aids_results_all_mice_min_max.Rdata")
)


# metrics
metrics <- paste0("T", 1:9)

#target vector
T_star <- c(
  T1 = 0,
  T2 = 0,
  T3 = 0.5,
  T4 = 0,
  T5 = 0,
  T6 = 0,
  T7 = 1,
  T8 = 0.5,
  T9 = 0
)

AIDS_Metrics_T1_T9_mice <- aids_results_all_mice|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(AIDS_Metrics_T1_T9_mice[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(AIDS_Metrics_T1_T9_mice[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

AIDS_Metrics_T1_T9_mice$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
AIDS_Metrics_T1_T9_mice$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
AIDS_Metrics_T1_T9_mice$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)

# Poids (modifiable)
w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

AIDS_Metrics_T1_T9_mice$L_total <-
  w_fid  * AIDS_Metrics_T1_T9_mice$L_fid +
  w_util * AIDS_Metrics_T1_T9_mice$L_util +
  w_conf * AIDS_Metrics_T1_T9_mice$L_conf

AIDS_Metrics_T1_T9_mice$dataset <- paste0("syn_", seq_len(n_synth))

# Save synthetic datasets
save(AIDS_Metrics_T1_T9_mice,
     file = here("results", "AIDS_Metrics_T1_T9_mice.Rdata")
)


# chimera no matching
load(here("results", "aids_synthetic_list_data_chimera_no_matching.Rdata"))
results_all <- lapply(seq_along(aids_synthetic_list_data_chimera_no_matching), function(i) {
  Assessment_function_unified_Ablation_Analysis(
    X_real = X_real,
    X_syn  = aids_synthetic_list_data_chimera_no_matching[[i]],
    S = 5,
    seed = seeds[i],
    formula = Surv(times, Censored) ~ .,
    modeling = "cox",
    timevar = "times",
    statusvar = "Censored",
    missing_rate = 0.15,
    iterative_masking=TRUE,
    target_col = "Homosexuality",
    feature_associate = feature_associate,
    is_survival = TRUE,
    trn = trn,
    val = val
  )
})

aids_results_all_chimera_no_matching <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(aids_results_all_chimera_no_matching,
     file = here("results", "aids_results_all_chimera_no_matching_min_max.Rdata")
)


# metrics
metrics <- paste0("T", 1:9)

#target vector
T_star <- c(
  T1 = 0,
  T2 = 0,
  T3 = 0.5,
  T4 = 0,
  T5 = 0,
  T6 = 0,
  T7 = 1,
  T8 = 0.5,
  T9 = 0
)

AIDS_Metrics_T1_T9_chimera_no_matching <- aids_results_all_chimera_no_matching|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(AIDS_Metrics_T1_T9_chimera_no_matching[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(AIDS_Metrics_T1_T9_chimera_no_matching[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

AIDS_Metrics_T1_T9_chimera_no_matching$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
AIDS_Metrics_T1_T9_chimera_no_matching$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
AIDS_Metrics_T1_T9_chimera_no_matching$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)


w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

AIDS_Metrics_T1_T9_chimera_no_matching$L_total <-
  w_fid  * AIDS_Metrics_T1_T9_chimera_no_matching$L_fid +
  w_util * AIDS_Metrics_T1_T9_chimera_no_matching$L_util +
  w_conf * AIDS_Metrics_T1_T9_chimera_no_matching$L_conf


AIDS_Metrics_T1_T9_chimera_no_matching$dataset <- paste0("syn_", seq_len(n_synth))

# Save synthetic datasets
save(AIDS_Metrics_T1_T9_chimera_no_matching,
     file = here("results", "AIDS_Metrics_T1_T9_chimera_no_matching.Rdata")
)


#================================================================================================
# Select the best dataset use for downstream analyses
#================================================================================================
X_real <- read.csv(here("data", "df_rein_without_na.csv"))
X_real <- X_real|> mutate_if(is.character,as.factor)

trn <- read.csv(here("data", "REIN_trn.csv"))
trn <- trn|> mutate_if(is.character,as.factor)
val <- read.csv(here("data", "REIN_hol.csv"))
val <- val|> mutate_if(is.character,as.factor)

# Fit logistic regression
logit_model <- glm(Diabetes ~ ., data = X_real|>dplyr::select(-c(times,Censored)), family = binomial)

coef_table <- summary(logit_model)$coefficients[-1, ] %>%
  as.data.frame() %>%
  rownames_to_column("Variable")

# Select top variables based on |z-value|
top_features <- coef_table %>%
  arrange(desc(abs(`z value`))) %>%
  slice(1:5) %>%
  pull(Variable)

top_features <- c("Body_mass_index","Age","Peripheral_artery_disease","Peripheral_artery_disease","Malignancy")

# GENERATE FEATURE COMBINATIONS
feature_associate <- lapply(seq_along(top_features), function(i) {
  top_features[1:i]
})

# Add empty set explicitly (useful for AIR baseline)
feature_associate <- c(list(character(0)), feature_associate)



# mice
load(here("results", "rein_synthetic_list_data_mice.Rdata"))
results_all <- lapply(seq_along(rein_synthetic_list_data_mice), function(i) {
  Assessment_function_unified_Ablation_Analysis(
    X_real = X_real,
    X_syn  = rein_synthetic_list_data_mice[[i]],
    S = 5,
    seed = seeds[i],
    formula = Surv(times, Censored) ~ .,
    modeling = "cox",
    timevar = "times",
    statusvar = "Censored",
    missing_rate = NULL,
    iterative_masking=FALSE,
    target_col = "Diabetes",
    feature_associate = feature_associate,
    is_survival = TRUE,
    trn = trn,
    val = val
  )
})

rein_results_all_mice <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(rein_results_all_mice,
     file = here("results", "rein_results_all_mice_min_max.Rdata")
)


# metrics
metrics <- paste0("T", 1:9)

#target vector
T_star <- c(
  T1 = 0,
  T2 = 0,
  T3 = 0.5,
  T4 = 0,
  T5 = 0,
  T6 = 0,
  T7 = 1,
  T8 = 0.5,
  T9 = 0
)

REIN_Metrics_T1_T9_mice <- rein_results_all_mice|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(REIN_Metrics_T1_T9_mice[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(REIN_Metrics_T1_T9_mice[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

REIN_Metrics_T1_T9_mice$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
REIN_Metrics_T1_T9_mice$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
REIN_Metrics_T1_T9_mice$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)


w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

REIN_Metrics_T1_T9_mice$L_total <-
  w_fid  * REIN_Metrics_T1_T9_mice$L_fid +
  w_util * REIN_Metrics_T1_T9_mice$L_util +
  w_conf * REIN_Metrics_T1_T9_mice$L_conf

REIN_Metrics_T1_T9_mice$dataset <- paste0("syn_", seq_len(n_synth))

# Save synthetic datasets
save(REIN_Metrics_T1_T9_mice,
     file = here("results", "REIN_Metrics_T1_T9_mice.Rdata")
)


# chimera no matching
load(here("results", "rein_synthetic_list_data_chimera_no_matching.Rdata"))
results_all <- lapply(seq_along(rein_synthetic_list_data_chimera_no_matching), function(i) {
  Assessment_function_unified_Ablation_Analysis(
    X_real = X_real,
    X_syn  = rein_synthetic_list_data_chimera_no_matching[[i]],
    S = 5,
    seed = seeds[i],
    formula = Surv(times, Censored) ~ .,
    modeling = "cox",
    timevar = "times",
    statusvar = "Censored",
    missing_rate = 0.30,
    iterative_masking=TRUE,
    target_col = "Diabetes",
    feature_associate = feature_associate,
    is_survival = TRUE,
    trn = trn,
    val = val
  )
})

rein_results_all_chimera_no_matching <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(rein_results_all_chimera_no_matching,
     file = here("results", "rein_results_all_chimera_no_matching_min_max.Rdata")
)


# metrics
metrics <- paste0("T", 1:9)

#target vector
T_star <- c(
  T1 = 0,
  T2 = 0,
  T3 = 0.5,
  T4 = 0,
  T5 = 0,
  T6 = 0,
  T7 = 1,
  T8 = 0.5,
  T9 = 0
)

REIN_Metrics_T1_T9_chimera_no_matching <- rein_results_all_chimera_no_matching|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(REIN_Metrics_T1_T9_chimera_no_matching[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(REIN_Metrics_T1_T9_chimera_no_matching[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

REIN_Metrics_T1_T9_chimera_no_matching$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
REIN_Metrics_T1_T9_chimera_no_matching$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
REIN_Metrics_T1_T9_chimera_no_matching$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)

w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

REIN_Metrics_T1_T9_chimera_no_matching$L_total <-
  w_fid  * REIN_Metrics_T1_T9_chimera_no_matching$L_fid +
  w_util * REIN_Metrics_T1_T9_chimera_no_matching$L_util +
  w_conf * REIN_Metrics_T1_T9_chimera_no_matching$L_conf

REIN_Metrics_T1_T9_chimera_no_matching$dataset <- paste0("syn_", seq_len(n_synth))

# Save synthetic datasets
save(REIN_Metrics_T1_T9_chimera_no_matching,
     file = here("results", "REIN_Metrics_T1_T9_chimera_no_matching.Rdata")
)
