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

#CHIMERA ========================================================================================
# ===============================================================================================
# DATA LOADING : PIMA
# ===============================================================================================
PIMA <- read.csv(here("data", "PIMA.csv"))
PIMA <- PIMA|> mutate_if(is.character,as.factor)

n_synth <- 50
set.seed(123)
seeds <- sample(1:1000, n_synth)

pima_synthetic_list_data_chimera <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- CHIMERA_generate(
      df=PIMA,
      missing_rate = 0.25,
      is_survival = FALSE,
      timevar = NULL,
      statusvar = NULL,
      nb_imput = 1,
      nb_max_it = 5,
      seed = seeds[i],
      method_conti = "pmm",
      verbose = TRUE
    )
    
  })["elapsed"]
  
  pima_synthetic_list_data_chimera[[i]] <- synthetic_data[["df_syn"]]
}

# -----------------------------------------------------------------------------------------------
# EXECUTION TIME SUMMARY
# -----------------------------------------------------------------------------------------------
execution_time_mean  <- mean(execution_times) # 2.988 s
execution_time_upper <- execution_time_mean + sd(execution_times) # 2.499 s
execution_time_lower <- execution_time_mean - sd(execution_times) # 3.478

# Save synthetic datasets
save(pima_synthetic_list_data_chimera,
     file = here("results", "pima_synthetic_list_data_chimera.Rdata")
)



# ===============================================================================================
# DATA LOADING : AIDS
# ===============================================================================================
AIDS <- read.csv(here("data", "AIDS.csv"))
AIDS <- AIDS|> mutate_if(is.character,as.factor)

aids_synthetic_list_data_chimera <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- CHIMERA_generate(
      df=AIDS,
      missing_rate = 0.15,
      is_survival = TRUE,
      timevar = "times",
      statusvar = "Censored",
      nb_imput = 1,
      nb_max_it = 5,
      seed = seeds[i],
      method_conti = "pmm",
      verbose = TRUE
    )
  })["elapsed"]
  
  aids_synthetic_list_data_chimera[[i]] <- synthetic_data[["df_syn"]]
}

# -----------------------------------------------------------------------------------------------
# EXECUTION TIME SUMMARY
# -----------------------------------------------------------------------------------------------
execution_time_mean  <- mean(execution_times) # 19.842
execution_time_upper <- execution_time_mean + sd(execution_times) # 7.024
execution_time_lower <- execution_time_mean - sd(execution_times) # 32.661

# Save synthetic datasets
save(aids_synthetic_list_data_chimera,
     file = here("results", "aids_synthetic_list_data_chimera.Rdata")
)


# ===============================================================================================
# DATA LOADING : REIN
# ===============================================================================================
df_rein <- read.csv(here("data", "df_rein_with_na.csv"))
df_rein <- df_rein|> mutate_if(is.character,as.factor)

rein_synthetic_list_data_chimera <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- CHIMERA_generate(
      df=df_rein,
      missing_rate = 0.30,
      is_survival = TRUE,
      timevar = "times",
      statusvar = "Censored",
      nb_imput = 1,
      nb_max_it = 5,
      seed = seeds[i],
      method_conti = "pmm",
      verbose = TRUE
    )
  })["elapsed"]
  
  rein_synthetic_list_data_chimera[[i]] <- synthetic_data[["df_syn"]]
}

# -----------------------------------------------------------------------------------------------
# EXECUTION TIME SUMMARY
# -----------------------------------------------------------------------------------------------
execution_time_mean  <- mean(execution_times) # 500.354
execution_time_upper <- execution_time_mean + sd(execution_times) # 492.944
execution_time_lower <- execution_time_mean - sd(execution_times) # 507.765

# Save synthetic datasets
save(rein_synthetic_list_data_chimera,
     file = here("results", "rein_synthetic_list_data_chimera.Rdata")
)


#SYNTHPOP ========================================================================================
# ===============================================================================================
# PIMA
# ===============================================================================================
pima_synthetic_list_data_synthpop <- vector("list", n_synth)

execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  
  execution_times[i] <- system.time({
    synthetic_data <- syn(PIMA,seed = seeds[i])
  })["elapsed"]
  
  pima_synthetic_list_data_synthpop[[i]] <- synthetic_data[["syn"]]
}

execution_time_mean  <- mean(execution_times) # 0.104
execution_time_upper <- execution_time_mean + sd(execution_times) # 0.137
execution_time_lower <- execution_time_mean - sd(execution_times) # 0.070s

# Save synthetic datasets
save(pima_synthetic_list_data_synthpop,
     file = here("results", "pima_synthetic_list_data_synthpop.Rdata")
)

# ===============================================================================================
# AIDS
# ===============================================================================================
aids_synthetic_list_data_synthpop <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  
  execution_times[i] <- system.time({
    synthetic_data <- syn(AIDS,seed = seeds[i])
  })["elapsed"]
  
  aids_synthetic_list_data_synthpop[[i]] <- synthetic_data[["syn"]]
}

execution_time_mean  <- mean(execution_times) # 1.105
execution_time_upper <- execution_time_mean + sd(execution_times) # 1.180
execution_time_lower <- execution_time_mean - sd(execution_times) # 1.029

# Save synthetic datasets
save(aids_synthetic_list_data_synthpop,
     file = here("results", "aids_synthetic_list_data_synthpop.Rdata")
)


# ===============================================================================================
# REIN
# ===============================================================================================
rein_synthetic_list_data_synthpop <- vector("list", n_synth)

execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  
  execution_times[i] <- system.time({
    synthetic_data <- syn(df_rein,seed = seeds[i])
  })["elapsed"]
  
  rein_synthetic_list_data_synthpop[[i]] <- synthetic_data[["syn"]]
}

execution_time_mean  <- mean(execution_times) # 57.312
execution_time_upper <- execution_time_mean + sd(execution_times) # 58.094
execution_time_lower <- execution_time_mean - sd(execution_times) # 56.530

# Save synthetic datasets
save(rein_synthetic_list_data_synthpop,
     file = here("reuslts", "rein_synthetic_list_data_synthpop.Rdata")
)


# CTGAN =========================================================================================
# PIMA
files <- sprintf(
  "%s/PIMA_CTGAN_%02d.csv",
  here("data","pima_synthetic_datasets_ctgan"),
  1:n_synth
)

pima_synthetic_list_data_ctgan <- vector("list", length(files))
pima_synthetic_list_data_ctgan <- lapply(files, read_csv, show_col_types = FALSE)

pima_synthetic_list_data_ctgan <- lapply(
  pima_synthetic_list_data_ctgan,
  \(df) df |> as.data.frame() |> dplyr::mutate(across(where(is.character), as.factor))
)

# Save synthetic datasets
save(pima_synthetic_list_data_ctgan,
     file = here("results", "pima_synthetic_list_data_ctgan.Rdata")
)

#execution_time_mean  = 8.34 sec
#execution_time_upper = 8.87 sec
#execution_time_lower = 7.81 sec

# AIDS
files <- sprintf(
  "%s/AIDS_CTGAN_%02d.csv",
  here("data","aids_synthetic_datasets_ctgan"),
  1:n_synth
)

aids_synthetic_list_data_ctgan <- vector("list", length(files))
aids_synthetic_list_data_ctgan <- lapply(files, read_csv, show_col_types = FALSE)

aids_synthetic_list_data_ctgan <- lapply(
  aids_synthetic_list_data_ctgan,
  \(df) df |> as.data.frame() |> dplyr::mutate(across(where(is.character), as.factor))
)

# Save synthetic datasets
save(aids_synthetic_list_data_ctgan,
     file = here("results", "aids_synthetic_list_data_ctgan.Rdata")
)

#execution_time_mean  = 39.79 sec
#execution_time_upper = 40.21 sec
#execution_time_lower = 39.37 sec

# REIN
files <- sprintf(
  "%s/REIN_CTGAN_%02d.csv",
  here("data","rein_synthetic_datasets_ctgan"),
  1:n_synth
)

rein_synthetic_list_data_ctgan <- vector("list", length(files))
rein_synthetic_list_data_ctgan <- lapply(files, read_csv, show_col_types = FALSE)

rein_synthetic_list_data_ctgan <- lapply(
  rein_synthetic_list_data_ctgan,
  \(df) df |> as.data.frame() |> dplyr::mutate(across(where(is.character), as.factor))
)



# Save synthetic datasets
save(rein_synthetic_list_data_ctgan,
     file = here("results", "rein_synthetic_list_data_ctgan.Rdata")
)

#execution_time_mean  = 786.21 sec
#execution_time_upper = 792.26 sec
#execution_time_lower = 780.16 sec


# CTGAN tain for holdout set =========================================================================================
# PIMA
files <- sprintf(
  "%s/PIMA_CTGAN_%02d.csv",
  here("data","pima_synthetic_datasets_ctgan_trn"),
  1:n_synth
)

pima_synthetic_list_data_ctgan_trn <- vector("list", length(files))
pima_synthetic_list_data_ctgan_trn <- lapply(files, read_csv, show_col_types = FALSE)

pima_synthetic_list_data_ctgan_trn <- lapply(
  pima_synthetic_list_data_ctgan,
  \(df) df |> as.data.frame() |> dplyr::mutate(across(where(is.character), as.factor))
)

# Save synthetic datasets
save(pima_synthetic_list_data_ctgan_trn,
     file = here("results", "pima_synthetic_list_data_ctgan_trn.Rdata")
)

# AIDS
files <- sprintf(
  "%s/AIDS_CTGAN_%02d.csv",
  here("data","aids_synthetic_datasets_ctgan_trn"),
  1:n_synth
)

aids_synthetic_list_data_ctgan_trn <- vector("list", length(files))
aids_synthetic_list_data_ctgan_trn <- lapply(files, read_csv, show_col_types = FALSE)

aids_synthetic_list_data_ctgan_trn <- lapply(
  aids_synthetic_list_data_ctgan,
  \(df) df |> as.data.frame() |> dplyr::mutate(across(where(is.character), as.factor))
)

# Save synthetic datasets
save(aids_synthetic_list_data_ctgan_trn,
     file = here("results", "aids_synthetic_list_data_ctgan_trn.Rdata")
)

# REIN
files <- sprintf(
  "%s/REIN_CTGAN_%02d.csv",
  here("data","rein_synthetic_datasets_ctgan_trn"),
  1:n_synth
)

rein_synthetic_list_data_ctgan_trn <- vector("list", length(files))
rein_synthetic_list_data_ctgan_trn <- lapply(files, read_csv, show_col_types = FALSE)

rein_synthetic_list_data_ctgan_trn <- lapply(
  rein_synthetic_list_data_ctgan_trn,
  \(df) df |> as.data.frame() |> dplyr::mutate(across(where(is.character), as.factor))
)


# Save synthetic datasets
save(rein_synthetic_list_data_ctgan_trn,
     file = here("results", "rein_synthetic_list_data_ctgan_trn.Rdata")
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


# chimera
load(here("results", "pima_synthetic_list_data_chimera.Rdata"))
results_all <- lapply(seq_along(pima_synthetic_list_data_chimera), function(i) {
  Assessment_function_unified(
    X_real = X_real,
    X_syn  = pima_synthetic_list_data_chimera[[i]],
    S = 5,
    seed = seeds[i],
    formula = Diabetes_diagnosis ~ .,
    modeling = "logistic",
    timevar = NULL,
    statusvar = NULL,
    missing_rate = 0.25,
    target_col = "Diabetes_diagnosis",
    feature_associate = feature_associate,
    is_survival = FALSE,
    generator= "chimera",
    trn = trn,
    val = val,
    syn_ctgan = NULL
  )
})

pima_results_all_chimera <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(pima_results_all_chimera,
     file = here("results", "pima_results_all_chimera_min_max.Rdata")
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

PIMA_Metrics_T1_T9_chimera <- pima_results_all_chimera|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(PIMA_Metrics_T1_T9_chimera[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(PIMA_Metrics_T1_T9_chimera[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

PIMA_Metrics_T1_T9_chimera$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
PIMA_Metrics_T1_T9_chimera$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
PIMA_Metrics_T1_T9_chimera$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)

# Poids (modifiable)
w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

PIMA_Metrics_T1_T9_chimera$L_total <-
  w_fid  * PIMA_Metrics_T1_T9_chimera$L_fid +
  w_util * PIMA_Metrics_T1_T9_chimera$L_util +
  w_conf * PIMA_Metrics_T1_T9_chimera$L_conf

PIMA_Metrics_T1_T9_chimera$dataset <- paste0("syn_", seq_len(n_synth))

# Save synthetic datasets
save(PIMA_Metrics_T1_T9_chimera,
     file = here("results", "PIMA_Metrics_T1_T9_chimera.Rdata")
)


# synthpop
load(here("results", "pima_synthetic_list_data_synthpop.Rdata"))
results_all <- lapply(seq_along(pima_synthetic_list_data_synthpop), function(i) {
  Assessment_function_unified(
    X_real = X_real,
    X_syn  = pima_synthetic_list_data_synthpop[[i]],
    S = 5,
    seed = seeds[i],
    formula = Diabetes_diagnosis ~ .,
    modeling = "logistic",
    timevar = NULL,
    statusvar = NULL,
    missing_rate = 0.25,
    target_col = "Diabetes_diagnosis",
    feature_associate = feature_associate,
    is_survival = FALSE,
    generator= "synthpop",
    trn = trn,
    val = val,
    syn_ctgan = NULL
  )
})

pima_results_all_synthpop <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(pima_results_all_synthpop,
     file = here("results", "pima_results_all_synthpop_min_max.Rdata")
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

PIMA_Metrics_T1_T9_synthpop <- pima_results_all_synthpop|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(PIMA_Metrics_T1_T9_synthpop[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(PIMA_Metrics_T1_T9_synthpop[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

PIMA_Metrics_T1_T9_synthpop$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
PIMA_Metrics_T1_T9_synthpop$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
PIMA_Metrics_T1_T9_synthpop$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)

# Poids (modifiable)
w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

PIMA_Metrics_T1_T9_synthpop$L_total <-
  w_fid  * PIMA_Metrics_T1_T9_synthpop$L_fid +
  w_util * PIMA_Metrics_T1_T9_synthpop$L_util +
  w_conf * PIMA_Metrics_T1_T9_synthpop$L_conf

PIMA_Metrics_T1_T9_synthpop$dataset <- paste0("syn_", seq_len(n_synth))

# Save synthetic datasets
save(PIMA_Metrics_T1_T9_synthpop,
     file = here("results", "PIMA_Metrics_T1_T9_synthpop.Rdata")
)


# CTGAN =========================================================================================
load(here("results", "pima_synthetic_list_data_ctgan.Rdata"))
load(here("results", "pima_synthetic_list_data_ctgan_trn.Rdata"))

results_all <- lapply(seq_along(pima_synthetic_list_data_ctgan), function(i) {
  Assessment_function_unified(
    X_real = X_real,
    X_syn  = pima_synthetic_list_data_ctgan[[i]],
    S = 5,
    seed = seeds[i],
    formula = Diabetes_diagnosis ~ .,
    modeling = "logistic",
    timevar = NULL,
    statusvar = NULL,
    missing_rate = 0.25,
    target_col = "Diabetes_diagnosis",
    feature_associate = feature_associate,
    is_survival = FALSE,
    generator= "ctgan",
    trn = trn,
    val = val,
    syn_ctgan = pima_synthetic_list_data_ctgan_trn[[i]]
  )
})

pima_results_all_ctgan <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(pima_results_all_ctgan,
     file = here("results", "pima_results_all_ctgan_min_max.Rdata")
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

PIMA_Metrics_T1_T9_ctgan <- pima_results_all_ctgan|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(PIMA_Metrics_T1_T9_ctgan[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(PIMA_Metrics_T1_T9_ctgan[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

PIMA_Metrics_T1_T9_ctgan$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
PIMA_Metrics_T1_T9_ctgan$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
PIMA_Metrics_T1_T9_ctgan$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)

# Poids (modifiable)
w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

PIMA_Metrics_T1_T9_ctgan$L_total <-
  w_fid  * PIMA_Metrics_T1_T9_ctgan$L_fid +
  w_util * PIMA_Metrics_T1_T9_ctgan$L_util +
  w_conf * PIMA_Metrics_T1_T9_ctgan$L_conf

PIMA_Metrics_T1_T9_ctgan$dataset <- paste0("syn_", seq_len(n_synth))

# Save synthetic datasets
save(PIMA_Metrics_T1_T9_ctgan,
     file = here("results", "PIMA_Metrics_T1_T9_ctgan.Rdata")
)


#================================================================================================
# Select the best dataset use for downstream analyses
#================================================================================================
X_real <- read.csv(here("data", "AIDS.csv"))
X_real <- X_real|> mutate_if(is.character,as.factor)
#X_real$Censored <- as.factor(X_real$Censored)

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

# chimera
load(here("results", "aids_synthetic_list_data_chimera.Rdata"))
results_all <- lapply(seq_along(aids_synthetic_list_data_chimera), function(i) {
  Assessment_function_unified(
    X_real = X_real,
    X_syn  = aids_synthetic_list_data_chimera[[i]],
    S = 5,
    seed = seeds[i],
    formula = Surv(times, Censored) ~ .,
    modeling = "cox",
    timevar = "times",
    statusvar = "Censored",
    missing_rate = 0.15,
    target_col = "Homosexuality",
    feature_associate = feature_associate,
    is_survival = TRUE,
    generator= "chimera",
    trn = trn,
    val = val,
    syn_ctgan = NULL
  )
})

aids_results_all_chimera <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(aids_results_all_chimera,
     file = here("results", "aids_results_all_chimera_min_max.Rdata")
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

AIDS_Metrics_T1_T9_chimera <- aids_results_all_chimera|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(AIDS_Metrics_T1_T9_chimera[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(AIDS_Metrics_T1_T9_chimera[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

AIDS_Metrics_T1_T9_chimera$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
AIDS_Metrics_T1_T9_chimera$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
AIDS_Metrics_T1_T9_chimera$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)

# Poids (modifiable)
w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

AIDS_Metrics_T1_T9_chimera$L_total <-
  w_fid  * AIDS_Metrics_T1_T9_chimera$L_fid +
  w_util * AIDS_Metrics_T1_T9_chimera$L_util +
  w_conf * AIDS_Metrics_T1_T9_chimera$L_conf

AIDS_Metrics_T1_T9_chimera$dataset <- paste0("syn_", seq_len(n_synth))

# Save synthetic datasets
save(AIDS_Metrics_T1_T9_chimera,
     file = here("results", "AIDS_Metrics_T1_T9_chimera.Rdata")
)


# synthpop
load(here("results", "aids_synthetic_list_data_synthpop.Rdata"))
results_all <- lapply(seq_along(aids_synthetic_list_data_synthpop), function(i) {
  Assessment_function_unified(
    X_real = X_real,
    X_syn  = aids_synthetic_list_data_synthpop[[i]],
    S = 5,
    seed = seeds[i],
    formula = Surv(times, Censored) ~ .,
    modeling = "cox",
    timevar = "times",
    statusvar = "Censored",
    missing_rate = 0.15,
    target_col = "Homosexuality",
    feature_associate = feature_associate,
    is_survival = TRUE,
    generator= "synthpop",
    trn = trn,
    val = val,
    syn_ctgan = NULL
  )
})

aids_results_all_synthpop <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(aids_results_all_synthpop,
     file = here("results", "aids_results_all_synthpop_min_max.Rdata")
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

AIDS_Metrics_T1_T9_synthpop <- aids_results_all_synthpop|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(AIDS_Metrics_T1_T9_synthpop[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(AIDS_Metrics_T1_T9_synthpop[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

AIDS_Metrics_T1_T9_synthpop$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
AIDS_Metrics_T1_T9_synthpop$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
AIDS_Metrics_T1_T9_synthpop$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)

# Poids (modifiable)
w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

AIDS_Metrics_T1_T9_synthpop$L_total <-
  w_fid  * AIDS_Metrics_T1_T9_synthpop$L_fid +
  w_util * AIDS_Metrics_T1_T9_synthpop$L_util +
  w_conf * AIDS_Metrics_T1_T9_synthpop$L_conf

AIDS_Metrics_T1_T9_synthpop$dataset <- paste0("syn_", seq_len(n_synth))

# Save synthetic datasets
save(AIDS_Metrics_T1_T9_synthpop,
     file = here("results", "AIDS_Metrics_T1_T9_synthpop.Rdata")
)


# ctgan
load(here("results", "aids_synthetic_list_data_ctgan.Rdata"))
load(here("results", "aids_synthetic_list_data_ctgan_trn.Rdata"))

results_all <- lapply(seq_along(aids_synthetic_list_data_ctgan), function(i) {
  Assessment_function_unified(
    X_real = X_real,
    X_syn  = aids_synthetic_list_data_ctgan[[i]],
    S = 5,
    seed = seeds[i],
    formula = Surv(times, Censored) ~ .,
    modeling = "cox",
    timevar = "times",
    statusvar = "Censored",
    missing_rate = 0.15,
    target_col = "Homosexuality",
    feature_associate = feature_associate,
    is_survival = TRUE,
    generator= "ctgan",
    trn = trn,
    val = val,
    syn_ctgan = aids_synthetic_list_data_ctgan_trn[[i]]
  )
})

aids_results_all_ctgan <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(aids_results_all_ctgan,
     file = here("results", "aids_results_all_ctgan_min_max.Rdata")
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

AIDS_Metrics_T1_T9_ctgan <- aids_results_all_ctgan|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(AIDS_Metrics_T1_T9_ctgan[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(AIDS_Metrics_T1_T9_ctgan[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

AIDS_Metrics_T1_T9_ctgan$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
AIDS_Metrics_T1_T9_ctgan$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
AIDS_Metrics_T1_T9_ctgan$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)

# Poids (modifiable)
w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

AIDS_Metrics_T1_T9_ctgan$L_total <-
  w_fid  * AIDS_Metrics_T1_T9_ctgan$L_fid +
  w_util * AIDS_Metrics_T1_T9_ctgan$L_util +
  w_conf * AIDS_Metrics_T1_T9_ctgan$L_conf

AIDS_Metrics_T1_T9_ctgan$dataset <- paste0("syn_", seq_len(n_synth))

# Save synthetic datasets
save(AIDS_Metrics_T1_T9_ctgan,
     file = here("results", "AIDS_Metrics_T1_T9_ctgan.Rdata")
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


# chimera
load(here("results", "rein_synthetic_list_data_chimera.Rdata"))
results_all <- lapply(seq_along(rein_synthetic_list_data_chimera), function(i) {
  Assessment_function_unified(
    X_real = X_real,
    X_syn  = rein_synthetic_list_data_chimera[[i]],
    S = 5,
    seed = seeds[i],
    formula = Surv(times, Censored) ~ .,
    modeling = "cox",
    timevar = "times",
    statusvar = "Censored",
    missing_rate = 0.30,
    target_col = "Diabetes",
    feature_associate = feature_associate,
    is_survival = TRUE,
    generator= "chimera",
    trn = trn,
    val = val,
    syn_ctgan = NULL
  )
})

rein_results_all_chimera <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(rein_results_all_chimera,
     file = here("results", "rein_results_all_chimera_min_max.Rdata")
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

REIN_Metrics_T1_T9_chimera <- rein_results_all_chimera|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(REIN_Metrics_T1_T9_chimera[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(REIN_Metrics_T1_T9_chimera[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

REIN_Metrics_T1_T9_chimera$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
REIN_Metrics_T1_T9_chimera$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
REIN_Metrics_T1_T9_chimera$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)

# Poids (modifiable)
w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

REIN_Metrics_T1_T9_chimera$L_total <-
  w_fid  * REIN_Metrics_T1_T9_chimera$L_fid +
  w_util * REIN_Metrics_T1_T9_chimera$L_util +
  w_conf * REIN_Metrics_T1_T9_chimera$L_conf

REIN_Metrics_T1_T9_chimera$dataset <- paste0("syn_", seq_len(n_synth))

# Save synthetic datasets
save(REIN_Metrics_T1_T9_chimera,
     file = here("results", "REIN_Metrics_T1_T9_chimera.Rdata")
)


# synthpop
load(here("results", "rein_synthetic_list_data_synthpop.Rdata"))
results_all <- lapply(seq_along(rein_synthetic_list_data_synthpop), function(i) {
  Assessment_function_unified(
    X_real = X_real,
    X_syn  = rein_synthetic_list_data_synthpop[[i]],
    S = 5,
    seed = seeds[i],
    formula = Surv(times, Censored) ~ .,
    modeling = "cox",
    timevar = "times",
    statusvar = "Censored",
    missing_rate = 0.30,
    target_col = "Diabetes",
    feature_associate = feature_associate,
    is_survival = TRUE,
    generator= "synthpop",
    trn = trn,
    val = val,
    syn_ctgan = NULL
  )
})

rein_results_all_synthpop <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(rein_results_all_synthpop,
     file = here("results", "rein_results_all_synthpop_min_max.Rdata")
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

REIN_Metrics_T1_T9_synthpop <- rein_results_all_synthpop|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(REIN_Metrics_T1_T9_synthpop[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(REIN_Metrics_T1_T9_synthpop[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

REIN_Metrics_T1_T9_synthpop$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
REIN_Metrics_T1_T9_synthpop$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
REIN_Metrics_T1_T9_synthpop$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)

# Poids (modifiable)
w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

REIN_Metrics_T1_T9_synthpop$L_total <-
  w_fid  * REIN_Metrics_T1_T9_synthpop$L_fid +
  w_util * REIN_Metrics_T1_T9_synthpop$L_util +
  w_conf * REIN_Metrics_T1_T9_synthpop$L_conf

REIN_Metrics_T1_T9_synthpop$dataset <- paste0("syn_", seq_len(n_synth))

# Save synthetic datasets
save(REIN_Metrics_T1_T9_synthpop,
     file = here("results", "REIN_Metrics_T1_T9_synthpop.Rdata")
)


# ctgan
load(here("results", "rein_synthetic_list_data_ctgan.Rdata"))
load(here("results", "rein_synthetic_list_data_ctgan_trn.Rdata"))


results_all <- lapply(seq_along(rein_synthetic_list_data_ctgan), function(i) {
  Assessment_function_unified(
    X_real = X_real,
    X_syn  = rein_synthetic_list_data_ctgan[[i]],
    S = 5,
    seed = seeds[i],
    formula = Surv(times, Censored) ~ .,
    modeling = "cox",
    timevar = "times",
    statusvar = "Censored",
    missing_rate = 0.30,
    target_col = "Diabetes",
    feature_associate = feature_associate,
    is_survival = TRUE,
    generator= "ctgan",
    trn = trn,
    val = val,
    syn_ctgan = rein_synthetic_list_data_ctgan_trn[[i]]
  )
})

rein_results_all_ctgan <- results_all |> 
  map_dfr(~ as.data.frame(t(.x)))

# Save synthetic datasets
save(rein_results_all_ctgan,
     file = here("results", "rein_results_all_ctgan_min_max.Rdata")
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

REIN_Metrics_T1_T9_ctgan <- rein_results_all_ctgan|>dplyr::select(-c(T10,T11,T12,T13))

# family
F_fid  <- c("T1", "T2", "T3")
F_util <- c("T4", "T5", "T6")
F_conf <- c("T7", "T8", "T9")


s_j <- sapply(REIN_Metrics_T1_T9_ctgan[metrics], sd, na.rm = TRUE)
loss_matrix <- sweep(REIN_Metrics_T1_T9_ctgan[metrics], 2, T_star, FUN = "-")
loss_matrix <- sweep(loss_matrix^2, 2, s_j^2, FUN = "/")

REIN_Metrics_T1_T9_ctgan$L_fid  <- rowMeans(loss_matrix[, F_fid],  na.rm = TRUE)
REIN_Metrics_T1_T9_ctgan$L_util <- rowMeans(loss_matrix[, F_util], na.rm = TRUE)
REIN_Metrics_T1_T9_ctgan$L_conf <- rowMeans(loss_matrix[, F_conf], na.rm = TRUE)

# Poids (modifiable)
w_fid  <- 1/3
w_util <- 1/3
w_conf <- 1/3

REIN_Metrics_T1_T9_ctgan$L_total <-
  w_fid  * REIN_Metrics_T1_T9_ctgan$L_fid +
  w_util * REIN_Metrics_T1_T9_ctgan$L_util +
  w_conf * REIN_Metrics_T1_T9_ctgan$L_conf

REIN_Metrics_T1_T9_ctgan$dataset <- paste0("syn_", seq_len(n_synth))

# Save synthetic datasets
save(REIN_Metrics_T1_T9_ctgan,
     file = here("results", "REIN_Metrics_T1_T9_ctgan.Rdata")
)




















