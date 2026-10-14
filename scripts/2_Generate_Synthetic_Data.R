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

# ===============================================================================================
# Directory Setup
# ===============================================================================================
args <- c("data", "results")

data_dir <- normalizePath(args[1], mustWork = TRUE)

if (!dir.exists(args[2])) {
  dir.create(args[2], recursive = TRUE)
}

results_dir <- normalizePath(args[2], mustWork = TRUE)


# ===============================================================================================
# Loading data
# ===============================================================================================

PIMA_train <- read.csv(here(data_dir, "PIMA_train.csv")) |> mutate_if(is.character, as.factor)
AIDS_train <- read.csv(here(data_dir, "AIDS_train.csv")) |> mutate_if(is.character, as.factor)
REIN_train <- read.csv(here(data_dir, "REIN_train.csv")) |> mutate_if(is.character, as.factor)

PIMA_trn <- read.csv(here(data_dir, "PIMA_trn.csv")) |> mutate_if(is.character, as.factor)
AIDS_trn <- read.csv(here(data_dir, "AIDS_trn.csv")) |> mutate_if(is.character, as.factor)
REIN_trn <- read.csv(here(data_dir, "REIN_trn.csv")) |> mutate_if(is.character, as.factor)


# CHIMERA Generator ============================================================
n_synth <- 50
set.seed(123)
seeds <- sample(1:1000, n_synth)

# PIMA
pima_synthetic_list_data_chimera <- vector("list", n_synth)
pima_nb_iteration_chimera <- as.numeric(n_synth)
execution_times <- numeric(n_synth)
missing_rate <- runif(1, min = 0.10, max = 0.30)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    results <- CHIMERA_generate(
      df=PIMA_train,
      missing_rate = missing_rate,
      is_survival = FALSE,
      timevar = NULL,
      statusvar = NULL,
      nb_imput = 1,
      nb_max_it = 5,
      seed = seeds[i],
      verbose = TRUE
    )

  })["elapsed"]

  pima_synthetic_list_data_chimera[[i]] <- results$synthetic_data
  pima_nb_iteration_chimera[i] <- results$nb_iter
}

summary(pima_nb_iteration_chimera)

# Execution time summary

execution_time_mean  <- mean(execution_times)
execution_time_upper <- execution_time_mean + sd(execution_times)
execution_time_lower <- execution_time_mean - sd(execution_times)

execition_time_results <- data.frame(
  execution_time_mean = numeric(),
  execution_time_upper = numeric(),
  execution_time_lower = numeric(),
  dataset = character(),
  method = character()
)

execition_time_results <- rbind(
  execition_time_results,
  data.frame(
    execution_time_mean = execution_time_mean,
    execution_time_upper = execution_time_upper,
    execution_time_lower = execution_time_lower,
    dataset = "pima",
    method = "CHIMERA"
  )
)

# Save synthetic datasets
save(pima_synthetic_list_data_chimera,
     file = here(results_dir, "pima_synthetic_list_data_chimera.Rdata")
)

# Generate data dor Holdout and linkability attack
pima_synthetic_list_data_chimera_trn <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    results <- CHIMERA_generate(
      df=PIMA_trn,
      missing_rate = missing_rate,
      is_survival = FALSE,
      timevar = NULL,
      statusvar = NULL,
      nb_imput = 1,
      nb_max_it = 5,
      seed = seeds[i],
      verbose = TRUE
    )

  })["elapsed"]

  pima_synthetic_list_data_chimera_trn[[i]] <- results$synthetic_data
}


# Save synthetic datasets
save(pima_synthetic_list_data_chimera_trn,
     file = here(results_dir, "pima_synthetic_list_data_chimera_trn.Rdata")
)


# AIDS
aids_synthetic_list_data_chimera <- vector("list", n_synth)
aids_nb_iteration_chimera <- as.numeric(n_synth)
execution_times <- numeric(n_synth)
missing_rate <- runif(1, min = 0.10, max = 0.30)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    results <- CHIMERA_generate(
      df=AIDS_train,
      missing_rate = missing_rate,
      is_survival = TRUE,
      timevar = "times",
      statusvar = "Censored",
      nb_imput = 1,
      nb_max_it = 5,
      seed = seeds[i],
      verbose = TRUE
    )
  })["elapsed"]

  aids_synthetic_list_data_chimera[[i]] <- results$synthetic_data
  aids_nb_iteration_chimera[i]  <- results$nb_iter
}

summary(aids_nb_iteration_chimera)

# -----------------------------------------------------------------------------------------------
# EXECUTION TIME SUMMARY
# -----------------------------------------------------------------------------------------------
execution_time_mean  <- mean(execution_times) # 19.842
execution_time_upper <- execution_time_mean + sd(execution_times) # 7.024
execution_time_lower <- execution_time_mean - sd(execution_times) # 32.661

execition_time_results <- rbind(
  execition_time_results,
  data.frame(
    execution_time_mean = execution_time_mean,
    execution_time_upper = execution_time_upper,
    execution_time_lower = execution_time_lower,
    dataset = "aids",
    method = "CHIMERA"
  )
)

# Save synthetic datasets
save(aids_synthetic_list_data_chimera,
     file = here(results_dir, "aids_synthetic_list_data_chimera.Rdata")
)



# Generate data dor Holdout and linkability attack
aids_synthetic_list_data_chimera_trn <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    results <- CHIMERA_generate(
      df=AIDS_trn,
      missing_rate = 0.15,
      is_survival = TRUE,
      timevar = "times",
      statusvar = "Censored",
      nb_imput = 1,
      nb_max_it = 5,
      seed = seeds[i],
      verbose = TRUE
    )
  })["elapsed"]

  aids_synthetic_list_data_chimera_trn[[i]] <- results$synthetic_data
}


#Save synthetic datasets
save(aids_synthetic_list_data_chimera_trn,
     file = here(results_dir, "aids_synthetic_list_data_chimera_trn.Rdata")
)


# REIN
rein_synthetic_list_data_chimera <- vector("list", n_synth)
rein_nb_iteration_chimera <- as.numeric(n_synth)
execution_times <- numeric(n_synth)
missing_rate <- runif(1, min = 0.10, max = 0.30)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    results <- CHIMERA_generate(
      df=REIN_train,
      missing_rate =missing_rate,
      is_survival = TRUE,
      timevar = "times",
      statusvar = "Censored",
      nb_imput = 1,
      nb_max_it = 5,
      seed = seeds[i],
      verbose = TRUE
    )
  })["elapsed"]

  rein_synthetic_list_data_chimera[[i]] <- results$synthetic_data
  rein_nb_iteration_chimera[i]  <- results$nb_iter
}

summary(rein_nb_iteration_chimera)

# -----------------------------------------------------------------------------------------------
# EXECUTION TIME SUMMARY
# -----------------------------------------------------------------------------------------------
execution_time_mean  <- mean(execution_times) # 500.354
execution_time_upper <- execution_time_mean + sd(execution_times) # 492.944
execution_time_lower <- execution_time_mean - sd(execution_times) # 507.765

execition_time_results <- rbind(
  execition_time_results,
  data.frame(
    execution_time_mean = execution_time_mean,
    execution_time_upper = execution_time_upper,
    execution_time_lower = execution_time_lower,
    dataset = "REIN",
    method = "CHIMERA"
  )
)

# Save synthetic datasets
save(rein_synthetic_list_data_chimera,
     file = here(results_dir, "rein_synthetic_list_data_chimera.Rdata")
)



# Generate data dor Holdout and linkability attack
rein_synthetic_list_data_chimera_trn <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    results <- CHIMERA_generate(
      df=REIN_trn,
      missing_rate = missing_rate,
      is_survival = TRUE,
      timevar = "times",
      statusvar = "Censored",
      nb_imput = 1,
      nb_max_it = 5,
      seed = seeds[i],
      verbose = TRUE
    )
  })["elapsed"]

  rein_synthetic_list_data_chimera_trn[[i]] <- results$synthetic_data
}



# Save synthetic datasets
save(rein_synthetic_list_data_chimera_trn,
     file = here(results_dir, "rein_synthetic_list_data_chimera_trn.Rdata")
)


# SYNTHPOP Generator ============================================================

# PIMA
pima_synthetic_list_data_synthpop <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {

  execution_times[i] <- system.time({
    synthetic_data <- syn(PIMA_train,seed = seeds[i])
  })["elapsed"]

  pima_synthetic_list_data_synthpop[[i]] <- synthetic_data[["syn"]]
}

execution_time_mean  <- mean(execution_times) # 0.104
execution_time_upper <- execution_time_mean + sd(execution_times) # 0.137
execution_time_lower <- execution_time_mean - sd(execution_times) # 0.070s

execition_time_results <- rbind(
  execition_time_results,
  data.frame(
    execution_time_mean = execution_time_mean,
    execution_time_upper = execution_time_upper,
    execution_time_lower = execution_time_lower,
    dataset = "pima",
    method = "SYNTHPOP"
  )
)

# Save synthetic datasets
save(pima_synthetic_list_data_synthpop,
     file = here(results_dir, "pima_synthetic_list_data_synthpop.Rdata")
)



# Generate data dor Holdout and linkability attack
pima_synthetic_list_data_synthpop_trn <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {

  execution_times[i] <- system.time({
    synthetic_data <- syn(PIMA_trn,seed = seeds[i])
  })["elapsed"]

  pima_synthetic_list_data_synthpop_trn[[i]] <- synthetic_data[["syn"]]
}


# Save synthetic datasets
save(pima_synthetic_list_data_synthpop_trn,
     file = here(results_dir, "pima_synthetic_list_data_synthpop_trn.Rdata")
)

# AIDS
aids_synthetic_list_data_synthpop <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {

  execution_times[i] <- system.time({
    synthetic_data <- syn(AIDS_train,seed = seeds[i])
  })["elapsed"]

  aids_synthetic_list_data_synthpop[[i]] <- synthetic_data[["syn"]]
}

execution_time_mean  <- mean(execution_times) # 1.105
execution_time_upper <- execution_time_mean + sd(execution_times) # 1.180
execution_time_lower <- execution_time_mean - sd(execution_times) # 1.029

execition_time_results <- rbind(
  execition_time_results,
  data.frame(
    execution_time_mean = execution_time_mean,
    execution_time_upper = execution_time_upper,
    execution_time_lower = execution_time_lower,
    dataset = "aids",
    method = "SYNTHPOP"
  )
)

# Save synthetic datasets
save(aids_synthetic_list_data_synthpop,
     file = here(results_dir, "aids_synthetic_list_data_synthpop.Rdata")
)



# Generate data dor Holdout and linkability attack
aids_synthetic_list_data_synthpop_trn <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {

  execution_times[i] <- system.time({
    synthetic_data <- syn(AIDS_trn,seed = seeds[i])
  })["elapsed"]

  aids_synthetic_list_data_synthpop_trn[[i]] <- synthetic_data[["syn"]]
}

# Save synthetic datasets
save(aids_synthetic_list_data_synthpop_trn,
     file = here(results_dir, "aids_synthetic_list_data_synthpop_trn.Rdata")
)



# REIN
rein_synthetic_list_data_synthpop <- vector("list", n_synth)

execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {

  execution_times[i] <- system.time({
    synthetic_data <- syn(REIN_train,seed = seeds[i])
  })["elapsed"]

  rein_synthetic_list_data_synthpop[[i]] <- synthetic_data[["syn"]]
}

execution_time_mean  <- mean(execution_times) # 57.312
execution_time_upper <- execution_time_mean + sd(execution_times) # 58.094
execution_time_lower <- execution_time_mean - sd(execution_times) # 56.530

execition_time_results <- rbind(
  execition_time_results,
  data.frame(
    execution_time_mean = execution_time_mean,
    execution_time_upper = execution_time_upper,
    execution_time_lower = execution_time_lower,
    dataset = "REIN",
    method = "SYNTHPOP"
  )
)

# Save synthetic datasets
save(rein_synthetic_list_data_synthpop,
     file = here(results_dir, "rein_synthetic_list_data_synthpop.Rdata")
)


# Save dataset
write.csv(execition_time_results, here(data_dir, "execution_time_results.csv"), row.names = FALSE)


# Generate data dor Holdout and linkability attack
rein_synthetic_list_data_synthpop_trn <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {

  execution_times[i] <- system.time({
    synthetic_data <- syn(REIN_trn,seed = seeds[i])
  })["elapsed"]

  rein_synthetic_list_data_synthpop_trn[[i]] <- synthetic_data[["syn"]]
}


# Save synthetic datasets
save(rein_synthetic_list_data_synthpop_trn,
     file = here(results_dir, "rein_synthetic_list_data_synthpop_trn.Rdata")
)


#####################################################################################
# Ablation study
####################################################################################

pima_synthetic_list_data_mice <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = PIMA_train,
      missing_rate = NULL,
      iterative_masking = FALSE,
      survival_data = FALSE,
      time_variable = NULL,
      event_variable = NULL,
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      verbose = TRUE
    )
    
  })["elapsed"]
  
  pima_synthetic_list_data_mice[[i]] <- synthetic_data[["synthetic_data"]]
}


# Save synthetic datasets
save(pima_synthetic_list_data_mice,
     file = here(results_dir, "pima_synthetic_list_data_mice.Rdata")
)

# For Holout set 
pima_synthetic_list_data_mice_trn <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = PIMA_trn,
      missing_rate = NULL,
      iterative_masking = FALSE,
      survival_data = FALSE,
      time_variable = NULL,
      event_variable = NULL,
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      verbose = TRUE
    )
    
  })["elapsed"]
  
  pima_synthetic_list_data_mice_trn[[i]] <- synthetic_data[["synthetic_data"]]
}


# Save synthetic datasets
save(pima_synthetic_list_data_mice_trn,
     file = here(results_dir, "pima_synthetic_list_data_mice_trn.Rdata")
)


# CHIMERA NO MATCHING
pima_synthetic_list_data_chimera_no_matching <- vector("list", n_synth)
execution_times <- numeric(n_synth)
missing_rate <- runif(1, min = 0.10, max = 0.30)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = PIMA_train,
      missing_rate = missing_rate,
      iterative_masking = TRUE,
      survival_data = FALSE,
      time_variable = NULL,
      event_variable = NULL,
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      verbose = TRUE
    )
    
  })["elapsed"]
  
  pima_synthetic_list_data_chimera_no_matching[[i]] <- synthetic_data[["synthetic_data"]]
}

# Save synthetic datasets
save(pima_synthetic_list_data_chimera_no_matching,
     file = here(results_dir, "pima_synthetic_list_data_chimera_no_matching.Rdata")
)


# for Holdout set
pima_synthetic_list_data_chimera_no_matching_trn <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = PIMA_trn,
      missing_rate = missing_rate,
      iterative_masking = TRUE,
      survival_data = FALSE,
      time_variable = NULL,
      event_variable = NULL,
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      verbose = TRUE
    )
    
  })["elapsed"]
  
  pima_synthetic_list_data_chimera_no_matching_trn[[i]] <- synthetic_data[["synthetic_data"]]
}

# Save synthetic datasets
save(pima_synthetic_list_data_chimera_no_matching_trn,
     file = here(results_dir, "pima_synthetic_list_data_chimera_no_matching_trn.Rdata")
)


# ===============================================================================================
# AIDS
# ===============================================================================================

# MICE =========================================================================
aids_synthetic_list_data_mice <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = AIDS_train,
      missing_rate = NULL,
      iterative_masking = FALSE,
      survival_data = TRUE,
      time_variable = "times",
      event_variable = "Censored",
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      verbose = TRUE
    )
    
  })["elapsed"]
  
  aids_synthetic_list_data_mice[[i]] <- synthetic_data[["synthetic_data"]]
}


# Save synthetic datasets
save(aids_synthetic_list_data_mice,
     file = here(results_dir, "aids_synthetic_list_data_mice.Rdata")
)


# for Holoutset
aids_synthetic_list_data_mice_trn <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = AIDS_trn,
      missing_rate = NULL,
      iterative_masking = FALSE,
      survival_data = TRUE,
      time_variable = "times",
      event_variable = "Censored",
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      verbose = TRUE
    )
    
  })["elapsed"]
  
  aids_synthetic_list_data_mice_trn[[i]] <- synthetic_data[["synthetic_data"]]
}

# Save synthetic datasets
save(aids_synthetic_list_data_mice_trn,
     file = here(results_dir, "aids_synthetic_list_data_mice_trn.Rdata")
)


# CHIMERA NO MATCHING
aids_synthetic_list_data_chimera_no_matching <- vector("list", n_synth)
execution_times <- numeric(n_synth)
missing_rate <- runif(1, min = 0.10, max = 0.30)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = AIDS_train,
      missing_rate = missing_rate,
      iterative_masking = TRUE,
      survival_data = TRUE,
      time_variable = "times",
      event_variable = "Censored",
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      verbose = TRUE
    )
    
  })["elapsed"]
  
  aids_synthetic_list_data_chimera_no_matching[[i]] <- synthetic_data[["synthetic_data"]]
}


# Save synthetic datasets
save(aids_synthetic_list_data_chimera_no_matching,
     file = here(results_dir, "aids_synthetic_list_data_chimera_no_matching.Rdata")
)

# Holdout set
aids_synthetic_list_data_chimera_no_matching_trn <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = AIDS_trn,
      missing_rate = missing_rate,
      iterative_masking = TRUE,
      survival_data = TRUE,
      time_variable = "times",
      event_variable = "Censored",
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      verbose = TRUE
    )
    
  })["elapsed"]
  
  aids_synthetic_list_data_chimera_no_matching_trn[[i]] <- synthetic_data[["synthetic_data"]]
}


# Save synthetic datasets
save(aids_synthetic_list_data_chimera_no_matching_trn,
     file = here(results_dir, "aids_synthetic_list_data_chimera_no_matching_trn.Rdata")
)

# ===============================================================================================
# DATA LOADING : REIN
# ===============================================================================================

rein_synthetic_list_data_mice <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = REIN_train,
      missing_rate = NULL,
      iterative_masking = FALSE,
      survival_data = TRUE,
      time_variable = "times",
      event_variable = "Censored",
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      verbose = TRUE
    )
    
  })["elapsed"]
  
  rein_synthetic_list_data_mice[[i]] <- synthetic_data[["synthetic_data"]]
}

# Save synthetic datasets
save(rein_synthetic_list_data_mice,
     file = here(results_dir, "rein_synthetic_list_data_mice.Rdata")
)

# for Holdout set 
rein_synthetic_list_data_mice_trn <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = REIN_trn,
      missing_rate = NULL,
      iterative_masking = FALSE,
      survival_data = TRUE,
      time_variable = "times",
      event_variable = "Censored",
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      verbose = TRUE
    )
    
  })["elapsed"]
  
  rein_synthetic_list_data_mice_trn[[i]] <- synthetic_data[["synthetic_data"]]
}


# Save synthetic datasets
save(rein_synthetic_list_data_mice_trn,
     file = here(results_dir, "rein_synthetic_list_data_mice_trn.Rdata")
)


# CHIMERA NO MATCHING
rein_synthetic_list_data_chimera_no_matching <- vector("list", n_synth)
execution_times <- numeric(n_synth)
missing_rate <- runif(1, min = 0.10, max = 0.30)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = REIN_train,
      missing_rate = missing_rate,
      iterative_masking = TRUE,
      survival_data = TRUE,
      time_variable = "times",
      event_variable = "Censored",
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      verbose = TRUE
    )
    
  })["elapsed"]
  
  rein_synthetic_list_data_chimera_no_matching[[i]] <- synthetic_data[["synthetic_data"]]
}


# Save synthetic datasets
save(rein_synthetic_list_data_chimera_no_matching,
     file = here(results_dir, "rein_synthetic_list_data_chimera_no_matching.Rdata")
)

# for Holdout set
rein_synthetic_list_data_chimera_no_matching_trn <- vector("list", n_synth)
execution_times <- numeric(n_synth)

for (i in seq_len(n_synth)) {
  execution_times[i] <- system.time({
    message(sprintf("▶ Generating dataset %d / %d (seed = %d)", i, n_synth, seeds[i]))
    synthetic_data <- Ablation_Analysis_generate_synthetic_data(
      data = REIN_trn,
      missing_rate = missing_rate,
      iterative_masking = TRUE,
      survival_data = TRUE,
      time_variable = "times",
      event_variable = "Censored",
      n_imputations = 1,
      max_iterations = 5,
      random_seed = seeds[i],
      verbose = TRUE
    )
    
  })["elapsed"]
  
  rein_synthetic_list_data_chimera_no_matching_trn[[i]] <- synthetic_data[["synthetic_data"]]
}


# Save synthetic datasets
save(rein_synthetic_list_data_chimera_no_matching_trn,
     file = here(results_dir, "rein_synthetic_list_data_chimera_no_matching_trn.Rdata")
)



# ============================================================
# Sensitivity analysis: CHIMERA across different p values
# ============================================================

# CHIMERA — Sensitivity analysis for missing_rate dataset: PIMA

n_synth <- 50

# Values of missing_rate to investigate
missing_rates <- c(0.05, 0.10, 0.15, 0.20, 0.25)

# Reproducibility
set.seed(123)
seeds <- sample(1:1000, n_synth * length(missing_rates))

# Store the seeds in a matrix: one row per missing_rate and one column per generation
seeds_matrix <- matrix(
  seeds,
  nrow = length(missing_rates),
  ncol = n_synth,
  byrow = TRUE
)

rownames(seeds_matrix) <- paste0(
  "p_",
  sprintf("%03d", round(missing_rates * 100))
)

colnames(seeds_matrix) <- paste0("generation_", seq_len(n_synth))


# ------------------------------------------------------------
# Initialize result lists
# ------------------------------------------------------------

# Synthetic datasets grouped by missing_rate
pima_synthetic_chimera_sensitivity_by_p <- setNames(
  vector("list", length(missing_rates)),
  rownames(seeds_matrix)
)


# ------------------------------------------------------------
# Generate 50 synthetic datasets per missing_rate
# ------------------------------------------------------------

for (j in seq_along(missing_rates)) {
  
  missing_rate <- missing_rates[j]
  
  # Name of the current parameter
  p_name <- paste0("p_", sprintf("%03d", round(missing_rate * 100)))
  
  message("\n")
  message("============================================================")
  message("CHIMERA sensitivity analysis: ", p_name)
  message("missing_rate = ", missing_rate)
  message("============================================================")
  
  # Initialize storage for the current missing_rate
  synthetic_list <- vector("list", n_synth)
  execution_times <- numeric(n_synth)
  
  for (i in seq_len(n_synth)) {
    
    current_seed <- seeds_matrix[j, i]
    
    message(
      sprintf(
        "Generating dataset %d/%d | missing_rate = %.2f | seed = %d",
        i,
        n_synth,
        missing_rate,
        current_seed
      )
    )
    
    # Measure execution time
    start_time <- Sys.time()
    
    results <- CHIMERA_generate(
      df = PIMA_train,
      missing_rate = missing_rate,
      is_survival = FALSE,
      timevar = NULL,
      statusvar = NULL,
      nb_imput = 1,
      nb_max_it = 5,
      seed = current_seed,
      verbose = TRUE
    )
    
    execution_times[i] <- as.numeric(
      difftime(Sys.time(), start_time, units = "secs")
    )
    
    # Store the synthetic dataset
    synthetic_list[[i]] <- results$synthetic_data
    
    message(
      sprintf(
        "Completed dataset %d/%d in %.2f seconds",
        i,
        n_synth,
        execution_times[i]
      )
    )
  }
  
  # Store the 50 datasets for this missing_rate
  pima_synthetic_chimera_sensitivity_by_p[[p_name]] <- synthetic_list
  message("Completed all generations for ", p_name)
}


# ------------------------------------------------------------
# Save all synthetic datasets
# ------------------------------------------------------------
# Save synthetic datasets
save(pima_synthetic_chimera_sensitivity_by_p,
     file = here(results_dir, "pima_synthetic_chimera_sensitivity_by_p.Rdata")
)


# for holdout Set and linkability attack
pima_synthetic_chimera_sensitivity_by_p_trn <- setNames(
  vector("list", length(missing_rates)),
  rownames(seeds_matrix)
)


# ------------------------------------------------------------
# Generate 50 synthetic datasets per missing_rate
# ------------------------------------------------------------

for (j in seq_along(missing_rates)) {
  
  missing_rate <- missing_rates[j]
  
  # Name of the current parameter
  p_name <- paste0("p_", sprintf("%03d", round(missing_rate * 100)))
  
  message("\n")
  message("============================================================")
  message("CHIMERA sensitivity analysis: ", p_name)
  message("missing_rate = ", missing_rate)
  message("============================================================")
  
  # Initialize storage for the current missing_rate
  synthetic_list <- vector("list", n_synth)
  execution_times <- numeric(n_synth)
  
  for (i in seq_len(n_synth)) {
    
    current_seed <- seeds_matrix[j, i]
    
    message(
      sprintf(
        "Generating dataset %d/%d | missing_rate = %.2f | seed = %d",
        i,
        n_synth,
        missing_rate,
        current_seed
      )
    )
    
    # Measure execution time
    start_time <- Sys.time()
    
    results <- CHIMERA_generate(
      df = PIMA_trn,
      missing_rate = missing_rate,
      is_survival = FALSE,
      timevar = NULL,
      statusvar = NULL,
      nb_imput = 1,
      nb_max_it = 5,
      seed = current_seed,
      verbose = TRUE
    )
    
    execution_times[i] <- as.numeric(
      difftime(Sys.time(), start_time, units = "secs")
    )
    
    # Store the synthetic dataset
    synthetic_list[[i]] <- results$synthetic_data
    
    message(
      sprintf(
        "Completed dataset %d/%d in %.2f seconds",
        i,
        n_synth,
        execution_times[i]
      )
    )
  }
  
  # Store the 50 datasets for this missing_rate
  pima_synthetic_chimera_sensitivity_by_p_trn[[p_name]] <- synthetic_list
  message("Completed all generations for ", p_name)
}

# Save synthetic datasets
save(pima_synthetic_chimera_sensitivity_by_p_trn,
     file = here(results_dir, "pima_synthetic_chimera_sensitivity_by_p_trn.Rdata")
)


# CHIMERA — Sensitivity analysis for missing_rate dataset: AIDS

# Synthetic datasets grouped by missing_rate
aids_synthetic_chimera_sensitivity_by_p <- setNames(
  vector("list", length(missing_rates)),
  rownames(seeds_matrix)
)


# ------------------------------------------------------------
# Generate 50 synthetic datasets per missing_rate
# ------------------------------------------------------------

for (j in seq_along(missing_rates)) {
  
  missing_rate <- missing_rates[j]
  
  # Name of the current parameter
  p_name <- paste0(
    "p_",
    sprintf("%03d", round(missing_rate * 100))
  )
  
  message("\n")
  message("============================================================")
  message("CHIMERA sensitivity analysis: AIDS — ", p_name)
  message("missing_rate = ", missing_rate)
  message("============================================================")
  
  # Initialize storage for the current missing_rate
  synthetic_list <- vector("list", n_synth)
  execution_times <- numeric(n_synth)
  
  for (i in seq_len(n_synth)) {
    
    current_seed <- seeds_matrix[j, i]
    
    message(
      sprintf(
        "Generating dataset %d/%d | missing_rate = %.2f | seed = %d",
        i,
        n_synth,
        missing_rate,
        current_seed
      )
    )
    
    # Measure execution time
    start_time <- Sys.time()
    
    results <- CHIMERA_generate(
      df = AIDS_train,
      missing_rate = missing_rate,
      is_survival = TRUE,
      timevar = "times",
      statusvar = "Censored",
      nb_imput = 1,
      nb_max_it = 5,
      seed = current_seed,
      verbose = TRUE
    )
    
    
    # Store the synthetic dataset
    synthetic_list[[i]] <- results$synthetic_data
    
    message(
      sprintf(
        "Completed dataset %d/%d in %.2f seconds",
        i,
        n_synth,
        execution_times[i]
      )
    )
  }
  
  # Store the 50 datasets for this missing_rate
  aids_synthetic_chimera_sensitivity_by_p[[p_name]] <- synthetic_list
  
  # Store execution times
  aids_execution_times_by_p[[p_name]] <- execution_times
  
  message("Completed all generations for ", p_name)
}


# ------------------------------------------------------------
# Save all synthetic datasets
# ------------------------------------------------------------

save(
  aids_synthetic_chimera_sensitivity_by_p,
  file = here(
    results_dir,
    "aids_synthetic_chimera_sensitivity_by_p.Rdata"
  )
)


# ------------------------------------------------------------
# Initialize result lists
# ------------------------------------------------------------

# Synthetic datasets grouped by missing_rate
aids_synthetic_chimera_sensitivity_by_p_trn <- setNames(
  vector("list", length(missing_rates)),
  rownames(seeds_matrix)
)


# ------------------------------------------------------------
# Generate 50 synthetic datasets per missing_rate
# ------------------------------------------------------------

for (j in seq_along(missing_rates)) {
  
  missing_rate <- missing_rates[j]
  
  # Name of the current parameter
  p_name <- paste0(
    "p_",
    sprintf("%03d", round(missing_rate * 100))
  )
  
  message("\n")
  message("============================================================")
  message("CHIMERA sensitivity analysis: AIDS — ", p_name)
  message("missing_rate = ", missing_rate)
  message("============================================================")
  
  # Initialize storage for the current missing_rate
  synthetic_list <- vector("list", n_synth)
  execution_times <- numeric(n_synth)
  
  for (i in seq_len(n_synth)) {
    
    current_seed <- seeds_matrix[j, i]
    
    message(
      sprintf(
        "Generating dataset %d/%d | missing_rate = %.2f | seed = %d",
        i,
        n_synth,
        missing_rate,
        current_seed
      )
    )
    
    # Measure execution time
    start_time <- Sys.time()
    
    synthetic_data <- CHIMERA_generate(
      df = AIDS_train,
      missing_rate = missing_rate,
      is_survival = TRUE,
      timevar = "times",
      statusvar = "Censored",
      nb_imput = 1,
      nb_max_it = 5,
      seed = current_seed,
      verbose = TRUE
    )
    
    
    # Store the synthetic dataset
    synthetic_list[[i]] <- synthetic_data
    
    message(
      sprintf(
        "Completed dataset %d/%d in %.2f seconds",
        i,
        n_synth,
        execution_times[i]
      )
    )
  }
  
  # Store the 50 datasets for this missing_rate
  aids_synthetic_chimera_sensitivity_by_p_trn[[p_name]] <- synthetic_list
  
  message("Completed all generations for ", p_name)
}


# ------------------------------------------------------------
# Save all synthetic datasets
# ------------------------------------------------------------

save(
  aids_synthetic_chimera_sensitivity_by_p_trn,
  file = here(
    results_dir,
    "aids_synthetic_chimera_sensitivity_by_p_trn.Rdata"
  )
)
