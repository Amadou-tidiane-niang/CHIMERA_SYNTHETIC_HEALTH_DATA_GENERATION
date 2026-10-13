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
  library(parallel)
  library(future)
  library(future.apply)
})

# Load custom functions (prepare_survival_data_time, compute_calibration,
# compute_calibration_synth, etc.)
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


# -----------------------------------------------------------------------------------------------
# GLOBAL PARALLELIZATION SETTINGS
# -----------------------------------------------------------------------------------------------

N_CORES <- max(1, parallel::detectCores() - 1)

# boot() requires a platform-appropriate cluster type
BOOT_PARALLEL_TYPE <- if (.Platform$OS.type == "windows") "snow" else "multicore"

# future_lapply() plan, used only where no internal parallelism is already at play
# (see the preprocessing_auc loops below)
future::plan(future::multisession, workers = N_CORES)

cat(sprintf("Parallelization enabled on %d cores (boot: %s)\n", N_CORES, BOOT_PARALLEL_TYPE))

#------------------------------------------------------------------------------------------------
## REIN DATASET
#------------------------------------------------------------------------------------------------

args <- c("data", "results")

data_dir <- normalizePath(args[1], mustWork = TRUE)

if (!dir.exists(args[2])) {
  dir.create(args[2], recursive = TRUE)
}

results_dir <- normalizePath(args[2], mustWork = TRUE)

REIN_train <- read.csv(here(data_dir, "REIN_train.csv")) |> mutate_if(is.character, as.factor)
REIN_test <- read.csv(here(data_dir, "REIN_test.csv")) |> mutate_if(is.character, as.factor)


n_synth <- 50
set.seed(123)
seeds <- sample(1:1000, n_synth)


# REIN synthetic data
load(here(results_dir, "rein_synthetic_list_data_chimera.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_synthpop.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_ctgan.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_gaussiancopula.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_tvae.Rdata"))
load(here(results_dir, "rein_synthetic_list_data_tabdiff.Rdata"))


synthetic_raw <- list(
  chimera        = rein_synthetic_list_data_chimera,
  synthpop       = rein_synthetic_list_data_synthpop,
  ctgan          = rein_synthetic_list_data_ctgan,
  gaussiancopula = rein_synthetic_list_data_gaussiancopula,
  tvae           = rein_synthetic_list_data_tvae,
  tabdiff        = rein_synthetic_list_data_tabdiff
)

method_labels <- c(
  chimera        = "CHIMERA",
  synthpop       = "SYNTHPOP",
  ctgan          = "CTGAN",
  gaussiancopula = "GAUSSIAN COPULA",
  tvae           = "TVAE",
  tabdiff        = "TABDIFF"
)


list2env(setNames(synthetic_raw, paste0("rein_synthetic_list_data_", names(synthetic_raw))), envir = environment())


###############################################################################################
## MULTIPLE IMPUTATION PREPROCESSING
###############################################################################################

nb_imputations <- 50

# Original dataset
results.original.imp <- preprocessing(REIN_train, nb_imputations = nb_imputations)
save(results.original.imp, file = here(results_dir, "results.original.imp.Rdata"))

# OPTIMIZATION: single loop over all synthetic methods.
# Stays sequential at the outer level (see note at the top of the file) because
# Imputed_data_function() already parallelizes imputation internally via parlmice().
results.imp <- list()

for (m in names(synthetic_raw)) {
  cat(sprintf("== Multiple imputation: %s ==\n", method_labels[[m]]))
  
  results.imp[[m]] <- lapply(
    synthetic_raw[[m]],
    preprocessing,
    nb_imputations = nb_imputations
  )
  
  save(
    list = "results.imp",
    file = here(results_dir, paste0("results.", m, ".imp.Rdata"))
  )
}


results.chimera.imp        <- results.imp$chimera
results.synthpop.imp       <- results.imp$synthpop
results.ctgan.imp          <- results.imp$ctgan
results.gaussiancopula.imp <- results.imp$gaussiancopula
results.tvae.imp           <- results.imp$tvae
results.tabdiff.imp        <- results.imp$tabdiff


###############################################################################################
## RISK SCORE ANALYSIS
## DEATH WITHIN 3 MONTHS AFTER INITIATION OF RENAL REPLACEMENT THERAPY
###############################################################################################
load(here(results_dir, "results.original.imp.Rdata"))
load(here(results_dir, "results.original.Rdata"))

# Original dataset
results.original.score <- Score(results.original.imp)
save(results.original.score, file = here(results_dir, "results.original.score.Rdata"))

# OPTIMIZATION: single loop over the 6 methods. Stays sequential at the outer level because
# Score() already parallelizes the bootstrap internally (boot(..., parallel = ...)).

results.score <- list()

for (m in names(synthetic_raw)) {
  
  results.score[[m]] <- vector("list", length(results.imp[[m]]))
  
  for (i in seq_along(results.imp[[m]])) {
    results.score[[m]][[i]] <- Score(results.imp[[m]][[i]])
    message(sprintf("%s iteration %d completed", method_labels[[m]], i))
  }
  
  assign(paste0("results.", m), results.score[[m]])
  save(list = paste0("results.", m), file = here(results_dir, paste0("results.", m, ".Rdata")))
}

results.chimera        <- results.score$chimera
results.synthpop       <- results.score$synthpop
results.ctgan          <- results.score$ctgan
results.gaussiancopula <- results.score$gaussiancopula
results.tvae           <- results.score$tvae
results.tabdiff        <- results.score$tabdiff


###############################################################################################
## VARIABLE SELECTION PERFORMANCE METRICS
## Sensitivity, Specificity, and Cohen's Kappa
###############################################################################################

tbl_df_real <- results.original
threshold   <- 70

selection_metrics <- list()

for (m in names(synthetic_raw)) {
  
  res_list_m <- lapply(results.score[[m]], function(df_synth) {
    compute_selection_metrics(
      tbl_df_real, df_synth,
      real_col = "nb_significative_coef.x", synth_col = "nb_significative_coef.y",
      threshold = threshold
    )
  })
  
  selection_metrics[[m]] <- do.call(rbind, lapply(res_list_m, as.data.frame))
  
  assign(paste0(m, "_selection_metrics_df"), selection_metrics[[m]])
  save(
    list = paste0(m, "_selection_metrics_df"),
    file = here(results_dir, paste0(m, "_selection_metrics_df.Rdata"))
  )
}


###############################################################################################
## BEST SYNTHETIC DATASET SELECTION
###############################################################################################
load(here(results_dir, "best_all_by_M50.Rdata"))
best_rein <- best_all_by_M50 |>
  dplyr::filter(dataset == "REIN") |>
  dplyr::mutate(method_key = tolower(gsub("[^A-Za-z]", "", method)))

best_row_chimera        <- best_rein |> dplyr::filter(method_key == "chimera") |> dplyr::pull(run)
best_row_synthpop       <- best_rein |> dplyr::filter(method_key == "synthpop") |> dplyr::pull(run)
best_row_ctgan          <- best_rein |> dplyr::filter(method_key == "ctgan") |> dplyr::pull(run)
best_row_gaussiancopula <- best_rein |> dplyr::filter(method_key == "gaussiancopula") |> dplyr::pull(run)
best_row_tvae           <- best_rein |> dplyr::filter(method_key == "tvae") |> dplyr::pull(run)
best_row_tabdiff        <- best_rein |> dplyr::filter(method_key == "tabdiff") |> dplyr::pull(run)

best_row <- c(
  chimera        = best_row_chimera,
  synthpop       = best_row_synthpop,
  ctgan          = best_row_ctgan,
  gaussiancopula = best_row_gaussiancopula,
  tvae           = best_row_tvae,
  tabdiff        = best_row_tabdiff
)

###############################################################################################
## LOAD BEST IMPUTED DATASETS
###############################################################################################

df.imp.original       <- results.original.imp
df.imp.list           <- Map(function(imp, i) imp[[i]], results.imp, best_row)

df.imp.chimera        <- df.imp.list$chimera
df.imp.synthpop       <- df.imp.list$synthpop
df.imp.ctgan          <- df.imp.list$ctgan
df.imp.gaussiancopula <- df.imp.list$gaussiancopula
df.imp.tvae           <- df.imp.list$tvae
df.imp.tabdiff        <- df.imp.list$tabdiff


###############################################################################################
## VARIABLE SELECTION + AUC (imputed test set) + POOLED ROC
## OPTIMIZATION: compute_auc_ci() computes predictions ONCE; the pooled ROC curve is then
## obtained via rowMeans() on those same predictions, without refitting any glm().
###############################################################################################

df_test <- REIN_test
ini  <- mice(df_test, maxit = 0)
meth <- ini$method
pred <- ini$predictorMatrix
meth[names(df_test)[sapply(df_test, is.numeric)]] <- "pmm"

df_test <- complete(mice(df_test, m = 1, maxit = 5, seed = 123, method = meth, predictorMatrix = pred, printFlag = TRUE))
df_test <- df_test |>
  rename(Censorship = Censored, Follow_up_time = times)

df_test$death <- as.factor(ifelse(df_test$Follow_up_time < 90 & df_test$Censorship == 1, "Yes", "No"))

# -----------------------------------------------------------------------------------------------
# ORIGINAL
# -----------------------------------------------------------------------------------------------

df.nbSignCoef.original <- results.original |> filter(nb_significative_coef > threshold)
name_variable_original <- clean_variable_names(df.nbSignCoef.original$facteur_risque, df.imp.original[[1]])

auc_ci_original <- compute_auc_ci(
  df.imp.original, df_test, outcome = "death",
  method_name = "REAL", name_variable = name_variable_original
)

auc_original <- auc_ci_original$metrics

# -----------------------------------------------------------------------------------------------
# SYNTHETIC METHODS
# -----------------------------------------------------------------------------------------------

name_variable_list <- list()

for (m in names(synthetic_raw)) {
  
  df.nbSignCoef <- lapply(
    results.score[[m]],
    function(x) {
      x |>
        dplyr::filter(nb_significative_coef > threshold)
    }
  )
  
  
  name_variable <- lapply(
    seq_along(df.nbSignCoef),
    function(i) {
      
      clean_variable_names(
        df.nbSignCoef[[i]]$facteur_risque,
        df.imp.original[[1]]
      )
    }
  )
  
  name_variable_list[[m]] <- name_variable
  
}


name_variable_chimera        <- name_variable_list$chimera
name_variable_synthpop       <- name_variable_list$synthpop
name_variable_ctgan          <- name_variable_list$ctgan
name_variable_gaussiancopula <- name_variable_list$gaussiancopula
name_variable_tvae           <- name_variable_list$tvae
name_variable_tabdiff        <- name_variable_list$tabdiff


# CHIMERA
auc_ci_chimera <- lapply(seq_len(50), function(i) {
  compute_auc_ci(
    results.chimera.imp[[i]],
    df_test,
    outcome = "death",
    method_name = "CHIMERA",
    name_variable = name_variable_chimera[[i]]
  )
})


# SYNTHPOP
auc_ci_synthpop <- lapply(seq_len(50), function(i) {
  compute_auc_ci(
    results.synthpop.imp[[i]],
    df_test,
    outcome = "death",
    method_name = "SYNTHPOP",
    name_variable = name_variable_synthpop[[i]]
  )
})


# CTGAN
auc_ci_ctgan <- lapply(seq_len(50), function(i) {
  compute_auc_ci(
    results.ctgan.imp[[i]],
    df_test,
    outcome = "death",
    method_name = "CTGAN",
    name_variable = name_variable_ctgan[[i]]
  )
})


# GAUSSIAN COPULA
auc_ci_gaussiancopula <- lapply(seq_len(50), function(i) {
  compute_auc_ci(
    results.gaussiancopula.imp[[i]],
    df_test,
    outcome = "death",
    method_name = "GAUSSIAN COPULA",
    name_variable = name_variable_gaussiancopula[[i]]
  )
})

# TVAE
auc_ci_tvae <- lapply(seq_len(50), function(i) {
  compute_auc_ci(
    results.tvae.imp[[i]],
    df_test,
    outcome = "death",
    method_name = "TVAE",
    name_variable = name_variable_tvae[[i]]
  )
})


# TABDIFF
auc_ci_tabdiff <- lapply(seq_len(50), function(i) {
  compute_auc_ci(
    results.tabdiff.imp[[i]],
    df_test,
    outcome = "death",
    method_name = "CTGAN",
    name_variable = name_variable_tabdiff[[i]]
  )
})


# Combine results from the six synthetic data generation methods
auc_ci_all <- list(
  "CHIMERA"          = auc_ci_chimera,
  "SYNTHPOP"         = auc_ci_synthpop,
  "CTGAN"            = auc_ci_ctgan,
  "GAUSSIAN COPULA"  = auc_ci_gaussiancopula,
  "TVAE"             = auc_ci_tvae,
  "TABDIFF"          = auc_ci_tabdiff
)

# Convert synthetic results into a single dataframe
auc_ci_df_synthetic <- imap_dfr(auc_ci_all, function(results, method) {
  
  map2_dfr(results, seq_along(results), function(dataset_results, dataset_id) {
    
    as.data.frame(dataset_results) %>%
      mutate(
        method = method,
        synthetic_dataset = dataset_id,
        .before = 1
      )
  })
})


# Convert real-data results into a dataframe
auc_ci_df_real <- auc_ci_original %>%
  mutate(
    method = "REAL",
    synthetic_dataset = 0,
    .before = 1
  )

# Combine real and synthetic results
auc_ci_df <- bind_rows(
  auc_ci_df_real,
  auc_ci_df_synthetic
)

auc_final_model_rein <- auc_ci_df|>select(-c("method"))
save(auc_final_model_rein, file = here(results_dir, "auc_final_model_rein.Rdata"))



# -------------------------------------------------------------------------
# Compute ROC curves for all methods
# -------------------------------------------------------------------------
roc_original <- compute_pooled_roc(
  imputed_data = df.imp.original,
  predictor_names = name_variable_original,
  test_data = df_test,
  outcome = "death",
  n_runs = 50
)


roc_chimera <- compute_method_rocs(
  imputed_datasets = results.chimera.imp,
  variable_list = name_variable_chimera,
  test_data = df_test,
  outcome = "death",
  n_runs = 50
)


roc_synthpop <- compute_method_rocs(
  imputed_datasets = results.synthpop.imp,
  variable_list = name_variable_synthpop,
  test_data = df_test,
  outcome = "death",
  n_runs = 50
)


roc_ctgan <- compute_method_rocs(
  imputed_datasets = results.ctgan.imp,
  variable_list = name_variable_ctgan,
  test_data = df_test,
  outcome = "death",
  n_runs = 50
)


roc_gaussiancopula <- compute_method_rocs(
  imputed_datasets = results.gaussiancopula.imp,
  variable_list = name_variable_gaussiancopula,
  test_data = df_test,
  outcome = "death",
  n_runs = 50
)


roc_tvae <- compute_method_rocs(
  imputed_datasets = results.tvae.imp,
  variable_list = name_variable_tvae,
  test_data = df_test,
  outcome = "death",
  n_runs = 50
)


roc_tabdiff <- compute_method_rocs(
  imputed_datasets = results.tabdiff.imp,
  variable_list = name_variable_tabdiff,
  test_data = df_test,
  outcome = "death",
  n_runs = 50
)

roc_curves <- list(
  REAL = roc_original,
  CHIMERA = roc_chimera,
  SYNTHPOP = roc_synthpop,
  CTGAN = roc_ctgan,
  GAUSSIAN_COPULA = roc_gaussiancopula,
  TVAE = roc_tvae,
  TABDIFF = roc_tabdiff
)


save(roc_curves, file = here(results_dir, "roc_curve_score.Rdata"))

# ============================================================
# AUC pour les 6 méthodes de données synthétiques
# ============================================================

# Liste des méthodes
methods <- c(
  "chimera",
  "synthpop",
  "ctgan",
  "gaussian_copula",
  "tvae",
  "tabdiff"
)

# Noms à afficher dans les résultats
method_names <- c(
  "CHIMERA",
  "SYNTHPOP",
  "CTGAN",
  "GAUSSIAN COPULA",
  "TVAE",
  "TABDIFF"
)

# ============================================================
# CHIMERA
# ============================================================

synthetic_data_list_chimera <- lapply(
  rein_synthetic_list_data_chimera,
  function(df) preprocessing(df, nb_imputations = 1)
)

df.nbSignCoef.chimera <- results.chimera[[best_row_chimera]] |>
  filter(nb_significative_coef > 70)

name_variable_chimera <- clean_variable_names(
  df.nbSignCoef.chimera$facteur_risque,
  df.imp.chimera[[1]]
)

synthetic_data_list_chimera <- lapply(
  synthetic_data_list_chimera,
  function(x) x[[1]]
)

auc_chimera <- compute_auc_ci(
  synthetic_data_list_chimera,
  df_test,
  outcome = "death",
  method_name = "CHIMERA",
  name_variable = name_variable_chimera
)


# ============================================================
# SYNTHPOP
# ============================================================

synthetic_data_list_synthpop <- lapply(
  rein_synthetic_list_data_synthpop,
  function(df) preprocessing(df, nb_imputations = 1)
)

df.nbSignCoef.synthpop <- results.synthpop[[best_row_synthpop]] |>
  filter(nb_significative_coef > 70)

name_variable_synthpop <- clean_variable_names(
  df.nbSignCoef.synthpop$facteur_risque,
  df.imp.synthpop[[1]]
)

synthetic_data_list_synthpop <- lapply(
  synthetic_data_list_synthpop,
  function(x) x[[1]]
)

auc_synthpop <- compute_auc_ci(
  synthetic_data_list_synthpop,
  df_test,
  outcome = "death",
  method_name = "SYNTHPOP",
  name_variable = name_variable_synthpop
)


# ============================================================
# CTGAN
# ============================================================

synthetic_data_list_ctgan <- lapply(
  rein_synthetic_list_data_ctgan,
  function(df) preprocessing(df, nb_imputations = 1)
)

df.nbSignCoef.ctgan <- results.ctgan[[best_row_ctgan]] |>
  filter(nb_significative_coef > 70)

name_variable_ctgan <- clean_variable_names(
  df.nbSignCoef.ctgan$facteur_risque,
  df.imp.ctgan[[1]]
)

synthetic_data_list_ctgan <- lapply(
  synthetic_data_list_ctgan,
  function(x) x[[1]]
)

auc_ctgan <- compute_auc_ci(
  synthetic_data_list_ctgan,
  df_test,
  outcome = "death",
  method_name = "CTGAN",
  name_variable = name_variable_ctgan
)


# ============================================================
# GAUSSIAN COPULA
# ============================================================

synthetic_data_list_gaussiancopula <- lapply(
  rein_synthetic_list_data_gaussiancopula,
  function(df) preprocessing(df, nb_imputations = 1)
)

df.nbSignCoef.gaussian_copula <- results.gaussiancopula[[best_row_gaussiancopula]] |>
  filter(nb_significative_coef > 70)

name_variable_gaussian_copula <- clean_variable_names(
  df.nbSignCoef.gaussian_copula$facteur_risque,
  df.imp.gaussiancopula[[1]]
)

synthetic_data_list_gaussiancopula <- lapply(
  synthetic_data_list_gaussiancopula,
  function(x) x[[1]]
)

auc_gaussian_copula <- compute_auc_ci(
  synthetic_data_list_gaussiancopula,
  df_test,
  outcome = "death",
  method_name = "GAUSSIAN COPULA",
  name_variable = name_variable_gaussian_copula
)


# ============================================================
# TVAE
# ============================================================

synthetic_data_list_tvae <- lapply(
  rein_synthetic_list_data_tvae,
  function(df) preprocessing(df, nb_imputations = 1)
)

df.nbSignCoef.tvae <- results.tvae[[best_row_tvae]] |>
  filter(nb_significative_coef > 70)

name_variable_tvae <- clean_variable_names(
  df.nbSignCoef.tvae$facteur_risque,
  df.imp.tvae[[1]]
)

synthetic_data_list_tvae <- lapply(
  synthetic_data_list_tvae,
  function(x) x[[1]]
)

auc_tvae <- compute_auc_ci(
  synthetic_data_list_tvae,
  df_test,
  outcome = "death",
  method_name = "TVAE",
  name_variable = name_variable_tvae
)


# ============================================================
# TABDIFF
# ============================================================

synthetic_data_list_tabdiff <- lapply(
  rein_synthetic_list_data_tabdiff,
  function(df) preprocessing(df, nb_imputations = 1)
)

df.nbSignCoef.tabdiff <- results.tabdiff[[best_row_tabdiff]] |>
  filter(nb_significative_coef > 70)

name_variable_tabdiff <- clean_variable_names(
  df.nbSignCoef.tabdiff$facteur_risque,
  df.imp.tabdiff[[1]]
)

synthetic_data_list_tabdiff <- lapply(
  synthetic_data_list_tabdiff,
  function(x) x[[1]]
)

auc_tabdiff <- compute_auc_ci(
  synthetic_data_list_tabdiff,
  df_test,
  outcome = "death",
  method_name = "TABDIFF",
  name_variable = name_variable_tabdiff
)


# ============================================================
# REGROUPER LES 6 MÉTHODES
# ============================================================

auc_all <- bind_rows(
  auc_chimera,
  auc_synthpop,
  auc_ctgan,
  auc_gaussian_copula,
  auc_tvae,
  auc_tabdiff
)

save(auc_all, file = here(results_dir, "auc_all_syn.Rdata"))

###############################################################################################
## CALIBRATION
###############################################################################################

synthetic_raw <- list(
  chimera        = synthetic_data_list_chimera,
  synthpop       = synthetic_data_list_synthpop,
  ctgan          = synthetic_data_list_ctgan,
  gaussiancopula = synthetic_data_list_gaussiancopula,
  tvae           = synthetic_data_list_tvae,
  tabdiff        = synthetic_data_list_tabdiff
)

calibration_results <- list()

for (m in names(synthetic_raw)) {
  calibration_results[[m]] <- compute_calibration_synth(
    results.score[[m]], df_test, "death", synthetic_raw[[m]], threshold = threshold
  )
  assign(paste0("res.", m, ".cal"), calibration_results[[m]])
  save(list = paste0("res.", m, ".cal"), file = here(results_dir, paste0("res.", m, ".cal.Rdata")))
}


# Original
res.original.cal <- compute_calibration(results.original, df_test, "death", df.imp.original, threshold = threshold)
save(res.original.cal, file = here(results_dir, "res.original.cal.Rdata"))


###############################################################################################
## END
###############################################################################################

