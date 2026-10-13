#===============================================================================================
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
  library(mice)         # Multiple imputation
  library(readxl)       # Excel import
  
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


# DATA LOADING | Pima  ------------------------------------------------------------
pima_indians_diabetes <- read.csv(here(data_dir, "pima-indians-diabetes.csv"))

# Rename columns
names(pima_indians_diabetes) <- c("Number_of_pregnancies",
                                  "Plasma_glucose_concentration",
                                  "Diastolic_blood_pressure",
                                  "Triceps_skin_fold_thickness",
                                  "Two_hour_serum_insulin",
                                  "Body_mass_index",
                                  "Diabetes_pedigree_function",
                                  "Age",
                                  "Diabetes_diagnosis")


# Display structure
str(pima_indians_diabetes) 

# Convert Outcome to factor
pima_indians_diabetes$Diabetes_diagnosis <- as.factor(pima_indians_diabetes$Diabetes_diagnosis) 
levels(pima_indians_diabetes$Diabetes_diagnosis)<-c("No","Yes")
PIMA <- pima_indians_diabetes

# Save original dataset
write.csv(PIMA, here(data_dir, "PIMA.csv"), row.names = FALSE)

# Train and test set 
set.seed(123)
n <- nrow(PIMA)
idx_train <- sample(seq_len(n), size = 0.7 * n)

PIMA_train <- PIMA[idx_train, ]
PIMA_test  <- PIMA[-idx_train, ]

row.names(PIMA_train) <- NULL
row.names(PIMA_test) <- NULL

write.csv(PIMA_train, here(data_dir, "PIMA_train.csv"), row.names = FALSE)
write.csv(PIMA_test, here(data_dir, "PIMA_test.csv"), row.names = FALSE)


# HoldoutSet : train and holdout set 
set.seed(123)
n <- nrow(PIMA)
idx_train <- sample(seq_len(n), size = 0.5 * n)

PIMA_trn <- PIMA[idx_train, ]
PIMA_hol  <- PIMA[-idx_train, ]

row.names(PIMA_trn) <- NULL
row.names(PIMA_hol) <- NULL

write.csv(PIMA_trn, here(data_dir, "PIMA_trn.csv"), row.names = FALSE)
write.csv(PIMA_hol, here(data_dir, "PIMA_hol.csv"), row.names = FALSE)

# ============================================================
# DATA LOADING | AIDS 
# ============================================================
AIDS <- read_delim(here(data_dir, "aids_original_data.csv"), delim = ";", escape_double = FALSE, trim_ws = TRUE)|> as.data.frame()

# Drop identifiers and variables that will not be used
AIDS <- AIDS|> dplyr::select(-c(pidnum,zprior,cd420,cd496,cd820,arms,strat))

# Variable renaming  (for readability)
AIDS <- AIDS %>%
  rename(
    Age=age,
    Weight = wtkg,
    Hemophilia = hemo,
    Homosexuality = homo,
    Intravenous_drug_use = drugs,
    Karnofsky_score = karnof,
    Prior_opportunistic_infections = oprior,
    Zidovudine_30mg = z30,
    Previous_antiretroviral_days = preanti,
    Race = race,
    Sex = gender,
    Antiretroviral_history = str2,
    Symptomatic_indicator = symptom,
    Treatment_indicator = treat,
    Treatment_discontinuation = offtrt,
    CD4_baseline = cd40,
    Clinical_progression_observed = r,
    CD8_baseline = cd80,
    Censored = cens,
    times = days
  )


# factor conversion
AIDS <- AIDS %>%
  mutate(
    Hemophilia = factor(Hemophilia, levels = c(0, 1), labels = c("No", "Yes")),
    Homosexuality = factor(Homosexuality, levels = c(0, 1), labels = c("No", "Yes")),
    Intravenous_drug_use = factor(Intravenous_drug_use, levels = c(0, 1), labels = c("No", "Yes")),
    Prior_opportunistic_infections = factor(Prior_opportunistic_infections, levels = c(0, 1), labels = c("No", "Yes")),
    Zidovudine_30mg = factor(Zidovudine_30mg, levels = c(0, 1), labels = c("No", "Yes")),
    Race = factor(Race, levels = c(0, 1), labels = c("White", "Others")),
    Sex = factor(Sex, levels = c(0, 1), labels = c("Female", "Male")),
    Antiretroviral_history = factor(Antiretroviral_history, levels = c(0, 1), labels = c("Naive", "Experienced")),
    Symptomatic_indicator = factor(Symptomatic_indicator, levels = c(0, 1), labels = c("No", "Yes")),
    Treatment_indicator = factor(Treatment_indicator, levels = c(0, 1), labels = c("ZDV only", "Others")),
    Treatment_discontinuation = factor(Treatment_discontinuation, levels = c(0, 1), labels = c("No", "Yes")),
    Clinical_progression_observed = factor(Clinical_progression_observed, levels = c(0, 1), labels = c("No", "Yes"))
  )


# Save dataset
write.csv(AIDS, here(data_dir, "AIDS.csv"), row.names = FALSE)


# Train and holdout set 
n <- nrow(AIDS)
set.seed(123)
idx_train <- sample(seq_len(n), size = 0.7 * n)

AIDS_train <- AIDS[idx_train, ]
AIDS_test  <- AIDS[-idx_train, ]

row.names(AIDS_train) <- NULL
row.names(AIDS_test) <- NULL

write.csv(AIDS_train, here(data_dir, "AIDS_train.csv"), row.names = FALSE)
write.csv(AIDS_test, here(data_dir, "AIDS_test.csv"), row.names = FALSE)

# HoldoutSet : train and holdout set 
set.seed(123)
idx_train <- sample(seq_len(n), size = 0.5 * n)

AIDS_trn <- AIDS[idx_train, ]
AIDS_hol  <- AIDS[-idx_train, ]

row.names(AIDS_trn) <- NULL
row.names(AIDS_hol) <- NULL

write.csv(AIDS_trn, here(data_dir, "AIDS_trn.csv"), row.names = FALSE)
write.csv(AIDS_hol, here(data_dir, "AIDS_hol.csv"), row.names = FALSE)


# ============================================================
# DATA LOADING | REIN
# ============================================================
REIN <- read.csv(here(data_dir, "REIN.csv")) |> mutate_if(is.character, as.factor)

# Train and test set 
n <- nrow(REIN)
set.seed(123)
idx_train <- sample(seq_len(n), size = 0.7 * n)

REIN_train <- REIN[idx_train, ]
REIN_test  <- REIN[-idx_train, ]

row.names(REIN_train) <- NULL
row.names(REIN_test) <- NULL

write.csv(REIN_train, here(data_dir, "REIN_train.csv"), row.names = FALSE)
write.csv(REIN_test, here(data_dir, "REIN_test.csv"), row.names = FALSE)


# HoldoutSet : train and holdout set 
set.seed(123)
idx_train <- sample(seq_len(n), size = 0.5 * n)

REIN_trn <- REIN[idx_train, ]
REIN_hol  <- REIN[-idx_train, ]

row.names(REIN_trn) <- NULL
row.names(REIN_hol) <- NULL

write.csv(REIN_trn, here(data_dir, "REIN_trn.csv"), row.names = FALSE)
write.csv(REIN_hol, here(data_dir, "REIN_hol.csv"), row.names = FALSE)




