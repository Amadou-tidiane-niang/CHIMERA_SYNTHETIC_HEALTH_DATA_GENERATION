#===============================================================================================
# INITIAL SETUP
# ===============================================================================================

# Clear environment
rm(list = ls())

# -----------------------------------------------------------------------------------------------
# LOAD LIBRARIES
# -----------------------------------------------------------------------------------------------
library(here)         # File path management
library(tidyverse)    # Data manipulation
library(mice)         # Multiple imputation
library(readxl)       # Excel import


# Load custom functions
source(here("scripts", "0_functions.R"))

# DATA LOADING | Pima  ------------------------------------------------------------
pima_indians_diabetes <- read.csv(here("data", "pima-indians-diabetes.csv"))

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
write.csv(PIMA, here("data", "PIMA.csv"), row.names = FALSE)

# HoldoutSet : train and holdout set 
set.seed(123)
n <- nrow(PIMA)
idx_train <- sample(seq_len(n), size = 0.5 * n)

PIMA_trn <- PIMA[idx_train, ]
PIMA_hol  <- PIMA[-idx_train, ]

row.names(PIMA_trn) <- NULL
row.names(PIMA_hol) <- NULL

write.csv(PIMA_trn, here("data", "PIMA_trn.csv"), row.names = FALSE)
write.csv(PIMA_hol, here("data", "PIMA_hol.csv"), row.names = FALSE)


# ============================================================
# DATA LOADING | AIDS 
# ============================================================
AIDS <- read_delim(here("data", "aids_original_data.csv"), delim = ";", escape_double = FALSE, trim_ws = TRUE)|> as.data.frame()

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
write.csv(AIDS, here("data", "AIDS.csv"), row.names = FALSE)


# HoldoutSet : train and holdout set 
n <- nrow(AIDS)
idx_train <- sample(seq_len(n), size = 0.5 * n)

AIDS_trn <- AIDS[idx_train, ]
AIDS_hol  <- AIDS[-idx_train, ]

row.names(AIDS_trn) <- NULL
row.names(AIDS_hol) <- NULL

write.csv(AIDS_trn, here("data", "AIDS_trn.csv"), row.names = FALSE)
write.csv(AIDS_hol, here("data", "AIDS_hol.csv"), row.names = FALSE)

# ============================================================
# DATA LOADING | AIDS 
# ============================================================
df_rein_with_na <- read.csv(here("data", "df_rein_with_na.csv"))


# --- Initial MICE setup
ini <- mice(df_rein_with_na, maxit = 0)
meth <- ini$method
pred <- ini$predictorMatrix
meth[names(df_rein_with)[sapply(df_rein_with_na, is.numeric)]] <- "pmm"

# --- Impute missing data
df_rein_without_na <- complete(mice(df_rein_with_na, m = 1, maxit = 5, seed = 123,method = meth, predictorMatrix = pred, printFlag = TRUE))

# Save dataset
write.csv(df_rein_without_na, here("data", "df_rein_without_na.csv"), row.names = FALSE)

# HoldoutSet : train and holdout set 
n <- nrow(df_rein_without_na)
idx_train <- sample(seq_len(n), size = 0.5 * n)

REIN_trn <- df_rein_without_na[idx_train, ]
REIN_hol  <- df_rein_without_na[-idx_train, ]

row.names(REIN_trn) <- NULL
row.names(REIN_hol) <- NULL

write.csv(REIN_trn, here("data", "REIN_trn.csv"), row.names = FALSE)
write.csv(REIN_hol, here("data", "REIN_hol.csv"), row.names = FALSE)


# create train-test set
set.seed(123)
n <- nrow(df_rein_with_na)
idx_train <- sample(seq_len(n), size = 0.7 * n)

df_rein_train_with_na_score <- df_rein_with_na[idx_train, ]
df_rein_test_with_na_score  <- df_rein_with_na[-idx_train, ]

row.names(df_rein_train_with_na_score) <- NULL
row.names(df_rein_test_with_na_score) <- NULL

write.csv(df_rein_train_with_na_score, here("data", "df_rein_train_with_na_score.csv"), row.names = FALSE)
write.csv(df_rein_test_with_na_score, here("data", "df_rein_test_with_na_score.csv"), row.names = FALSE)

# --- Initial MICE setup
ini <- mice(df_rein_train_with_na_score, maxit = 0)
meth <- ini$method
pred <- ini$predictorMatrix
meth[names(df_rein_train_with_na_score)[sapply(df_rein_train_with_na_score, is.numeric)]] <- "pmm"

# --- Impute missing data
df_rein_train_without_na_score <- complete(mice(df_rein_train_with_na_score, m = 1, maxit = 5, seed = 123,method = meth, predictorMatrix = pred, printFlag = TRUE))

# Save dataset
write.csv(df_rein_train_without_na_score, here("data", "df_rein_train_without_na_score.csv"), row.names = FALSE)

# HoldoutSet : train and holdout set 
n <- nrow(df_rein_train_without_na_score)
idx_train <- sample(seq_len(n), size = 0.5 * n)

REIN_trn_score <- df_rein_train_without_na_score[idx_train, ]
REIN_hol_score  <- df_rein_train_without_na_score[-idx_train, ]

row.names(REIN_trn_score) <- NULL
row.names(REIN_hol_score) <- NULL

write.csv(REIN_trn_score, here("data", "REIN_trn_score.csv"), row.names = FALSE)
write.csv(REIN_hol_score, here("data", "REIN_hol_score.csv"), row.names = FALSE)
