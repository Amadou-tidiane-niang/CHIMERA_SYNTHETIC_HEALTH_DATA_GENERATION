# ===============================================================================================
# INITIAL SETUP
# ===============================================================================================

# Clear environment
rm(list = ls())

# -----------------------------------------------------------------------------------------------
# LOAD LIBRARIES
# -----------------------------------------------------------------------------------------------
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
library(gtsummary)

# Load custom functions
source(here("scripts", "0_functions.R"))                   

# =========================================================
# Load datasets
# =========================================================
PIMA <- read.csv(here("data", "PIMA.csv"))
PIMA <- PIMA|> mutate_if(is.character,as.factor)

load(here("results", "PIMA_Metrics_T1_T9_chimera.Rdata"))
load(here("results", "PIMA_Metrics_T1_T9_synthpop.Rdata"))
load(here("results", "PIMA_Metrics_T1_T9_ctgan.Rdata"))

# Select the dataset minimizing the global distance score
# (lower distance = better overall performance)

best_row_chimera  <- which.min(PIMA_Metrics_T1_T9_chimera$L_total)
best_row_synthpop <- which.min(PIMA_Metrics_T1_T9_synthpop$L_total)
best_row_ctgan    <- which.min(PIMA_Metrics_T1_T9_ctgan$L_total)

load(here("results", "pima_synthetic_list_data_chimera.Rdata"))
load(here("results", "pima_synthetic_list_data_synthpop.Rdata"))
load(here("results", "pima_synthetic_list_data_ctgan.Rdata"))


# Best synthetic data
PIMA_CHIMERA <- pima_synthetic_list_data_chimera[[best_row_chimera]]
PIMA_SYNTHPOP <- pima_synthetic_list_data_synthpop[[best_row_synthpop]]
PIMA_CTGAN <- pima_synthetic_list_data_ctgan[[best_row_ctgan]]

# Save synthetic datasets
save(PIMA,file = here("data", "PIMA.Rdata"))
save(PIMA_CHIMERA,file = here("data", "PIMA_CHIMERA.Rdata"))
save(PIMA_SYNTHPOP,file = here("data", "PIMA_SYNTHPOP.Rdata"))
save(PIMA_CTGAN,file = here("data", "PIMA_CTGAN.Rdata"))

# Convert Outcome to factor
PIMA$Diabetes_diagnosis <- factor(PIMA$Diabetes_diagnosis,levels=c("No","Yes"),labels=c(0,1))
PIMA_CHIMERA$Diabetes_diagnosis <- factor(PIMA_CHIMERA$Diabetes_diagnosis,levels=c("No","Yes"),labels=c(0,1))
PIMA_SYNTHPOP$Diabetes_diagnosis <- factor(PIMA_SYNTHPOP$Diabetes_diagnosis,levels=c("No","Yes"),labels=c(0,1))
PIMA_CTGAN$Diabetes_diagnosis <- factor(PIMA_CTGAN$Diabetes_diagnosis,levels=c("No","Yes"),labels=c(0,1))


# Analytical utility  ---------------------------------------------------------------------
# ORIGINAL
PIMA$Age <- PIMA$Age/10
PIMA$Plasma_glucose_concentration <- PIMA$Plasma_glucose_concentration/5
PIMA$Diastolic_blood_pressure <- PIMA$Diastolic_blood_pressure/5
PIMA$Two_hour_serum_insulin <- PIMA$Two_hour_serum_insulin/5
PIMA$Diabetes_pedigree_function <- as.numeric(scale(PIMA$Diabetes_pedigree_function))

# CHIMERA
PIMA_CHIMERA$Age <- PIMA_CHIMERA$Age/10
PIMA_CHIMERA$Plasma_glucose_concentration <- PIMA_CHIMERA$Plasma_glucose_concentration/5
PIMA_CHIMERA$Diastolic_blood_pressure <- PIMA_CHIMERA$Diastolic_blood_pressure/5
PIMA_CHIMERA$Two_hour_serum_insulin <- PIMA_CHIMERA$Two_hour_serum_insulin/5
PIMA_CHIMERA$Diabetes_pedigree_function <- as.numeric(scale(PIMA_CHIMERA$Diabetes_pedigree_function))


# SYNTHPOP
PIMA_SYNTHPOP$Age <- PIMA_SYNTHPOP$Age/10
PIMA_SYNTHPOP$Plasma_glucose_concentration <- PIMA_SYNTHPOP$Plasma_glucose_concentration/5
PIMA_SYNTHPOP$Diastolic_blood_pressure <- PIMA_SYNTHPOP$Diastolic_blood_pressure/5
PIMA_SYNTHPOP$Two_hour_serum_insulin <- PIMA_SYNTHPOP$Two_hour_serum_insulin/5
PIMA_SYNTHPOP$Diabetes_pedigree_function <- as.numeric(scale(PIMA_SYNTHPOP$Diabetes_pedigree_function))

# CTGAN
PIMA_CTGAN$Age <- PIMA_CTGAN$Age/10
PIMA_CTGAN$Plasma_glucose_concentration <- PIMA_CTGAN$Plasma_glucose_concentration/5
PIMA_CTGAN$Diastolic_blood_pressure <- PIMA_CTGAN$Diastolic_blood_pressure/5
PIMA_CTGAN$Two_hour_serum_insulin <- PIMA_CTGAN$Two_hour_serum_insulin/5
PIMA_CTGAN$Diabetes_pedigree_function <- as.numeric(scale(PIMA_CTGAN$Diabetes_pedigree_function))

# Logistic Regression ---------------------------------------------------------------------
Variables_num <- PIMA |> dplyr::select(where(is.numeric)) |> names()

# --- Logistic regression on original data ---
# Create formula for logistic regression with Outcome as response and all numeric variables as predictors
formula_reg.log <- as.formula(paste("Diabetes_diagnosis ~", paste0(Variables_num, collapse = "+"), sep = ""))

# Fit logistic regression model on the original data
reg.log.original <- glm(formula = formula_reg.log, data = PIMA, family = "binomial")

# Summary of the model
summary_original <- summary(reg.log.original)

# Extract coefficients table from the summary
coef_table <- summary_original$coefficients

# Create a regression table with exponentiated coefficients (Odds Ratios)
tab_reg_original <- tbl_regression(reg.log.original, exponentiate = TRUE)

# Get predicted probabilities for the original data
predictions.original <- predict(reg.log.original, type = "response")

# Compute ROC curve for the original data predictions
pima.roc.original <- roc(PIMA$Diabetes_diagnosis, predictions.original)
save(pima.roc.original,file=here("results","pima.roc.original.Rdata"))

# Create a table summarizing Odds Ratios (OR), 95% confidence intervals (CI), and standard errors for original data
tbl_original <- data.frame(
  Covariates = Variables_num,
  OR_original = round(na.omit(tab_reg_original$table_body$estimate), 3),
  IC95_low_original = round(na.omit(tab_reg_original$table_body$conf.low), 3),
  IC95_high_original = round(na.omit(tab_reg_original$table_body$conf.high), 3),
  Std_Error_original = round(na.omit(tab_reg_original[["table_body"]][["std.error"]]), 3)  # Exclude intercept
)


# --- Logistic regression on synthetic data generated by chimera ---
# Same formula as before
formula_reg.log <- as.formula(paste("Diabetes_diagnosis ~", paste0(Variables_num, collapse = "+"), sep = ""))

# Fit logistic regression on MICE synthetic dataset
reg.log.chimera <- glm(formula = formula_reg.log, data = PIMA_CHIMERA, family = "binomial")

# Model summary
summary_chimera <- summary(reg.log.chimera)

# Extract coefficient table
coef_table <- summary_chimera$coefficients

# Regression table with ORs
tab_reg_chimera <- tbl_regression(reg.log.chimera, exponentiate = TRUE)

# Predicted probabilities on synthetic MICE data
predictions.chimera <- predict(reg.log.chimera, type = "response")

# ROC curve for synthetic MICE data
pima.roc.chimera <- roc(PIMA_CHIMERA$Diabetes_diagnosis, predictions.chimera)
save(pima.roc.chimera,file=here("results","pima.roc.chimera.Rdata"))

# Summary table for MICE synthetic data
tbl_chimera <- data.frame(
  Covariates = Variables_num,
  OR_chimera = round(na.omit(tab_reg_chimera$table_body$estimate), 3),
  IC95_low_chimera = round(na.omit(tab_reg_chimera$table_body$conf.low), 3),
  IC95_high_chimera = round(na.omit(tab_reg_chimera$table_body$conf.high),3),
  Std_Error_chimera = round(coef_table[-1, "Std. Error"], 3)
)



# --- Logistic regression on synthetic data generated by Synthpop ---
# Same formula
formula_reg.log <- as.formula(paste("Diabetes_diagnosis ~", paste0(Variables_num, collapse = "+"), sep = ""))

# Fit logistic regression on Synthpop synthetic dataset
reg.log.synthpop <- glm(formula = formula_reg.log, data = PIMA_SYNTHPOP, family = "binomial")

# Model summary
summary_synthpop <- summary(reg.log.synthpop)

# Extract coefficients
coef_table <- summary_synthpop$coefficients

# Regression table with ORs
tab_reg_synthpop <- tbl_regression(reg.log.synthpop, exponentiate = TRUE)

# Predictions on Synthpop synthetic data
predictions.synthpop <- predict(reg.log.synthpop, type = "response")

# ROC curve for Synthpop synthetic data
pima.roc.synthpop <- roc(PIMA_SYNTHPOP$Diabetes_diagnosis, predictions.synthpop)
save(pima.roc.synthpop,file=here("results","pima.roc.synthpop.Rdata"))

# Summary table for Synthpop synthetic data
tbl_synthpop <- data.frame(
  Covariates = Variables_num,
  OR_synthpop = round(na.omit(tab_reg_synthpop$table_body$estimate), 3),
  IC95_low_synthpop = round(na.omit(tab_reg_synthpop$table_body$conf.low), 3),
  IC95_high_synthpop = round(na.omit(tab_reg_synthpop$table_body$conf.high), 3),
  Std_Error_Synthpop = round(coef_table[-1, "Std. Error"], 3)
)


# --- Logistic regression on synthetic data generated by CTGAN ---
# Same formula
formula_reg.log <- as.formula(paste("Diabetes_diagnosis ~", paste0(Variables_num, collapse = "+"), sep = ""))

# Fit logistic regression on CTGAN synthetic dataset
reg.log.ctgan <- glm(formula = formula_reg.log, data = PIMA_CTGAN, family = "binomial")

# Model summary
summary_ctgan <- summary(reg.log.ctgan)

# Extract coefficients
coef_table <- summary_ctgan$coefficients

# Regression table with ORs
tab_reg_ctgan <- tbl_regression(reg.log.ctgan, exponentiate = TRUE)

# Predictions on CTGAN synthetic data
predictions.ctgan <- predict(reg.log.ctgan, type = "response")

# ROC curve for CTGAN synthetic data
pima.roc.ctgan <- roc(PIMA_CTGAN$Diabetes_diagnosis, predictions.ctgan)
save(pima.roc.ctgan,file=here("results","pima.roc.ctgan.Rdata"))

# Summary table for ctgan synthetic data
tbl_ctgan <- data.frame(
  Covariates = Variables_num,
  OR_ctgan = round(na.omit(tab_reg_ctgan$table_body$estimate), 3),
  IC95_low_ctgan = round(na.omit(tab_reg_ctgan$table_body$conf.low), 3),
  IC95_high_ctgan = round(na.omit(tab_reg_ctgan$table_body$conf.high), 3),
  Std_Error_ctgan = round(coef_table[-1, "Std. Error"], 3)
)


# --- Merge all summary tables by covariables ---
tbl_df <- merge(tbl_original, tbl_chimera, by = "Covariates")
tbl_df <- merge(tbl_df, tbl_synthpop, by = "Covariates")
tbl_df <- merge(tbl_df, tbl_ctgan, by = "Covariates")


# Calculate standardized mean differences (SMD) between log odds ratios of synthetic and original datasets
# SMD formula: |(log(OR_synthetic) - log(OR_original)) / sqrt(Std_Error_synthetic^2 + Std_Error_original^2)|
tbl_df$SMD_log1 <- abs(round(
  (log(tbl_df$OR_chimera) - log(tbl_df$OR_original)) / 
    sqrt(tbl_df$Std_Error_chimera^2 + tbl_df$Std_Error_original^2), 3))

tbl_df$SMD_log2 <- abs(round(
  (log(tbl_df$OR_synthpop) - log(tbl_df$OR_original)) / 
    sqrt(tbl_df$Std_Error_Synthpop^2 + tbl_df$Std_Error_original^2), 3))

tbl_df$SMD_log3 <- abs(round(
  (log(tbl_df$OR_ctgan) - log(tbl_df$OR_original)) / 
    sqrt(tbl_df$Std_Error_ctgan^2 + tbl_df$Std_Error_original^2), 3))

PIMA_tbl_OR <- tbl_df
save(PIMA_tbl_OR,file=here("results","PIMA_tbl_OR.Rdata"))


# =========================================================
# Load datasets
# =========================================================
AIDS <- read.csv(here("data", "AIDS.csv"))
AIDS <- AIDS|> mutate_if(is.character,as.factor)

load(here("results", "AIDS_Metrics_T1_T9_chimera.Rdata"))
load(here("results", "AIDS_Metrics_T1_T9_synthpop.Rdata"))
load(here("results", "AIDS_Metrics_T1_T9_ctgan.Rdata"))

# Select the dataset minimizing the global distance score
# (lower distance = better overall performance)

best_row_chimera  <- which.min(AIDS_Metrics_T1_T9_chimera$L_total)
best_row_synthpop <- which.min(AIDS_Metrics_T1_T9_synthpop$L_total)
best_row_ctgan    <- which.min(AIDS_Metrics_T1_T9_ctgan$L_total)

load(here("results", "aids_synthetic_list_data_chimera.Rdata"))
load(here("results", "aids_synthetic_list_data_synthpop.Rdata"))
load(here("results", "aids_synthetic_list_data_ctgan.Rdata"))


# Best synthetic data
AIDS_CHIMERA <- aids_synthetic_list_data_chimera[[best_row_chimera]]
AIDS_SYNTHPOP <- aids_synthetic_list_data_synthpop[[best_row_synthpop]]
AIDS_CTGAN <- aids_synthetic_list_data_ctgan[[best_row_ctgan]]

# Save synthetic datasets
save(AIDS,file = here("data", "AIDS.Rdata"))
save(AIDS_CHIMERA,file = here("data", "AIDS_CHIMERA.Rdata"))
save(AIDS_SYNTHPOP,file = here("data", "AIDS_SYNTHPOP.Rdata"))
save(AIDS_CTGAN,file = here("data", "AIDS_CTGAN.Rdata"))

# -----------------------------------------------------------------------------------------
# utility Evaluation ---------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
AIDS$group <- "REAL"
AIDS_CHIMERA$group <- "CHIMERA"
AIDS_SYNTHPOP$group <- "SYNTHPOP"
AIDS_CTGAN$group <- "CTGAN"

# Combine all datasets into one: Original, synthpop, Synthpop, and CTGAN
data <- rbind(AIDS, AIDS_CHIMERA, AIDS_SYNTHPOP, AIDS_CTGAN)

# Ensure 'cens' is numeric and convert 'days' into weeks
data$Censored <- as.numeric(data$Censored)

# Convert 'group' to a factor to allow group-wise survival analysis
data$group <- factor(data$group, levels = c("REAL", "CHIMERA", "SYNTHPOP", "CTGAN"))

# Fit Kaplan-Meier survival curves by group
km_fit <- survfit(Surv(times, Censored) ~ group, data = data)
km_aids <- survminer::surv_summary(km_fit, data = data)
save(km_aids,file = here("results", "km_aids.Rdata"))

# cox model --------------------------------------------------------------------------
# filter
df_original <- AIDS|> select(-group)
# CHIMRA
df_chimera <- AIDS_CHIMERA|> select(-group)
# Synthpop
df_synthpop <- AIDS_SYNTHPOP|> select(-group)
# CTGAN
df_ctgan <- AIDS_CTGAN|> select(-group)


# Define time points
time_points <- seq(0, 1100, by = 100)
set.seed(123)
idx <- sample(seq_len(nrow(df_original)), size = 0.7 * nrow(df_original))
train <- df_original[idx, ]
test  <- df_original[-idx, ]

# Compute auc curve
auc_original <- predictive_utility_surv_curve(train, test, time_points,n_boot = 200,conf_level = 0.95) |>
  mutate(source = "REAL")

set.seed(123)
idx <- sample(seq_len(nrow(df_chimera)), size = 0.7 * nrow(df_chimera))
train <- df_chimera[idx, ]
test  <- df_chimera[-idx, ]

# Compute auc curve
auc_chimera <- predictive_utility_surv_curve(train, test, time_points,n_boot = 200,conf_level = 0.95) |>
  mutate(source = "CHIMERA")

set.seed(123)
idx <- sample(seq_len(nrow(df_synthpop)), size = 0.7 * nrow(df_synthpop))
train <- df_synthpop[idx, ]
test  <- df_synthpop[-idx, ]

auc_synthpop <- predictive_utility_surv_curve(train, test, time_points,n_boot = 200,conf_level = 0.95) |>
  mutate(source = "SYNTHPOP")

set.seed(123)
idx <- sample(seq_len(nrow(df_ctgan)), size = 0.7 * nrow(df_ctgan))
train <- df_ctgan[idx, ]
test  <- df_ctgan[-idx, ]

auc_ctgan <- predictive_utility_surv_curve(train, test, time_points,n_boot = 200,conf_level = 0.95) |>
  mutate(source = "CTGAN")

# Regrouper les résultats
auc_all_aids <- bind_rows(auc_original, auc_chimera, auc_synthpop, auc_ctgan)
save(auc_all_aids,file=here("results","auc_all_aids.Rdata"))


# standardized --------------------------------------------------------------------------------
#Original
df_original$Age <- df_original$Age/10
df_original$Karnofsky_score <- df_original$Karnofsky_score/5
df_original$CD4_baseline <- df_original$CD4_baseline/100
df_original$CD8_baseline <- df_original$CD8_baseline/100

#Chimera
df_chimera$Age <- df_chimera$Age/10
df_chimera$Karnofsky_score <- df_chimera$Karnofsky_score/5
df_chimera$CD4_baseline <- df_chimera$CD4_baseline/100
df_chimera$CD8_baseline <- df_chimera$CD8_baseline/100

#Synthpop
df_synthpop$Age <- df_synthpop$Age/10
df_synthpop$Karnofsky_score <- df_synthpop$Karnofsky_score/5
df_synthpop$CD4_baseline <- df_synthpop$CD4_baseline/100
df_synthpop$CD8_baseline <- df_synthpop$CD8_baseline/100

#Ctgan
df_ctgan$Age <- df_ctgan$Age/10
df_ctgan$Karnofsky_score <- df_ctgan$Karnofsky_score/5
df_ctgan$CD4_baseline <- df_ctgan$CD4_baseline/100
df_ctgan$CD8_baseline <- df_ctgan$CD8_baseline/100


#Original
cox.model.original <- coxph(Surv(times, Censored) ~ ., data = df_original)
summary_original <- summary(cox.model.original)
coef_table <- summary_original$coefficients

# IC 95% et HR
conf_int <- summary_original$conf.int

# Créer un tableau résumé
tbl_original <- data.frame(
  Covariates = rownames(coef_table),
  HR_original = round(conf_int[, "exp(coef)"], 3),
  IC95_low = round(conf_int[, "lower .95"], 3),
  IC95_high = round(conf_int[, "upper .95"], 3),
  SE_original = round(coef_table[, "se(coef)"], 6)
)


# chimera
cox.model.chimera <- coxph(Surv(times, Censored) ~ ., data = df_chimera)
summary_chimera <- summary(cox.model.chimera)
coef_table <- summary_chimera$coefficients
# IC 95% et HR
conf_int <- summary_chimera$conf.int

# Créer un tableau résumé
tbl_chimera <- data.frame(
  Covariates = rownames(coef_table),
  HR_chimera = round(conf_int[, "exp(coef)"], 3),
  IC95_low_chimera = round(conf_int[, "lower .95"], 3),
  IC95_high_chimera = round(conf_int[, "upper .95"], 3),
  SE_chimera = round(coef_table[, "se(coef)"], 6)
)



# synthpop
cox.model.synthpop <- coxph(Surv(times,Censored) ~ ., data = df_synthpop)
summary_synthpop <- summary(cox.model.synthpop)
coef_table <- summary_synthpop$coefficients
# IC 95% et HR
conf_int <- summary_synthpop$conf.int

# Créer un tableau résumé
tbl_synthpop <- data.frame(
  Covariates = rownames(coef_table),
  HR_synthpop = round(conf_int[, "exp(coef)"], 3),
  IC95_low_synthpop = round(conf_int[, "lower .95"], 3),
  IC95_high_synthpop = round(conf_int[, "upper .95"], 3),
  SE_synthpop = round(coef_table[, "se(coef)"], 6)
)


#Avatar CTGAN
cox.model.ctgan <- coxph(Surv(times, Censored) ~ ., data = df_ctgan)
summary_ctgan <- summary(cox.model.ctgan)
coef_table <- summary_ctgan$coefficients
# IC 95% et HR
conf_int <- summary_ctgan$conf.int

# Créer un tableau résumé
tbl_ctgan <- data.frame(
  Covariates = rownames(coef_table),
  HR_ctgan = round(conf_int[, "exp(coef)"], 3),
  IC95_low_ctgan = round(conf_int[, "lower .95"], 3),
  IC95_high_ctgan = round(conf_int[, "upper .95"], 3),
  SE_ctgan = round(coef_table[, "se(coef)"], 6)
)

tbl_df <- merge(tbl_original,tbl_chimera,by="Covariates",all.x = TRUE)
tbl_df <- merge(tbl_df,tbl_synthpop,by="Covariates",all.x = TRUE)
tbl_df <- merge(tbl_df,tbl_ctgan,by="Covariates",all.x = TRUE)

tbl_df$SMD_log1 <- abs(round((log(tbl_df$HR_chimera)-log(tbl_df$HR_original))/(sqrt(tbl_df$SE_original**2+tbl_df$SE_chimera**2)),2))
tbl_df$SMD_log2 <- abs(round((log(tbl_df$HR_synthpop)-log(tbl_df$HR_original))/(sqrt(tbl_df$SE_original**2+tbl_df$SE_synthpop**2)),2))
tbl_df$SMD_log3 <- abs(round((log(tbl_df$HR_ctgan)-log(tbl_df$HR_original))/(sqrt(tbl_df$SE_original**2+tbl_df$SE_ctgan**2)),2))

tbl_AIDS_HR <- tbl_df

# Forest plot
save(tbl_AIDS_HR,file=here("results","tbl_AIDS_HR.Rdata"))


# =========================================================
# Load datasets
# =========================================================
REIN <- read.csv(here("data", "df_rein_without_na.csv"))
REIN <- REIN|> mutate_if(is.character,as.factor)

load(here("results", "REIN_Metrics_T1_T9_chimera.Rdata"))
load(here("results", "REIN_Metrics_T1_T9_synthpop.Rdata"))
load(here("results", "REIN_Metrics_T1_T9_ctgan.Rdata"))

# Select the dataset minimizing the global distance score
# (lower distance = better overall performance)

best_row_chimera  <- which.min(REIN_Metrics_T1_T9_chimera$L_total)
best_row_synthpop <- which.min(REIN_Metrics_T1_T9_synthpop$L_total)
best_row_ctgan    <- which.min(REIN_Metrics_T1_T9_ctgan$L_total)

load(here("results", "rein_synthetic_list_data_chimera.Rdata"))
load(here("results", "rein_synthetic_list_data_synthpop.Rdata"))
load(here("results", "rein_synthetic_list_data_ctgan.Rdata"))


# Best synthetic data
REIN_CHIMERA <- rein_synthetic_list_data_chimera[[best_row_chimera]]
REIN_SYNTHPOP <- rein_synthetic_list_data_synthpop[[best_row_synthpop]]
REIN_CTGAN <- rein_synthetic_list_data_ctgan[[best_row_ctgan]]

# Save synthetic datasets
save(REIN,file = here("data", "REIN.Rdata"))
save(REIN_CHIMERA,file = here("data", "REIN_CHIMERA.Rdata"))
save(REIN_SYNTHPOP,file = here("data", "REIN_SYNTHPOP.Rdata"))
save(REIN_CTGAN,file = here("data", "REIN_CTGAN.Rdata"))

# filter
df.original <- REIN|>
  select(c(Age,Sex,Serum_albumin_level,Body_mass_index,Diabetes,Heart_failure,Peripheral_artery_disease,Coronary_artery_disease,Myocardial_infarction,Stroke,Cardiac_arrhythmia,Chronic_respiratory_failure,Malignancy,Cirrhosis,Severe_behavioral_disorders,Walking_autonomy,times,Censored))



df.chimera <- REIN_CHIMERA|>
  select(c(Age,Sex,Serum_albumin_level,Body_mass_index,Diabetes,Heart_failure,Peripheral_artery_disease,Coronary_artery_disease,Myocardial_infarction,Stroke,Cardiac_arrhythmia,Chronic_respiratory_failure,Malignancy,Cirrhosis,Severe_behavioral_disorders,Walking_autonomy,times,Censored))


df.synthpop <- REIN_SYNTHPOP|>
  select(c(Age,Sex,Serum_albumin_level,Body_mass_index,Diabetes,Heart_failure,Peripheral_artery_disease,Coronary_artery_disease,Myocardial_infarction,Stroke,Cardiac_arrhythmia,Chronic_respiratory_failure,Malignancy,Cirrhosis,Severe_behavioral_disorders,Walking_autonomy,times,Censored))


df.ctgan <- REIN_CTGAN|>
  select(c(Age,Sex,Serum_albumin_level,Body_mass_index,Diabetes,Heart_failure,Peripheral_artery_disease,Coronary_artery_disease,Myocardial_infarction,Stroke,Cardiac_arrhythmia,Chronic_respiratory_failure,Malignancy,Cirrhosis,Severe_behavioral_disorders,Walking_autonomy,times,Censored))

df.original$group <- "REAL"
df.chimera$group <- "CHIMERA"
df.synthpop$group <- "SYNTHPOP"
df.ctgan$group <- "CTGAN"

# Combine all datasets into one: Original, synthpop, Synthpop, and CTGAN
data <- rbind(df.original,df.chimera,df.synthpop,df.ctgan)

# Convert 'group' to a factor to allow group-wise survival analysis
data$group <- factor(data$group, levels = c("REAL", "CHIMERA", "SYNTHPOP", "CTGAN"))

# Fit Kaplan-Meier survival curves by group
km_fit <- survfit(Surv(times, Censored) ~ group, data = data)
km_rein <- survminer::surv_summary(km_fit, data = data)
save(km_rein,file = here("results", "km_rein.Rdata"))

df.original <- df.original|> select(-group)
df.chimera <- df.chimera|> select(-group)
df.synthpop <- df.synthpop|> select(-group)
df.ctgan <- df.ctgan|> select(-group)

# -----------------------------------------------------------------------------------------
# utility Evaluation ---------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# cox model --------------------------------------------------------------------------
# Define time points
time_points <- seq(0, 3900, by = 100)
set.seed(123)
idx <- sample(seq_len(nrow(df.original)), size = 0.7 * nrow(df.original))
train <- df.original[idx, ]
test  <- df.original[-idx, ]


# Compute auc curve
auc_original <- predictive_utility_surv_curve(train, test, time_points,n_boot = 200,conf_level = 0.95) |>
  mutate(source = "REAL")

set.seed(123)
idx <- sample(seq_len(nrow(df.chimera)), size = 0.7 * nrow(df.chimera))
train <- df.chimera[idx, ]
test  <- df.chimera[-idx, ]

# Compute auc curve
auc_chimera <- predictive_utility_surv_curve(train, test, time_points,n_boot = 200,conf_level = 0.95) |>
  mutate(source = "CHIMERA")

set.seed(123)
idx <- sample(seq_len(nrow(df.synthpop)), size = 0.7 * nrow(df.synthpop))
train <- df.synthpop[idx, ]
test  <- df.synthpop[-idx, ]

auc_synthpop <- predictive_utility_surv_curve(train, test, time_points,n_boot = 200,conf_level = 0.95) |>
  mutate(source = "SYNTHPOP")

set.seed(123)
idx <- sample(seq_len(nrow(df.ctgan)), size = 0.7 * nrow(df.ctgan))
train <- df.ctgan[idx, ]
test  <- df.ctgan[-idx, ]

auc_ctgan <- predictive_utility_surv_curve(train, test, time_points,n_boot = 200,conf_level = 0.95) |>
  mutate(source = "CTGAN")

# Regrouper les résultats
auc_all_rein <- bind_rows(auc_original, auc_chimera, auc_synthpop, auc_ctgan)
save(auc_all_rein,file=here("results","auc_all_rein.Rdata"))


# standardized --------------------------------------------------------------------------------
#Original
df.original$Age <- df.original$Age/10

#Chimera
df.chimera$Age <- df.chimera$Age/10

#Synthpop
df.synthpop$Age <- df.synthpop$Age/10

#Ctgan
df.ctgan$Age <- df.ctgan$Age/10

#Original
cox.model.original <- coxph(Surv(times, Censored) ~ ., data = df.original)
summary_original <- summary(cox.model.original)
coef_table <- summary_original$coefficients

# IC 95% et HR
conf_int <- summary_original$conf.int

# Créer un tableau résumé
tbl_original <- data.frame(
  Covariates = rownames(coef_table),
  HR_original = round(conf_int[, "exp(coef)"], 3),
  IC95_low = round(conf_int[, "lower .95"], 3),
  IC95_high = round(conf_int[, "upper .95"], 3),
  SE_original = round(coef_table[, "se(coef)"], 6)
)


# chimera
cox.model.chimera <- coxph(Surv(times, Censored) ~ ., data = df.chimera)
summary_chimera <- summary(cox.model.chimera)
coef_table <- summary_chimera$coefficients
# IC 95% et HR
conf_int <- summary_chimera$conf.int

# Créer un tableau résumé
tbl_chimera <- data.frame(
  Covariates = rownames(coef_table),
  HR_chimera = round(conf_int[, "exp(coef)"], 3),
  IC95_low_chimera = round(conf_int[, "lower .95"], 3),
  IC95_high_chimera = round(conf_int[, "upper .95"], 3),
  SE_chimera = round(coef_table[, "se(coef)"], 6)
)



# synthpop
cox.model.synthpop <- coxph(Surv(times,Censored) ~ ., data = df.synthpop)
summary_synthpop <- summary(cox.model.synthpop)
coef_table <- summary_synthpop$coefficients
# IC 95% et HR
conf_int <- summary_synthpop$conf.int

# Créer un tableau résumé
tbl_synthpop <- data.frame(
  Covariates = rownames(coef_table),
  HR_synthpop = round(conf_int[, "exp(coef)"], 3),
  IC95_low_synthpop = round(conf_int[, "lower .95"], 3),
  IC95_high_synthpop = round(conf_int[, "upper .95"], 3),
  SE_synthpop = round(coef_table[, "se(coef)"], 6)
)


#Avatar CTGAN
cox.model.ctgan <- coxph(Surv(times, Censored) ~ ., data = df.ctgan)
summary_ctgan <- summary(cox.model.ctgan)
coef_table <- summary_ctgan$coefficients
# IC 95% et HR
conf_int <- summary_ctgan$conf.int

# Créer un tableau résumé
tbl_ctgan <- data.frame(
  Covariates = rownames(coef_table),
  HR_ctgan = round(conf_int[, "exp(coef)"], 3),
  IC95_low_ctgan = round(conf_int[, "lower .95"], 3),
  IC95_high_ctgan = round(conf_int[, "upper .95"], 3),
  SE_ctgan = round(coef_table[, "se(coef)"], 6)
)

tbl_df <- merge(tbl_original,tbl_chimera,by="Covariates",all.x = TRUE)
tbl_df <- merge(tbl_df,tbl_synthpop,by="Covariates",all.x = TRUE)
tbl_df <- merge(tbl_df,tbl_ctgan,by="Covariates",all.x = TRUE)

tbl_df$SMD_log1 <- abs(round((log(tbl_df$HR_chimera)-log(tbl_df$HR_original))/(sqrt(tbl_df$SE_original**2+tbl_df$SE_chimera**2)),2))
tbl_df$SMD_log2 <- abs(round((log(tbl_df$HR_synthpop)-log(tbl_df$HR_original))/(sqrt(tbl_df$SE_original**2+tbl_df$SE_synthpop**2)),2))
tbl_df$SMD_log3 <- abs(round((log(tbl_df$HR_ctgan)-log(tbl_df$HR_original))/(sqrt(tbl_df$SE_original**2+tbl_df$SE_ctgan**2)),2))

tbl_REIN_HR <- tbl_df
save(tbl_REIN_HR,file=here("results","tbl_REIN_HR.Rdata"))











