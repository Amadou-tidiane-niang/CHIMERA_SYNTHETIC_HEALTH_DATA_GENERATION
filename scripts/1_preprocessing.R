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
# DATA LOADING | REIN 
# ============================================================
data <- read_excel(here("data", "cohorte_incident_2012-2021.xlsx"))
data_comp <- read_excel(here("data", "complement_janv2025.xlsx"))

# merge
data <- merge(data,data_comp,by="id_ano",all.x="TRUE")
str(data)

# filter
df <- data|>
  dplyr::select(-c(EQD_REG_LIB,res_dep_cod,NEPH_COD,id_ano,liste_longue))

# Transform all character variables into factor
df <- df|>mutate_if(is.character,as.factor)

# Data Management
df <- df %>%
  mutate(sex = as.factor(df$sex),
         URGn = as.factor(df$URGn),
         KTTINIn = as.factor(df$KTTINIn),
         PBRn = as.factor(df$PBRn),
         IDEn = as.factor(df$IDEn),
         DIABn = as.factor(df$DIABn),
         IRCn = as.factor(df$IRCn),
         O2n = as.factor(df$O2n),
         SASn = as.factor(df$SASn),
         ICn = as.factor(df$ICn),
         IDMn = as.factor(df$IDMn),
         RYTHMn = as.factor(df$RYTHMn),
         ANEVn = as.factor(df$ANEVn), 
         AMIn = as.factor(df$AMIn),
         AVCAITn = as.factor(df$AVCAITn),
         KCn = as.factor(df$KCn),
         VHBn = as.factor(df$VHBn),
         VHCn= as.factor(df$VHCn),
         CIRHn = as.factor(df$CIRHn),
         HANDn = as.factor(df$HANDn),
         AMPn = as.factor(df$AMPn),
         PLEGn = as.factor(df$PLEGn),
         CECITEn = as.factor(df$CECITEn),
         COMPORTn = as.factor(df$COMPORTn),
         sero = as.factor(df$sero),
         coro = as.factor(df$coro),
         ICOROn= as.factor(df$ICOROn))




# the end date of the study
date_cens_adminis <- "2022-12-31"

# compute survival time and censored time
df <- prepare_survival_data_time(df,date_cens_adminis,"DDIRT",c("dgrf","dsvr","dpdv","ddc"),censored_compute = TRUE)


# remove all variables related to survival time and censoring
df <- df|>
  dplyr::select(-c(DDIRT,dgrf,dsvr,ddc,dpdv,DINSCMED))

# Variable with % NAs less than 50
variables_names <- names(df)[colMeans(is.na(df)) < 0.5]
df <- df[,variables_names]

# Rename variables
df <- df %>%
  rename(
    Age=age,
    Sex=sex,
    Serum_albumin_level=ALBINI,
    Body_mass_index=bmi,
    eGFR=CKEPI,
    Diabetes=DIABn,
    Heart_failure=STADICn,
    Peripheral_artery_disease=STDAMIn,
    Coronary_artery_disease=ICOROn,
    Myocardial_infarction=IDMn,
    Stroke=AVCAITn,
    Cardiac_arrhythmia=RYTHMn,
    Chronic_respiratory_failure=IRCn,
    Malignancy=KCn,
    Cirrhosis=CIRHn,
    Severe_behavioral_disorders=COMPORTn,
    Walking_autonomy=MARCHn,
    Censored=censored
  )

# Select and rename relevant variables from the original dataset
df <- df |>
  dplyr::select(c(Age, Sex, Serum_albumin_level,Body_mass_index,eGFR,Diabetes, Heart_failure,
                  Peripheral_artery_disease, Coronary_artery_disease, Myocardial_infarction, Stroke,
                  Cardiac_arrhythmia, Chronic_respiratory_failure, Malignancy, Cirrhosis,
                  Severe_behavioral_disorders, Walking_autonomy, times, Censored))

# convert factor
df <- df|>
  mutate(
    Diabetes = factor(Diabetes, levels = c(0, 1), labels = c("No", "Yes")),
    Sex = factor(Sex, levels = c("Femme", "Homme"), labels = c("Female", "Male")),
    Heart_failure=factor(Heart_failure,levels=c("Pas d'insuffisance cardiaque","Insuffisance cardiaque Stade I-II","Insuffisance cardiaque Stade III-IV"),labels=c("None","NYHA I–II","NYHA III–IV")),
    Peripheral_artery_disease=factor(Peripheral_artery_disease,levels=c("Pas d'artérite des membres inférieurs","Artérite des membres inférieurs stade I-II","Artérite des membres inférieurs stade III-IV"),labels=c("None", "stage I-II", "stage III-IV")),
    Coronary_artery_disease= factor(Coronary_artery_disease , levels = c(0, 1), labels = c("No", "Yes")),
    Myocardial_infarction=factor(Myocardial_infarction, levels = c(0, 1), labels = c("No", "Yes")),
    Stroke=factor(Stroke, levels = c(0, 1), labels = c("No", "Yes")),
    Cardiac_arrhythmia=factor(Cardiac_arrhythmia, levels = c(0, 1), labels = c("No", "Yes")),
    Chronic_respiratory_failure=factor(Chronic_respiratory_failure, levels = c(0, 1), labels = c("No", "Yes")),
    Malignancy=factor(Malignancy, levels = c(0, 1), labels = c("No", "Yes")),
    Cirrhosis=factor(Cirrhosis, levels = c(0, 1), labels = c("No", "Yes")),
    Severe_behavioral_disorders=factor(Severe_behavioral_disorders, levels = c(0, 1), labels = c("No", "Yes")),
    Walking_autonomy=factor(Walking_autonomy,levels=c("Marche autonome","Nécessité d'une tierce personne","Incapacité totale"),labels=c("Independent","Assisted","Unable"))
)
    
df <- df[which(df$Age>=75),]
row.names(df) <- NULL

df_rein_with_na <- df

# Save dataset
write.csv(df_rein_with_na, here("data", "df_rein_with_na.csv"), row.names = FALSE)

# --- Initial MICE setup
ini <- mice(df, maxit = 0)
meth <- ini$method
pred <- ini$predictorMatrix
meth[names(df)[sapply(df, is.numeric)]] <- "pmm"

# --- Impute missing data
df_rein_without_na <- complete(mice(df, m = 1, maxit = 5, seed = 123,method = meth, predictorMatrix = pred, printFlag = TRUE))

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
