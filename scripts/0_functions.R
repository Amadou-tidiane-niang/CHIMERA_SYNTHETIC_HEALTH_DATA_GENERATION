# ============================================================
# 0. Useful function 
# ============================================================
prepare_survival_data_time <- function(data_time,
                                       date_cens_adminis,
                                       ddir_col,
                                       event_cols,
                                       censored_compute = FALSE) {
  
  # Convert dates
  data_time[[ddir_col]] <- as.POSIXct(data_time[[ddir_col]], tz = "UTC")
  data_time[event_cols] <- lapply(data_time[event_cols], as.POSIXct, tz = "UTC")
  
  # Administrative censoring
  censor_date <- as.POSIXct(date_cens_adminis, tz = "UTC")
  
  data_time <- data_time %>%
    mutate(across(all_of(event_cols), ~pmin(.x, censor_date)))
  
  # Earliest event
  min_event_time <- do.call(
    pmin,
    c(data_time[event_cols], na.rm = TRUE)
  )
  
  data_time$times <- as.POSIXct(min_event_time, tz = "UTC")
  data_time$times[is.na(data_time$times)] <- censor_date
  
  # Convert to time (days)
  data_time$times <- as.numeric(
    difftime(data_time$times, data_time[[ddir_col]], units = "days")
  )
  
  # Optional censoring indicator
  if (censored_compute) {
    data_time$censored <- with(data_time, ifelse(
      !is.na(ddc) &
        (is.na(dgrf) | ddc < dgrf) &
        (is.na(dsvr) | ddc < dsvr) &
        (is.na(dpdv) | ddc < dpdv),
      1, 0
    ))
    
    data_time$censored <- as.integer(data_time$censored)
  }
  
  return(data_time)
}


# Freedman- Draconis
fd_bins <- function(x) {
  N <- length(x)
  h <- 2 * (quantile(x, 0.75) - quantile(x, 0.25)) / (N^(1/3))
  if (h == 0) h <- 1e-8
  breaks <-seq(min(x), max(x) + h, by = h)
  return(breaks)
}

# Jensen Shanon divergence
JSD <- function(p, q) {
  m <- 0.5 * (p + q)
  kl <- function(a, b) sum(ifelse(a == 0, 0, a * log(a / b)))
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

# ============================================================
# 1. CHIMERA: Synthetic Data Generation Function
# ============================================================
#missing_rate = 0.10;is_survival = TRUE;timevar = NULL; statusvar = NULL;nb_imput = 1;nb_max_it = 10;seed = 123;method_conti = "pmm";verbose = TRUE
CHIMERA_generate <-function(df,missing_rate = 0.10,is_survival = TRUE,timevar = NULL, statusvar = NULL,nb_imput = 1,nb_max_it = 10,seed = 123,method_conti = "pmm",verbose = TRUE){
  
  # ============================================================
  # DESCRIPTION
  # ============================================================
  # Generate synthetic data using:
  # 1. MCAR missingness simulation
  # 2. MICE imputation
  # 3. Mahalanobis matching
  # Optional: survival data handling with Nelson-Aalen estimator
  
  set.seed(seed)
  
  # --- Checks
  if (is_survival) {
    if (is.null(timevar) || is.null(statusvar)) {
      stop("timevar and statusvar must be provided for survival data")
    }
    if (!(timevar %in% names(df)) || !(statusvar %in% names(df))) {
      stop("timevar or statusvar not found in df")
    }
  }
  
  # Backup of Real missing values
  if (any(is.na(df))) {
    na_indices_initial <- which(is.na(df), arr.ind = TRUE)
    
    # --- Initial MICE setup
    ini <- mice(df, maxit = 0)
    meth <- ini$method
    pred <- ini$predictorMatrix
    meth[names(df)[sapply(df, is.numeric)]] <- method_conti
    
    # --- Impute missing data
    df <- complete(mice(df, m = nb_imput, maxit = nb_max_it, seed = seed,
                        method = meth, predictorMatrix = pred, printFlag = verbose))
  }else{na_indices_initial <- NULL}
  
  # --- If survival data, compute Nelson-Aalen estimate
  if (is_survival) {
    df <- df |> rename(Censored = all_of(statusvar), times = all_of(timevar))
    df$H0 <- nelsonaalen(df, timevar = "times", statusvar = "Censored")
    df_ref <- df[, !(names(df) %in% "times")]
  } else {
    df_ref <- df
  }
  
  # --- Create an empty synthetic dataset
  df_syn <- as.data.frame(matrix(NA, nrow = nrow(df_ref), ncol = ncol(df_ref)))
  names(df_syn) <- names(df_ref)
  
  # --- Introduce MCAR missingness
  data_with_na <- delete_MCAR(df_ref, p = missing_rate, cols_mis = seq_len(ncol(df_ref)))
  df_switched <- df_ref
  
  # === MAIN LOOP: iterate until synthetic dataset is complete
  nb_iter <- 0
  while (any(is.na(df_syn))) {
    na_indices <- which(is.na(data_with_na), arr.ind = TRUE)
    
    # --- Apply NA mask to switched dataset
    df_switched[na_indices] <- NA
    
    # Create a binary matrix indicating the positions of missing values in df.switched
    # Each element is 1 if the value is NA (missing), and 0 otherwise.
    mat_na_before <- +is.na(df_switched)
    
    # --- MICE setup and imputation
    ini <- mice(data_with_na, maxit = 0)
    meth <- ini$method
    pred <- ini$predictorMatrix
    meth[names(df_ref)[sapply(df_ref, is.numeric)]] <- method_conti
    
    df_imputed <- complete(mice(data_with_na, m = nb_imput, maxit = nb_max_it, seed = seed,
                                method = meth, predictorMatrix = pred, printFlag = verbose))
    
    # --- Mahalanobis matching
    df_ref$Matching <- "Real"
    df_imputed$Matching <- "Avatar"
    match_data <- rbind(df_ref, df_imputed)
    match_data$Matching <- as.factor(match_data$Matching)
    
    match <- matchit(Matching ~ ., data = match_data, method = "nearest", distance = "mahalanobis")
    
    # -----Extract the matched indices and convert the matrix of matched row indices to a numeric vector
    matched_indices <- as.numeric(match[["match.matrix"]])
    
    # Subset the matched dataset using the matched indices and also remove the "Matching" column which is not needed for further analysis
    df_imputed <- match_data[matched_indices, ][, !(names(match_data) %in% "Matching")]
    
    # Remove the "Matching" column from the Real copy as well, for consistency
    df_ref <- df_ref |> dplyr::select(-Matching)
    
    # Reset the row names of the imputed dataframe to be sequential from 1 to n
    # This ensures alignment with df_ref (useful for merging or comparing)
    rownames(df_imputed) <- seq_len(nrow(df_ref))
    
    # --- Track changed NA positions
    df_switched_after <- df_switched
    df_switched_after$ind_after <- (nrow(df_ref)+1):nrow(match_data)
    df_switched_after <- df_switched_after[match(matched_indices, df_switched_after$ind_after), ]
    df_switched_after <- df_switched_after[, !(names(df_switched_after) %in% "ind_after")]
    
    #---- Pattern of missing data after matching, using the Real dataset
    mat_na_after <- +is.na(df_switched_after)
    
    # Compare the missing data patterns before and after matching (TRUE if unchanged, FALSE if altered)
    mat_compare <- mat_na_before == mat_na_after
    
    # --- Fill in avatarized data
    for(col in names(df_syn)){
      df_syn[,col] <- ifelse(is.na(df_switched[,col]),df_imputed[,col], df_syn[,col])
    }
    
    
    # Identify the positions to be transferred to NA: where NA has changed
    changed_na <- which(!mat_compare, arr.ind = TRUE)
    
    # Convert the two index matrices into single strings to compare them
    initial_na_keys <- paste(na_indices[,1], na_indices[,2], sep = "-")
    changed_keys <- paste(changed_na[,1], changed_na[,2], sep = "-")
    
    # Find the positions in changed_na that correspond to initial NAs
    to_na_keys <- intersect(initial_na_keys, changed_keys)
    
    # Transform into row/column indices
    rows <- as.integer(sub("-.*", "", to_na_keys))
    cols <- as.integer(sub(".*-", "", to_na_keys))
    
    # Apply the NA update in df_syn
    for (k in seq_along(rows)) {
      df_syn[rows[k], cols[k]] <- NA
    }
    
    # --- Prepare for next iteration
    df_switched <- df_ref
    
    # Restore factor levels
    for (col in names(df_ref)) {
      if (is.factor(df_ref[[col]])) {
        df_switched[[col]] <- factor(df_switched[[col]], levels = levels(df_ref[[col]]))
      }
    }
    
    # --- Reapply NA for next round if necessary
    p_na <- min(colMeans(is.na(df_syn)))
    na_indices <- which(is.na(df_syn), arr.ind = TRUE)
    
    if(missing_rate<p_na){
      # Compute the number of elements to be replaced by NA
      nb_na <- ceiling(missing_rate * nrow(df_ref))
      na_indices <- na_indices[sample(nrow(na_indices), nb_na * ncol(df_ref)), , drop = FALSE]
    }else{na_indices <- na_indices}
    
    
    df_switched[na_indices] <- NA
    data_with_na <- df_switched
    
    # Transform categorical variables into factors              
    for (i in names(df_syn)) {
      if (is.factor(df_ref[, i])) {
        df_syn[, i] <- as.factor(df_syn[, i])
        levels(df_syn[, i]) <- levels(df_ref[, i])
      }
    }
    nb_iter <- nb_iter+1
    # Loop exit condition
    if (all(!is.na(df_syn)) | p_na<missing_rate) {
      break
    }
  }
  
  # ----Last step
  
  # --- MICE setup and imputation
  ini <- mice(df_switched, maxit = 0)
  meth <- ini$method
  pred <- ini$predictorMatrix
  meth[names(df_ref)[sapply(df_ref, is.numeric)]] <- method_conti
  
  # impute dataset
  df_imputed <- complete(mice(df_switched, m = nb_imput, maxit = nb_max_it, seed = seed,
                              method = meth, predictorMatrix = pred, printFlag = verbose))
  
  # --- Fill in avatarized data
  for(col in names(df_syn)){
    df_syn[,col] <- ifelse(is.na(df_switched[,col]),df_imputed[,col], df_syn[,col])
  }
  
  # Transform categorical variables into factors              
  for (i in names(df_syn)) {
    if (is.factor(df_ref[, i])) {
      df_syn[, i] <- as.factor(df_syn[, i])
      levels(df_syn[, i]) <- levels(df_ref[, i])
    }
  }
  
  # --- Mahalanobis matching
  df_ref$Matching <- "Real"
  df_syn$Matching <- "Avatar"
  match_data <- rbind(df_ref, df_syn)
  match_data$Matching <- as.factor(match_data$Matching)
  match <- matchit(Matching ~ ., data = match_data, method = "nearest", distance = "mahalanobis")
  matched_indices <- as.numeric(match[["match.matrix"]])
  
  df_syn <- match_data[matched_indices, ][, !(names(match_data) %in% "Matching")]
  row.names(df_syn) <- 1:nrow(df)
  
  # --- Predict back survival times using Nelson-Aalen + spline
  if (is_survival) {
    spline_fit <- smooth.spline(df$H0, df$times)
    predicted_times <- predict(spline_fit, df_syn$H0)$y
    df_syn$times <- round(predicted_times)
    
    df_syn <- df_syn[, names(df)]
    df_syn <- df_syn[, !(names(df_syn) %in% "H0")]
    df <- df[, !(names(df) %in% "H0")]
  }
  
  # Restore Real names
  if (is_survival) {
    # Administrative censorship management
    max_time_event <- max(df[which(df$Censored==1),"times"])
    df_syn[which(df_syn$times> max_time_event),"Censored"] <- 0
    
    df_syn <- df_syn |> rename(
      !!statusvar := Censored,
      !!timevar := times
    )
    
    df <- df |> rename(
      !!statusvar := Censored,
      !!timevar := times
    )
    
  }
  
  # --- Reapply Real missingness mask
  if (exists("na_indices_initial")) {
    return(list(df=df,df_syn=df_syn,na_indices_initial=na_indices_initial,nb_iter=nb_iter))
  }else{return(list(df=df,df_syn=df_syn,nb_iter=nb_iter))}
  
  
}



# ============================================================
# Privacy Metrics : AIR
# ============================================================
# 8. Attribute Inference Risk (AIR)

privacy_AIR <- function(
    real_data,
    syn_data,
    target_col,
    feature_associate
) {
  
  #' Compute Attribute Inference Risk (AIR)
  #'
  #' @param real_data Original dataset
  #' @param syn_data  Synthetic dataset
  #' @param target_col Target variable to infer
  #' @param feature_associate Variables excluded from predictors
  #' @return F1-score
  
  colnames(syn_data) <- colnames(real_data)
  
  # ------------------------------------------------------------
  # Separate predictors and targets
  # ------------------------------------------------------------
  real_X <- real_data |> dplyr::select(-all_of(target_col))
  syn_X  <- syn_data  |> dplyr::select(-all_of(target_col))
  
  y_true <- real_data |> dplyr::select(all_of(target_col))
  syn_y  <- syn_data  |> dplyr::select(all_of(target_col))
  
  stopifnot(identical(colnames(real_X), colnames(syn_X)))
  
  # ------------------------------------------------------------
  # Harmonize categorical variables
  # ------------------------------------------------------------
  cat_cols <- names(real_X)[
    sapply(real_X, is.character) | sapply(real_X, is.factor)
  ]
  
  for (col in cat_cols) {
    levels_union <- union(unique(real_X[[col]]), unique(syn_X[[col]]))
    
    real_X[[col]] <- factor(real_X[[col]], levels = levels_union)
    syn_X[[col]]  <- factor(syn_X[[col]],  levels = levels_union)
  }
  
  # Remove associated features
  real_X <- real_X |> dplyr::select(-all_of(feature_associate))
  syn_X  <- syn_X  |> dplyr::select(-all_of(feature_associate))
  
  # ------------------------------------------------------------
  # Nearest neighbor matching
  # ------------------------------------------------------------
  nearest <- gower_topn(x = real_X, y = syn_X, n = 1)
  nearest_idx <- nearest$index
  
  # Predicted targets
  y_pred <- syn_y[nearest_idx, , drop = FALSE]
  rownames(y_pred) <- NULL
  
  # ------------------------------------------------------------
  # Compute F1 score
  # ------------------------------------------------------------
  y_true <- as.factor(y_true[[1]])
  y_pred <- factor(y_pred[[1]], levels = levels(y_true))
  
  f1 <- MLmetrics::F1_Score(y_pred = y_pred, y_true = y_true)
  
  return(f1)
}



# ==============================================================
# Assessment_function_unified
# ==============================================================
# PURPOSE:
# Comprehensive evaluation of synthetic data quality
# covering fidelity, utility, and privacy risks.
# ==============================================================

Assessment_function_unified <- function(
    X_real,
    X_syn,
    S = 5,
    seed = 123,
    formula,
    modeling = c("logistic", "cox"),
    timevar = NULL,
    statusvar = NULL,
    missing_rate = 0.25,
    target_col,
    feature_associate,
    is_survival=FALSE,
    generator=NULL,
    trn   = trn,
    val   = hol,
    syn_ctgan_trn = NULL
) {
  
  # ============================================================
  # 0. PREPROCESSING
  # ============================================================
  
  X_real <- X_real |> dplyr::mutate(across(where(is.character), as.factor))
  X_syn  <- X_syn  |> dplyr::mutate(across(where(is.character), as.factor))
  
  # ============================================================
  # 1. VARIABLE TYPES
  # ============================================================
  
  continuous_vars  <- names(X_real)[sapply(X_real, is.numeric)]
  categorical_vars <- names(X_real)[sapply(X_real, is.factor)]
  
  # ============================================================
  # 2. PREPROCESSING (BINNING + ENCODING)
  # ============================================================
  
  X_real_proc <- X_real
  X_syn_proc  <- X_syn
  
  for (v in continuous_vars) {
    brks <- fd_bins(X_real[[v]])
    X_real_proc[[v]] <- cut(X_real[[v]], brks, labels = FALSE, include.lowest = TRUE)
    X_syn_proc[[v]]  <- cut(X_syn[[v]],  brks, labels = FALSE, include.lowest = TRUE)
  }
  
  for (v in categorical_vars) {
    lvls <- union(levels(X_real[[v]]), levels(X_syn[[v]]))
    X_real_proc[[v]] <- as.integer(factor(X_real[[v]], levels = lvls))
    X_syn_proc[[v]]  <- as.integer(factor(X_syn[[v]],  levels = lvls))
  }
  
  # ============================================================
  # T1 — UNIVARIATE FIDELITY (JSD)
  # ============================================================
  
  T1 <- mean(sapply(names(X_real_proc), function(v) {
    p <- table(X_real_proc[[v]]) / nrow(X_real_proc)
    q <- table(X_syn_proc[[v]])  / nrow(X_syn_proc)
    lv <- union(names(p), names(q))
    p <- p[lv]; p[is.na(p)] <- 0
    q <- q[lv]; q[is.na(q)] <- 0
    JSD(p, q)
  }))
  
  # ============================================================
  # T2 — BIVARIATE FIDELITY
  # ============================================================
  
  T2_vals <- c()
  
  if (length(continuous_vars) > 1) {
    T2_vals <- c(T2_vals,
                 sqrt(sum((cor(X_real[, continuous_vars]) - cor(X_syn[, continuous_vars]))^2))
    )
  }
  
  if (length(categorical_vars) > 1) {
    cramer <- function(df) {
      n <- ncol(df)
      M <- matrix(1, n, n)
      for (i in 1:n) for (j in i:n) {
        if (i != j) {
          v <- vcd::assocstats(table(df[[i]], df[[j]]))$cramer
          M[i, j] <- M[j, i] <- v
        }
      }
      M
    }
    
    T2_vals <- c(T2_vals,
                 sqrt(sum((cramer(X_real[, categorical_vars]) - cramer(X_syn[, categorical_vars]))^2))
    )
  }
  
  T2 <- mean(T2_vals)
  
  # ============================================================
  # T3 — DISCRIMINATION
  # ============================================================
  
  X <- rbind(X_real, X_syn)
  
  if (modeling == "cox") {
    X <- X |> dplyr::select(-all_of(c(timevar, statusvar)))
  }
  
  y <- factor(c(rep(0, nrow(X_real)), rep(1, nrow(X_syn))))
  folds <- caret::createFolds(y, k = S)
  
  T3 <- mean(sapply(folds, function(idx) {
    train <- X[-idx, ]
    test  <- X[idx, ]
    y_tr  <- y[-idx]
    y_te  <- y[idx]
    
    fit <- glm(y_tr ~ ., data = data.frame(y_tr, train), family = binomial)
    prob <- predict(fit, test, type = "response")
    pROC::auc(y_te, prob)
  }))
  
  # ============================================================
  # T4 — STANDARDIZED DIFFERENCE
  # ============================================================
  
  if (modeling == "logistic") {
    fit_real <- glm(formula, data = X_real, family = binomial)
    fit_syn  <- glm(formula, data = X_syn,  family = binomial)
    
    coef_real <- coef(summary(fit_real))[-1, ]
    coef_syn  <- coef(summary(fit_syn))[-1, ]
    
    effect_real <- coef_real[, "Estimate"]
    se_real     <- coef_real[, "Std. Error"]
    effect_syn  <- coef_syn[, "Estimate"]
    se_syn      <- coef_syn[, "Std. Error"]
    
  } else {
    fit_real <- survival::coxph(formula, data = X_real)
    fit_syn  <- survival::coxph(formula, data = X_syn)
    
    sr <- summary(fit_real)
    ss <- summary(fit_syn)
    
    effect_real <- sr$coefficients[, "coef"]
    se_real     <- sr$coefficients[, "se(coef)"]
    effect_syn  <- ss$coefficients[, "coef"]
    se_syn      <- ss$coefficients[, "se(coef)"]
  }
  
  SDiff <- abs((effect_syn - effect_real) / sqrt(se_syn^2 + se_real^2))
  SD_q05 <- as.numeric(quantile(SDiff, probs = 0.05))
  SD_q95 <- as.numeric(quantile(SDiff, probs = 0.95))
  T4 <- as.numeric(mean(SDiff))
  
  # ============================================================
  # T5 / T6 — PREDICTIVE UTILITY
  # ============================================================
  
  target <- target_col
  folds <- caret::createFolds(X_real[[target]], k = S)
  
  auc_trtr <- auc_tstr <- auc_tsrtr <- numeric(S)
  
  for (s in seq_len(S)) {
    test  <- X_real[folds[[s]], ]
    train <- X_real[-folds[[s]], ]
    
    if (modeling == "logistic") {
      fit_trtr <- glm(formula, train, family = binomial)
      fit_tstr <- glm(formula, X_syn, family = binomial)
      fit_mix  <- glm(formula, rbind(train, X_syn), family = binomial)
      
      auc_trtr[s] <- pROC::auc(test[[target]], predict(fit_trtr, test, "response"))
      auc_tstr[s] <- pROC::auc(test[[target]], predict(fit_tstr, test, "response"))
      auc_tsrtr[s]<- pROC::auc(test[[target]], predict(fit_mix,  test, "response"))
      
    } else {
      Surv_test <- survival::Surv(test[[timevar]], test[[statusvar]])
      times <- quantile(train[[timevar]],probs = seq(0.1, 0.9, by=0.1), na.rm=TRUE)
      cox_auc <- function(fit, train_data) {
        lp_tr <- predict(fit, train_data)
        lp_te <- predict(fit, test)
        
        mean(survAUC::AUC.cd(
          survival::Surv(train_data[[timevar]], train_data[[statusvar]]),
          Surv_test, lp_tr, lp_te,times
        )$auc, na.rm = TRUE)
      }
      
      auc_trtr[s]  <- cox_auc(coxph(formula, train), train)
      auc_tstr[s]  <- cox_auc(coxph(formula, X_syn), X_syn)
      auc_tsrtr[s] <- cox_auc(coxph(formula, rbind(train, X_syn)), rbind(train, X_syn))
    }
  }
  
  T5 <- mean(auc_tstr)  - mean(auc_trtr)
  T6 <- mean(auc_tsrtr) - mean(auc_trtr)
  
  # ============================================================
  # T7 — RECORD MATCHING
  # ============================================================
  
  T7 <- mean(gower::gower_topn(x = X_real, y = X_syn, n = 1)$distance)
  
  # ============================================================
  # T8 — HOLDOUT PRIVACY
  # ============================================================
  row.names(trn) <-NULL
  
  # ------------------------------------------------------------
  # GENERATE SYNTHETIC DATA (depending on generator)
  # ------------------------------------------------------------
  if (generator == "chimera") {
    
    synthetic_data <- CHIMERA_generate(
      df = trn,
      missing_rate = missing_rate,
      is_survival = is_survival,
      timevar = timevar,
      statusvar = statusvar,
      nb_imput = 1,
      nb_max_it = 5,
      seed = seed,
      method_conti = "pmm",
      verbose = TRUE
    )
    
    syn <- synthetic_data[["df_syn"]]
    
  } else if (generator == "synthpop") {
    
    syn_obj <- synthpop::syn(trn)
    syn <- syn_obj[["syn"]]
    
  }else if (generator == "ctgan"){
    syn <- as.data.frame(syn_ctgan_trn)
  }
  
  all <- bind_rows(trn, val, syn)
  all_hot <- dummyVars(" ~ .", data = all, fullRank = TRUE) |> predict(newdata = all) |> as.matrix()
  
  n_trn <- nrow(trn)
  n_val <- nrow(val)
  n_syn <- nrow(syn)
  
  trn_hot <- all_hot[1:n_trn, ]
  val_hot <- all_hot[(n_trn+1):(n_trn+n_val), ]
  syn_hot <- all_hot[(n_trn+n_val+1):(n_trn+n_val+n_syn), ]
  
  d_trn <- get.knnx(trn_hot, syn_hot, k = 1)$nn.dist[,1]
  d_val <- get.knnx(val_hot, syn_hot, k = 1)$nn.dist[,1]
  
  T8 <- mean(d_trn < d_val) + (n_trn/(n_trn+n_val))*mean(d_trn == d_val)
  
  # ============================================================
  # T9 — ATTRIBUTE INFERENCE
  # ============================================================
  
  AIR <- unlist(lapply(feature_associate, function(vars) {
    privacy_AIR(X_real, X_syn, target_col, vars)
  }))
  
  T9 <- mean(AIR)
  AIR_min <- min(AIR)
  AIR_max <- max(AIR)
  
  # ============================================================
  # OUTPUT
  # ============================================================
  
  return(c(
    T1=T1, T2=T2, T3=T3, T4=T4,
    T5=T5, T6=T6, T7=T7, T8=T8, T9=T9,
    T10=SD_q05, T11=SD_q95,
    T12=AIR_min, T13=AIR_max
  ))
}




Ablation_Analysis_generate_synthetic_data <- function(
    data,
    missing_rate = 0.10,
    iterative_masking = FALSE,
    survival_data = FALSE,
    time_variable = NULL,
    event_variable = NULL,
    n_imputations = 1,
    max_iterations = 5,
    random_seed = 123,
    continuous_method = "pmm",
    verbose = TRUE
){
  
  #=============================================================================
  # Abaltion Analysis
  #=============================================================================
  # Description:
  #   Generates synthetic datasets using MICE-based data synthesis.
  #   Two synthesis strategies are available:
  #
  #   1. Standard synthesis:
  #      Direct full-data synthesis using MICE.
  #
  #   2. Iterative masking synthesis:
  #      Progressive MCAR masking + iterative MICE reconstruction.
  #
  #   The function also supports survival data through:
  #      - Nelson-Aalen cumulative hazard estimation
  #      - spline-based reconstruction of survival times
  
  # ===========================================================================
  # 1. Initialization
  # ===========================================================================
  
  set.seed(random_seed)
  
  # ===========================================================================
  # 2. Validate survival arguments
  # ===========================================================================
  
  if (survival_data) {
    
    if (is.null(time_variable) || is.null(event_variable)) {
      
      stop(
        "time variable and event variable must be provided ",
        "when survival_data = TRUE."
      )
    }
    
    if (!(time_variable %in% names(data)) ||
        !(event_variable %in% names(data))) {
      
      stop(
        "Specified survival columns do not exist in the dataset."
      )
    }
  }
  
  # ===========================================================================
  # 3. Handle pre-existing missing values
  # ===========================================================================
  
  if (any(is.na(data))) {
    
    original_missing_positions <- which(is.na(data),arr.ind = TRUE)
    
    # Initialize MICE configuration
    mice_init <- mice(data, maxit = 0)
    
    imputation_methods <- mice_init$method
    
    predictor_matrix <- mice_init$predictorMatrix
    
    # Define continuous variable imputation method
    imputation_methods[names(data)[sapply(data, is.numeric)]] <- continuous_method
    
    # Initial imputation of the original dataset
    data <- complete(
      
      mice(
        data,
        m = n_imputations,
        maxit = max_iterations,
        seed = random_seed,
        method = imputation_methods,
        predictorMatrix = predictor_matrix,
        printFlag = verbose
      )
    )
    
  } else {
    
    original_missing_positions <- NULL
  }
  
  # ===========================================================================
  # 4. Survival preprocessing
  # ===========================================================================
  
  if (survival_data) {
    
    # Rename variables internally
    data <- data |>
      rename(event_status = all_of(event_variable),
             survival_time = all_of(time_variable)
             )
    
    # Compute Nelson-Aalen cumulative hazard estimator
    data$cumulative_hazard <- nelsonaalen(data,timevar = "survival_time",statusvar = "event_status")
    
    # Remove survival time before synthesis
    working_data <- data[,!(names(data) %in% "survival_time")]
    
  } else {
    working_data <- data
  }
  
  # ===========================================================================
  # 5. STANDARD SYNTHESIS MODE
  # ===========================================================================
  
  if (!iterative_masking) {
    
    # Create synthesis mask
    synthesis_mask <- make.where(data, "all")
    
    # Define synthesis methods
    synthesis_methods <- make.method(data,where = synthesis_mask)
    
    # Generate synthetic dataset
    synthetic_data <- complete(
      mice(
        data,
        m = n_imputations,
        maxit = max_iterations,
        seed = random_seed,
        method = synthesis_methods,
        where = synthesis_mask,
        printFlag = verbose
      )
    )
    
    iteration_count <- 1
    
  } else {
    
    # =========================================================================
    # 6. ITERATIVE MASKING SYNTHESIS MODE
    # =========================================================================
    
    # Initialize empty synthetic dataset
    synthetic_data <- as.data.frame(matrix(NA,nrow = nrow(working_data),ncol = ncol(working_data)))
    names(synthetic_data) <- names(working_data)
    
    # Introduce MCAR missing values
    masked_data <- delete_MCAR(working_data,p = missing_rate,cols_mis = seq_len(ncol(working_data)))
    temporary_data <- working_data
    iteration_count <- 0
    
    # -------------------------------------------------------------------------
    # Main iterative synthesis loop
    # -------------------------------------------------------------------------
    
    while (any(is.na(synthetic_data))) {
      
      current_missing_positions <- which(is.na(masked_data),arr.ind = TRUE)
      
      # Apply missing mask
      temporary_data[current_missing_positions] <- NA
      
      # Initialize MICE
      mice_init <- mice(masked_data, maxit = 0)
      
      imputation_methods <- mice_init$method
      
      predictor_matrix <- mice_init$predictorMatrix
      
      # Continuous variable imputation method
      imputation_methods[names(working_data)[sapply(working_data, is.numeric)]] <- continuous_method
      
      # Impute masked dataset
      imputed_data <- complete(
        
        mice(
          masked_data,
          m = n_imputations,
          maxit = max_iterations,
          seed = random_seed,
          method = imputation_methods,
          predictorMatrix = predictor_matrix,
          printFlag = verbose
        )
      )
      
      # -----------------------------------------------------------------------
      # Fill synthetic dataset progressively
      # -----------------------------------------------------------------------
      
      for (variable_name in names(synthetic_data)) {
        synthetic_data[, variable_name] <- ifelse(is.na(temporary_data[, variable_name]),
                                                  imputed_data[, variable_name],
                                                  synthetic_data[, variable_name])
        }
      
      # -----------------------------------------------------------------------
      # Restore original structure
      # -----------------------------------------------------------------------
      
      temporary_data <- working_data
      
      # Restore factor levels
      for (variable_name in names(working_data)) {
        
        if (is.factor(working_data[[variable_name]])) {
          temporary_data[[variable_name]] <- factor(
          temporary_data[[variable_name]],
          levels = levels(working_data[[variable_name]])
          )
        }
      }
      
      # -----------------------------------------------------------------------
      # Prepare next masking iteration
      # -----------------------------------------------------------------------
      
      minimum_missing_fraction <- min(colMeans(is.na(synthetic_data)))
      remaining_missing_positions <- which(is.na(synthetic_data),arr.ind = TRUE)
      
      # Control masking intensity
      if (missing_rate < minimum_missing_fraction) {
        number_missing <- ceiling(missing_rate * nrow(working_data))
        remaining_missing_positions <- remaining_missing_positions[
          sample(nrow(remaining_missing_positions),number_missing * ncol(working_data)),
          ,
          drop = FALSE
        ]
      }
      
      temporary_data[remaining_missing_positions] <- NA
      masked_data <- temporary_data
      
      # -----------------------------------------------------------------------
      # Restore factor levels in synthetic dataset
      # -----------------------------------------------------------------------
      
      for (variable_name in names(synthetic_data)) {
        if (is.factor(working_data[, variable_name])) {
          synthetic_data[, variable_name] <- as.factor(synthetic_data[, variable_name])
          levels(synthetic_data[, variable_name]) <- levels(working_data[, variable_name])
        }
      }
      
      iteration_count <- iteration_count + 1
      
      # -----------------------------------------------------------------------
      # Stop condition
      # -----------------------------------------------------------------------
      
      if (
        all(!is.na(synthetic_data)) || minimum_missing_fraction < missing_rate
      ) {
        break
      }
    }
    
    # =========================================================================
    # 7. Final imputation step
    # =========================================================================
    
    mice_init <- mice(masked_data, maxit = 0)
    imputation_methods <- mice_init$method
    predictor_matrix <- mice_init$predictorMatrix
    imputation_methods[names(working_data)[sapply(working_data, is.numeric)]] <- continuous_method
    final_imputed_data <- complete(
      
      mice(
        masked_data,
        m = n_imputations,
        maxit = max_iterations,
        seed = random_seed,
        method = imputation_methods,
        predictorMatrix = predictor_matrix,
        printFlag = verbose
      )
    )
    
    # Fill remaining missing values
    for (variable_name in names(synthetic_data)) {
      synthetic_data[, variable_name] <- ifelse(is.na(masked_data[, variable_name]),final_imputed_data[, variable_name],synthetic_data[, variable_name])
    }
    
    row.names(synthetic_data) <- NULL
  }
  
  # Transform categorical variables into factors              
  for (i in names(synthetic_data)) {
    if (is.factor(data[, i])) {
      synthetic_data[, i] <- as.factor(synthetic_data[, i])
      levels(synthetic_data[, i]) <- levels(data[, i])
    }
  }
  
  # ===========================================================================
  # 8. Survival time reconstruction
  # ===========================================================================
  
  if (survival_data) {
    
    # Fit spline between cumulative hazard and survival time
    spline_model <- smooth.spline(
      data$cumulative_hazard,
      data$survival_time
    )
    
    # Predict synthetic survival times
    predicted_times <- predict(
      spline_model,
      synthetic_data$cumulative_hazard
    )$y
    
    synthetic_data$survival_time <- round(predicted_times)
    
    # Reorder columns
    synthetic_data <- synthetic_data[, names(data)]
    
    # Remove cumulative hazard variable
    synthetic_data <- synthetic_data[,!(names(synthetic_data) %in% "cumulative_hazard")]
    
    data <- data[,!(names(data) %in% "cumulative_hazard")]
  }
  
  # ===========================================================================
  # 9. Restore original survival variable names
  # ===========================================================================
  
  if (survival_data) {
    
    # Administrative censoring correction
    maximum_event_time <- max(
      data[data$event_status == 1, "survival_time"]
    )
    
    synthetic_data[
      which(synthetic_data$survival_time > maximum_event_time),
      "event_status"
    ] <- 0
    
    # Restore original variable names
    synthetic_data <- synthetic_data |>
      rename(
        !!event_variable := event_status,
        !!time_variable  := survival_time
      )
    
    data <- data |>
      rename(
        !!event_variable := event_status,
        !!time_variable  := survival_time
      )
  }
  
  # ===========================================================================
  # 10. Output
  # ===========================================================================
  
  return(list(original_data = data,synthetic_data = synthetic_data,original_missing_positions = original_missing_positions))
}


Assessment_function_unified_Ablation_Analysis <- function(
    X_real,
    X_syn,
    S = 5,
    seed = 123,
    formula,
    modeling = c("logistic", "cox"),
    timevar = NULL,
    statusvar = NULL,
    missing_rate = 0.25,
    iterative_masking=NULL,
    target_col,
    feature_associate,
    is_survival=FALSE,
    trn   = trn,
    val   = hol
) {
  
  # ============================================================
  # 0. PREPROCESSING
  # ============================================================
  
  X_real <- X_real |> dplyr::mutate(across(where(is.character), as.factor))
  X_syn  <- X_syn  |> dplyr::mutate(across(where(is.character), as.factor))
  
  # ============================================================
  # 1. VARIABLE TYPES
  # ============================================================
  
  continuous_vars  <- names(X_real)[sapply(X_real, is.numeric)]
  categorical_vars <- names(X_real)[sapply(X_real, is.factor)]
  
  # ============================================================
  # 2. PREPROCESSING (BINNING + ENCODING)
  # ============================================================
  
  X_real_proc <- X_real
  X_syn_proc  <- X_syn
  
  for (v in continuous_vars) {
    brks <- fd_bins(X_real[[v]])
    X_real_proc[[v]] <- cut(X_real[[v]], brks, labels = FALSE, include.lowest = TRUE)
    X_syn_proc[[v]]  <- cut(X_syn[[v]],  brks, labels = FALSE, include.lowest = TRUE)
  }
  
  for (v in categorical_vars) {
    lvls <- union(levels(X_real[[v]]), levels(X_syn[[v]]))
    X_real_proc[[v]] <- as.integer(factor(X_real[[v]], levels = lvls))
    X_syn_proc[[v]]  <- as.integer(factor(X_syn[[v]],  levels = lvls))
  }
  
  # ============================================================
  # T1 — UNIVARIATE FIDELITY (JSD)
  # ============================================================
  
  T1 <- mean(sapply(names(X_real_proc), function(v) {
    p <- table(X_real_proc[[v]]) / nrow(X_real_proc)
    q <- table(X_syn_proc[[v]])  / nrow(X_syn_proc)
    lv <- union(names(p), names(q))
    p <- p[lv]; p[is.na(p)] <- 0
    q <- q[lv]; q[is.na(q)] <- 0
    JSD(p, q)
  }))
  
  # ============================================================
  # T2 — BIVARIATE FIDELITY
  # ============================================================
  
  T2_vals <- c()
  
  if (length(continuous_vars) > 1) {
    T2_vals <- c(T2_vals,
                 sqrt(sum((cor(X_real[, continuous_vars]) - cor(X_syn[, continuous_vars]))^2))
    )
  }
  
  if (length(categorical_vars) > 1) {
    cramer <- function(df) {
      n <- ncol(df)
      M <- matrix(1, n, n)
      for (i in 1:n) for (j in i:n) {
        if (i != j) {
          v <- vcd::assocstats(table(df[[i]], df[[j]]))$cramer
          M[i, j] <- M[j, i] <- v
        }
      }
      M
    }
    
    T2_vals <- c(T2_vals,
                 sqrt(sum((cramer(X_real[, categorical_vars]) - cramer(X_syn[, categorical_vars]))^2))
    )
  }
  
  T2 <- mean(T2_vals)
  
  # ============================================================
  # T3 — DISCRIMINATION
  # ============================================================
  
  X <- rbind(X_real, X_syn)
  
  if (modeling == "cox") {
    X <- X |> dplyr::select(-all_of(c(timevar, statusvar)))
  }
  
  y <- factor(c(rep(0, nrow(X_real)), rep(1, nrow(X_syn))))
  folds <- caret::createFolds(y, k = S)
  
  T3 <- mean(sapply(folds, function(idx) {
    train <- X[-idx, ]
    test  <- X[idx, ]
    y_tr  <- y[-idx]
    y_te  <- y[idx]
    
    fit <- glm(y_tr ~ ., data = data.frame(y_tr, train), family = binomial)
    prob <- predict(fit, test, type = "response")
    pROC::auc(y_te, prob)
  }))
  
  # ============================================================
  # T4 — STANDARDIZED DIFFERENCE
  # ============================================================
  
  if (modeling == "logistic") {
    fit_real <- glm(formula, data = X_real, family = binomial)
    fit_syn  <- glm(formula, data = X_syn,  family = binomial)
    
    coef_real <- coef(summary(fit_real))[-1, ]
    coef_syn  <- coef(summary(fit_syn))[-1, ]
    
    effect_real <- coef_real[, "Estimate"]
    se_real     <- coef_real[, "Std. Error"]
    effect_syn  <- coef_syn[, "Estimate"]
    se_syn      <- coef_syn[, "Std. Error"]
    
  } else {
    fit_real <- survival::coxph(formula, data = X_real)
    fit_syn  <- survival::coxph(formula, data = X_syn)
    
    sr <- summary(fit_real)
    ss <- summary(fit_syn)
    
    effect_real <- sr$coefficients[, "coef"]
    se_real     <- sr$coefficients[, "se(coef)"]
    effect_syn  <- ss$coefficients[, "coef"]
    se_syn      <- ss$coefficients[, "se(coef)"]
  }
  
  SDiff <- abs((effect_syn - effect_real) / sqrt(se_syn^2 + se_real^2))
  SD_q05 <- as.numeric(quantile(SDiff, probs = 0.05))
  SD_q95 <- as.numeric(quantile(SDiff, probs = 0.95))
  T4 <- as.numeric(mean(SDiff))
  
  # ============================================================
  # T5 / T6 — PREDICTIVE UTILITY
  # ============================================================
  
  target <- target_col
  folds <- caret::createFolds(X_real[[target]], k = S)
  
  auc_trtr <- auc_tstr <- auc_tsrtr <- numeric(S)
  
  for (s in seq_len(S)) {
    test  <- X_real[folds[[s]], ]
    train <- X_real[-folds[[s]], ]
    
    if (modeling == "logistic") {
      fit_trtr <- glm(formula, train, family = binomial)
      fit_tstr <- glm(formula, X_syn, family = binomial)
      fit_mix  <- glm(formula, rbind(train, X_syn), family = binomial)
      
      auc_trtr[s] <- pROC::auc(test[[target]], predict(fit_trtr, test, "response"))
      auc_tstr[s] <- pROC::auc(test[[target]], predict(fit_tstr, test, "response"))
      auc_tsrtr[s]<- pROC::auc(test[[target]], predict(fit_mix,  test, "response"))
      
    } else {
      Surv_test <- survival::Surv(test[[timevar]], test[[statusvar]])
      times <- quantile(train[[timevar]],probs = seq(0.1, 0.9, by=0.1), na.rm=TRUE)
      cox_auc <- function(fit, train_data) {
        lp_tr <- predict(fit, train_data)
        lp_te <- predict(fit, test)
        
        mean(survAUC::AUC.cd(
          survival::Surv(train_data[[timevar]], train_data[[statusvar]]),
          Surv_test, lp_tr, lp_te,times
        )$auc, na.rm = TRUE)
      }
      
      auc_trtr[s]  <- cox_auc(coxph(formula, train), train)
      auc_tstr[s]  <- cox_auc(coxph(formula, X_syn), X_syn)
      auc_tsrtr[s] <- cox_auc(coxph(formula, rbind(train, X_syn)), rbind(train, X_syn))
    }
  }
  
  T5 <- mean(auc_tstr)  - mean(auc_trtr)
  T6 <- mean(auc_tsrtr) - mean(auc_trtr)
  
  # ============================================================
  # T7 — RECORD MATCHING
  # ============================================================
  
  T7 <- mean(gower::gower_topn(x = X_real, y = X_syn, n = 1)$distance)
  
  # ============================================================
  # T8 — HOLDOUT PRIVACY
  # ============================================================
  row.names(trn) <-NULL
  
  # ------------------------------------------------------------
  # GENERATE SYNTHETIC DATA (depending on generator)
  # ------------------------------------------------------------
  synthetic_data <- Ablation_Analysis_generate_synthetic_data(
    data = trn,
    missing_rate = missing_rate,
    iterative_masking = iterative_masking,
    survival_data = is_survival,
    time_variable = timevar,
    event_variable = statusvar,
    n_imputations = 1,
    max_iterations = 5,
    random_seed = seed,
    continuous_method = "pmm",
    verbose = TRUE
    )
    
  syn <- synthetic_data[["synthetic_data"]]
  
  all <- bind_rows(trn, val, syn)
  all_hot <- dummyVars(" ~ .", data = all, fullRank = TRUE) |> predict(newdata = all) |> as.matrix()
  
  n_trn <- nrow(trn)
  n_val <- nrow(val)
  n_syn <- nrow(syn)
  
  trn_hot <- all_hot[1:n_trn, ]
  val_hot <- all_hot[(n_trn+1):(n_trn+n_val), ]
  syn_hot <- all_hot[(n_trn+n_val+1):(n_trn+n_val+n_syn), ]
  
  d_trn <- get.knnx(trn_hot, syn_hot, k = 1)$nn.dist[,1]
  d_val <- get.knnx(val_hot, syn_hot, k = 1)$nn.dist[,1]
  
  T8 <- mean(d_trn < d_val) + (n_trn/(n_trn+n_val))*mean(d_trn == d_val)
  
  # ============================================================
  # T9 — ATTRIBUTE INFERENCE
  # ============================================================
  
  AIR <- unlist(lapply(feature_associate, function(vars) {
    privacy_AIR(X_real, X_syn, target_col, vars)
  }))
  
  T9 <- mean(AIR)
  AIR_min <- min(AIR)
  AIR_max <- max(AIR)
  
  # ============================================================
  # OUTPUT
  # ============================================================
  
  return(c(
    T1=T1, T2=T2, T3=T3, T4=T4,
    T5=T5, T6=T6, T7=T7, T8=T8, T9=T9,
    T10=SD_q05, T11=SD_q95,
    T12=AIR_min, T13=AIR_max
  ))
}


# SCORE

# ===============================================================================================
# Multiple imputation function for survival datasets
# ===============================================================================================
#
# This function performs multiple imputation using the MICE framework
# on datasets containing censored survival information.
#
# Workflow
# --------
# 1. Estimate the cumulative baseline hazard using the Nelson–Aalen estimator.
# 2. Include the cumulative hazard as an auxiliary predictor during imputation.
# 3. Remove survival outcome variables before imputation to avoid leakage.
# 4. Perform multiple imputation using MICE.
# 5. Reattach survival outcome variables to each completed dataset.
#
# Parameters
# ----------
# df : data.frame
#     Input dataset containing missing values and survival information.
#
# nb_imputations : integer
#     Number of multiple imputations to generate.
#
# seed : numeric, default = 123
#     Random seed for reproducibility.
#
# Returns
# -------
# A list of completed imputed datasets.
#
# Notes
# -----
# - The cumulative hazard estimated by the Nelson–Aalen estimator
#   is incorporated as an auxiliary variable to preserve the
#   underlying survival-risk structure during imputation.
#
# - Variables `Follow_up_time` and `death` are excluded from the
#   imputation process and reintroduced afterward.
#
# - The function assumes:
#       * Follow_up_time = survival/censoring time
#       * Censorship     = event indicator (1 = observed event)
#
# ===============================================================================================

Imputed_data_function <- function(
    df,
    nb_imputations = NULL,
    seed = 123
) {
  
  # ---------------------------------------------------------------------------------------------
  # Copy input data
  # ---------------------------------------------------------------------------------------------
  
  data <- df
  
  
  # ---------------------------------------------------------------------------------------------
  # Estimate cumulative baseline hazard using Nelson–Aalen estimator
  # ---------------------------------------------------------------------------------------------
  
  H0 <- nelsonaalen(
    data,
    timevar   = "Follow_up_time",
    statusvar = "Censorship"
  )
  
  data$H0 <- H0
  
  
  # ---------------------------------------------------------------------------------------------
  # Remove survival outcome variables before imputation
  # ---------------------------------------------------------------------------------------------
  #
  # This avoids information leakage during the imputation process.
  #
  # ---------------------------------------------------------------------------------------------
  
  data <- data |>
    select(-c("Follow_up_time", "death"))
  
  
  # ---------------------------------------------------------------------------------------------
  # Initialize MICE methods and predictor matrix
  # ---------------------------------------------------------------------------------------------
  
  ini <- mice(data, maxit = 0)
  
  meth <- ini$method
  
  pred <- ini$predictorMatrix
  
  
  # ---------------------------------------------------------------------------------------------
  # Perform multiple imputation if missing values are present
  # ---------------------------------------------------------------------------------------------
  
  if (any(is.na(data))) {
    
    imputed_data <- mice(
      data,
      m               = nb_imputations,
      maxit           = 5,
      seed            = seed,
      method          = meth,
      predictorMatrix = pred,
      verbose         = FALSE
    )
    
    
    # -------------------------------------------------------------------------------------------
    # Extract all completed datasets
    # -------------------------------------------------------------------------------------------
    
    df.imp <- complete(
      imputed_data,
      action  = "all",
      include = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------------------------
  # Reattach survival variables removed before imputation
  # ---------------------------------------------------------------------------------------------
  
  for (i in 1:nb_imputations) {
    
    cat(sprintf(
      "── Imputation %d / %d\n",
      i,
      nb_imputations
    ))
    
    # Remove auxiliary cumulative hazard
    df.imp[[i]] <- df.imp[[i]] |>
      select(-c(H0))
    
    # Reintroduce survival variables
    df.imp[[i]][["Follow_up_time"]] <- df$Follow_up_time
    
    df.imp[[i]][["death"]] <- df$death
  }
  
  
  # ---------------------------------------------------------------------------------------------
  # Return completed datasets
  # ---------------------------------------------------------------------------------------------
  
  return(df.imp)
}



# ===============================================================================================
# Preprocessing pipeline for survival-data imputation
# ===============================================================================================
#
# This function prepares a survival dataset before multiple imputation.
#
# Workflow
# --------
# 1. Convert character variables into factors.
# 2. Rename survival-related variables.
# 3. Select variables used for modeling.
# 4. Create a derived binary mortality outcome.
# 5. Generate multiple imputed datasets using the
#    `Imputed_data_function()`.
#
# Parameters
# ----------
# df : data.frame
#     Input survival dataset.
#
# nb_imputations : integer, default = 50
#     Number of multiple imputations to generate.
#
# seed : numeric, default = 123
#     Random seed for reproducibility.
#
# Returns
# -------
# A list of imputed datasets.
#
# Notes
# -----
# - The binary outcome `death` is defined as:
#
#       death = "Yes"
#           if Follow_up_time < 90 days
#           and Censorship == 1
#
#       otherwise "No"
#
# - This preprocessing pipeline is specifically designed
#   for the REIN survival dataset used in the CHIMERA framework.
#
# ===============================================================================================

preprocessing <- function(
    df,
    nb_imputations = 50,
    seed = 123
) {
  
  # ---------------------------------------------------------------------------------------------
  # Convert character variables into factors
  # ---------------------------------------------------------------------------------------------
  
  df <- df |>
    mutate(
      across(
        where(is.character),
        as.factor
      )
    )
  
  
  # ---------------------------------------------------------------------------------------------
  # Rename survival variables
  # ---------------------------------------------------------------------------------------------
  
  df <- df %>%
    rename(
      Censorship    = Censored,
      Follow_up_time = times
    )
  
  
  # ---------------------------------------------------------------------------------------------
  # Select variables used in the analysis
  # ---------------------------------------------------------------------------------------------
  
  df <- df |>
    select(
      c(
        Age,
        Sex,
        Serum_albumin_level,
        Body_mass_index,
        Diabetes,
        Heart_failure,
        Peripheral_artery_disease,
        Coronary_artery_disease,
        Myocardial_infarction,
        Stroke,
        Cardiac_arrhythmia,
        Chronic_respiratory_failure,
        Malignancy,
        Cirrhosis,
        Severe_behavioral_disorders,
        Walking_autonomy,
        Follow_up_time,
        Censorship
      )
    )
  
  
  # ---------------------------------------------------------------------------------------------
  # Create derived mortality outcome
  # ---------------------------------------------------------------------------------------------
  #
  # Death within 90 days:
  #   - observed event
  #   - follow-up time < 90 days
  #
  # ---------------------------------------------------------------------------------------------
  
  df$death <- as.factor(
    ifelse(
      df$Follow_up_time < 90 &
        df$Censorship == 1,
      "Yes",
      "No"
    )
  )
  
  
  # ---------------------------------------------------------------------------------------------
  # Generate multiple imputed datasets
  # ---------------------------------------------------------------------------------------------
  
  results.imp <- Imputed_data_function(
    df,
    nb_imputations = nb_imputations,
    seed = seed
  )
  
  
  # ---------------------------------------------------------------------------------------------
  # Return imputed datasets
  # ---------------------------------------------------------------------------------------------
  
  return(results.imp)
}

# ===============================================================================================
# Bootstrap function for logistic regression coefficient significance
# ===============================================================================================
#
# This function is designed to be used within the `boot()` framework.
# For each bootstrap resample:
#   1. A logistic regression model is fitted.
#   2. Wald-test p-values for all coefficients are extracted.
#
# Parameters
# ----------
# data : data.frame
#     Input dataset containing predictors and outcome variable.
#
# indices : vector
#     Bootstrap resampling indices generated by `boot()`.
#
# Returns
# -------
# Numeric vector containing Wald-test p-values for all model coefficients.
#
# Notes
# -----
# - The outcome variable is assumed to be `death`.
# - Variables `Censorship` and `Follow_up_time` are excluded from the model.
#
# ===============================================================================================

bootstrap_model <- function(data, indices) {
  
  # Generate bootstrap sample
  d <- data[indices, ]
  
  # Fit logistic regression model
  model <- glm(
    death ~ .,
    data = d |> select(-c(Censorship, Follow_up_time)),
    family = binomial
  )
  
  # Extract Wald-test p-values
  wald_test <- summary(model)$coefficients[, 4]
  
  return(wald_test)
}



# ===============================================================================================
# Stability score based on bootstrap variable selection frequency
# ===============================================================================================
#
# This function evaluates the stability of variable significance across
# multiple imputed datasets using bootstrap resampling.
#
# Workflow
# --------
# 1. Apply bootstrap logistic regression independently to each imputed dataset.
# 2. Extract coefficient p-values from each bootstrap replicate.
# 3. Combine p-values for multi-level categorical variables using Fisher's method.
# 4. Compute the percentage of bootstrap replicates in which each variable
#    is statistically significant (p < 0.05).
#
# Parameters
# ----------
# df.imp : list
#     List of imputed datasets.
#
# Returns
# -------
# data.frame with:
#   - facteur_risque :
#       Variable name.
#
#   - nb_significative_coef :
#       Percentage of bootstrap replicates where the variable
#       is statistically significant.
#
# Notes
# -----
# - Fisher's method is used to combine p-values from dummy variables
#   corresponding to multi-level categorical predictors.
#
# - The percentage is computed over:
#       number_of_imputed_datasets × bootstrap_replicates
#
# - The current implementation assumes:
#       R = 100 bootstrap replicates
#       50 imputed datasets
#       => total = 5000 evaluations
#
# ===============================================================================================

Score <- function(df.imp) {
  
  # ---------------------------------------------------------------------------------------------
  # Apply bootstrap procedure to each imputed dataset
  # ---------------------------------------------------------------------------------------------
  
  res <- lapply(
    df.imp,
    function(data)
      boot(
        data = data,
        statistic = bootstrap_model,
        R = 100
      )
  )
  
  
  # ---------------------------------------------------------------------------------------------
  # Extract coefficient names from logistic regression model
  # ---------------------------------------------------------------------------------------------
  
  variable_names <- rownames(
    summary(
      glm(
        death ~ .,
        data = df.imp[[1]] |> select(-c(Censorship, Follow_up_time)),
        family = binomial
      )
    )$coefficients
  )
  
  
  # ---------------------------------------------------------------------------------------------
  # Combine bootstrap p-values across all imputed datasets
  # ---------------------------------------------------------------------------------------------
  
  colnames(res[[1]][["t"]]) <- variable_names
  
  df.pval <- res[[1]][["t"]]
  
  for(i in 2:length(df.imp)) {
    
    colnames(res[[i]][["t"]]) <- variable_names
    
    df.pval <- rbind(
      df.pval,
      res[[i]][["t"]]
    )
  }
  
  
  # ---------------------------------------------------------------------------------------------
  # Remove intercept term
  # ---------------------------------------------------------------------------------------------
  
  df.pval <- df.pval[, -1]
  df.pval <- as.data.frame(df.pval)
  
  
  # =============================================================================================
  # Fisher combination for multi-level categorical variables
  # =============================================================================================
  
  
  # ---------------------------------------------------------------------------------------------
  # Heart failure
  # ---------------------------------------------------------------------------------------------
  
  pvalues <- df.pval[, 6:7]
  
  X2_stat <- -2 * rowSums(log(pvalues))
  
  p_combined_fisher <- as.numeric(
    pchisq(X2_stat, df = 4, lower.tail = FALSE)
  )
  
  df.pval$Heart_failure <- p_combined_fisher
  
  
  # ---------------------------------------------------------------------------------------------
  # Peripheral artery disease
  # ---------------------------------------------------------------------------------------------
  
  pvalues <- df.pval[, 8:9]
  
  X2_stat <- -2 * rowSums(log(pvalues))
  
  p_combined_fisher <- as.numeric(
    pchisq(X2_stat, df = 4, lower.tail = FALSE)
  )
  
  df.pval$Peripheral_artery_disease <- p_combined_fisher
  
  
  # ---------------------------------------------------------------------------------------------
  # Walking autonomy
  # ---------------------------------------------------------------------------------------------
  
  pvalues <- df.pval[, 18:19]
  
  X2_stat <- -2 * rowSums(log(pvalues))
  
  p_combined_fisher <- as.numeric(
    pchisq(X2_stat, df = 4, lower.tail = FALSE)
  )
  
  df.pval$Walking_autonomy <- p_combined_fisher
  
  
  # ---------------------------------------------------------------------------------------------
  # Remove original dummy-variable p-values
  # ---------------------------------------------------------------------------------------------
  
  df.pval <- df.pval[, -c(6,7,8,9,18,19)]
  
  
  # =============================================================================================
  # Compute selection frequency
  # =============================================================================================
  
  
  # Count the number of significant bootstrap replicates
  proportions <- colSums(
    df.pval < 0.05,
    na.rm = TRUE
  )
  
  
  # Convert to output dataframe
  df.nbSignCoef <- data.frame(
    facteur_risque = colnames(df.pval)
  )
  
  df.nbSignCoef$nb_significative_coef <- as.numeric(proportions)
  
  
  # Convert counts into percentages
  df.nbSignCoef$nb_significative_coef <- round(
    (df.nbSignCoef$nb_significative_coef * 100) / 5000
  )
  
  return(df.nbSignCoef)
}



# ===============================================================================================
# Compute variable-selection agreement metrics
# ===============================================================================================
#
# This function compares the set of significant variables identified
# in real versus synthetic datasets.
#
# Variables are classified as "selected" if their stability score
# exceeds a predefined threshold.
#
# The following metrics are computed:
#
#   - Sensitivity (T1):
#       Ability of synthetic data to recover variables selected
#       in the real dataset.
#
#   - Specificity (T2):
#       Ability of synthetic data to avoid selecting variables
#       not selected in the real dataset.
#
#   - Cohen's Kappa (T3):
#       Agreement between real and synthetic selections
#       beyond chance.
#
# Parameters
# ----------
# tbl_df_real : data.frame
#     Stability scores computed on the real dataset.
#
# tbl_df_synth : data.frame
#     Stability scores computed on the synthetic dataset.
#
# real_col : character
#     Column name containing stability scores for the real dataset.
#
# synth_col : character
#     Column name containing stability scores for the synthetic dataset.
#
# threshold : numeric, default = 70
#     Selection threshold expressed as a percentage.
#
# Returns
# -------
# List containing:
#
#   T1 : sensitivity
#   T2 : specificity
#   T3 : Cohen's kappa
#
# ===============================================================================================

compute_selection_metrics <- function(
    tbl_df_real,
    tbl_df_synth,
    real_col,
    synth_col,
    threshold = 70
) {
  
  # ---------------------------------------------------------------------------------------------
  # Merge real and synthetic results
  # ---------------------------------------------------------------------------------------------
  
  df <- dplyr::left_join(
    tbl_df_real,
    tbl_df_synth,
    by = "facteur_risque"
  )
  
  
  # ---------------------------------------------------------------------------------------------
  # Identify selected variables
  # ---------------------------------------------------------------------------------------------
  
  real_vec <- df |>
    dplyr::filter(.data[[real_col]] > threshold) |>
    dplyr::pull(facteur_risque)
  
  synth_vec <- df |>
    dplyr::filter(.data[[synth_col]] > threshold) |>
    dplyr::pull(facteur_risque)
  
  
  # ---------------------------------------------------------------------------------------------
  # Convert selections to binary vectors
  # ---------------------------------------------------------------------------------------------
  
  all_vars <- df$facteur_risque
  
  real_bin  <- all_vars %in% real_vec
  synth_bin <- all_vars %in% synth_vec
  
  
  # ---------------------------------------------------------------------------------------------
  # Confusion matrix components
  # ---------------------------------------------------------------------------------------------
  
  TP <- sum(real_bin & synth_bin)
  TN <- sum(!real_bin & !synth_bin)
  FP <- sum(!real_bin & synth_bin)
  FN <- sum(real_bin & !synth_bin)
  
  
  # ---------------------------------------------------------------------------------------------
  # Sensitivity and specificity
  # ---------------------------------------------------------------------------------------------
  
  sensitivity <- TP / (TP + FN)
  
  specificity <- TN / (TN + FP)
  
  
  # ---------------------------------------------------------------------------------------------
  # Cohen's kappa
  # ---------------------------------------------------------------------------------------------
  
  total <- TP + TN + FP + FN
  
  po <- (TP + TN) / total
  
  pe <- (
    ((TP + FN) / total) * ((TP + FP) / total)
  ) +
    (
      ((TN + FP) / total) * ((TN + FN) / total)
    )
  
  kappa <- (po - pe) / (1 - pe)
  
  
  # ---------------------------------------------------------------------------------------------
  # Return metrics
  # ---------------------------------------------------------------------------------------------
  
  list(
    T1 = sensitivity,
    T2 = specificity,
    T3 = kappa
  )
}

clean_variable_names <- function(variable_names, data) {
  
  column_name <- colnames(data)  # Extract original column names from the reference dataset
  
  # For each variable name in the input list, try to find the original column
  clean_names <- sapply(variable_names, function(var) {
    
    # Find any column name that matches the beginning of the variable name
    matching <- column_name[sapply(column_name, function(col) startsWith(var, col))]
    
    # If there's at least one match, return the first one found
    if (length(matching) > 0) {
      return(matching[1])
    } else {
      return(NA)  # If no match is found, return NA
    }
  })
  
  # Remove NAs and ensure the result is unique
  unique(clean_names[!is.na(clean_names)])
}



########################################################################################
# Calibration curve functions
########################################################################################

# --------------------------------------------------------------------------------------
# Function: compute_calibration_synth
# --------------------------------------------------------------------------------------
# Description:
# This function evaluates calibration performance for multiple synthetic datasets.
# For each synthetic dataset:
#   1. Variables are selected based on bootstrap stability results.
#   2. A logistic regression model is fitted on the synthetic dataset.
#   3. Predictions are generated on the real test dataset.
#   4. Calibration curves and calibration statistics are computed.
#
# Calibration is assessed using:
#   - Restricted cubic spline (RCS)-based calibration curves
#   - Brier score
#   - Calibration intercept
#   - Calibration slope
#
# The function returns both:
#   - All calibration curves (long format)
#   - Calibration statistics for each synthetic dataset
#
# Arguments:
#   results.synth      : List of variable-selection results obtained with Score()
#   df_test            : External test dataset
#   outcome            : Name of the binary outcome variable
#   synthetic_data_list: List of synthetic datasets
#   grid               : Probability grid (unused but kept for compatibility)
#   seed               : Random seed for reproducibility
#
# Returns:
#   A list containing:
#     - calibration_curves:
#         Combined calibration curves for all synthetic datasets
#     - stats_per_imputation:
#         Calibration metrics for each synthetic dataset
# --------------------------------------------------------------------------------------

compute_calibration_synth <- function(
    results.synth,
    df_test,
    outcome,
    synthetic_data_list,
    grid = seq(0, 1, by = 0.01),
    seed = 123
) {
  
  M <- length(synthetic_data_list)
  
  all_curves <- vector("list", M)
  all_stats  <- vector("list", M)
  
  for (m in seq_len(M)) {
    
    cat(sprintf("── Synthetic dataset %d / %d\n", m, M))
    
    # ------------------------------------------------------------------------------
    # Variable selection based on bootstrap stability
    # ------------------------------------------------------------------------------
    df.nbSignCoef.synth <- results.synth[[m]] |>
      filter(nb_significative_coef > 70)
    
    name_variable <- clean_variable_names(
      df.nbSignCoef.synth$facteur_risque,
      df.imp.original[[m]]
    )
    
    name_variable <- append(name_variable, "death")
    
    # ------------------------------------------------------------------------------
    # Prepare test dataset
    # ------------------------------------------------------------------------------
    df_test_m <- df_test |>
      select(all_of(name_variable))
    
    df_test_m[[outcome]] <- ifelse(
      df_test_m[[outcome]] == "Yes",
      1,
      0
    )
    
    # ------------------------------------------------------------------------------
    # Prepare synthetic dataset
    # ------------------------------------------------------------------------------
    df <- synthetic_data_list[[m]] |>
      mutate(across(where(is.character), as.factor)) |>
      rename(
        Censorship   = censored,
        Follow_up_time = times
      ) |>
      mutate(
        death = as.factor(
          ifelse(Follow_up_time < 90 & Censorship == 1, "Yes", "No")
        )
      )
    
    # Remove survival-related variables before modeling
    df <- df |>
      select(-c(Follow_up_time, Censorship))
    
    # ------------------------------------------------------------------------------
    # Multiple imputation using MICE
    # ------------------------------------------------------------------------------
    ini  <- mice(df, maxit = 0)
    
    meth <- ini$method
    pred <- ini$predictorMatrix
    
    # Predictive mean matching for numeric variables
    meth[names(df)[sapply(df, is.numeric)]] <- "pmm"
    
    df <- complete(
      mice(
        df,
        m = 1,
        maxit = 5,
        seed = seed,
        method = meth,
        predictorMatrix = pred,
        printFlag = FALSE
      )
    )
    
    # Keep only selected variables
    df <- df |>
      select(all_of(name_variable))
    
    df[[outcome]] <- ifelse(df[[outcome]] == "Yes", 1, 0)
    
    # ------------------------------------------------------------------------------
    # Logistic regression model
    # ------------------------------------------------------------------------------
    formula <- as.formula(paste(outcome, "~ ."))
    
    glmFit <- glm(
      formula,
      data = df,
      family = binomial
    )
    
    # ------------------------------------------------------------------------------
    # Predictions on external test dataset
    # ------------------------------------------------------------------------------
    pHat  <- predict(glmFit, newdata = df_test_m, type = "response")
    yTest <- df_test_m[[outcome]]
    
    # ------------------------------------------------------------------------------
    # Calibration curve computation
    # ------------------------------------------------------------------------------
    calPerf <- val.prob.ci.2(
      pHat,
      yTest,
      smooth = "rcs"
    )
    
    calData <- calPerf$CalibrationCurves$RCS |>
      select(-knots) |>
      setNames(c("pHat", "obs", "ymin", "ymax"))
    
    # ------------------------------------------------------------------------------
    # Store calibration curve
    # ------------------------------------------------------------------------------
    all_curves[[m]] <- calData |>
      mutate(imputation = m)
    
    # ------------------------------------------------------------------------------
    # Store calibration statistics
    # ------------------------------------------------------------------------------
    all_stats[[m]] <- data.frame(
      m         = m,
      Brier     = calPerf[["stats"]][["Brier"]],
      intercept = calPerf[["stats"]][["Intercept"]],
      slope     = calPerf[["stats"]][["Slope"]]
    )
  }
  
  # ------------------------------------------------------------------------------
  # Return results
  # ------------------------------------------------------------------------------
  return(list(
    calibration_curves   = bind_rows(all_curves),
    stats_per_imputation = do.call(rbind, all_stats)
  ))
}


# --------------------------------------------------------------------------------------
# Function: compute_calibration
# --------------------------------------------------------------------------------------
# Description:
# This function evaluates calibration performance using multiply imputed
# real datasets.
#
# For each imputed dataset:
#   1. Variables are selected using bootstrap stability analysis.
#   2. A logistic regression model is fitted.
#   3. Predictions are generated on an external test dataset.
#   4. Calibration curves and calibration metrics are computed.
#
# Calibration metrics include:
#   - Brier score
#   - Calibration intercept
#   - Calibration slope
#
# Arguments:
#   results.original : Variable-selection results obtained with Score()
#   df_test          : External test dataset
#   outcome          : Name of the binary outcome variable
#   df.imp.original  : List of imputed datasets
#   grid             : Probability grid (unused but retained for compatibility)
#   seed             : Random seed for reproducibility
#
# Returns:
#   A list containing:
#     - calibration_curves:
#         Combined calibration curves for all imputations
#     - stats_per_imputation:
#         Calibration statistics for each imputed dataset
# --------------------------------------------------------------------------------------

compute_calibration <- function(
    results.original,
    df_test,
    outcome,
    df.imp.original,
    grid = seq(0, 1, by = 0.01),
    seed = 123
) {
  
  M <- length(df.imp.original)
  
  all_curves <- vector("list", M)
  all_stats  <- vector("list", M)
  
  for (m in seq_len(M)) {
    
    cat(sprintf("── Imputation %d / %d\n", m, M))
    
    # ------------------------------------------------------------------------------
    # Variable selection based on bootstrap stability
    # ------------------------------------------------------------------------------
    df.nbSignCoef <- results.original |>
      filter(nb_significative_coef > 70)
    
    name_variable <- clean_variable_names(
      df.nbSignCoef$facteur_risque,
      df.imp.original[[m]]
    )
    
    name_variable <- append(name_variable, "death")
    
    # ------------------------------------------------------------------------------
    # Prepare external test dataset
    # ------------------------------------------------------------------------------
    df_test_m <- df_test |>
      select(all_of(name_variable))
    
    df_test_m[[outcome]] <- ifelse(
      df_test_m[[outcome]] == "Yes",
      1,
      0
    )
    
    # ------------------------------------------------------------------------------
    # Prepare imputed dataset
    # ------------------------------------------------------------------------------
    df <- df.imp.original[[m]]
    
    # Remove survival-related variables before modeling
    df <- df |>
      select(-c(Follow_up_time, Censorship))
    
    # ------------------------------------------------------------------------------
    # Multiple imputation using MICE
    # ------------------------------------------------------------------------------
    ini  <- mice(df, maxit = 0)
    
    meth <- ini$method
    pred <- ini$predictorMatrix
    
    # Predictive mean matching for numeric variables
    meth[names(df)[sapply(df, is.numeric)]] <- "pmm"
    
    df <- complete(
      mice(
        df,
        m = 1,
        maxit = 5,
        seed = seed,
        method = meth,
        predictorMatrix = pred,
        printFlag = FALSE
      )
    )
    
    # Keep selected variables only
    df <- df |>
      select(all_of(name_variable))
    
    df[[outcome]] <- ifelse(df[[outcome]] == "Yes", 1, 0)
    
    # ------------------------------------------------------------------------------
    # Logistic regression model
    # ------------------------------------------------------------------------------
    formula <- as.formula(paste(outcome, "~ ."))
    
    glmFit <- glm(
      formula,
      data = df,
      family = binomial
    )
    
    # ------------------------------------------------------------------------------
    # Predictions on external test dataset
    # ------------------------------------------------------------------------------
    pHat  <- predict(glmFit, newdata = df_test_m, type = "response")
    yTest <- df_test_m[[outcome]]
    
    # ------------------------------------------------------------------------------
    # Calibration curve computation
    # ------------------------------------------------------------------------------
    calPerf <- val.prob.ci.2(
      pHat,
      yTest,
      smooth = "rcs"
    )
    
    calData <- calPerf$CalibrationCurves$RCS |>
      select(-knots) |>
      setNames(c("pHat", "obs", "ymin", "ymax"))
    
    # ------------------------------------------------------------------------------
    # Store calibration curve
    # ------------------------------------------------------------------------------
    all_curves[[m]] <- calData |>
      mutate(imputation = m)
    
    # ------------------------------------------------------------------------------
    # Store calibration statistics
    # ------------------------------------------------------------------------------
    all_stats[[m]] <- data.frame(
      m         = m,
      Brier     = calPerf[["stats"]][["Brier"]],
      intercept = calPerf[["stats"]][["Intercept"]],
      slope     = calPerf[["stats"]][["Slope"]]
    )
  }
  
  # ------------------------------------------------------------------------------
  # Return results
  # ------------------------------------------------------------------------------
  return(list(
    calibration_curves   = bind_rows(all_curves),
    stats_per_imputation = do.call(rbind, all_stats)
  ))
}