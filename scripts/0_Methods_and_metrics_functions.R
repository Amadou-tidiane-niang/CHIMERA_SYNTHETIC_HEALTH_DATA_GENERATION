# ============================================================
# Useful function 
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


# Build file paths for synthetic datasets stored as individual CSV files.

default_file_paths <- function(dataset, synthesizer, suffix, n_synth) {
  folder <- here(
    data_dir,
    sprintf("%s_synthetic_datasets_%s%s", toupper(dataset), synthesizer, suffix)
  )
  sprintf("%s/%s_%s_%02d.csv", folder, dataset, synthesizer, seq_len(n_synth))
}

# Build file paths for TabDiff datasets stored in separate run folders

tabdiff_file_paths <- function(dataset, synthesizer, suffix, n_synth) {
  folder <- here(
    data_dir,
    sprintf("%s_synthetic_datasets_%s%s", toupper(dataset), synthesizer, suffix)
  )
  sprintf("%s/run_%d/samples.csv", folder, seq_len(n_synth))
}


# Select the appropriate file-path builder for a synthesizer.

build_synthetic_file_paths <- function(dataset, synthesizer, suffix, n_synth) {
  builder <- file_path_builders[[synthesizer]]
  if (is.null(builder)) {
    stop("No file-path builder registered for synthesizer: ", synthesizer)
  }
  builder(dataset, synthesizer, suffix, n_synth)
}

# Convert character columns to factors and recode "nothing" as "None".

default_postprocess <- function(df) {
  df |>
    as.data.frame() |>
    dplyr::mutate(
      dplyr::across(
        dplyr::where(is.character),
        as.factor
      )
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::where(is.factor),
        ~ dplyr::recode(.x, "nothing" = "None")
      )
    )
}

# Return the post-processing function for a dataset and synthesizer.
get_postprocess_fn <- function(dataset, synthesizer) {
  default_postprocess
}

# Read and post-process all synthetic datasets for one combination.

read_synthetic_datasets <- function(dataset, synthesizer, suffix, n_synth) {
  files <- build_synthetic_file_paths(dataset, synthesizer, suffix, n_synth)
  
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    stop(sprintf(
      "[%s - %s%s] %d/%d file(s) not found, e.g.: %s",
      dataset, synthesizer, suffix, length(missing_files), n_synth, missing_files[1]
    ))
  }
  
  postprocess_fn <- get_postprocess_fn(dataset, synthesizer)
  
  synth_list <- lapply(files, readr::read_csv, show_col_types = FALSE)
  synth_list <- lapply(synth_list, postprocess_fn)
  
  synth_list
}

# Generate the standardised object name for a dataset and synthesizer.

synthetic_object_name <- function(dataset, synthesizer, suffix) {
  sprintf("%s_synthetic_list_data_%s%s", tolower(dataset), tolower(synthesizer), suffix)
}

# Save a list of synthetic datasets as an RData file.

save_synthetic_datasets <- function(synth_list, dataset, synthesizer, suffix) {
  obj_name <- synthetic_object_name(dataset, synthesizer, suffix)
  out_path <- here(results_dir, paste0(obj_name, ".Rdata"))
  
  assign(obj_name, synth_list)
  save(list = obj_name, file = out_path)
  
  message(sprintf(
    "Saved %-45s (%d datasets) -> %s", obj_name, length(synth_list), out_path
  ))
  
  invisible(out_path)
}

# Read, post-process, and save one synthetic dataset combination.

process_synthetic_combination <- function(dataset, synthesizer, suffix, n_synth) {
  synth_list <- read_synthetic_datasets(dataset, synthesizer, suffix, n_synth)
  out_path <- save_synthetic_datasets(synth_list, dataset, synthesizer, suffix)
  rm(synth_list)
  out_path
}


# Import and save all dataset, synthesizer, and suffix combinations.

import_all_synthetic_datasets <- function(datasets, synthesizers, suffixes, n_synth) {
  saved_paths <- character(0)
  
  for (dataset in datasets) {
    for (synthesizer in synthesizers) {
      for (suffix_label in names(suffixes)) {
        suffix <- suffixes[[suffix_label]]
        
        message(sprintf("Processing: %s | %s | %s", dataset, synthesizer, suffix_label))
        
        out_path <- process_synthetic_combination(
          dataset     = dataset,
          synthesizer = synthesizer,
          suffix      = suffix,
          n_synth     = n_synth
        )
        
        saved_paths <- c(saved_paths, out_path)
      }
    }
  }
  
  invisible(saved_paths)
}

# load data
load_first <- function(file) {
  
  e <- new.env()
  
  object_name <- load(
    here::here(results_dir, file),
    envir = e
  )
  
  get(
    object_name[1],
    envir = e
  )
}

# impute data with mice
impute_mice <- function(
    data,
    m = 1,
    maxit = 5,
    seed = 123,
    method = NULL,
    predictorMatrix = NULL,
    parallel = TRUE,
    workers = NULL,
    printFlag = TRUE
) {
  
  #  Check input
  
  
  is_list <- is.list(data) && !is.data.frame(data)
  
  if (!is_list && !is.data.frame(data)) {
    stop("'data' must be a data.frame or a list of data.frames.")
  }
  
  if (is_list && !all(vapply(data, is.data.frame, logical(1)))) {
    stop("All elements of 'data' must be data.frames.")
  }
  
  data_list <- if (is_list) data else list(data)
  
  
  # Determine MICE settings automatically

  
  if (is.null(method) || is.null(predictorMatrix)) {
    
    ini <- mice(
      data_list[[1]],
      maxit = 0,
      printFlag = FALSE
    )
    
    if (is.null(method)) {
      method <- ini$method
    }
    
    if (is.null(predictorMatrix)) {
      predictorMatrix <- ini$predictorMatrix
    }
  }
  
  
  # Define imputation function
  
  impute_one <- function(dat) {
    
    imp <- mice(
      dat,
      m = m,
      maxit = maxit,
      seed = seed,
      method = method,
      predictorMatrix = predictorMatrix,
      printFlag = printFlag
    )
    
    complete(imp, 1)
  }
  
  
  # Parallel or sequential execution
  
  if (parallel && length(data_list) > 1) {
    
    if (is.null(workers)) {
      workers <- max(1, future::availableCores() - 1)
    }
    
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    
    future::plan(
      future::multisession,
      workers = workers
    )
    
    result <- future.apply::future_lapply(
      data_list,
      impute_one,
      future.seed = TRUE
    )
    
  } else {
    
    set.seed(seed)
    
    result <- lapply(
      data_list,
      impute_one
    )
  }
  
  
  # Return same structure as input
  
  if (!is_list) {
    return(result[[1]])
  }
  
  names(result) <- names(data)
  
  return(result)
}

# ============================================================
# 1. CHIMERA: Synthetic Data Generation Function
# ============================================================

CHIMERA_generate <-function(df,missing_rate = 0.10,is_survival = TRUE,timevar = NULL, statusvar = NULL,nb_imput = 1,nb_max_it = 10,seed = NULL,verbose = TRUE){
  
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
    #meth[names(df)[sapply(df, is.numeric)]] <- method_conti
    
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
  #meth[names(df_ref)[sapply(df_ref, is.numeric)]] <- method_conti
  
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
  
  # NA mask
  # Backup of Real missing values
  if (any(is.na(df))) {
    df_syn[na_indices_initial] <- NA
  }
  
  # --- If survival data, compute Nelson-Aalen estimate
  if (is_survival) {
    df <- df |> rename(Censored = all_of(statusvar), times = all_of(timevar))
    df$H0 <- nelsonaalen(df, timevar = "times", statusvar = "Censored")
    df_ref <- df[, !(names(df) %in% "times")]
  } else {
    df_ref <- df
  }
  
  # --- Return
  return(df_syn)
  
  
  
}



compute_metrics_loss <- function(results_all,method_name,n_synth = NULL,w_fid = 1/3,w_util = 1/3,w_conf = 1/3) {
  
  
  # Convert results to data frame
  results_df <- as.data.frame(results_all)
  
  
  # Metrics used for the composite loss
  metrics <- paste0("T", 1:9)
  
  # Preserve T10-T13
  additional_metrics <- paste0("T", 10:13)
  
  # Target vector
  T_star <- c(T1 = 0,T2 = 0,T3 = 0.5,T4 = 0,T5 = 0,T6 = 0,T7 = 1,T8 = 0.5,T9 = 0)
  
  # Metric families
  F_fid <- c("T1","T2","T3")
  F_util <- c("T4","T5","T6")
  F_conf <- c("T7","T8","T9")
  
  # Extract T1-T9
  metrics_df <- results_df |>
    dplyr::select(
      dplyr::all_of(metrics)
    )
  
  # Standard deviations
  s_j <- sapply(metrics_df,sd,na.rm = TRUE)
  
  # Check zero or undefined SD
  if (any(s_j == 0 | is.na(s_j))) {
    problematic_metrics <- names(s_j)[s_j == 0 | is.na(s_j)]
    
    stop(
      paste0(
        "Zero or undefined standard deviation for: ",
        paste(
          problematic_metrics,
          collapse = ", "
        ),
        ". Standardized loss cannot be computed."
      )
    )
  }
  
  # Standardized squared loss
  loss_matrix <- sweep(metrics_df,2,T_star,FUN = "-")
  loss_matrix <- sweep(loss_matrix^2,2,s_j^2,FUN = "/")
  
  # Family-specific losses
  results_df$L_fid <- rowMeans(loss_matrix[,F_fid,drop = FALSE],na.rm = TRUE)
  results_df$L_util <- rowMeans(loss_matrix[,F_util,drop = FALSE],na.rm = TRUE)
  results_df$L_conf <- rowMeans(loss_matrix[,F_conf,drop = FALSE],na.rm = TRUE)
  
  # Check weights
  if (!isTRUE(all.equal(w_fid + w_util + w_conf,1))) {
    stop(
      "w_fid, w_util and w_conf must sum to 1."
    )
  }
  
  # Total loss
  results_df$L_total <- w_fid  * results_df$L_fid + w_util * results_df$L_util + w_conf * results_df$L_conf
  
  # Dataset identifiers
  results_df$dataset <- paste0("syn_",seq_len(n_synth))
  
  # Method identifier
  results_df$method <- method_name
  
  
  
  # Return final data frame
  return(results_df)
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


# Function to compute time-dependent AUC with 95% bootstrap CI
predictive_utility_surv_curve <- function(train, test, time_points = seq(0, 1000, by = 100), n_boot = 200, conf_level = 0.95) {
  
  # Fit Cox proportional hazards model on training and testing data
  cox_model <- coxph(Surv(times, Censored) ~ ., data = train)
  
  # Predict linear predictors on train/test
  lp_train <- predict(cox_model, newdata = train)
  lp_test  <- predict(cox_model, newdata = test)
  
  # Create survival objects
  Surv_train <- Surv(train$times, train$Censored)
  Surv_test  <- Surv(test$times, test$Censored)
  
  # Original AUC values
  auc_cd <- AUC.cd(Surv.rsp = Surv_train,
                   Surv.rsp.new = Surv_test,
                   lp = lp_train,
                   lpnew = lp_test,
                   times = time_points)$auc
  
  # Bootstrap to estimate CI
  auc_boot <- matrix(NA, nrow = n_boot, ncol = length(time_points))
  
  for (b in seq_len(n_boot)) {
    # Resample test set (with replacement)
    boot_idx <- sample(seq_len(nrow(test)), replace = TRUE)
    test_boot <- test[boot_idx, ]
    
    # Predict on bootstrap sample
    lp_test_boot <- predict(cox_model, newdata = test_boot)
    Surv_test_boot <- Surv(test_boot$times, test_boot$Censored)
    
    # Compute AUC for bootstrap sample
    auc_boot[b, ] <- AUC.cd(Surv.rsp = Surv_train,
                            Surv.rsp.new = Surv_test_boot,
                            lp = lp_train,
                            lpnew = lp_test_boot,
                            times = time_points)$auc
  }
  
  # Confidence intervals from bootstrap distribution
  alpha <- 1 - conf_level
  ci_lower <- apply(auc_boot, 2, quantile, probs = alpha / 2, na.rm = TRUE)
  ci_upper <- apply(auc_boot, 2, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  
  # Return results as a data frame
  return(data.frame(
    time  = time_points,
    auc   = auc_cd,
    lower = ci_lower,
    upper = ci_upper
  ))
}


##==============================================================================
# Evaluation metrics for synthetics data
##==============================================================================

# Attribute Inference Risk (AIR)

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


## -----------------------------------------------------------------------
## Core distance computation : DCR and nndr
## -----------------------------------------------------------------------
compute_dcr <- function(synthetic, reference, eps = 1e-8, weights = NULL) {
  stopifnot(
    nrow(reference) >= 2,
    all(names(synthetic) %in% names(reference))
  )
  
  reference <- reference[, names(synthetic), drop = FALSE]
  
  res <- gower::gower_topn(
    x = synthetic,
    y = reference,
    n = 2,
    weights = weights
  )
  
  dcr1 <- as.numeric(res$distance[1, ])
  dcr2 <- as.numeric(res$distance[2, ])
  
  # Use eps only when the second-nearest-neighbor distance is zero
  nndr <- ifelse(
    dcr2 == 0,
    dcr1 / eps,
    dcr1 / dcr2
  )
  
  tibble::tibble(
    dcr1 = dcr1,
    dcr2 = dcr2,
    nndr = nndr
  )
}





# ------------------------------------------------------------
# linkability_attack : vectorisee + sous-echantillonnage.
# Fournie mais DESACTIVEE par defaut dans evaluate_synthetic(),
# comme dans votre pipeline de production actuel.
# ------------------------------------------------------------

# Sous-echantillonnage utilise par les metriques O(n*m)
subsample_df <- function(df, max_n, seed = NULL) {
  if (nrow(df) <= max_n) return(df)
  if (!is.null(seed)) set.seed(seed)
  df[sample(seq_len(nrow(df)), size = max_n, replace = FALSE), , drop = FALSE]
}

batch_top_k <- function(targets_df, synthetic_df, k, weights = NULL) {
  res <- gower::gower_topn(x = targets_df, y = synthetic_df, n = k, weights = weights)
  list(index = res$index, distance = res$distance)
}

linkability_attack <- function(train_data, synthetic_data, A, B, control_data,
                               n_targets = 100, k = 1, weights = NULL, seed = 1234,
                               max_syn_n = 5000) {
  
  set.seed(seed)
  if (length(intersect(A, B)) > 0) stop("A and B must be disjoint.")
  if (k > nrow(synthetic_data)) k <- nrow(synthetic_data)
  
  synthetic_data <- subsample_df(synthetic_data, max_syn_n, seed + 2)
  
  n_targets_real    <- min(n_targets, nrow(train_data))
  n_targets_control <- min(n_targets, nrow(control_data))
  
  idx_real    <- sample(seq_len(nrow(train_data)),   n_targets_real)
  idx_control <- sample(seq_len(nrow(control_data)), n_targets_control)
  
  targets_real    <- train_data[idx_real, , drop = FALSE]
  targets_control <- control_data[idx_control, , drop = FALSE]
  
  synthetic_A <- synthetic_data[, A, drop = FALSE]
  synthetic_B <- synthetic_data[, B, drop = FALSE]
  
  run_attack <- function(targets, target_idx, name) {
    n_t <- nrow(targets)
    nb_A <- batch_top_k(targets[, A, drop = FALSE], synthetic_A, k, weights)
    nb_B <- batch_top_k(targets[, B, drop = FALSE], synthetic_B, k, weights)
    idx_A <- matrix(as.integer(nb_A$index), nrow = k)
    idx_B <- matrix(as.integer(nb_B$index), nrow = k)
    link <- vapply(seq_len(n_t), function(i) length(intersect(idx_A[, i], idx_B[, i])) > 0, logical(1))
    tibble::tibble(attack = name, target_id = seq_len(n_t), original_id = target_idx, k = k, successful_link = link)
  }
  
  training_results <- run_attack(targets_real,    idx_real,    "training")
  control_results  <- run_attack(targets_control, idx_control, "control")
  
  n_syn <- nrow(synthetic_data)
  idx_A_n <- t(replicate(n_targets_real, sample(n_syn, k)))
  idx_B_n <- t(replicate(n_targets_real, sample(n_syn, k)))
  if (k == 1) { idx_A_n <- matrix(idx_A_n, ncol = 1); idx_B_n <- matrix(idx_B_n, ncol = 1) }
  link_naive <- vapply(seq_len(n_targets_real), function(i) length(intersect(idx_A_n[i, ], idx_B_n[i, ])) > 0, logical(1))
  naive_results <- tibble::tibble(attack = "naive", target_id = seq_len(n_targets_real),
                                  original_id = idx_real, k = k, successful_link = link_naive)
  
  attack_results <- dplyr::bind_rows(training_results, control_results, naive_results)
  
  evaluation <- attack_results %>%
    dplyr::group_by(attack, k) %>%
    dplyr::summarise(n_targets = dplyr::n(), n_successful_links = sum(successful_link),
                     linkability_rate = mean(successful_link), .groups = "drop")
  
  comparison <- evaluation %>%
    dplyr::select(attack, k, linkability_rate) %>%
    tidyr::pivot_wider(names_from = attack, values_from = linkability_rate) %>%
    dplyr::mutate(training_excess = training - naive, control_excess = control - naive)
  
  list(evaluation = evaluation, comparison = comparison)
}

compute_rrd_thresholds <- function(X_train, X_test,alpha = 0.02, seed = 0) {
  
  if (alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be between 0 and 1.")
  }
  
  # Compute Real-to-Real distances
  
  rrd <- compute_dcr(X_train, X_test)
  
  # Compute the alpha-th percentile thresholds
  
  rrd_alpha_dcr <- quantile(
    rrd$dcr1,
    probs = alpha,
    na.rm = TRUE
  )
  
  rrd_alpha_nndr <- quantile(
    rrd$nndr,
    probs = alpha,
    na.rm = TRUE
  )
  
  n_rrd <- nrow(X_train)
  
  # Return the thresholds
 
  
  return(list(
    rrd_alpha_dcr = as.numeric(rrd_alpha_dcr),
    rrd_alpha_nndr = as.numeric(rrd_alpha_nndr),
    n_rrd          = as.numeric(n_rrd)
  ))
}

# ==============================================================
# Assessment_function_unified
# ==============================================================
# PURPOSE:
# Comprehensive evaluation of synthetic data quality
# covering fidelity, utility, and privacy risks.
# ==============================================================

Assessment_function_unified <- function(
    X_real,               # real dataset
    X_train,              # train dataset
    X_syn,                # test dataset
    X_test,               # test dataset
    S = 5,                # 5-fold cross validation
    formula,
    modeling = c("logistic", "cox"),
    timevar = NULL,
    statusvar = NULL,
    target_col,
    feature_associate,
    is_survival=FALSE,
    trn   = trn,
    val   = hol,
    syn_trn = NULL,
    seed,
    alpha=NULL,
    beta=NULL,
    n_rrd=NULL,
    k=1,               # k-neighboors
    n_workers = 8, 
    run_linkability = TRUE
) {
  
  # PREPROCESSING
  
  X_train <- X_train |> dplyr::mutate(across(where(is.character), as.factor))
  X_syn  <- X_syn  |> dplyr::mutate(across(where(is.character), as.factor))
  
  # VARIABLE TYPES
  
  continuous_vars  <- names(X_train)[sapply(X_train, is.numeric)]
  categorical_vars <- names(X_train)[sapply(X_train, is.factor)]
  
  
  # PREPROCESSING (BINNING + ENCODING)

  
  X_train_proc <- X_train
  X_syn_proc  <- X_syn
  
  for (v in continuous_vars) {
    brks <- fd_bins(X_train[[v]])
    X_train_proc[[v]] <- cut(X_train[[v]], brks, labels = FALSE, include.lowest = TRUE)
    X_syn_proc[[v]]  <- cut(X_syn[[v]],  brks, labels = FALSE, include.lowest = TRUE)
  }
  
  for (v in categorical_vars) {
    lvls <- union(levels(X_train[[v]]), levels(X_syn[[v]]))
    X_train_proc[[v]] <- as.integer(factor(X_train[[v]], levels = lvls))
    X_syn_proc[[v]]  <- as.integer(factor(X_syn[[v]],  levels = lvls))
  }
  
  # ============================================================
  # T1 — UNIVARIATE FIDELITY (JSD)
  # ============================================================
  
  JSD <- mean(sapply(names(X_train_proc), function(v) {
    p <- table(X_train_proc[[v]]) / nrow(X_train_proc)
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
                 sqrt(sum((cor(X_train[, continuous_vars]) - cor(X_syn[, continuous_vars]))^2))
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
                 sqrt(sum((cramer(X_train[, categorical_vars]) - cramer(X_syn[, categorical_vars]))^2))
    )
  }
  
  CFD <- mean(T2_vals)
  
  # ============================================================
  # T3 — DISCRIMINATION
  # ============================================================
  
  X <- rbind(X_train, X_syn)
  
  if (modeling == "cox") {
    X <- X |> dplyr::select(-all_of(c(timevar, statusvar)))
  }
  
  y <- factor(c(rep(0, nrow(X_train)), rep(1, nrow(X_syn))))
  folds <- caret::createFolds(y, k = S)
  
  auc_mean <- mean(sapply(folds, function(idx) {
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
    fit_real <- glm(formula, data = X_train, family = binomial)
    fit_syn  <- glm(formula, data = X_syn,  family = binomial)
    
    coef_real <- coef(summary(fit_real))[-1, ]
    coef_syn  <- coef(summary(fit_syn))[-1, ]
    
    effect_real <- coef_real[, "Estimate"]
    se_real     <- coef_real[, "Std. Error"]
    effect_syn  <- coef_syn[, "Estimate"]
    se_syn      <- coef_syn[, "Std. Error"]
    
  } else {
    fit_real <- survival::coxph(formula, data = X_train)
    fit_syn  <- survival::coxph(formula, data = X_syn)
    
    sr <- summary(fit_real)
    ss <- summary(fit_syn)
    
    effect_real <- sr$coefficients[, "coef"]
    se_real     <- sr$coefficients[, "se(coef)"]
    effect_syn  <- ss$coefficients[, "coef"]
    se_syn      <- ss$coefficients[, "se(coef)"]
  }
  
  SDiff <- abs((effect_syn - effect_real) / sqrt(se_syn^2 + se_real^2))
  SDiff_q05 <- as.numeric(quantile(SDiff, probs = 0.05))
  SDiff_q25 <- as.numeric(quantile(SDiff, probs = 0.25))
  SDiff_median <- as.numeric(quantile(SDiff, probs = 0.50))
  SDiff_q75 <- as.numeric(quantile(SDiff, probs = 0.75))
  SDiff_q95 <- as.numeric(quantile(SDiff, probs = 0.95))
  SDiff_mean <- as.numeric(mean(SDiff))
  
  # ============================================================
  # T5 / T6 — PREDICTIVE UTILITY
  # ============================================================
  
  target <- target_col
  
  
  # ============================================================
  # LOGISTIC REGRESSION
  # ============================================================
  
  if (modeling == "logistic") {
    
    # ----------------------------------------------------------
    # TRTR: Train on Real, Test on Real
    # ----------------------------------------------------------
    fit_trtr <- glm(formula,data = X_train,family = binomial)
    pred_trtr <- predict(fit_trtr,newdata = X_test,type = "response")
    auc_trtr <- as.numeric(pROC::auc(response = X_test[[target]],predictor = pred_trtr))
    
    
    # ----------------------------------------------------------
    # TSTR: Train on Synthetic, Test on Real
    # ----------------------------------------------------------
    fit_tstr <- glm(formula,data = X_syn,family = binomial)
    pred_tstr <- predict(fit_tstr,newdata = X_test,type = "response")
    auc_tstr <- as.numeric(pROC::auc(response = X_test[[target]],predictor = pred_tstr))
    
    
    # ----------------------------------------------------------
    # TSRTR: Train on Synthetic + Real, Test on Real
    # ----------------------------------------------------------
    train_mix <- rbind(X_train,X_syn)
    fit_tsrtr <- glm(formula,data = train_mix,family = binomial)
    pred_tsrtr <- predict(fit_tsrtr,newdata = X_test,type = "response")
    auc_tsrtr <- as.numeric(pROC::auc(response = X_test[[target]],predictor = pred_tsrtr))
    
    
    # ============================================================
    # SURVIVAL ANALYSIS
    # ============================================================
    
  } else {
    
    # ----------------------------------------------------------
    # Independent real test set
    # ----------------------------------------------------------
    Surv_test <- survival::Surv(X_test[[timevar]],X_test[[statusvar]])
    
    # Evaluation times are determined exclusively from the
    # real training data to avoid information leakage.
    times <- quantile(X_train[[timevar]],probs = seq(0.1, 0.9, by = 0.1),na.rm = TRUE)
    
    
    # ----------------------------------------------------------
    # Time-dependent AUC
    # ----------------------------------------------------------
    cox_auc <- function(fit, train_data) {
      lp_tr <- predict(fit,newdata = train_data)
      lp_te <- predict(fit,newdata = X_test)
      mean(survAUC::AUC.cd(
        survival::Surv(
          train_data[[timevar]],
          train_data[[statusvar]]
        ),
        Surv_test,
        lp_tr,
        lp_te,
        times
      )$auc,
      na.rm = TRUE
      )
    }
    
    
    # ----------------------------------------------------------
    # TRTR: Train on Real, Test on Real
    # ----------------------------------------------------------
    fit_trtr <- survival::coxph(formula,data = X_train)
    auc_trtr <- cox_auc(fit_trtr,X_train)
    
    
    # ----------------------------------------------------------
    # TSTR: Train on Synthetic, Test on Real
    # ----------------------------------------------------------
    fit_tstr <- survival::coxph(formula,data = X_syn)
    auc_tstr <- cox_auc(fit_tstr,X_syn)
    
    
    # ----------------------------------------------------------
    # TSRTR: Train on Synthetic + Real, Test on Real
    # ----------------------------------------------------------
    train_mix <- rbind(X_train,X_syn)
    fit_tsrtr <- survival::coxph(formula,data = train_mix)
    auc_tsrtr <- cox_auc(fit_tsrtr,train_mix)
  }
  
  
  # T5: difference between TSTR and the real-data baseline
  delta_tstr <- abs(auc_tstr - auc_trtr)
  
  # T6: difference between TSRTR and the real-data baseline
  delta_tsrtr <- abs(auc_tsrtr - auc_trtr)
  
  
  # ============================================================
  # DCR
  # ============================================================
  srd <- compute_dcr(synthetic = X_syn,reference = X_real)
  srd_n_below_dcr <- sum(srd$dcr1 < alpha)
  rrd_n_below_dcr <- (0.02 * n_rrd)
  dcr <- srd_n_below_dcr /rrd_n_below_dcr 
  
  res <- compute_dcr(synthetic = X_syn,reference = X_train)
  rm=mean(res$dcr1)
  nndr <- as.numeric(quantile(res$nndr, probs = 0.05, na.rm = TRUE))
  
  # ============================================================
  # T8 — HOLDOUT PRIVACY
  # ============================================================
  all <- bind_rows(trn, val, syn_trn)
  all_hot <- dummyVars(" ~ .", data = all, fullRank = TRUE) |> predict(newdata = all) |> as.matrix()
  
  n_trn <- nrow(trn)
  n_val <- nrow(val)
  n_syn <- nrow(syn_trn)
  
  trn_hot <- all_hot[1:n_trn, ]
  val_hot <- all_hot[(n_trn+1):(n_trn+n_val), ]
  syn_hot <- all_hot[(n_trn+n_val+1):(n_trn+n_val+n_syn), ]
  
  d_trn <- get.knnx(trn_hot, syn_hot, k = 1)$nn.dist[,1]
  d_val <- get.knnx(val_hot, syn_hot, k = 1)$nn.dist[,1]
  
  shape <- mean(d_trn < d_val) + (n_trn/(n_trn+n_val))*mean(d_trn == d_val)
  
  # ============================================================
  # T9 — ATTRIBUTE INFERENCE
  # ============================================================
  
  air <- unlist(lapply(feature_associate, function(vars) {
    privacy_AIR(X_train, X_syn, target_col, vars)
  }))
  
  air_mean <- mean(air)
  air_min <- min(air)
  air_max <- max(air)
  
  if (run_linkability) {
    
    message("Running linkability attack...")
    
    # ------------------------------------------------------------
    # Prepare training data
    # ------------------------------------------------------------
    
    tr <- trn
    
    if (is_survival) {
      
      # Rename survival variables to the names expected by linkability_attack()
      if (statusvar != "Censored") {
        names(tr)[names(tr) == statusvar] <- "Censored"
      }
      
      if (timevar != "times") {
        names(tr)[names(tr) == timevar] <- "times"
      }
      
      # Exclude survival-specific variables
      vars <- names(
        tr[, !(names(tr) %in% c("Censored", "times"))]
      )
      
    } else {
      
      # All variables are available for the attack
      vars <- names(tr)
    }
    
    # ------------------------------------------------------------
    # Randomly split variables into A and B
    # ------------------------------------------------------------
    
    vr <- sample(vars)
    
    nA <- floor(length(vr) / 2)
    
    A <- vr[1:nA]
    B <- vr[(nA + 1):length(vr)]
    
    # ------------------------------------------------------------
    # Number of target records
    # ------------------------------------------------------------
    
    n_targets <- floor(nrow(tr) * 0.5)
    
    # ------------------------------------------------------------
    # Run linkability attack
    # ------------------------------------------------------------
    
    link_out <- linkability_attack(
      train_data    = tr,
      synthetic_data = syn_trn,
      A             = A,
      B             = B,
      control_data  = val,
      n_targets     = n_targets,
      k             = k,
      seed           = seed
    )
    
    # ------------------------------------------------------------
    # Extract linkability rates
    # ------------------------------------------------------------
    
    lr <- link_out$evaluation$linkability_rate
    n_succes_link <- link_out$evalution$n_successful_links
    
    
  }
  
  
  R=(lr[3]-lr[1])/(1-lr[1])
  delta_linkability_training_naive=lr[3]-lr[2]
  
  # ============================================================
  # OUTPUT
  # ============================================================
  
  return(c(
    T1=JSD, 
    T2=CFD,
    T3=auc_mean, 
    T4=SDiff_mean,
    T5=delta_tstr, 
    T6=delta_tsrtr, 
    T7=rm, 
    T8=shape, 
    T9=air_mean,
    T10=dcr,
    T11=nndr,
    T12=delta_linkability_training_naive,
    T13=R,
    linkability_rate_training=lr[3],
    linkability_rate_control=lr[1],
    linkability_rate_naive=lr[2],
    n_succes_link=n_succes_link,
    air_min=air_min,
    air_max=air_max,
    SDiff_q05=SDiff_q05,
    SDiff_median=SDiff_median,
    SDiff_q25=SDiff_q25,
    SDiff_q75=SDiff_q75,
    SDiff_q95=SDiff_q95,
    srd_n_below_dcr=srd_n_below_dcr,
    rrd_n_below_dcr=rrd_n_below_dcr
    
  ))
}



## ------------------------------------------------------------------------
## Core function: evaluate one (dataset, method) combination
## ------------------------------------------------------------------------

get_method_display_name <- function(method) {
  method_display_names[[method]] %||% toupper(method)
}

run_method_assessment <- function(
    dataset,
    method,
    X_real,
    X_train,
    X_test,
    trn,
    val,
    feature_associate,
    seeds,
    formula,
    target_col,
    modeling    = "logistic",
    timevar     = NULL,
    statusvar   = NULL,
    is_survival = FALSE,
    S           = 5,
    k           = 1,
    n_synth     = 50,
    n_workers   = 8,
    alpha=NULL,
    beta=NULL,
    n_rrd=NULL,
    run_linkability = TRUE
) {
  dataset_lc <- tolower(dataset)
  dataset_uc <- toupper(dataset)
  
  main_obj_name <- sprintf("%s_synthetic_list_data_%s", dataset_lc, method)
  trn_obj_name  <- sprintf("%s_synthetic_list_data_%s_trn", dataset_lc, method)
  
  main_path <- here(results_dir, paste0(main_obj_name, ".Rdata"))
  trn_path  <- here(results_dir, paste0(trn_obj_name, ".Rdata"))
  
  if (!file.exists(main_path)) stop("Missing file: ", main_path)
  if (!file.exists(trn_path))  stop("Missing file: ", trn_path)
  
  # Load into a private environment so successive methods never clash and
  # nothing leaks into .GlobalEnv.
  env <- new.env()
  load(main_path, envir = env)
  load(trn_path, envir = env)
  
  synth_list     <- env[[main_obj_name]]
  synth_list_trn <- env[[trn_obj_name]]
  
  method_display <- get_method_display_name(method)
  
  message(sprintf(
    "\n=== %s | %s : %d synthetic replicates ===",
    dataset_uc, method_display, length(synth_list)
  ))
  
  results_all <- lapply(seq_along(synth_list), function(i) {
    start_time <- Sys.time()
    
    cat(sprintf(
      "\n[%s - %s] Iteration %d of %d...\n",
      dataset_uc, method_display, i, length(synth_list)
    ))
    
    result <- Assessment_function_unified(
      X_real                = X_real,
      X_train            = X_train,
      X_syn             = synth_list[[i]],
      X_test            = X_test,
      S                 = S,
      formula           = formula,
      modeling          = modeling,
      timevar           = timevar,
      statusvar         = statusvar,
      target_col        = target_col,
      feature_associate = feature_associate,
      is_survival       = is_survival,
      trn               = trn,
      val               = val,
      syn_trn           = synth_list_trn[[i]],
      seed              = seeds[i],
      k                 = k,
      n_workers         = n_workers,
      alpha             = alpha,
      beta              = beta,
      n_rrd             = n_rrd,
      run_linkability   = run_linkability
    )
    
    elapsed <- difftime(Sys.time(), start_time, units = "secs")
    cat(sprintf(
      "[%s - %s] Iteration %d completed in %.2f seconds.\n",
      dataset_uc, method_display, i, as.numeric(elapsed)
    ))
    
    result
  })
  
  # ---- Per-replicate metrics -----------------------------------------
  fup_metrics <- results_all |> purrr::map_dfr(~ as.data.frame(t(.x)))
  
  # ---- compute target function ---------------------------------------
  metrics_loss <- compute_metrics_loss(
    results_all = fup_metrics,
    method_name = method_display,
    n_synth     = n_synth
  )
  
  fup_obj_name <- sprintf("%s_FUP_Metrics_%s", dataset_uc, method)
  assign(fup_obj_name, metrics_loss)
  save(list = fup_obj_name, file = here(results_dir, paste0(fup_obj_name, ".Rdata")))
  
  message(sprintf("Saved %s ", fup_obj_name))
  
  invisible(list(
    fup_metrics  = fup_metrics,
    metrics_loss = metrics_loss,
    results_all  = results_all
  ))
}


## ------------------------------------------------------------------------
## Driver: evaluate every method for one dataset
## ------------------------------------------------------------------------

run_dataset_evaluation <- function(
    dataset,
    methods,
    X_real,
    X_train,
    X_test,
    trn,
    val,
    feature_associate,
    seeds,
    formula,
    target_col,
    modeling    = "logistic",
    timevar     = NULL,
    statusvar   = NULL,
    is_survival = FALSE,
    S           = 5,
    k           = 1,
    n_synth     = 50,
    alpha       = NULL,
    beta        = NULL,
    n_rrd       = NULL,
    run_linkability = TRUE
) {
  results <- vector("list", length(methods))
  names(results) <- methods
  
  for (method in methods) {
    results[[method]] <- run_method_assessment(
      dataset           = dataset,
      method            = method,
      X_real            = X_real,
      X_train           = X_train,
      X_test            = X_test,
      trn               = trn,
      val               = val,
      feature_associate = feature_associate,
      seeds             = seeds,
      formula           = formula,
      target_col        = target_col,
      modeling          = modeling,
      timevar           = timevar,
      statusvar         = statusvar,
      is_survival       = is_survival,
      S                 = S,
      k                 = k,
      n_synth           = n_synth,
      alpha             = alpha,
      beta              = beta,
      n_rrd             = n_rrd,
      run_linkability   = run_linkability
    )
  }
  
  invisible(results)
}



run_all_evaluations <- function(
    dataset_configs,
    S       = 5,
    k       = 1,
    n_synth = 50,
    run_linkability = TRUE
) {
  all_results <- vector("list", length(dataset_configs))
  names(all_results) <- names(dataset_configs)
  
  for (dataset in names(dataset_configs)) {
    cfg <- dataset_configs[[dataset]]
    
    if (is.null(cfg$methods)) {
      stop(
        "No methods specified for dataset '", dataset,
        "'. Please add methods = c(...) to dataset_configs."
      )
    }
    
    all_results[[dataset]] <- run_dataset_evaluation(
      dataset           = dataset,
      methods           = cfg$methods,
      X_real            = cfg$X_real,
      X_train           = cfg$X_train,
      X_test            = cfg$X_test,
      trn               = cfg$trn,
      val               = cfg$val,
      feature_associate = cfg$feature_associate,
      seeds             = cfg$seeds,
      formula           = cfg$formula,
      target_col        = cfg$target_col,
      modeling          = cfg$modeling %||% "logistic",
      timevar           = cfg$timevar,
      statusvar         = cfg$statusvar,
      is_survival       = cfg$is_survival %||% FALSE,
      S                 = cfg$S %||% S,
      k                 = cfg$k %||% k,
      n_synth           = cfg$n_synth %||% n_synth,
      alpha             = cfg$alpha %||% NULL,
      beta              = cfg$beta %||% NULL,
      n_rrd             = cfg$n_rrd %||% NULL,
      run_linkability   = cfg$run_linkability %||% run_linkability
    )
  }
  
  invisible(all_results)
}


analyse_M <- function(
    data,
    method,
    objective_col = "L_total",
    datasets = c("PIMA", "AIDS", "REIN"),
    M_values = c(1, 5, 10, 25, 50),
    seed = 1234,
    minimize = TRUE,
    threshold_gain = 5
) {
  
  # =======================================================
  # 1. Input validation
  # =======================================================
  
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  
  if (!"method" %in% names(data)) {
    stop("The 'method' column is missing from the input data.")
  }
  
  if (!"dataset" %in% names(data)) {
    stop("The 'dataset' column is missing from the input data.")
  }
  
  if (!"run" %in% names(data)) {
    stop("The 'run' column is missing from the input data.")
  }
  
  if (!objective_col %in% names(data)) {
    stop(
      paste0(
        "The objective column '",
        objective_col,
        "' is missing from the input data."
      )
    )
  }
  
  if (length(M_values) == 0) {
    stop("'M_values' must contain at least one value.")
  }
  
  if (any(M_values < 1)) {
    stop("All values in 'M_values' must be >= 1.")
  }
  
  if (any(M_values != as.integer(M_values))) {
    stop("All values in 'M_values' must be integers.")
  }
  
  if (threshold_gain < 0) {
    stop("'threshold_gain' must be >= 0.")
  }
  
  
  # =======================================================
  # 2. Select requested method and datasets
  # =======================================================
  
  data_method <- data |>
    dplyr::filter(
      .data$method == method,
      .data$dataset %in% datasets
    )
  
  
  # =======================================================
  # 3. Define evaluation metrics
  # =======================================================
  
  metrics <- paste0("T", 1:13)
  
  available_runs <- sort(unique(data_method$run))
  
  n_runs <- length(available_runs)
  
  
  # =======================================================
  # 4. Check number of available runs
  # =======================================================
  
  if (n_runs == 0) {
    stop(
      "No runs are available for the selected method and datasets."
    )
  }
  
  if (max(M_values) > n_runs) {
    stop(
      paste0(
        "The maximum requested M (",
        max(M_values),
        ") is larger than the number of available runs (",
        n_runs,
        ")."
      )
    )
  }
  
  
  # =======================================================
  # 5. Initialize result containers
  # =======================================================
  
  results <- list()
  
  all_best <- list()
  
  all_intrinsic <- list()
  
  all_median_objective <- list()
  
  
  # =======================================================
  # 6. Generate a common random ordering of runs
  # =======================================================
  
  set.seed(seed)
  
  run_order <- sample(
    available_runs,
    size = n_runs,
    replace = FALSE
  )
  
  
  # =======================================================
  # 7. Evaluate each value of M
  # =======================================================
  
  for (M in M_values) {
    
    # -------------------------------------------------------
    # Nested candidate set
    # -------------------------------------------------------
    
    sampled_runs <- run_order[1:M]
    
    tmp <- data_method |>
      dplyr::filter(
        .data$run %in% sampled_runs
      )
    
    
    # -------------------------------------------------------
    # Descriptive statistics for T1-T13
    # -------------------------------------------------------
    
    intrinsic <- tmp |>
      dplyr::select(
        dataset,
        method,
        run,
        dplyr::all_of(metrics)
      ) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(metrics),
        names_to = "metric",
        values_to = "value"
      ) |>
      dplyr::group_by(
        method,
        dataset,
        metric
      ) |>
      dplyr::summarise(
        
        mean = mean(
          value,
          na.rm = TRUE
        ),
        
        sd = sd(
          value,
          na.rm = TRUE
        ),
        
        median = median(
          value,
          na.rm = TRUE
        ),
        
        Q25 = quantile(
          value,
          0.25,
          na.rm = TRUE,
          names = FALSE
        ),
        
        Q75 = quantile(
          value,
          0.75,
          na.rm = TRUE,
          names = FALSE
        ),
        
        min = min(
          value,
          na.rm = TRUE
        ),
        
        max = max(
          value,
          na.rm = TRUE
        ),
        
        .groups = "drop"
      ) |>
      dplyr::mutate(
        M = M
      )
    
    
    # -------------------------------------------------------
    # Objective summary across the M runs
    # -------------------------------------------------------
    
    objective_summary <- tmp |>
      dplyr::group_by(
        dataset
      ) |>
      dplyr::summarise(
        
        mean_objective = mean(
          .data[[objective_col]],
          na.rm = TRUE
        ),
        
        sd_objective = sd(
          .data[[objective_col]],
          na.rm = TRUE
        ),
        
        median_objective = median(
          .data[[objective_col]],
          na.rm = TRUE
        ),
        
        min_objective = min(
          .data[[objective_col]],
          na.rm = TRUE
        ),
        
        max_objective = max(
          .data[[objective_col]],
          na.rm = TRUE
        ),
        
        n = sum(
          !is.na(.data[[objective_col]])
        ),
        
        .groups = "drop"
      ) |>
      dplyr::mutate(
        method = method,
        M = M
      ) |>
      dplyr::select(
        dataset,
        method,
        M,
        n,
        mean_objective,
        sd_objective,
        median_objective,
        min_objective,
        max_objective
      )
    
    
    # -------------------------------------------------------
    # Identify the best run according to objective
    # -------------------------------------------------------
    
    if (minimize) {
      
      best_rows <- tmp |>
        dplyr::group_by(
          dataset
        ) |>
        dplyr::slice_min(
          order_by = .data[[objective_col]],
          n = 1,
          with_ties = FALSE
        ) |>
        dplyr::ungroup()
      
    } else {
      
      best_rows <- tmp |>
        dplyr::group_by(
          dataset
        ) |>
        dplyr::slice_max(
          order_by = .data[[objective_col]],
          n = 1,
          with_ties = FALSE
        ) |>
        dplyr::ungroup()
    }
    
    
    # -------------------------------------------------------
    # Retain the complete best run
    # -------------------------------------------------------
    
    best_rows <- best_rows |>
      dplyr::select(
        dataset,
        method,
        run,
        dplyr::all_of(metrics),
        dplyr::all_of(objective_col)
      ) |>
      dplyr::mutate(
        M = M
      )
    
    
    # -------------------------------------------------------
    # Store results for current M
    # -------------------------------------------------------
    
    all_best[[paste0("M_", M)]] <- best_rows
    
    all_intrinsic[[paste0("M_", M)]] <- intrinsic
    
    all_median_objective[[paste0("M_", M)]] <- objective_summary
    
    results[[paste0("M_", M)]] <- list(
      M = M,
      sampled_runs = sampled_runs,
      intrinsic = intrinsic,
      objective_summary = objective_summary,
      best_rows = best_rows
    )
  }
  
  
  # =======================================================
  # 8. Combine best runs across all M
  # =======================================================
  
  best_all <- dplyr::bind_rows(
    all_best
  ) |>
    dplyr::arrange(
      dataset,
      M
    )
  
  
  # =======================================================
  # 9. Combine median objective across all M
  # =======================================================
  
  median_objective_all <- dplyr::bind_rows(
    all_median_objective
  ) |>
    dplyr::arrange(
      dataset,
      M
    )
  
  
  # =======================================================
  # 10. Calculate incremental gain based on median objective
  # =======================================================
  
  gain <- median_objective_all |>
    dplyr::select(
      dataset,
      M,
      median_objective
    ) |>
    dplyr::group_by(
      dataset
    ) |>
    dplyr::arrange(
      M,
      .by_group = TRUE
    ) |>
    dplyr::mutate(
      
      # ---------------------------------------------------
      # Median objective at previous M
      # ---------------------------------------------------
      
      previous_median_objective = dplyr::lag(
        median_objective
      ),
      
      # ---------------------------------------------------
      # Absolute improvement
      # ---------------------------------------------------
      
      gain_absolute = dplyr::case_when(
        
        is.na(previous_median_objective) ~ NA_real_,
        
        minimize ~
          previous_median_objective -
          median_objective,
        
        TRUE ~
          median_objective -
          previous_median_objective
      ),
      
      # ---------------------------------------------------
      # Relative improvement (%)
      # ---------------------------------------------------
      
      gain_relative = dplyr::case_when(
        
        is.na(previous_median_objective) ~ NA_real_,
        
        previous_median_objective == 0 ~ NA_real_,
        
        TRUE ~
          100 *
          abs(gain_absolute) /
          abs(previous_median_objective)
      ),
      
      # ---------------------------------------------------
      # Meaningful improvement
      # ---------------------------------------------------
      
      meaningful_gain = dplyr::case_when(
        
        is.na(gain_relative) ~ NA,
        
        gain_relative > threshold_gain ~ TRUE,
        
        TRUE ~ FALSE
      ),
      
      gain_threshold = threshold_gain
    ) |>
    dplyr::ungroup()
  
  
  # =======================================================
  # Identify first M with non-meaningful gain
  # =======================================================
  
  first_negligible_M <- gain |>
    dplyr::filter(
      meaningful_gain == FALSE
    ) |>
    dplyr::group_by(
      dataset
    ) |>
    dplyr::slice_min(
      order_by = M,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      dataset,
      M,
      gain_relative,
      gain_threshold
    ) |>
    dplyr::rename(
      first_negligible_M = M
    )
  
  
  # =======================================================
  # Return all results
  # =======================================================
  
  return(
    list(
      
      # ---------------------------------------------------
      # Results for each M
      # ---------------------------------------------------
      
      results = results,
      
      # ---------------------------------------------------
      # Best run for each dataset and M
      # ---------------------------------------------------
      
      best_all = best_all,
      
      # ---------------------------------------------------
      # Median objective for each dataset and M
      # ---------------------------------------------------
      
      median_objective = median_objective_all,
      
      # ---------------------------------------------------
      # Incremental gains based on median objective
      # ---------------------------------------------------
      
      gain = gain,
      
      # ---------------------------------------------------
      # First M where gain <= threshold
      # ---------------------------------------------------
      
      first_negligible_M = first_negligible_M,
      
      # ---------------------------------------------------
      # Common random permutation of runs
      # ---------------------------------------------------
      
      run_order = run_order,
      
      # ---------------------------------------------------
      # Analysis parameters
      # ---------------------------------------------------
      
      parameters = list(
        method = method,
        objective = objective_col,
        datasets = datasets,
        M_values = M_values,
        seed = seed,
        minimize = minimize,
        threshold_gain = threshold_gain,
        n_available_runs = n_runs
      )
    )
  )
}




# ============================================================
# Bootstrap confidence intervals and comparison of methods
#
# Input data format:
#   T1, T2, ..., T13, dataset, method, run
#
# Target-based metrics:
#   median(|value - target|)
#
# Non-target metrics:
#   median(value)
#
# Bootstrap comparison is based on the difference between
# the two method-specific statistics.
# ============================================================


bootstrap_compare_methods <- function(
    data_chimera,
    data_comparator,
    comparator_name,
    metric_direction,
    datasets = c("PIMA", "AIDS", "REIN"),
    metrics = paste0("T", 1:13),
    B = 10000,
    confidence_level = 0.95,
    seed = 1234
) {
  
  # ============================================================
  # 1. Input checks
  # ============================================================
  
  required_columns <- c("dataset", "run")
  
  if (!all(required_columns %in% names(data_chimera))) {
    stop(
      "data_chimera must contain at least: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  if (!all(required_columns %in% names(data_comparator))) {
    stop(
      "data_comparator must contain at least: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  missing_chimera_metrics <- setdiff(
    metrics,
    names(data_chimera)
  )
  
  missing_comparator_metrics <- setdiff(
    metrics,
    names(data_comparator)
  )
  
  if (length(missing_chimera_metrics) > 0) {
    stop(
      "Missing metrics in data_chimera: ",
      paste(missing_chimera_metrics, collapse = ", ")
    )
  }
  
  if (length(missing_comparator_metrics) > 0) {
    stop(
      "Missing metrics in data_comparator: ",
      paste(missing_comparator_metrics, collapse = ", ")
    )
  }
  
  # ============================================================
  # 2. Check metric_direction
  # ============================================================
  
  required_direction_columns <- c(
    "metric",
    "direction",
    "target"
  )
  
  if (!all(required_direction_columns %in%
           names(metric_direction))) {
    
    stop(
      "metric_direction must contain the columns: ",
      paste(required_direction_columns, collapse = ", ")
    )
  }
  
  # Keep only requested metrics
  metric_direction <- metric_direction[
    metric_direction$metric %in% metrics,
    ,
    drop = FALSE
  ]
  
  # Check that all requested metrics are present
  missing_direction <- setdiff(
    metrics,
    metric_direction$metric
  )
  
  if (length(missing_direction) > 0) {
    stop(
      "Missing metric definitions for: ",
      paste(missing_direction, collapse = ", ")
    )
  }
  
  # Reorder according to metrics
  metric_direction <- metric_direction[
    match(metrics, metric_direction$metric),
    ,
    drop = FALSE
  ]
  
  # ============================================================
  # 3. Validate directions
  # ============================================================
  
  allowed_directions <- c(
    "lower",
    "higher",
    "target"
  )
  
  invalid_directions <- setdiff(
    unique(metric_direction$direction),
    allowed_directions
  )
  
  if (length(invalid_directions) > 0) {
    stop(
      "Invalid direction(s): ",
      paste(invalid_directions, collapse = ", "),
      ". Allowed values are: lower, higher, target."
    )
  }
  
  # ============================================================
  # 4. Check target values
  #
  # Every metric with direction = "target" MUST have
  # a non-missing target.
  #
  # Metrics with lower/higher direction do not need a target.
  # ============================================================
  
  target_metrics <- metric_direction$metric[
    metric_direction$direction == "target"
  ]
  
  missing_target <- metric_direction$metric[
    metric_direction$direction == "target" &
      is.na(metric_direction$target)
  ]
  
  if (length(missing_target) > 0) {
    stop(
      "Target value missing for target metric(s): ",
      paste(missing_target, collapse = ", ")
    )
  }
  
  # ============================================================
  # 5. Convert WIDE -> LONG
  # ============================================================
  
  chimera_long <- do.call(
    rbind,
    lapply(metrics, function(m) {
      
      data.frame(
        Dataset = data_chimera$dataset,
        Run = data_chimera$run,
        Metric = m,
        Value = data_chimera[[m]],
        stringsAsFactors = FALSE
      )
      
    })
  )
  
  comparator_long <- do.call(
    rbind,
    lapply(metrics, function(m) {
      
      data.frame(
        Dataset = data_comparator$dataset,
        Run = data_comparator$run,
        Metric = m,
        Value = data_comparator[[m]],
        stringsAsFactors = FALSE
      )
      
    })
  )
  
  # ============================================================
  # 6. Restrict to requested datasets and metrics
  # ============================================================
  
  chimera_long <- chimera_long[
    chimera_long$Dataset %in% datasets &
      chimera_long$Metric %in% metrics,
    ,
    drop = FALSE
  ]
  
  comparator_long <- comparator_long[
    comparator_long$Dataset %in% datasets &
      comparator_long$Metric %in% metrics,
    ,
    drop = FALSE
  ]
  
  # ============================================================
  # 7. Remove non-finite values
  # ============================================================
  
  chimera_long <- chimera_long[
    is.finite(chimera_long$Value),
    ,
    drop = FALSE
  ]
  
  comparator_long <- comparator_long[
    is.finite(comparator_long$Value),
    ,
    drop = FALSE
  ]
  
  # ============================================================
  # 8. Bootstrap median for lower/higher metrics
  # ============================================================
  
  bootstrap_median <- function(
    x,
    B,
    confidence_level
  ) {
    
    x <- x[is.finite(x)]
    
    if (length(x) < 2) {
      
      return(
        list(
          estimate = NA_real_,
          boot = numeric(0),
          lower = NA_real_,
          upper = NA_real_,
          n = length(x)
        )
      )
    }
    
    # Observed statistic
    estimate <- median(x)
    
    # Bootstrap distribution
    boot <- replicate(
      B,
      median(
        sample(
          x,
          size = length(x),
          replace = TRUE
        )
      )
    )
    
    alpha <- 1 - confidence_level
    
    ci <- quantile(
      boot,
      probs = c(
        alpha / 2,
        1 - alpha / 2
      ),
      na.rm = TRUE,
      names = FALSE
    )
    
    list(
      estimate = estimate,
      boot = boot,
      lower = ci[1],
      upper = ci[2],
      n = length(x)
    )
  }
  
  
  # ============================================================
  # 9. Bootstrap median distance to target
  #
  # USED FOR ALL target METRICS
  #
  # Statistic:
  #
  #   median(|x_i - target|)
  #
  # This is applied to:
  #   T3, T5, T6, T7, T8, T10, T11, T12, T13
  # according to the supplied metric_direction table.
  # ============================================================
  
  bootstrap_median_distance <- function(
    x,
    target,
    B,
    confidence_level
  ) {
    
    x <- x[is.finite(x)]
    
    if (length(x) < 2) {
      
      return(
        list(
          estimate = NA_real_,
          boot = numeric(0),
          lower = NA_real_,
          upper = NA_real_,
          n = length(x)
        )
      )
    }
    
    # ----------------------------------------------------------
    # Distance to target for EACH generation
    # ----------------------------------------------------------
    
    distances <- abs(
      x - target
    )
    
    # ----------------------------------------------------------
    # Point estimate:
    #
    # median of generation-specific distances
    # ----------------------------------------------------------
    
    estimate <- median(
      distances
    )
    
    # ----------------------------------------------------------
    # Bootstrap distribution
    # ----------------------------------------------------------
    
    boot <- replicate(
      B,
      {
        
        x_boot <- sample(
          x,
          size = length(x),
          replace = TRUE
        )
        
        distances_boot <- abs(
          x_boot - target
        )
        
        median(
          distances_boot
        )
      }
    )
    
    # ----------------------------------------------------------
    # Bootstrap percentile CI
    # ----------------------------------------------------------
    
    alpha <- 1 - confidence_level
    
    ci <- quantile(
      boot,
      probs = c(
        alpha / 2,
        1 - alpha / 2
      ),
      na.rm = TRUE,
      names = FALSE
    )
    
    list(
      estimate = estimate,
      boot = boot,
      lower = ci[1],
      upper = ci[2],
      n = length(x)
    )
  }
  
  
  # ============================================================
  # 10. Set random seed
  # ============================================================
  
  set.seed(seed)
  
  
  # ============================================================
  # 11. Initialize results
  # ============================================================
  
  results <- vector(
    "list",
    length(datasets) * length(metrics)
  )
  
  result_index <- 1
  
  
  # ============================================================
  # 12. Loop over datasets
  # ============================================================
  
  for (dataset in datasets) {
    
    # ----------------------------------------------------------
    # Loop over metrics
    # ----------------------------------------------------------
    
    for (metric in metrics) {
      
      # --------------------------------------------------------
      # Extract values for CHIMERA
      # --------------------------------------------------------
      
      x <- chimera_long$Value[
        chimera_long$Dataset == dataset &
          chimera_long$Metric == metric
      ]
      
      # --------------------------------------------------------
      # Extract values for comparator
      # --------------------------------------------------------
      
      y <- comparator_long$Value[
        comparator_long$Dataset == dataset &
          comparator_long$Metric == metric
      ]
      
      # Remove non-finite values
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
      
      # --------------------------------------------------------
      # Retrieve direction and target
      # --------------------------------------------------------
      
      metric_info <- metric_direction[
        metric_direction$metric == metric,
        ,
        drop = FALSE
      ]
      
      direction <- metric_info$direction[1]
      
      target <- metric_info$target[1]
      
      # --------------------------------------------------------
      # Check sample size
      # --------------------------------------------------------
      
      if (length(x) < 2 || length(y) < 2) {
        
        results[[result_index]] <- data.frame(
          
          Dataset = dataset,
          
          Metric = metric,
          
          Comparator = comparator_name,
          
          CHIMERA_Median = if (
            length(x) > 0
          ) {
            median(x)
          } else {
            NA_real_
          },
          
          Comparator_Median = if (
            length(y) > 0
          ) {
            median(y)
          } else {
            NA_real_
          },
          
          CHIMERA_Estimate = NA_real_,
          
          Comparator_Estimate = NA_real_,
          
          Difference = NA_real_,
          
          CI_Lower = NA_real_,
          
          CI_Upper = NA_real_,
          
          Significant = NA,
          
          Better_Method = NA_character_
          
          
        )
        
        result_index <- result_index + 1
        
        next
      }
      
      
      # ========================================================
      # 13. TARGET METRIC
      #
      # For ALL target metrics:
      #
      #   statistic = median(|x - target|)
      #
      # Smaller distance = better.
      # ========================================================
      
      if (direction == "target") {
        
        boot_chimera <- bootstrap_median_distance(
          x = x,
          target = target,
          B = B,
          confidence_level = confidence_level
        )
        
        boot_comparator <- bootstrap_median_distance(
          x = y,
          target = target,
          B = B,
          confidence_level = confidence_level
        )
        
        chimera_estimate <-
          boot_chimera$estimate
        
        comparator_estimate <-
          boot_comparator$estimate
        
        # ------------------------------------------------------
        # Difference
        #
        # CHIMERA distance - comparator distance
        #
        # Negative -> CHIMERA closer to target
        # Positive -> comparator closer to target
        # ------------------------------------------------------
        
        difference <-
          chimera_estimate -
          comparator_estimate
        
        # ------------------------------------------------------
        # Bootstrap difference
        # ------------------------------------------------------
        
        boot_difference <-
          boot_chimera$boot -
          boot_comparator$boot
        
        alpha <- 1 - confidence_level
        
        ci <- quantile(
          boot_difference,
          probs = c(
            alpha / 2,
            1 - alpha / 2
          ),
          na.rm = TRUE,
          names = FALSE
        )
        
        ci_lower <- ci[1]
        
        ci_upper <- ci[2]
        
        # ------------------------------------------------------
        # Significance
        #
        # CI excludes zero
        # ------------------------------------------------------
        
        significant <- (
          ci_lower > 0 ||
            ci_upper < 0
        )
        
        # ------------------------------------------------------
        # Better method
        # Smaller distance = better
        # ------------------------------------------------------
        
        if (significant) {
          
          if (difference < 0) {
            
            better_method <- "CHIMERA"
            
          } else {
            
            better_method <- comparator_name
          }
          
        } else {
          
          better_method <-
            "No significant difference"
        }
      }
      
      
      # ========================================================
      # 14. LOWER IS BETTER
      # ========================================================
      
      else if (direction == "lower") {
        
        boot_chimera <- bootstrap_median(
          x = x,
          B = B,
          confidence_level = confidence_level
        )
        
        boot_comparator <- bootstrap_median(
          x = y,
          B = B,
          confidence_level = confidence_level
        )
        
        chimera_estimate <-
          boot_chimera$estimate
        
        comparator_estimate <-
          boot_comparator$estimate
        
        # CHIMERA - comparator
        difference <-
          chimera_estimate -
          comparator_estimate
        
        boot_difference <-
          boot_chimera$boot -
          boot_comparator$boot
        
        alpha <- 1 - confidence_level
        
        ci <- quantile(
          boot_difference,
          probs = c(
            alpha / 2,
            1 - alpha / 2
          ),
          na.rm = TRUE,
          names = FALSE
        )
        
        ci_lower <- ci[1]
        
        ci_upper <- ci[2]
        
        significant <- (
          ci_lower > 0 ||
            ci_upper < 0
        )
        
        # Smaller value = better
        if (significant) {
          
          if (difference < 0) {
            
            better_method <- "CHIMERA"
            
          } else {
            
            better_method <- comparator_name
          }
          
        } else {
          
          better_method <-
            "No significant difference"
        }
      }
      
      
      # ========================================================
      # 15. HIGHER IS BETTER
      # ========================================================
      
      else if (direction == "higher") {
        
        boot_chimera <- bootstrap_median(
          x = x,
          B = B,
          confidence_level = confidence_level
        )
        
        boot_comparator <- bootstrap_median(
          x = y,
          B = B,
          confidence_level = confidence_level
        )
        
        chimera_estimate <-
          boot_chimera$estimate
        
        comparator_estimate <-
          boot_comparator$estimate
        
        # CHIMERA - comparator
        difference <-
          chimera_estimate -
          comparator_estimate
        
        boot_difference <-
          boot_chimera$boot -
          boot_comparator$boot
        
        alpha <- 1 - confidence_level
        
        ci <- quantile(
          boot_difference,
          probs = c(
            alpha / 2,
            1 - alpha / 2
          ),
          na.rm = TRUE,
          names = FALSE
        )
        
        ci_lower <- ci[1]
        
        ci_upper <- ci[2]
        
        significant <- (
          ci_lower > 0 ||
            ci_upper < 0
        )
        
        # Larger value = better
        if (significant) {
          
          if (difference > 0) {
            
            better_method <- "CHIMERA"
            
          } else {
            
            better_method <- comparator_name
          }
          
        } else {
          
          better_method <-
            "No significant difference"
        }
      }
      
      
      # ========================================================
      # 16. Store result
      # ========================================================
      
      results[[result_index]] <- data.frame(
        
        Dataset = dataset,
        
        Metric = metric,
        
        Comparator = comparator_name,
        
        # ------------------------------------------------------
        # Raw median of the original metric
        # ------------------------------------------------------
        
        CHIMERA_Median = median(x),
        
        Comparator_Median = median(y),
        
        # ------------------------------------------------------
        # Statistic actually used for comparison
        #
        # lower/higher:
        #   median(metric)
        #
        # target:
        #   median(abs(metric - target))
        # ------------------------------------------------------
        
        CHIMERA_Estimate =
          chimera_estimate,
        
        Comparator_Estimate =
          comparator_estimate,
        
        # ------------------------------------------------------
        # Difference
        # ------------------------------------------------------
        
        Difference =
          difference,
        
        # ------------------------------------------------------
        # Bootstrap CI for difference
        # ------------------------------------------------------
        
        CI_Lower =
          ci_lower,
        
        CI_Upper =
          ci_upper,
        
        # ------------------------------------------------------
        # Significance
        # ------------------------------------------------------
        
        Significant =
          significant,
        
        # ------------------------------------------------------
        # Better method
        # ------------------------------------------------------
        
        Better_Method =
          better_method,
  
        stringsAsFactors = FALSE
      )
      
      result_index <- result_index + 1
    }
  }
  
  
  # ============================================================
  # 17. Combine all results
  # ============================================================
  
  results <- do.call(
    rbind,
    results
  )
  
  rownames(results) <- NULL
  
  
  # ============================================================
  # 19. Return
  # ============================================================
  
  return(results)
}

# Analytic et predictive utility function
fit_and_summarise_classic <- function(
    df,
    label,
    run = NULL,
    reference_table = NULL
) {
  
  # ---------------------------------------------------------------------------
  # Logistic regression
  # ---------------------------------------------------------------------------
  
  model <- glm(
    formula_reg.log,
    data = df,
    family = binomial
  )
  
  
  # ---------------------------------------------------------------------------
  # Model coefficients
  # ---------------------------------------------------------------------------
  
  coef_tab <- summary(model)$coefficients
  
  
  # ---------------------------------------------------------------------------
  # Regression summary
  # ---------------------------------------------------------------------------
  
  tab_reg <- gtsummary::tbl_regression(
    model,
    exponentiate = TRUE
  )
  
  body <- tab_reg$table_body |>
    dplyr::filter(
      !is.na(estimate)
    )
  
  
  # ---------------------------------------------------------------------------
  # Predictions on the same dataset used to fit the model
  # ---------------------------------------------------------------------------
  
  preds <- predict(
    model,
    newdata = df,
    type = "response"
  )
  
  
  # ---------------------------------------------------------------------------
  # ROC curve and AUC
  # ---------------------------------------------------------------------------
  
  roc_obj <- pROC::roc(
    response = df$Diabetes_diagnosis,
    predictor = preds,
    quiet = TRUE
  )
  
  
  # ---------------------------------------------------------------------------
  # Odds ratios, 95% CI and standard errors
  # ---------------------------------------------------------------------------
  
  tbl <- data.frame(
    Covariates = Variables_num,
    OR = body$estimate,
    CI95_lower = body$conf.low,
    CI95_upper = body$conf.high,
    Std_Error = coef_tab[-1, "Std. Error"]
  )
  
  
  # ---------------------------------------------------------------------------
  # SMD of log odds ratios
  #
  # SMD is calculated ONLY for synthetic datasets.
  # The original dataset is the reference and therefore does not receive
  # an SMD value.
  # ---------------------------------------------------------------------------
  
  if (!is.null(reference_table)) {
    
    tbl <- tbl |>
      dplyr::left_join(
        reference_table |>
          dplyr::select(
            Covariates,
            OR_original,
            Std_Error_original
          ),
        by = "Covariates"
      ) |>
      dplyr::mutate(
        SMD_log = abs(
          (
            log(OR) -
              log(OR_original)
          ) /
            sqrt(
              Std_Error^2 +
                Std_Error_original^2
            )
        )
      ) |>
      dplyr::select(
        Covariates,
        OR,
        CI95_lower,
        CI95_upper,
        Std_Error,
        SMD_log
      )
  }
  
  
  # ---------------------------------------------------------------------------
  # Round numerical results
  # ---------------------------------------------------------------------------
  
  tbl <- tbl |>
    dplyr::mutate(
      dplyr::across(
        where(is.numeric),
        ~ round(.x, 3)
      )
    )
  
  
  # ---------------------------------------------------------------------------
  # Add method name to column names
  # ---------------------------------------------------------------------------
  
  names(tbl)[-1] <- paste0(
    names(tbl)[-1],
    "_",
    label
  )
  
  
  # ---------------------------------------------------------------------------
  # Return results
  # ---------------------------------------------------------------------------
  
  list(
    roc = roc_obj,
    table = tbl
  )
}


# ===============================================================================================
# FUNCTION FOR COX MODEL + HR + TIME-DEPENDENT AUC
# ===============================================================================================

fit_and_summarise_survival <- function(
    df,
    label,
    run = NULL,
    reference_table = NULL,
    test_data,
    time_points = seq(0, 1100, by = 100),
    n_boot = 200,
    conf_level = 0.95
) {
  
  
  # =============================================================================================
  # 10.1 Cox proportional hazards model
  # =============================================================================================
  
  model <- survival::coxph(
    survival::Surv(
      times,
      Censored
    ) ~ .,
    data = df
  )
  
  
  # =============================================================================================
  # 10.2 Model summary
  # =============================================================================================
  
  model_summary <- summary(model)
  
  coef_tab <- model_summary$coefficients
  
  conf_int <- model_summary$conf.int
  
  
  # =============================================================================================
  # 10.3 HR, 95% CI and standard error
  # =============================================================================================
  
  tbl <- data.frame(
    Covariates = rownames(coef_tab),
    
    HR = conf_int[
      ,
      "exp(coef)"
    ],
    
    CI95_lower = conf_int[
      ,
      "lower .95"
    ],
    
    CI95_upper = conf_int[
      ,
      "upper .95"
    ],
    
    SE = coef_tab[
      ,
      "se(coef)"
    ],
    
    stringsAsFactors = FALSE
  )
  
  
  # =============================================================================================
  # 10.4 SMD of log-HR
  #
  # The original dataset is the reference.
  # SMD is calculated only for synthetic datasets.
  # =============================================================================================
  
  if (!is.null(reference_table)) {
    
    tbl <- tbl |>
      dplyr::left_join(
        
        reference_table |>
          dplyr::select(
            Covariates,
            HR_REAL,
            SE_REAL
          ),
        
        by = "Covariates"
      ) |>
      
      dplyr::mutate(
        
        SMD_log = abs(
          
          (
            log(HR) -
              log(HR_REAL)
          ) /
            
            sqrt(
              SE^2 +
                SE_REAL^2
            )
        )
      ) |>
      
      dplyr::select(
        Covariates,
        HR,
        CI95_lower,
        CI95_upper,
        SE,
        SMD_log
      )
  }
  
  
  # =============================================================================================
  # 10.5 Round regression results
  # =============================================================================================
  
  tbl <- tbl |>
    dplyr::mutate(
      dplyr::across(
        where(is.numeric),
        ~ round(
          .x,
          3
        )
      )
    )
  
  
  # =============================================================================================
  # 10.6 Add method name to table columns
  # =============================================================================================
  
  names(tbl)[-1] <- paste0(
    names(tbl)[-1],
    "_",
    toupper(label)
  )
  
  
  # =============================================================================================
  # 10.7 Time-dependent AUC
  #
  # Model is trained on df.
  # Predictions are evaluated on the independent AIDS_test dataset.
  # =============================================================================================
  
  auc_td <- predictive_utility_surv_curve(
    
    train = df,
    
    test = test_data,
    
    time_points = time_points,
    
    n_boot = n_boot,
    
    conf_level = conf_level
  ) |>
    
    dplyr::mutate(
      source = toupper(label),
      run = run
    )
  
  
  # =============================================================================================
  # 10.8 Return all results
  # =============================================================================================
  
  list(
    table = tbl,
    
    auc_time_dependent = auc_td
  )
}




#################################################################################################"
# Ablation analysis functions
#################################################################################################
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
  # Survival time reconstruction
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
    nb_imputations = nb_imputations,
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


preprocessing_auc <- function(
    df,
    nb_imputations = 1,
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
  # --- Initial MICE setup
  ini <- mice(df, maxit = 0)
  meth <- ini$method
  pred <- ini$predictorMatrix
  meth[names(df)[sapply(df, is.numeric)]] <- "pmm"
  
  # --- Impute missing data
  df <- complete(mice(df, m = nb_imputations, maxit = 5, seed = 123,method = meth, predictorMatrix = pred, printFlag = TRUE))
  
  
  # ---------------------------------------------------------------------------------------------
  # Return imputed datasets
  # ---------------------------------------------------------------------------------------------
  
  return(df)
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


compute_auc_ci <- function(df_list,df_complete,outcome,drop_vars = c("Censorship", "Follow_up_time"),method_name = "METHOD",name_variable = NULL) {
  
  
  res_list <- lapply(seq_along(df_list), function(i) {
    
    # Fit a logistic regression model on the i-th dataset
    reg.log <- glm(
      formula = as.formula(
        paste(outcome, " ~ ", paste(name_variable, collapse = " + "))
      ),
      data = df_list[[i]] |> dplyr::select(-all_of(drop_vars)),
      family = "binomial"
    )
    
    # Generate predicted probabilities on the test dataset
    preds <- predict(reg.log, newdata = df_test, type = "response")
    
    # Compute ROC curve and Area Under the Curve (AUC)
    # using the outcome from the complete (reference) dataset
    roc_obj <- roc(df_complete[[outcome]], preds, quiet = TRUE)
    
    # Extract AUC value
    auc_val <- auc(roc_obj)
    
    # Compute confidence interval for the AUC
    ci_vals <- ci.auc(roc_obj)
    
    # Store results in a data frame
    data.frame(
      imputation = i,          # Index of the dataset
      source = method_name,    # Name of the method used
      auc = as.numeric(auc_val),
      ci_lower = ci_vals[1],   # Lower bound of AUC CI
      ci_upper = ci_vals[3]    # Upper bound of AUC CI
    )
  })
  
  # Combine results from all datasets into a single data frame
  bind_rows(res_list)
}


# -------------------------------------------------------------------------
# Function to fit logistic models across 50 synthetic datasets,
# pool predictions by averaging them, and compute the ROC curve
# -------------------------------------------------------------------------

compute_pooled_roc <- function(
    imputed_data,
    predictor_names,
    test_data,
    outcome = "death",
    n_runs = 50
) {
  
  # Check the number of imputations
  if (length(imputed_data) < n_runs) {
    stop(
      "Not enough imputations: ",
      length(imputed_data),
      " available, but ",
      n_runs,
      " required."
    )
  }
  
  # Check that predictors are provided
  if (length(predictor_names) == 0) {
    stop("No predictor variables provided.")
  }
  
  # Remove duplicated predictors
  predictor_names <- unique(predictor_names)
  
  # Clean variable names using the first imputed dataset
  predictor_names <- clean_variable_names(
    predictor_names,
    imputed_data[[1]]
  )
  
  # Check that all predictors exist in the test dataset
  missing_predictors <- setdiff(
    predictor_names,
    names(test_data)
  )
  
  if (length(missing_predictors) > 0) {
    stop(
      "Predictors missing from test_data: ",
      paste(missing_predictors, collapse = ", ")
    )
  }
  
  # Build the logistic regression formula
  model_formula <- reformulate(
    termlabels = predictor_names,
    response = outcome
  )
  
  # Initialize predictions
  pooled_predictions <- numeric(nrow(test_data))
  
  # Fit one model per imputation
  for (i in seq_len(n_runs)) {
    
    logistic_model <- glm(
      formula = model_formula,
      data = imputed_data[[i]],
      family = binomial()
    )
    
    # Predict mortality probabilities
    predictions <- predict(
      logistic_model,
      newdata = test_data,
      type = "response"
    )
    
    # Accumulate predictions
    pooled_predictions <- pooled_predictions + predictions
  }
  
  # Average predictions over the 50 imputations
  pooled_predictions <- pooled_predictions / n_runs
  
  # Compute the ROC curve
  roc_curve <- pROC::roc(
    response = test_data[[outcome]],
    predictor = pooled_predictions,
    quiet = TRUE
  )
  
  return(roc_curve)
}


# -------------------------------------------------------------------------
# Compute one ROC curve for each synthetic dataset
# -------------------------------------------------------------------------

compute_method_rocs <- function(
    imputed_datasets,
    variable_list,
    test_data,
    outcome = "death",
    n_runs = 50
) {
  
  # Check that each synthetic dataset has its own variable list
  if (length(imputed_datasets) != length(variable_list)) {
    stop(
      "The number of synthetic datasets does not match ",
      "the number of variable lists."
    )
  }
  
  # Compute one ROC curve per synthetic dataset
  roc_list <- lapply(seq_along(imputed_datasets), function(i) {
    
    message("Processing synthetic dataset ", i)
    
    compute_pooled_roc(
      imputed_data = imputed_datasets[[i]],
      predictor_names = variable_list[[i]],
      test_data = test_data,
      outcome = outcome,
      n_runs = n_runs
    )
  })
  
  # Name the ROC curves
  names(roc_list) <- paste0(
    "synthetic_dataset_",
    seq_along(roc_list)
  )
  
  return(roc_list)
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
    seed = 123,
    threshold=70
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
      filter(nb_significative_coef > threshold)
    
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
        Censorship   = Censored,
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
    seed = 123,
    threshold = 70
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
      filter(nb_significative_coef > threshold)
    
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






