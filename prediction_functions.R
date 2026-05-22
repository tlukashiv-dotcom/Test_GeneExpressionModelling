#' @title Prediction Functions for Dynamic Gene Expression Modelling
#' @description Core prediction wrapper for forward validation, backward validation,
#' and future forecasting.
#' @author Taras Lukashiv and collaborators
#' @date 2026-05-20


# =============================================================================
# Predict gene expression at a target timepoint
# =============================================================================

#' Predict Gene Expression at a Target Timepoint
#'
#' @description
#' Performs gene-level stochastic prediction using mixture parameter estimation,
#' Wasserstein-based transition matrices, weighted state parameters, and the
#' Semi-Markov / Brownian Bridge simulation engine.
#'
#' The function supports:
#' - forward validation / forecasting;
#' - backward validation / backcasting via time-axis reflection.
#'
#' @param df Data frame containing a numeric TP column and gene expression columns.
#' @param gene Character. Gene column name to model.
#' @param start_tp Numeric. First biological timepoint used for training.
#' @param end_tp Numeric. Last biological timepoint used for training.
#' @param prediction_tp Numeric. Target biological timepoint to predict.
#' @param n_opt Integer. Number of mixture states/components. Default is 2.
#' @param N_SIM Integer. Number of stochastic simulations. Default is 100.
#' @param N_Steps Integer. Number of Brownian Bridge interpolation steps.
#' @param memory_decay Numeric. Exponential memory decay parameter.
#' @param direction Character. Either "forward" or "backward".
#'
#' @return A list containing predicted log-values, predicted count-values,
#' simulation object, transition matrix, mixture parameters, state probabilities,
#' and metadata.
#'
# =============================================================================
# Forward-only prediction function
# =============================================================================

predict_gene_expression_forward <- function(df,
                                            gene,
                                            start_tp,
                                            end_tp,
                                            prediction_tp,
                                            n_opt = 2,
                                            N_SIM = 100,
                                            N_Steps = 200,
                                            memory_decay = 1) {
  
  if (!is.data.frame(df)) stop("df must be a data frame.")
  if (!"TP" %in% colnames(df)) stop("df must contain a numeric TP column.")
  if (!gene %in% colnames(df)) stop("Gene not found: ", gene)
  
  df$TP <- as.numeric(as.character(df$TP))
  
  if (any(is.na(df$TP))) {
    stop("TP column contains NA values after numeric conversion.")
  }
  
  gene_values_raw <- suppressWarnings(as.numeric(df[[gene]]))
  
  if (all(is.na(gene_values_raw))) {
    stop("Gene values are all NA for gene: ", gene)
  }
  
  df_temp <- data.frame(
    TP = df$TP,
    gene_value = log2(gene_values_raw + 1)
  )
  
  df_temp <- df_temp[
    is.finite(df_temp$TP) &
      is.finite(df_temp$gene_value),
  ]
  
  all_tps <- sort(unique(df_temp$TP))
  
  train_tps <- all_tps[
    all_tps >= start_tp &
      all_tps <= end_tp
  ]
  
  train_tps <- sort(train_tps)
  
  if (length(train_tps) < 2) {
    stop("Not enough training timepoints.")
  }
  
  if (prediction_tp <= max(train_tps)) {
    stop("For forward prediction, prediction_tp must be greater than max(train_tps).")
  }
  
  param_df_all <- estimate_mixture_params(
    gene_values = df_temp$gene_value,
    timepoints = df_temp$TP,
    n_opt = n_opt
  )
  
  if (!all(c("TP", "pi", "mu", "sigma") %in% colnames(param_df_all))) {
    stop("estimate_mixture_params() must return TP, pi, mu, sigma.")
  }
  
  param_train <- param_df_all[param_df_all$TP %in% train_tps, ]
  
  param_train$TP_order <- factor(
    param_train$TP,
    levels = train_tps
  )
  
  param_train <- param_train[order(param_train$TP_order), ]
  param_train$TP_order <- NULL
  param_train$TP <- as.numeric(as.character(param_train$TP))
  
  state_dim <- max(table(param_train$TP))
  
  Weights <- exp(-memory_decay * (max(train_tps) - train_tps))
  un_tr <- train_tps[train_tps < max(train_tps)]
  Weights_tr <- exp(-memory_decay * (max(un_tr) - un_tr))
  
  Weights <- as.numeric(Weights)
  Weights_tr <- as.numeric(Weights_tr)
  
  Weights <- Weights / sum(Weights)
  Weights_tr <- Weights_tr / sum(Weights_tr)
  
  Prob <- build_transition_matrix(
    param_df = param_train,
    weights_tr = Weights_tr,
    state_dim = state_dim
  )
  
  Prob <- as.matrix(Prob)
  Prob[!is.finite(Prob)] <- 0
  
  row_sums <- rowSums(Prob)
  if (any(row_sums == 0)) {
    for (zr in which(row_sums == 0)) {
      Prob[zr, ] <- rep(1 / state_dim, state_dim)
    }
  }
  
  Prob <- normalize_rows(Prob)
  
  par_check <- get_weighted_check_params(
    param_df = param_train,
    weights = Weights,
    n_opt = n_opt
  )
  
  par_check <- as.matrix(par_check)
  
  first_tp <- train_tps[1]
  first_param <- param_train[param_train$TP == first_tp, ]
  
  state_probs <- as.numeric(first_param$pi[seq_len(state_dim)])
  
  if (any(!is.finite(state_probs)) || sum(state_probs) <= 0) {
    state_probs <- rep(1 / state_dim, state_dim)
  } else {
    state_probs <- state_probs / sum(state_probs)
  }
  
  Times_sim <- train_tps
  Times_forec <- prediction_tp
  
  full_time_grid <- c(Times_sim, Times_forec)
  
  if (any(diff(full_time_grid) <= 0)) {
    stop("Simulation time grid must be strictly increasing.")
  }
  
  simulation_res <- SIMUL_PROC_rand_init(
    P = Prob,
    Par = par_check,
    Times_sim = Times_sim,
    Times_forec = Times_forec,
    N = N_SIM,
    N_Steps = N_Steps,
    state_probs = state_probs
  )
  
  target_idx <- which.min(abs(simulation_res$time - Times_forec))
  
  pred_log <- as.numeric(simulation_res$X[, target_idx])
  pred_counts <- 2^pred_log - 1
  pred_counts[pred_counts < 0] <- 0
  
  X_counts <- 2^simulation_res$X - 1
  X_counts[X_counts < 0] <- 0
  
  return(list(
    gene = gene,
    direction = "forward",
    start_tp = start_tp,
    end_tp = end_tp,
    prediction_tp = prediction_tp,
    train_tps = train_tps,
    simulation_times = Times_sim,
    forecast_times = Times_forec,
    target_time = prediction_tp,
    pred_log = pred_log,
    pred_counts = pred_counts,
    X_log = simulation_res$X,
    X_counts = X_counts,
    time = simulation_res$time,
    simulation = simulation_res,
    transition_matrix = Prob,
    parameters = param_train,
    mixture_parameters = param_train,
    weighted_parameters = par_check,
    state_probs = state_probs,
    memory_weights = Weights,
    transition_weights = Weights_tr
  ))
}

# =============================================================================
# Backward validation via time reflection
# =============================================================================

predict_gene_expression_backward <- function(df,
                                             gene,
                                             start_tp,
                                             end_tp,
                                             prediction_tp,
                                             n_opt = 2,
                                             N_SIM = 100,
                                             N_Steps = 200,
                                             memory_decay = 1) {
  
  if (!is.data.frame(df)) stop("df must be a data frame.")
  if (!"TP" %in% colnames(df)) stop("df must contain TP column.")
  if (!gene %in% colnames(df)) stop("Gene not found: ", gene)
  
  df_reflected <- df
  df_reflected$TP <- as.numeric(as.character(df_reflected$TP))
  
  if (any(is.na(df_reflected$TP))) {
    stop("TP column contains NA values after numeric conversion.")
  }
  
  max_tp <- max(df_reflected$TP, na.rm = TRUE)
  min_tp <- min(df_reflected$TP, na.rm = TRUE)
  
  reflect_tp <- function(tp) {
    max_tp - tp + min_tp
  }
  
  df_reflected$TP_original <- df_reflected$TP
  df_reflected$TP <- reflect_tp(df_reflected$TP)
  
  start_tp_reflected <- reflect_tp(start_tp)
  end_tp_reflected <- reflect_tp(end_tp)
  prediction_tp_reflected <- reflect_tp(prediction_tp)
  
  start_forward <- min(start_tp_reflected, end_tp_reflected)
  end_forward <- max(start_tp_reflected, end_tp_reflected)
  
  res <- predict_gene_expression_forward(
    df = df_reflected,
    gene = gene,
    start_tp = start_forward,
    end_tp = end_forward,
    prediction_tp = prediction_tp_reflected,
    n_opt = n_opt,
    N_SIM = N_SIM,
    N_Steps = N_Steps,
    memory_decay = memory_decay
  )
  
  res$direction <- "backward"
  
  res$original_start_tp <- start_tp
  res$original_end_tp <- end_tp
  res$original_prediction_tp <- prediction_tp
  
  res$reflected_start_tp <- start_tp_reflected
  res$reflected_end_tp <- end_tp_reflected
  res$reflected_prediction_tp <- prediction_tp_reflected
  
  res$start_tp <- start_tp
  res$end_tp <- end_tp
  res$prediction_tp <- prediction_tp
  res$target_time <- prediction_tp
  
  tp_original <- as.numeric(as.character(df$TP))
  
  res$train_tps <- sort(
    unique(tp_original[tp_original <= start_tp & tp_original >= end_tp]),
    decreasing = TRUE
  )
  
  res$simulation_times_reflected <- res$simulation_times
  res$forecast_times_reflected <- res$forecast_times
  res$time_reflected <- res$time
  
  res$simulation_times <- res$train_tps
  res$forecast_times <- prediction_tp
  res$time <- reflect_tp(res$time_reflected)
  
  
  if (!is.null(res$parameters) && "TP" %in% colnames(res$parameters)) {
    res$parameters_reflected <- res$parameters
    res$parameters$TP <- reflect_tp(res$parameters$TP)
    res$parameters <- res$parameters[
      order(res$parameters$TP, decreasing = TRUE),
    ]
  }
  
  if (!is.null(res$mixture_parameters) &&
      "TP" %in% colnames(res$mixture_parameters)) {
    
    res$mixture_parameters_reflected <- res$mixture_parameters
    res$mixture_parameters$TP <- reflect_tp(res$mixture_parameters$TP)
    res$mixture_parameters <- res$mixture_parameters[
      order(res$mixture_parameters$TP, decreasing = TRUE),
    ]
  }
  
  return(res)
}

# =============================================================================
# Unified wrapper
# =============================================================================

#' @export
predict_gene_expression <- function(df,
                                    gene,
                                    start_tp,
                                    end_tp,
                                    prediction_tp,
                                    n_opt = 2,
                                    N_SIM = 100,
                                    N_Steps = 200,
                                    memory_decay = 1,
                                    direction = c("forward", "backward")) {
  
  direction <- match.arg(direction)
  
  if (direction == "forward") {
    
    return(
      predict_gene_expression_forward(
        df = df,
        gene = gene,
        start_tp = start_tp,
        end_tp = end_tp,
        prediction_tp = prediction_tp,
        n_opt = n_opt,
        N_SIM = N_SIM,
        N_Steps = N_Steps,
        memory_decay = memory_decay
      )
    )
    
  } else {
    
    return(
      predict_gene_expression_backward(
        df = df,
        gene = gene,
        start_tp = start_tp,
        end_tp = end_tp,
        prediction_tp = prediction_tp,
        n_opt = n_opt,
        N_SIM = N_SIM,
        N_Steps = N_Steps,
        memory_decay = memory_decay
      )
    )
  }
}