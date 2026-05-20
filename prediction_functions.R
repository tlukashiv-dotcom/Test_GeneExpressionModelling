#' @title Prediction Functions for Dynamic Gene Expression Modelling
#' @description Core prediction wrapper for forward validation, backward validation,
#' and future forecasting.
#' @author Taras Lukashiv and collaborators
#' @date 2026-05-06


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
#' - backward validation / backcasting.
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
                                    direction = "forward" # "backward"
                                    ) {
  
  #direction <- match.arg(direction)
  
  # ---------------------------------------------------------------------------
  # 0. Basic checks
  # ---------------------------------------------------------------------------
  
  if (!is.data.frame(df)) {
    stop("df must be a data frame.")
  }
  
  if (!"TP" %in% colnames(df)) {
    stop("df must contain a numeric TP column.")
  }
  
  if (!gene %in% colnames(df)) {
    stop("Gene not found in expression matrix: ", gene)
  }
  
  if (!is.numeric(df$TP)) {
    df$TP <- as.numeric(as.character(df$TP))
  }
  
  if (any(is.na(df$TP))) {
    stop("TP column contains NA values after numeric conversion.")
  }
  
  # ---------------------------------------------------------------------------
  # 1. Prepare gene-level data
  # ---------------------------------------------------------------------------
  
  gene_values_raw <- suppressWarnings(as.numeric(df[[gene]]))
  
  if (all(is.na(gene_values_raw))) {
    stop("Gene values are all NA for gene: ", gene)
  }
  
  df_temp <- data.frame(
    TP = as.numeric(df$TP),
    gene_value = log2(gene_values_raw + 1)
  )
  
  df_temp <- df_temp[
    is.finite(df_temp$TP) &
      is.finite(df_temp$gene_value),
  ]
  
  # ---------------------------------------------------------------------------
  #  Change time points by directions
  # ---------------------------------------------------------------------------
  
  if(direction =="backward"){
    max_tp = max(df_temp$TP)
    df_temp$TP = max_tp - df_temp$TP
    start_tp = max_tp - start_tp
    end_tp = max_tp - end_tp
    prediction_tp = max_tp - prediction_tp
  }
  
  all_tps <- sort(unique(df_temp$TP))
  
  if (length(all_tps) < 3) {
    stop("Not enough timepoints for prediction.")
  }
  
  if (!prediction_tp %in% all_tps && prediction_tp <= max(all_tps)) {
    stop("prediction_tp is inside observed range but not present in df$TP.")
  }
 
  # ---------------------------------------------------------------------------
  # 2. Define training timepoints
  # ---------------------------------------------------------------------------
  train_tps <- all_tps[
    all_tps >= min(c(start_tp,end_tp)) &
      all_tps <= max(c(start_tp,end_tp))
    ]
  train_tps <- sort(train_tps)
  
  if (length(train_tps) < 2) {
    stop("Not enough training timepoints.")
  }
  
  # For validation, prediction_tp may be observed.
  # For forecasting, prediction_tp may be outside the observed range.
  if (prediction_tp <= max(train_tps)) {
    stop("For forward prediction, prediction_tp must be greater than max(train_tps).")
  }
  
  # ---------------------------------------------------------------------------
  # 3. Estimate mixture parameters
  # ---------------------------------------------------------------------------
  
  param_df_all <- estimate_mixture_params(
    gene_values = df_temp$gene_value,
    timepoints = df_temp$TP,
    n_opt = n_opt
  )
  
  if (!all(c("TP", "pi", "mu", "sigma") %in% colnames(param_df_all))) {
    stop("estimate_mixture_params() must return columns: TP, pi, mu, sigma.")
  }
  
  param_train <- param_df_all[param_df_all$TP %in% train_tps, ]
  
  if (nrow(param_train) == 0) {
    stop("No mixture parameters estimated for training timepoints.")
  }
  
  # ---------------------------------------------------------------------------
  # 4. Reorder parameter table according to temporal direction
  #    This is essential for backward transition estimation.
  # ---------------------------------------------------------------------------
  
  param_train$TP_order <- factor(
    param_train$TP,
    levels = train_tps
  )
  
  param_train <- param_train[order(param_train$TP_order), ]
  param_train$TP_order <- NULL
  
  param_train$TP <- as.numeric(as.character(param_train$TP))
  
  state_dim <- max(table(param_train$TP))
  
  if (!is.finite(state_dim) || state_dim < 1) {
    stop("Invalid state_dim estimated.")
  }
  
  if (state_dim != n_opt) {
    warning(
      "state_dim (", state_dim, ") differs from n_opt (", n_opt,
      "). Proceeding with state_dim."
    )
  }
  
  # ---------------------------------------------------------------------------
  # 5. Memory weights
  # ---------------------------------------------------------------------------
  
  Weights <- exp(-memory_decay * (max(train_tps) - train_tps))
  
  un_tr <- train_tps[train_tps < max(train_tps)]
  
  if (length(un_tr) == 0) {
    stop("Not enough transition timepoints for forward prediction.")
  }
  
  Weights_tr <- exp(-memory_decay * (max(un_tr) - un_tr))
  
  Weights <- as.numeric(Weights)
  Weights_tr <- as.numeric(Weights_tr)
  
  Weights <- Weights / sum(Weights)
  Weights_tr <- Weights_tr / sum(Weights_tr)
  
  if (any(!is.finite(Weights)) || any(!is.finite(Weights_tr))) {
    stop("Non-finite memory weights detected.")
  }
  
  # ---------------------------------------------------------------------------
  # 6. Transition matrix
  # ---------------------------------------------------------------------------
  
  Prob <- build_transition_matrix(
    param_df = param_train,
    weights_tr = Weights_tr,
    state_dim = state_dim
  )
  
  Prob <- as.matrix(Prob)
  
  if (nrow(Prob) != state_dim || ncol(Prob) != state_dim) {
    stop("Transition matrix has invalid dimensions.")
  }
  
  Prob[!is.finite(Prob)] <- 0
  
  row_sums <- rowSums(Prob)
  if (any(row_sums == 0)) {
    zero_rows <- which(row_sums == 0)
    for (zr in zero_rows) {
      Prob[zr, ] <- rep(1 / state_dim, state_dim)
    }
  }
  
  Prob <- normalize_rows(Prob)
  
  # ---------------------------------------------------------------------------
  # 7. Weighted state parameters
  # ---------------------------------------------------------------------------
  
  par_check <- get_weighted_check_params(
    param_df = param_train,
    weights = Weights,
    n_opt = n_opt
  )
  
  par_check <- as.matrix(par_check)
  
  # SIMUL_PROC_rand_init() accepts either [2 x state_dim] or [state_dim x 2].
  if (!(all(dim(par_check) == c(2, state_dim)) ||
        all(dim(par_check) == c(state_dim, 2)))) {
    stop("Weighted parameter matrix has invalid dimensions.")
  }
  
  # ---------------------------------------------------------------------------
  # 8. Initial state probabilities
  # ---------------------------------------------------------------------------
  
  first_tp <- train_tps[1]
  first_param <- param_train[param_train$TP == first_tp, ]
  
  if (nrow(first_param) < state_dim) {
    stop("Not enough mixture components for initial timepoint.")
  }
  
  state_probs <- as.numeric(first_param$pi[seq_len(state_dim)])
  
  if (any(!is.finite(state_probs)) || sum(state_probs) <= 0) {
    state_probs <- rep(1 / state_dim, state_dim)
  } else {
    state_probs <- state_probs / sum(state_probs)
  }
  
  # ---------------------------------------------------------------------------
  # 9. Simulation timeline
  # ---------------------------------------------------------------------------
  
  # Biological and simulation time are aligned.
  Times_sim <- seq(min(train_tps), max(train_tps))
  Times_forec <- c(prediction_tp)
  
  if (any(diff(c(Times_sim, Times_forec)) <= 0)) {
    stop("Simulation time grid must be strictly increasing.")
  }
  
  # ---------------------------------------------------------------------------
  # 10. Stochastic simulation
  # ---------------------------------------------------------------------------
  
  simulation_res <- SIMUL_PROC_rand_init(
    P = Prob,
    Par = par_check,
    Times_sim = Times_sim,
    Times_forec = Times_forec,
    N = N_SIM,
    N_Steps = N_Steps,
    state_probs = state_probs
  )
  
  if (!is.list(simulation_res) ||
      !"X" %in% names(simulation_res) ||
      !"time" %in% names(simulation_res)) {
    stop("SIMUL_PROC_rand_init() must return list with X and time.")
  }
  
  # ---------------------------------------------------------------------------
  # 11. Extract target prediction
  # ---------------------------------------------------------------------------
  
  target_idx <- which.min(abs(simulation_res$time - prediction_tp))
  
  if (length(target_idx) == 0 || is.na(target_idx)) {
    stop("Could not locate target prediction time in simulation output.")
  }
  
  pred_log <- as.numeric(simulation_res$X[, target_idx])
  pred_counts <- 2^pred_log - 1
  
  pred_counts[pred_counts < 0] <- 0
  
  # ---------------------------------------------------------------------------
  # 12. Inverse time change for direction = 'backward'
  # ---------------------------------------------------------------------------
  
  if(direction =="backward"){
    df_temp$TP = max_tp - df_temp$TP
    start_tp = max_tp - start_tp
    end_tp = max_tp - end_tp
    prediction_tp = max_tp - prediction_tp
    train_tps = max_tp - train_tps
    Times_sim = max_tp - Times_sim
    Times_forec = max_tp - Times_forec
    simulation_res$time = max_tp - simulation_res$time
  }
  
  # ---------------------------------------------------------------------------
  # 13. Return output
  # ---------------------------------------------------------------------------
  
  return(list(
    gene = gene,
    direction = direction,
    
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
    X_counts = 2^simulation_res$X - 1,
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
