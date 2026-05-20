#' @title Trajectory Analysis Utilities
#' @description Functions for prediction decay, trajectory uncertainty,
#' transition stability, state entropy, and prediction interval calibration.
#' @author Taras Lukashiv and collaborators
#' @date 2026-05-06


# =============================================================================
# 1. Trajectory Summary Statistics
# =============================================================================

summarize_trajectories <- function(time,
                                   trajectory_matrix,
                                   probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975),
                                   save_path = NULL) {
  
  if (is.null(dim(trajectory_matrix))) {
    stop("trajectory_matrix must be a matrix.")
  }
  
  summary_df <- data.frame(
    Time = time,
    Mean = apply(trajectory_matrix, 2, mean, na.rm = TRUE),
    SD = apply(trajectory_matrix, 2, sd, na.rm = TRUE)
  )
  
  quant_mat <- t(apply(
    trajectory_matrix,
    2,
    stats::quantile,
    probs = probs,
    na.rm = TRUE
  ))
  
  colnames(quant_mat) <- paste0("Q_", probs)
  
  summary_df <- cbind(summary_df, quant_mat)
  
  if (!is.null(save_path)) {
    write.csv(summary_df, save_path, row.names = FALSE)
  }
  
  return(summary_df)
}


# =============================================================================
# 2. Final Prediction Interval Summary
# =============================================================================

summarize_final_prediction <- function(pred_counts,
                                       probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975),
                                       save_path = NULL) {
  
  out <- data.frame(
    Mean = mean(pred_counts, na.rm = TRUE),
    SD = sd(pred_counts, na.rm = TRUE)
  )
  
  quant <- stats::quantile(pred_counts, probs = probs, na.rm = TRUE)
  
  for (i in seq_along(probs)) {
    out[[paste0("Q_", probs[i])]] <- quant[i]
  }
  
  if (!is.null(save_path)) {
    write.csv(out, save_path, row.names = FALSE)
  }
  
  return(out)
}


# =============================================================================
# 3. Prediction Decay Across Horizons
# =============================================================================

analyze_prediction_decay <- function(horizon_results,
                                     horizon_labels = NULL,
                                     p_col = "t_p_val",
                                     fc_col = "log2FC",
                                     p_threshold = 0.05,
                                     fc_threshold = 1,
                                     save_path = NULL) {
  
  if (is.null(horizon_labels)) {
    if (!is.null(names(horizon_results))) {
      horizon_labels <- names(horizon_results)
    } else {
      horizon_labels <- paste0("Horizon_", seq_along(horizon_results))
    }
  }
  
  decay_df <- data.frame()
  
  for (i in seq_along(horizon_results)) {
    
    res <- horizon_results[[i]]
    
    if (!p_col %in% colnames(res)) {
      stop("p_col not found in horizon result: ", p_col)
    }
    
    if (!fc_col %in% colnames(res)) {
      stop("fc_col not found in horizon result: ", fc_col)
    }
    
    adj_p <- p.adjust(res[[p_col]], method = "BH")
    abs_fc <- abs(res[[fc_col]])
    
    temp <- data.frame(
      Horizon = horizon_labels[i],
      Success_Rate = mean(adj_p >= p_threshold & abs_fc < fc_threshold, na.rm = TRUE) * 100,
      Small_Diff_Signif_Rate = mean(adj_p < p_threshold & abs_fc < fc_threshold, na.rm = TRUE) * 100,
      Fail_Rate = mean(adj_p < p_threshold & abs_fc >= fc_threshold, na.rm = TRUE) * 100,
      Large_Diff_Insignif_Rate = mean(adj_p >= p_threshold & abs_fc >= fc_threshold, na.rm = TRUE) * 100,
      Mean_abs_FC = mean(abs_fc, na.rm = TRUE),
      Median_abs_FC = median(abs_fc, na.rm = TRUE),
      Mean_W_log = if ("W_log_dist" %in% colnames(res)) mean(res$W_log_dist, na.rm = TRUE) else NA_real_,
      Mean_W_orig = if ("W_orig_dist" %in% colnames(res)) mean(res$W_orig_dist, na.rm = TRUE) else NA_real_,
      N_genes = sum(!is.na(res[[p_col]])),
      stringsAsFactors = FALSE
    )
    
    decay_df <- rbind(decay_df, temp)
  }
  
  if (!is.null(save_path)) {
    write.csv(decay_df, save_path, row.names = FALSE)
  }
  
  return(decay_df)
}


# =============================================================================
# 4. Plot Prediction Decay
# =============================================================================

plot_prediction_decay <- function(decay_df,
                                  save_path = NULL,
                                  metric = c("Success_Rate", "Fail_Rate", "Mean_abs_FC", "Mean_W_log")) {
  
  metric <- match.arg(metric)
  
  p <- ggplot2::ggplot(
    decay_df,
    ggplot2::aes(x = factor(Horizon, levels = Horizon), y = .data[[metric]], group = 1)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      title = paste("Prediction Decay:", metric),
      x = "Prediction Horizon",
      y = metric
    ) +
    ggplot2::theme_minimal()
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = 8,
      height = 5,
      dpi = 300,
      bg = "white"
    )
  }
  
  return(p)
}


# =============================================================================
# 5. Transition Matrix Entropy
# =============================================================================

calculate_transition_entropy <- function(P,
                                         normalize = TRUE) {
  
  P <- as.matrix(P)
  
  P[P <= 0 | is.na(P)] <- NA
  
  row_entropy <- apply(P, 1, function(row) {
    probs <- row[!is.na(row)]
    -sum(probs * log(probs))
  })
  
  if (normalize) {
    max_entropy <- log(ncol(P))
    row_entropy <- row_entropy / max_entropy
  }
  
  data.frame(
    State = seq_len(nrow(P)),
    Entropy = row_entropy
  )
}


# =============================================================================
# 6. Transition Matrix Stability
# =============================================================================

compare_transition_matrices <- function(P1,
                                        P2,
                                        method = c("frobenius", "mean_abs", "max_abs")) {
  
  method <- match.arg(method)
  
  P1 <- as.matrix(P1)
  P2 <- as.matrix(P2)
  
  if (!all(dim(P1) == dim(P2))) {
    stop("Transition matrices must have the same dimensions.")
  }
  
  diff_mat <- P1 - P2
  
  score <- switch(
    method,
    frobenius = sqrt(sum(diff_mat^2, na.rm = TRUE)),
    mean_abs = mean(abs(diff_mat), na.rm = TRUE),
    max_abs = max(abs(diff_mat), na.rm = TRUE)
  )
  
  list(
    score = score,
    difference_matrix = diff_mat
  )
}


# =============================================================================
# 7. Dominant Transition Path
# =============================================================================

extract_dominant_transition_path <- function(P,
                                             start_state = NULL,
                                             n_steps = 5) {
  
  P <- as.matrix(P)
  
  if (is.null(start_state)) {
    start_state <- which.max(rowSums(P))
  }
  
  states <- numeric(n_steps + 1)
  states[1] <- start_state
  
  for (i in seq_len(n_steps)) {
    states[i + 1] <- which.max(P[states[i], ])
  }
  
  data.frame(
    Step = seq_len(n_steps + 1) - 1,
    State = states
  )
}


# =============================================================================
# 8. State Occupancy Summary
# =============================================================================

summarize_state_occupancy <- function(state_matrix) {
  
  if (is.null(dim(state_matrix))) {
    stop("state_matrix must be a matrix with simulations as rows and timepoints as columns.")
  }
  
  timepoints <- seq_len(ncol(state_matrix))
  
  out <- data.frame()
  
  for (t in timepoints) {
    
    tab <- table(state_matrix[, t])
    
    temp <- data.frame(
      Time_Index = t,
      State = as.integer(names(tab)),
      Count = as.integer(tab),
      Proportion = as.numeric(tab) / nrow(state_matrix)
    )
    
    out <- rbind(out, temp)
  }
  
  return(out)
}


# =============================================================================
# 9. Prediction Interval Calibration for One Gene
# =============================================================================

evaluate_gene_interval_calibration <- function(real_counts,
                                               pred_counts,
                                               lower_prob = 0.025,
                                               upper_prob = 0.975) {
  
  lower <- stats::quantile(pred_counts, probs = lower_prob, na.rm = TRUE)
  upper <- stats::quantile(pred_counts, probs = upper_prob, na.rm = TRUE)
  
  covered <- real_counts >= lower & real_counts <= upper
  
  data.frame(
    Lower = as.numeric(lower),
    Upper = as.numeric(upper),
    Coverage_Rate = mean(covered, na.rm = TRUE) * 100,
    Mean_Real = mean(real_counts, na.rm = TRUE),
    Mean_Pred = mean(pred_counts, na.rm = TRUE),
    Median_Real = median(real_counts, na.rm = TRUE),
    Median_Pred = median(pred_counts, na.rm = TRUE),
    Interval_Width = as.numeric(upper - lower)
  )
}


# =============================================================================
# 10. Genome-wide Calibration
# =============================================================================

evaluate_genomewide_calibration <- function(df,
                                            df_pred,
                                            target_tp,
                                            genes,
                                            lower_prob = 0.025,
                                            upper_prob = 0.975,
                                            save_path = NULL) {
  
  out <- data.frame()
  
  for (gene in genes) {
    
    if (!gene %in% colnames(df)) next
    if (!gene %in% df_pred$Gene) next
    
    real_counts <- as.numeric(df[df$TP == target_tp, gene])
    
    pred_counts <- as.numeric(
      df_pred[df_pred$Gene == gene, 2:ncol(df_pred)]
    )
    
    temp <- evaluate_gene_interval_calibration(
      real_counts = real_counts,
      pred_counts = pred_counts,
      lower_prob = lower_prob,
      upper_prob = upper_prob
    )
    
    temp$Gene <- gene
    
    out <- rbind(out, temp)
  }
  
  out <- out[, c(
    "Gene",
    "Lower",
    "Upper",
    "Coverage_Rate",
    "Mean_Real",
    "Mean_Pred",
    "Median_Real",
    "Median_Pred",
    "Interval_Width"
  )]
  
  if (!is.null(save_path)) {
    write.csv(out, save_path, row.names = FALSE)
  }
  
  return(out)
}


# =============================================================================
# 11. Plot Genome-wide Calibration
# =============================================================================

plot_genomewide_calibration <- function(calibration_df,
                                        save_path = NULL) {
  
  p <- ggplot2::ggplot(
    calibration_df,
    ggplot2::aes(x = Coverage_Rate)
  ) +
    ggplot2::geom_histogram(
      bins = 30,
      fill = "steelblue",
      color = "white",
      alpha = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = mean(calibration_df$Coverage_Rate, na.rm = TRUE),
      linetype = "dashed",
      linewidth = 1
    ) +
    ggplot2::labs(
      title = "Genome-wide Prediction Interval Calibration",
      subtitle = paste(
        "Mean coverage:",
        round(mean(calibration_df$Coverage_Rate, na.rm = TRUE), 2),
        "%"
      ),
      x = "Per-gene Coverage Rate (%)",
      y = "Number of Genes"
    ) +
    ggplot2::theme_minimal()
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = 8,
      height = 5,
      dpi = 300,
      bg = "white"
    )
  }
  
  return(p)
}


# =============================================================================
# 12. Trajectory Deviation from Observed Timepoints
# =============================================================================

calculate_trajectory_observed_deviation <- function(df,
                                                    gene,
                                                    trajectory_time,
                                                    trajectory_matrix,
                                                    count_scale = TRUE) {
  
  if (!gene %in% colnames(df)) {
    stop("Gene not found in df: ", gene)
  }
  
  if (!"TP" %in% colnames(df)) {
    stop("df must contain TP column.")
  }
  
  mean_traj <- apply(trajectory_matrix, 2, mean, na.rm = TRUE)
  
  observed_tps <- sort(unique(df$TP))
  common_tps <- intersect(observed_tps, trajectory_time)
  
  out <- data.frame()
  
  for (tp in common_tps) {
    
    obs <- as.numeric(df[df$TP == tp, gene])
    
    if (!count_scale) {
      obs <- log2(obs + 1)
    }
    
    pred_mean <- mean_traj[which(trajectory_time == tp)]
    
    temp <- data.frame(
      Gene = gene,
      TP = tp,
      Observed_Mean = mean(obs, na.rm = TRUE),
      Predicted_Mean = pred_mean,
      Mean_Difference = mean(obs, na.rm = TRUE) - pred_mean,
      Abs_Mean_Difference = abs(mean(obs, na.rm = TRUE) - pred_mean)
    )
    
    out <- rbind(out, temp)
  }
  
  return(out)
}