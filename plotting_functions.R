#' @title Plotting Functions for Dynamic Gene Expression Modelling
#' @description Publication-ready plotting utilities for validation, forecasting,
#' density comparisons, trajectory visualization, and failure analysis.
#' @author Taras Lukashiv and collaborators
#' @date 2026-05-20


# =============================================================================
# 1. Violin Plot for Gene Expression by Timepoint
# =============================================================================

plot_violin_expression <- function(df,
                                   gene,
                                   save_path = NULL,
                                   log_transform = TRUE) {
  
  if (!gene %in% colnames(df)) {
    stop("Gene not found in df: ", gene)
  }
  
  if (!"TP" %in% colnames(df)) {
    stop("df must contain TP column.")
  }
  
  gene_values <- as.numeric(df[[gene]])
  
  if (log_transform) {
    gene_values <- log2(gene_values + 1)
    y_label <- paste0("log2(", gene, " + 1)")
  } else {
    y_label <- paste0(gene, " expression")
  }
  
  plot_df <- data.frame(
    TP_factor = as.factor(df$TP),
    gene_val = gene_values
  )
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = TP_factor, y = gene_val)) +
    ggplot2::geom_violin(
      trim = FALSE,
      fill = "#A4A4A4",
      color = "darkred",
      alpha = 0.5
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 1),
      alpha = 0.4,
      size = 1
    ) +
    ggplot2::geom_boxplot(
      width = 0.1,
      outlier.shape = NA,
      fill = "white",
      alpha = 0.5
    ) +
    ggplot2::labs(
      title = paste("Violin Plot for Gene:", gene),
      x = "Timepoint",
      y = y_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = 7,
      height = 5,
      dpi = 300,
      bg = "white"
    )
  }
  
  return(p)
}


# =============================================================================
# 2. Mixture Density Plot
# =============================================================================

plot_mixture_density <- function(mix_model,
                                 gene,
                                 tp,
                                 save_path = NULL) {
  
  vals <- as.numeric(mix_model$data)
  
  if (length(vals) == 0 || all(is.na(vals))) {
    stop("No valid data found in mix_model.")
  }
  
  x_range <- seq(
    min(vals, na.rm = TRUE) - diff(range(vals, na.rm = TRUE)) * 0.1,
    max(vals, na.rm = TRUE) + diff(range(vals, na.rm = TRUE)) * 0.1,
    length.out = 400
  )
  
  density_df <- data.frame(x = x_range)
  
  total_density <- mclust::dens(
    modelName = mix_model$modelName,
    data = x_range,
    parameters = mix_model$parameters
  )
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_histogram(
      data = data.frame(value = vals),
      ggplot2::aes(x = value, y = ggplot2::after_stat(density)),
      bins = 30,
      fill = "grey85",
      color = "white",
      alpha = 0.8
    ) +
    ggplot2::geom_line(
      data = data.frame(x = x_range, y = total_density),
      ggplot2::aes(x = x, y = y),
      color = "black",
      linewidth = 1
    ) +
    ggplot2::labs(
      title = paste("Mixture Density for", gene, "at Week", tp),
      x = paste0("log2(", gene, " + 1)"),
      y = "Density"
    ) +
    ggplot2::theme_minimal()
  
  n_components <- mix_model$G
  
  vars <- mix_model$parameters$variance$sigmasq
  
  for (k in seq_len(n_components)) {
    
    mu <- mix_model$parameters$mean[k]
    sigma <- sqrt(if (length(vars) > 1) vars[k] else vars)
    pro <- mix_model$parameters$pro[k]
    
    comp_density <- pro * stats::dnorm(x_range, mean = mu, sd = sigma)
    
    p <- p +
      ggplot2::geom_line(
        data = data.frame(x = x_range, y = comp_density),
        ggplot2::aes(x = x, y = y),
        linetype = "dashed",
        linewidth = 0.8
      )
  }
  
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
# 3. Mclust Classification / Uncertainty Plot
# =============================================================================

plot_mixture_classification <- function(mix_model,
                                        gene,
                                        tp,
                                        save_path = NULL) {
  
  vals <- as.numeric(mix_model$data)
  
  uncertainty <- 1 - apply(mix_model$z, 1, max)
  
  plot_df <- data.frame(
    expression = vals,
    uncertainty = uncertainty,
    cluster = as.factor(mix_model$classification)
  )
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = expression, y = uncertainty, color = cluster)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::geom_vline(
      xintercept = mix_model$parameters$mean,
      linetype = "dotted",
      linewidth = 0.8
    ) +
    ggplot2::labs(
      title = paste("Classification Uncertainty for", gene, "at Week", tp),
      x = paste0("log2(", gene, " + 1)"),
      y = "Classification Uncertainty",
      color = "Cluster"
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
# 4. Simulated Trajectories Plot
# =============================================================================

plot_simulated_trajectories <- function(time,
                                        trajectory_matrix,
                                        gene,
                                        save_path = NULL,
                                        y_cap = 1000,
                                        interval = c(0.05, 0.95)) {
  
  if (is.null(dim(trajectory_matrix))) {
    stop("trajectory_matrix must be a matrix.")
  }
  
  avg_traj <- apply(trajectory_matrix, 2, mean, na.rm = TRUE)
  lower <- apply(trajectory_matrix, 2, stats::quantile, probs = interval[1], na.rm = TRUE)
  upper <- apply(trajectory_matrix, 2, stats::quantile, probs = interval[2], na.rm = TRUE)
  
  plot_df <- data.frame(
    time = time,
    mean = avg_traj,
    lower = lower,
    upper = upper
  )
  
  traj_df <- as.data.frame(t(trajectory_matrix))
  traj_df$time <- time
  
  traj_long <- tidyr::pivot_longer(
    traj_df,
    cols = -time,
    names_to = "simulation",
    values_to = "expression"
  )
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = traj_long,
      ggplot2::aes(x = time, y = pmin(expression, y_cap), group = simulation),
      alpha = 0.08
    ) +
    ggplot2::geom_ribbon(
      data = plot_df,
      ggplot2::aes(x = time, ymin = pmin(lower, y_cap), ymax = pmin(upper, y_cap)),
      alpha = 0.25
    ) +
    ggplot2::geom_line(
      data = plot_df,
      ggplot2::aes(x = time, y = pmin(mean, y_cap)),
      linewidth = 1.1
    ) +
    ggplot2::labs(
      title = paste("Simulated Trajectories for", gene),
      x = "Time",
      y = paste("Expression of", gene)
    ) +
    ggplot2::theme_minimal()
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = 9,
      height = 6,
      dpi = 300,
      bg = "white"
    )
  }
  
  return(p)
}


# =============================================================================
# 5. Average Integer Trajectory Plot
# =============================================================================

plot_average_integer_trajectory <- function(time,
                                            trajectory_matrix,
                                            gene,
                                            save_path = NULL,
                                            y_cap = 100) {
  
  avg_traj_int <- round(apply(trajectory_matrix, 2, mean, na.rm = TRUE), 0)
  
  plot_df <- data.frame(
    time = time,
    expression = pmin(avg_traj_int, y_cap)
  )
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time, y = expression)) +
    ggplot2::geom_area(alpha = 0.25) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::labs(
      title = paste("Average Integer Trajectory:", gene),
      x = "Time",
      y = gene
    ) +
    ggplot2::theme_minimal()
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = 9,
      height = 6,
      dpi = 300,
      bg = "white"
    )
  }
  
  return(p)
}


# =============================================================================
# 6. Mean Difference Histogram
# =============================================================================

plot_mean_difference_histogram <- function(results_df,
                                           save_path = NULL,
                                           lower_bound = -10) {
  
  plot_data <- subset(
    results_df,
    mean_diff > lower_bound & is.finite(mean_diff)
  )
  
  avg_val <- mean(plot_data$mean_diff, na.rm = TRUE)
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = mean_diff)) +
    ggplot2::geom_histogram(
      bins = 200,
      fill = "steelblue",
      color = "white",
      alpha = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = avg_val,
      color = "firebrick",
      linetype = "dashed",
      linewidth = 1
    ) +
    ggplot2::coord_cartesian(xlim = c(-3, 3)) +
    ggplot2::labs(
      title = "Distribution of Mean Differences",
      subtitle = paste("Global Mean difference:", round(avg_val, 3)),
      x = "Mean Difference (True - Simulated)",
      y = "Frequency"
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.major = ggplot2::element_line(color = "grey95"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black")
    )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = 10,
      height = 6,
      dpi = 300,
      bg = "white"
    )
  }
  
  return(p)
}


# =============================================================================
# 7. Scatter Plots: Adjusted P-values and Log2FC
# =============================================================================

plot_metric_scatter_panels <- function(results_df,
                                       save_path_combined = NULL,
                                       save_path_pvals = NULL,
                                       save_path_fc = NULL) {
  
  plot_df <- results_df
  plot_df$adj_p_t <- p.adjust(plot_df$t_p_val, method = "BH")
  plot_df$adj_p_wilcox <- p.adjust(plot_df$Wilcox_p_val, method = "BH")
  
  p_pvals <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = -log10(adj_p_t), y = -log10(adj_p_wilcox))
  ) +
    ggplot2::geom_point(alpha = 0.4, color = "steelblue") +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "red"
    ) +
    ggplot2::labs(
      title = "Comparison of Adjusted P-values",
      x = "-log10(Adj. P-value) T-test",
      y = "-log10(Adj. P-value) Wilcoxon"
    ) +
    ggplot2::theme_minimal()
  
  p_fc <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = log2FC, y = log2FC_counts)
  ) +
    ggplot2::geom_point(alpha = 0.4, color = "darkgreen") +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "red"
    ) +
    ggplot2::labs(
      title = "Comparison of Log2FC",
      x = "Log2FC (Log-transformed data)",
      y = "Log2FC (Counts data)"
    ) +
    ggplot2::theme_minimal()
  
  combined <- patchwork::wrap_plots(p_pvals, p_fc, ncol = 2)
  
  if (!is.null(save_path_combined)) {
    ggplot2::ggsave(
      filename = save_path_combined,
      plot = combined,
      width = 12,
      height = 6,
      dpi = 300,
      bg = "white"
    )
  }
  
  if (!is.null(save_path_pvals)) {
    ggplot2::ggsave(
      filename = save_path_pvals,
      plot = p_pvals,
      width = 6,
      height = 6,
      dpi = 300,
      bg = "white"
    )
  }
  
  if (!is.null(save_path_fc)) {
    ggplot2::ggsave(
      filename = save_path_fc,
      plot = p_fc,
      width = 6,
      height = 6,
      dpi = 300,
      bg = "white"
    )
  }
  
  return(list(
    p_values = p_pvals,
    log2fc = p_fc,
    combined = combined
  ))
}


# =============================================================================
# 8. Mean-Variance / CV Plot
# =============================================================================

plot_mean_variance_cv <- function(results_df,
                                  save_path = NULL) {
  
  plot_df <- results_df
  
  plot_df$cv_real <- plot_df$sd_real_count / plot_df$mean_real_count
  plot_df$cv_pred <- plot_df$sd_pred_count / plot_df$mean_pred_count
  
  p <- ggplot2::ggplot(plot_df) +
    ggplot2::geom_point(
      ggplot2::aes(x = log10(mean_pred_count + 1), y = cv_pred, color = "Predicted"),
      alpha = 0.3
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = log10(mean_real_count + 1), y = cv_real, color = "Real"),
      alpha = 0.3
    ) +
    ggplot2::scale_color_manual(
      values = c("Real" = "grey50", "Predicted" = "firebrick")
    ) +
    ggplot2::labs(
      title = "Mean-Variance Relationship",
      x = "Log10(Mean Expression + 1)",
      y = "Coefficient of Variation (CV)",
      color = "Data Type"
    ) +
    ggplot2::theme_minimal()
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = 10,
      height = 6,
      dpi = 300,
      bg = "white"
    )
  }
  
  return(p)
}


# =============================================================================
# 9. Volcano Failure Plot
# =============================================================================

plot_volcano_failures <- function(results_df,
                                  save_path = NULL,
                                  label_top_n = 10) {
  
  plot_df <- results_df
  plot_df$adj_p_t <- p.adjust(plot_df$t_p_val, method = "BH")
  
  fail_list <- plot_df %>%
    dplyr::filter(adj_p_t < 0.05 & abs(log2FC) >= 1) %>%
    dplyr::arrange(dplyr::desc(abs(log2FC)))
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = log2FC, y = -log10(adj_p_t))) +
    ggplot2::geom_point(alpha = 0.2, color = "grey70") +
    ggplot2::geom_point(
      data = fail_list,
      ggplot2::aes(color = "Critical Failures"),
      alpha = 0.6
    ) +
    ggrepel::geom_text_repel(
      data = head(fail_list, label_top_n),
      ggplot2::aes(label = Names),
      size = 3,
      max.overlaps = 10
    ) +
    ggplot2::scale_color_manual(values = c("Critical Failures" = "firebrick")) +
    ggplot2::labs(
      title = "Identification of Critical Failure Genes",
      x = "Log2 Fold Change (Observed vs Predicted)",
      y = "-log10(Adjusted P-value)"
    ) +
    ggplot2::theme_minimal()
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = 8,
      height = 6,
      dpi = 300,
      bg = "white"
    )
  }
  
  return(list(plot = p, fail_list = fail_list))
}


# =============================================================================
# 10. Density Comparison for Real vs Predicted Counts
# =============================================================================

plot_real_predicted_density <- function(gene,
                                        real_counts,
                                        pred_counts,
                                        save_path = NULL,
                                        title_prefix = "Distribution Comparison") {
  
  plot_df <- data.frame(
    expression = c(real_counts, pred_counts),
    source = factor(
      c(rep("Real", length(real_counts)), rep("Predicted", length(pred_counts))),
      levels = c("Real", "Predicted")
    )
  )
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = expression, fill = source, color = source)
  ) +
    ggplot2::geom_density(alpha = 0.5, adjust = 0.5) +
    ggplot2::geom_jitter(
      ggplot2::aes(y = -0.002),
      height = 0.001,
      alpha = 0.3,
      size = 1
    ) +
    ggplot2::labs(
      title = paste(title_prefix, ":", gene),
      subtitle = "Empirical distribution (Real) vs. Stochastic prediction (Predicted)",
      x = "Expression Counts",
      y = "Density",
      fill = "Data Source",
      color = "Data Source"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
    )
  
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
# 11. Stratified Accuracy Plot
# =============================================================================

plot_stratified_accuracy <- function(results_df,
                                     save_path = NULL,
                                     n_bins = 10) {
  
  plot_df <- results_df
  plot_df$adj_p_t <- p.adjust(plot_df$t_p_val, method = "BH")
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      decile = dplyr::ntile(mean_real_count, n_bins),
      is_success = ifelse(adj_p_t >= 0.05 & abs(log2FC) < 1, 1, 0)
    )
  
  strat_accuracy <- plot_df %>%
    dplyr::group_by(decile) %>%
    dplyr::summarise(
      mean_exp = mean(mean_real_count, na.rm = TRUE),
      accuracy = mean(is_success, na.rm = TRUE) * 100,
      n_genes = dplyr::n(),
      .groups = "drop"
    )
  
  p <- ggplot2::ggplot(strat_accuracy, ggplot2::aes(x = factor(decile), y = accuracy)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(round(accuracy, 1), "%")),
      vjust = -0.5,
      size = 3.5
    ) +
    ggplot2::labs(
      title = "Model Accuracy Stratified by Expression Magnitude",
      subtitle = "Failures are concentrated in the lowest expression deciles",
      x = "Expression Decile (1 = Lowest, 10 = Highest)",
      y = "Predictive Success Rate (%)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = 10,
      height = 6,
      dpi = 300,
      bg = "white"
    )
  }
  
  return(list(plot = p, stats = strat_accuracy))
}


# =============================================================================
# 12. Additional Validation Plotting Helpers
# =============================================================================

# -----------------------------------------------------------------------------
# 1. plot_gene_density
# -----------------------------------------------------------------------------

plot_gene_density <- function(gene_name, real_data, pred_data, title_prefix = "Validation") {
  plot_df <- data.frame(
    expression = c(real_data, pred_data),
    source = factor(
      c(rep("Real", length(real_data)), rep("Predicted", length(pred_data))),
      levels = c("Real", "Predicted")
    )
  )
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = expression, fill = source, color = source)) +
    ggplot2::geom_density(alpha = 0.5, adjust = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste(title_prefix, ":", gene_name),
      x = "Expression Counts",
      y = "Density",
      fill = "Data Source",
      color = "Data Source"
    )
}

# -----------------------------------------------------------------------------
# 2. plot_stratified_deviation_distribution
# -----------------------------------------------------------------------------

plot_stratified_deviation_distribution <- function(stratified_object,
                                                   y_var = c("abs_log2FC", "abs_log2FC_counts", "abs_mean_diff"),
                                                   title = "Stratified Deviation Distribution") {
  
  y_var <- match.arg(y_var)
  
  df <- stratified_object$data
  
  ggplot2::ggplot(df, ggplot2::aes(x = factor(stratum), y = .data[[y_var]])) +
    ggplot2::geom_violin(fill = "grey80", color = "grey40", alpha = 0.8) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.7) +
    ggplot2::labs(
      title = title,
      x = "Expression Stratum",
      y = y_var
    ) +
    ggplot2::theme_minimal()
}

# -----------------------------------------------------------------------------
# 3. plot_extended_density_comparison
# -----------------------------------------------------------------------------

plot_extended_density_comparison <- function(gene_name,
                                             real_counts,
                                             predicted_counts,
                                             title_prefix = "Density Comparison",
                                             save_path = NULL) {
  
  plot_df <- data.frame(
    expression = c(real_counts, predicted_counts),
    source = factor(
      c(
        rep("Real", length(real_counts)),
        rep("Predicted", length(predicted_counts))
      ),
      levels = c("Real", "Predicted")
    )
  )
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = expression, fill = source, color = source)
  ) +
    ggplot2::geom_density(alpha = 0.45, adjust = 0.5) +
    ggplot2::geom_jitter(
      ggplot2::aes(y = -0.002),
      height = 0.001,
      alpha = 0.25,
      size = 0.8
    ) +
    ggplot2::labs(
      title = paste(title_prefix, ":", gene_name),
      x = "Expression Counts",
      y = "Density",
      fill = "Source",
      color = "Source"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
    )
  
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
