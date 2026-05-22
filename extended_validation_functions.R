#' @title Extended Genome-scale Validation Functions
#' @description Additional manuscript-grade validation utilities:
#' intermediate-horizon accuracy tables, density comparisons,
#' stratified deviation distributions, and per-gene DE agreement.
#' @author Taras Lukashiv and collaborators
#' @date 2026-05-06


# =============================================================================
# 1. Accuracy Metrics Helper
# =============================================================================

calculate_accuracy_metrics <- function(p_val_vec, fc_vec, p_method = "BH") {
  
  adj_p <- p.adjust(p_val_vec, method = p_method)
  abs_fc <- abs(fc_vec)
  
  out <- c(
    Success = mean(adj_p >= 0.05 & abs_fc < 1, na.rm = TRUE) * 100,
    Small_Diff_Signif = mean(adj_p < 0.05 & abs_fc < 1, na.rm = TRUE) * 100,
    Fail = mean(adj_p < 0.05 & abs_fc >= 1, na.rm = TRUE) * 100,
    Large_Diff_Insignif = mean(adj_p >= 0.05 & abs_fc >= 1, na.rm = TRUE) * 100
  )
  
  round(out, 2)
}


# =============================================================================
# 2. Prediction Accuracy Table
# =============================================================================

generate_prediction_accuracy_table <- function(results_df,
                                               save_path = NULL,
                                               latex_path = NULL,
                                               caption = "Comparative Analysis of Prediction Accuracy Metrics",
                                               label = "tab:prediction-accuracy") {
  
  required_cols <- c(
    "t_p_val",
    "Wilcox_p_val",
    "Wilcox_p_val_rounded",
    "log2FC",
    "log2FC_counts"
  )
  
  missing_cols <- setdiff(required_cols, colnames(results_df))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  t_test_log <- calculate_accuracy_metrics(
    results_df$t_p_val,
    results_df$log2FC
  )
  
  wilcox_raw <- calculate_accuracy_metrics(
    results_df$Wilcox_p_val,
    results_df$log2FC_counts
  )
  
  wilcox_rounded <- calculate_accuracy_metrics(
    results_df$Wilcox_p_val_rounded,
    results_df$log2FC_counts
  )
  
  summary_table <- data.frame(
    Category_Criteria = c(
      "Success (p_adj >= 0.05, |log2FC| < 1)",
      "Small Diff, Signif (p_adj < 0.05, |log2FC| < 1)",
      "Fail (p_adj < 0.05, |log2FC| >= 1)",
      "Large Diff, Insignif (p_adj >= 0.05, |log2FC| >= 1)"
    ),
    T_test_Log = paste0(t_test_log, "%"),
    Wilcoxon_Raw = paste0(wilcox_raw, "%"),
    Wilcoxon_Rounded = paste0(wilcox_rounded, "%"),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(save_path)) {
    write.csv(summary_table, save_path, row.names = FALSE)
  }
  
  if (!is.null(latex_path)) {
    
    latex_table <- paste0(
      "\\begin{table}[H]\n",
      "\\centering\n",
      "\\caption{", caption, "}\n",
      "\\label{", label, "}\n",
      "\\begin{tabular}{|l|c|c|c|}\n",
      "\\hline\n",
      "\\textbf{Category (Criteria)} & ",
      "\\begin{tabular}[c]{@{}c@{}}\\textbf{T-test} \\\\ \\textbf{(Log-data)}\\end{tabular} & ",
      "\\begin{tabular}[c]{@{}c@{}}\\textbf{Wilcoxon} \\\\ \\textbf{(Raw)}\\end{tabular} & ",
      "\\begin{tabular}[c]{@{}c@{}}\\textbf{Wilcoxon} \\\\ \\textbf{(Rounded)}\\end{tabular} \\\\ \\hline \n",
      
      "Success ($p_{adj} \\\\ge 0.05, |\\\\log_2 \\\\mathrm{FC}| < 1$) & ",
      "\\textbf{", summary_table$T_test_Log[1], "} & ",
      summary_table$Wilcoxon_Raw[1], " & ",
      summary_table$Wilcoxon_Rounded[1], " \\\\ \\hline \n",
      
      "Small Diff, Signif ($p_{adj} < 0.05, |\\\\log_2 \\\\mathrm{FC}| < 1$) & ",
      summary_table$T_test_Log[2], " & ",
      summary_table$Wilcoxon_Raw[2], " & ",
      "\\textbf{", summary_table$Wilcoxon_Rounded[2], "} \\\\ \\hline \n",
      
      "Fail ($p_{adj} < 0.05, |\\\\log_2 \\\\mathrm{FC}| \\\\ge 1$) & ",
      summary_table$T_test_Log[3], " & ",
      summary_table$Wilcoxon_Raw[3], " & ",
      summary_table$Wilcoxon_Rounded[3], " \\\\ \\hline \n",
      
      "Large Diff, Insignif ($p_{adj} \\\\ge 0.05, |\\\\log_2 \\\\mathrm{FC}| \\\\ge 1$) & ",
      summary_table$T_test_Log[4], " & ",
      summary_table$Wilcoxon_Raw[4], " & ",
      summary_table$Wilcoxon_Rounded[4], " \\\\ \\hline\n",
      "\\end{tabular}\n",
      "\\end{table}\n"
    )
    
    writeLines(latex_table, con = latex_path)
  }
  
  return(summary_table)
}


# =============================================================================
# 3. Intermediate-horizon Accuracy Table
# =============================================================================

generate_horizon_accuracy_table <- function(results_list,
                                            horizon_labels = NULL,
                                            save_path = NULL) {
  
  if (is.null(names(results_list)) && is.null(horizon_labels)) {
    horizon_labels <- paste0("Horizon_", seq_along(results_list))
  }
  
  if (is.null(horizon_labels)) {
    horizon_labels <- names(results_list)
  }
  
  out <- data.frame()
  
  for (i in seq_along(results_list)) {
    
    res <- results_list[[i]]
    
    acc_t <- calculate_accuracy_metrics(res$t_p_val, res$log2FC)
    acc_w <- calculate_accuracy_metrics(res$Wilcox_p_val, res$log2FC_counts)
    acc_wr <- calculate_accuracy_metrics(res$Wilcox_p_val_rounded, res$log2FC_counts)
    
    temp <- data.frame(
      Horizon = horizon_labels[i],
      
      T_test_Success = acc_t["Success"],
      T_test_Small_Diff_Signif = acc_t["Small_Diff_Signif"],
      T_test_Fail = acc_t["Fail"],
      T_test_Large_Diff_Insignif = acc_t["Large_Diff_Insignif"],
      
      Wilcoxon_Raw_Success = acc_w["Success"],
      Wilcoxon_Raw_Small_Diff_Signif = acc_w["Small_Diff_Signif"],
      Wilcoxon_Raw_Fail = acc_w["Fail"],
      Wilcoxon_Raw_Large_Diff_Insignif = acc_w["Large_Diff_Insignif"],
      
      Wilcoxon_Rounded_Success = acc_wr["Success"],
      Wilcoxon_Rounded_Small_Diff_Signif = acc_wr["Small_Diff_Signif"],
      Wilcoxon_Rounded_Fail = acc_wr["Fail"],
      Wilcoxon_Rounded_Large_Diff_Insignif = acc_wr["Large_Diff_Insignif"],
      
      stringsAsFactors = FALSE
    )
    
    out <- rbind(out, temp)
  }
  
  if (!is.null(save_path)) {
    write.csv(out, save_path, row.names = FALSE)
  }
  
  return(out)
}


# =============================================================================
# 4. Stratified Deviation Distributions
# =============================================================================

analyze_stratified_deviations <- function(results_df,
                                          n_bins = 10,
                                          stratify_by = c("mean_real_count", "cv_real", "mean_pred_count"),
                                          save_path = NULL) {
  
  stratify_by <- match.arg(stratify_by)
  
  if (!stratify_by %in% colnames(results_df)) {
    stop("Column not found in results_df: ", stratify_by)
  }
  
  df <- results_df
  
  if (!"cv_real" %in% colnames(df)) {
    df$cv_real <- df$sd_real_count / df$mean_real_count
  }
  
  if (!"cv_pred" %in% colnames(df)) {
    df$cv_pred <- df$sd_pred_count / df$mean_pred_count
  }
  
  df$abs_log2FC <- abs(df$log2FC)
  df$abs_log2FC_counts <- abs(df$log2FC_counts)
  df$abs_mean_diff <- abs(df$mean_diff)
  df$adj_p_t <- p.adjust(df$t_p_val, method = "BH")
  df$adj_p_wilcox <- p.adjust(df$Wilcox_p_val, method = "BH")
  
  df <- df[is.finite(df[[stratify_by]]), ]
  
  df$stratum <- dplyr::ntile(df[[stratify_by]], n_bins)
  
  strat_summary <- df %>%
    dplyr::group_by(stratum) %>%
    dplyr::summarise(
      n_genes = dplyr::n(),
      mean_stratifier = mean(.data[[stratify_by]], na.rm = TRUE),
      median_stratifier = median(.data[[stratify_by]], na.rm = TRUE),
      mean_abs_log2FC = mean(abs_log2FC, na.rm = TRUE),
      median_abs_log2FC = median(abs_log2FC, na.rm = TRUE),
      mean_abs_log2FC_counts = mean(abs_log2FC_counts, na.rm = TRUE),
      median_abs_log2FC_counts = median(abs_log2FC_counts, na.rm = TRUE),
      mean_abs_mean_diff = mean(abs_mean_diff, na.rm = TRUE),
      success_rate = mean(adj_p_t >= 0.05 & abs_log2FC < 1, na.rm = TRUE) * 100,
      fail_rate = mean(adj_p_t < 0.05 & abs_log2FC >= 1, na.rm = TRUE) * 100,
      .groups = "drop"
    )
  
  if (!is.null(save_path)) {
    write.csv(strat_summary, save_path, row.names = FALSE)
  }
  
  return(list(
    data = df,
    summary = strat_summary
  ))
}


# =============================================================================
# 5. Plot Stratified Deviation Distributions
# =============================================================================

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


# =============================================================================
# 6. Extended Density Comparisons
# =============================================================================

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


# =============================================================================
# 7. Select Representative Genes
# =============================================================================

select_representative_genes <- function(results_df,
                                        n_success = 5,
                                        n_failure = 5,
                                        n_small_diff = 5) {
  
  df <- results_df
  df$adj_p_t <- p.adjust(df$t_p_val, method = "BH")
  df$abs_log2FC <- abs(df$log2FC)
  
  success_genes <- df %>%
    dplyr::filter(adj_p_t >= 0.05, abs_log2FC < 1) %>%
    dplyr::arrange(abs_log2FC) %>%
    dplyr::slice_head(n = n_success)
  
  failure_genes <- df %>%
    dplyr::filter(adj_p_t < 0.05, abs_log2FC >= 1) %>%
    dplyr::arrange(dplyr::desc(abs_log2FC)) %>%
    dplyr::slice_head(n = n_failure)
  
  small_diff_genes <- df %>%
    dplyr::filter(adj_p_t < 0.05, abs_log2FC < 1) %>%
    dplyr::arrange(abs_log2FC) %>%
    dplyr::slice_head(n = n_small_diff)
  
  return(list(
    success = success_genes,
    failure = failure_genes,
    small_diff_significant = small_diff_genes
  ))
}


# =============================================================================
# 8. Per-gene DE Counts and Agreement
# =============================================================================

calculate_per_gene_DE_counts <- function(results_df,
                                         p_col = "t_p_val",
                                         fc_col = "log2FC",
                                         p_threshold = 0.05,
                                         fc_threshold = 1,
                                         p_method = "BH",
                                         save_path = NULL) {
  
  if (!p_col %in% colnames(results_df)) {
    stop("p_col not found in results_df: ", p_col)
  }
  
  if (!fc_col %in% colnames(results_df)) {
    stop("fc_col not found in results_df: ", fc_col)
  }
  
  df <- results_df
  
  df$adj_p <- p.adjust(df[[p_col]], method = p_method)
  df$abs_fc <- abs(df[[fc_col]])
  
  df$DE_status <- dplyr::case_when(
    df$adj_p < p_threshold & df$abs_fc >= fc_threshold ~ "Strong_DE_disagreement",
    df$adj_p < p_threshold & df$abs_fc < fc_threshold ~ "Statistical_only_difference",
    df$adj_p >= p_threshold & df$abs_fc >= fc_threshold ~ "Biological_only_difference",
    df$adj_p >= p_threshold & df$abs_fc < fc_threshold ~ "Concordant_non_DE",
    TRUE ~ NA_character_
  )
  
  summary_counts <- df %>%
    dplyr::count(DE_status, name = "n_genes") %>%
    dplyr::mutate(percentage = round(n_genes / sum(n_genes, na.rm = TRUE) * 100, 2))
  
  if (!is.null(save_path)) {
    write.csv(summary_counts, save_path, row.names = FALSE)
  }
  
  return(list(
    per_gene = df,
    summary = summary_counts
  ))
}


# =============================================================================
# 9. Compare Observed vs Predicted DE Gene Sets
# =============================================================================

compare_DE_gene_sets <- function(observed_de_genes,
                                 predicted_de_genes) {
  
  observed_de_genes <- unique(observed_de_genes)
  predicted_de_genes <- unique(predicted_de_genes)
  
  tp <- length(intersect(observed_de_genes, predicted_de_genes))
  fp <- length(setdiff(predicted_de_genes, observed_de_genes))
  fn <- length(setdiff(observed_de_genes, predicted_de_genes))
  
  precision <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
  recall <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
  f1 <- ifelse(
    is.na(precision) | is.na(recall) | (precision + recall) == 0,
    NA_real_,
    2 * precision * recall / (precision + recall)
  )
  
  jaccard <- ifelse(
    length(union(observed_de_genes, predicted_de_genes)) == 0,
    NA_real_,
    length(intersect(observed_de_genes, predicted_de_genes)) /
      length(union(observed_de_genes, predicted_de_genes))
  )
  
  data.frame(
    Observed_DE = length(observed_de_genes),
    Predicted_DE = length(predicted_de_genes),
    True_Positive = tp,
    False_Positive = fp,
    False_Negative = fn,
    Precision = precision,
    Recall = recall,
    F1 = f1,
    Jaccard = jaccard
  )
}


# =============================================================================
# 10. Prediction Interval Calibration
# =============================================================================

evaluate_prediction_calibration <- function(real_counts,
                                            prediction_matrix,
                                            lower_prob = 0.025,
                                            upper_prob = 0.975) {
  
  lower <- apply(prediction_matrix, 1, quantile, probs = lower_prob, na.rm = TRUE)
  upper <- apply(prediction_matrix, 1, quantile, probs = upper_prob, na.rm = TRUE)
  pred_mean <- apply(prediction_matrix, 1, mean, na.rm = TRUE)
  
  covered <- real_counts >= lower & real_counts <= upper
  
  data.frame(
    coverage = mean(covered, na.rm = TRUE) * 100,
    mean_interval_width = mean(upper - lower, na.rm = TRUE),
    mean_prediction = mean(pred_mean, na.rm = TRUE),
    mean_real = mean(real_counts, na.rm = TRUE)
  )
}