#' @title Extended Validation Functions for Dynamic Gene Expression Modelling
#' @description Extended genome-scale validation utilities, stratified deviation analysis, DE agreement, and calibration checks.
#' @author Taras Lukashiv and collaborators
#' @date 2026-05-20

# =============================================================================
# Functions
# =============================================================================

# -----------------------------------------------------------------------------
# 1. generate_horizon_accuracy_table
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# 2. analyze_stratified_deviations
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# 3. select_representative_genes
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# 4. calculate_per_gene_DE_counts
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# 5. compare_DE_gene_sets
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# 6. evaluate_prediction_calibration
# -----------------------------------------------------------------------------

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
