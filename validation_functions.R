#' @title Validation Functions for Dynamic Gene Expression Modelling
#' @description Core validation metrics, accuracy summaries, marker integrity checks, and prediction accuracy tables.
#' @author Taras Lukashiv and collaborators
#' @date 2026-05-20

# =============================================================================
# Functions
# =============================================================================

# -----------------------------------------------------------------------------
# 1. calculate_validation_metrics
# -----------------------------------------------------------------------------

calculate_validation_metrics <- function(sample_data_log,
                                         sample_simul_log,
                                         sample_data_counts,
                                         sample_simul_counts) {
  
  Wilc <- tryCatch(
    wilcox.test(sample_data_counts, sample_simul_counts)$p.value,
    error = function(e) NA_real_
  )
  
  Wilc_round <- tryCatch(
    wilcox.test(sample_data_counts, round(sample_simul_counts))$p.value,
    error = function(e) NA_real_
  )
  
  T_t <- tryCatch(
    t.test(sample_data_log, sample_simul_log)$p.value,
    error = function(e) NA_real_
  )
  
  W_log <- tryCatch(
    transport::wasserstein1d(sample_data_log, sample_simul_log),
    error = function(e) NA_real_
  )
  
  W_orig <- tryCatch(
    transport::wasserstein1d(sample_data_counts, sample_simul_counts),
    error = function(e) NA_real_
  )
  
  mean_data <- mean(sample_data_log, na.rm = TRUE)
  mean_simul <- mean(sample_simul_log, na.rm = TRUE)
  
  median_data_counts <- median(sample_data_counts, na.rm = TRUE)
  median_simul_counts <- median(sample_simul_counts, na.rm = TRUE)
  
  data.frame(
    t_p_val = T_t,
    Wilcox_p_val = Wilc,
    W_log_dist = W_log,
    W_orig_dist = W_orig,
    Converg = 1,
    mean_real_count = mean(sample_data_counts, na.rm = TRUE),
    sd_real_count = sd(sample_data_counts, na.rm = TRUE),
    log2FC = log2((mean_data + 1) / (mean_simul + 1)),
    log2FC_counts = log2((median_data_counts + 1) / (median_simul_counts + 1)),
    mean_pred_count = mean(sample_simul_counts, na.rm = TRUE),
    sd_pred_count = sd(sample_simul_counts, na.rm = TRUE),
    Wilcox_p_val_rounded = Wilc_round,
    mean_diff = mean_data - mean_simul,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# 2. calculate_accuracy_metrics
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# 3. summarize_model_performance
# -----------------------------------------------------------------------------

summarize_model_performance <- function(results_df) {
  
  required_cols <- c(
    "t_p_val",
    "Wilcox_p_val",
    "Wilcox_p_val_rounded",
    "log2FC",
    "log2FC_counts"
  )
  
  missing_cols <- setdiff(required_cols, colnames(results_df))
  if (length(missing_cols) > 0) {
    stop("results_df is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  t_test_log <- calculate_accuracy_metrics(results_df$t_p_val, results_df$log2FC)
  wilcox_raw <- calculate_accuracy_metrics(results_df$Wilcox_p_val, results_df$log2FC_counts)
  wilcox_rounded <- calculate_accuracy_metrics(
    results_df$Wilcox_p_val_rounded,
    results_df$log2FC_counts
  )
  
  data.frame(
    Category = c(
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
}

# -----------------------------------------------------------------------------
# 4. check_marker_integrity
# -----------------------------------------------------------------------------

check_marker_integrity <- function(results_df,
                                   marker_list,
                                   p_adj_method = "BH",
                                   p_threshold = 0.05,
                                   fc_threshold = 1.0) {
  
  required_cols <- c("Names", "t_p_val", "log2FC", "mean_real_count")
  missing_cols <- setdiff(required_cols, colnames(results_df))
  
  if (length(missing_cols) > 0) {
    stop("results_df is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  marker_status <- results_df[results_df$Names %in% marker_list, ]
  
  if (nrow(marker_status) == 0) {
    warning("None of the provided markers were found in the results data frame.")
    return(NULL)
  }
  
  marker_status$adj_p_t <- p.adjust(marker_status$t_p_val, method = p_adj_method)
  marker_status$Status <- ifelse(
    marker_status$adj_p_t >= p_threshold & abs(marker_status$log2FC) < fc_threshold,
    "Success",
    "Fail"
  )
  
  marker_status <- marker_status[, c("Names", "Status", "adj_p_t", "log2FC", "mean_real_count")]
  marker_status[order(marker_status$Status, decreasing = TRUE), ]
}

# -----------------------------------------------------------------------------
# 5. generate_prediction_accuracy_table
# -----------------------------------------------------------------------------

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
