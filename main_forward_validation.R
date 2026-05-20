#' @title Main Forward Validation Pipeline
#' @author Taras Lukashiv, Igor Malyk, Mathias Galati, Ahmed Hemedan and Venkata Satagopam
#' @date 2026-05-06

# =============================================================================
# 00. Clean Environment
# =============================================================================

rm(list = ls())
gc()

# =============================================================================
# 0. Load Modules
# =============================================================================

source("requirements.R")

library(dplyr)

source("simulation_functions.R")
source("parameter_estimation.R")
source("validation_functions.R")
source("extended_validation_functions.R")
source("plotting_functions.R")
source("process_cef_data.R")
source("prediction_functions.R")

# =============================================================================
# 1. Initialization
# =============================================================================

set.seed(125)

output_dir <- "Results_Forward_Validation_New"
if (!dir.exists(output_dir)) dir.create(output_dir)

N_SIM <- 100
n_opt <- 2
VALIDATION_MODE <- "forward"

# =============================================================================
# 2. Data Loading
# =============================================================================

df <- process_cef_data("GSE76381_EmbryoMoleculeCounts.cef.txt")

genes <- read.csv("all_with_sig_identif_new.csv", header = TRUE)
genes <- genes[genes$significant == 1, ]

N_f <- nrow(genes)

# =============================================================================
# 3. Initialize Results Storage
# =============================================================================

df_pred <- data.frame(matrix(NA, nrow = N_f, ncol = N_SIM + 1))
colnames(df_pred) <- c("Gene", paste0("X_", seq_len(N_SIM)))

df_pred_log <- data.frame(matrix(NA, nrow = N_f, ncol = N_SIM + 1))
colnames(df_pred_log) <- c("Gene", paste0("X_", seq_len(N_SIM)))

RES_df <- data.frame(
  Names = character(),
  t_p_val = double(),
  Wilcox_p_val = double(),
  W_log_dist = double(),
  W_orig_dist = double(),
  Converg = double(),
  mean_real_count = double(),
  sd_real_count = double(),
  log2FC = double(),
  log2FC_counts = double(),
  mean_pred_count = double(),
  sd_pred_count = double(),
  Wilcox_p_val_rounded = double(),
  mean_diff = double(),
  stringsAsFactors = FALSE
)

# =============================================================================
# 4. Main Forward Validation Loop using prediction_functions.R
# =============================================================================

df$TP <- as.numeric(substr(df$Timepoint, 6, nchar(df$Timepoint)))
all_tps <- sort(unique(na.omit(df$TP)))

start_tp <- min(all_tps)
test_tp  <- max(all_tps)
end_tp   <- max(all_tps[all_tps < test_tp])

for (g in seq_len(N_f)) {
  
  Correctness <- 1
  gene <- genes[g, 1]
  
  message("Processing gene ", g, " : ", gene)
  
  if (!gene %in% colnames(df)) {
    warning("Gene not found in expression matrix: ", gene)
    next
  }
  
  pred_res <- tryCatch(
    predict_gene_expression(
      df = df,
      gene = gene,
      start_tp = start_tp,
      end_tp = end_tp,
      prediction_tp = test_tp,
      n_opt = n_opt,
      N_SIM = N_SIM,
      N_Steps = 200,
      memory_decay = 1,
      direction = VALIDATION_MODE
    ),
    error = function(e) {
      warning("Prediction failed for gene ", gene, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(pred_res)) next
  
  Sample_data <- log2(as.numeric(df[df$TP == test_tp, gene]) + 1)
  Sample_simul <- pred_res$pred_log
  
  Sample_data_transf <- 2^Sample_data - 1
  Sample_simul_transf <- pred_res$pred_counts
  
  metrics <- calculate_validation_metrics(
    sample_data_log = Sample_data,
    sample_simul_log = Sample_simul,
    sample_data_counts = Sample_data_transf,
    sample_simul_counts = Sample_simul_transf
  )
  
  RES_df[g, 1] <- gene
  RES_df[g, 2] <- metrics$t_p_val
  RES_df[g, 3] <- metrics$Wilcox_p_val
  RES_df[g, 4] <- metrics$W_log_dist
  RES_df[g, 5] <- metrics$W_orig_dist
  RES_df[g, 6] <- Correctness
  RES_df[g, 7] <- metrics$mean_real_count
  RES_df[g, 8] <- metrics$sd_real_count
  RES_df[g, 9] <- metrics$log2FC
  RES_df[g, 10] <- metrics$log2FC_counts
  RES_df[g, 11] <- metrics$mean_pred_count
  RES_df[g, 12] <- metrics$sd_pred_count
  RES_df[g, 13] <- metrics$Wilcox_p_val_rounded
  RES_df[g, 14] <- metrics$mean_diff
  
  df_pred_log[g, 1] <- gene
  df_pred_log[g, 2:(N_SIM + 1)] <- Sample_simul
  
  df_pred[g, 1] <- gene
  df_pred[g, 2:(N_SIM + 1)] <- round(Sample_simul_transf, digits = 0)
}

RES_df <- RES_df[!is.na(RES_df$Names), ]

message("Successfully processed genes: ", nrow(RES_df))

if (nrow(RES_df) == 0) {
  stop("RES_df is empty. Check warnings from predict_gene_expression().")
}

# =============================================================================
# 5. Save Core Results
# =============================================================================

write.csv(
  df_pred,
  file.path(output_dir, "Prediction11_int.csv"),
  row.names = FALSE
)

write.csv(
  df_pred_log,
  file.path(output_dir, "Prediction11_log.csv"),
  row.names = FALSE
)

write.csv(
  RES_df,
  file.path(output_dir, "Validation_results11.csv"),
  row.names = FALSE
)

message("Forward validation complete. Results saved to: ", output_dir)

# =============================================================================
# 6. Failure and Marker Analysis
# =============================================================================

da_markers <- c("TH", "NR4A2", "LMX1A", "PITX3", "FOXA2", "EN1", "LMX1B")

marker_report <- check_marker_integrity(RES_df, da_markers)

if (!is.null(marker_report)) {
  write.csv(
    marker_report,
    file.path(output_dir, "Marker_Integrity_Report.csv"),
    row.names = FALSE
  )
}

performance_summary <- summarize_model_performance(RES_df)
print(performance_summary)

write.csv(
  performance_summary,
  file.path(output_dir, "Performance_Summary.csv"),
  row.names = FALSE
)

# =============================================================================
# 7. Prediction Accuracy Summary Table
# =============================================================================

message("Generating prediction accuracy summary table...")

summary_table <- generate_prediction_accuracy_table(
  results_df = RES_df,
  save_path = file.path(output_dir, "Prediction_Accuracy_Summary.csv"),
  latex_path = file.path(output_dir, "Prediction_Accuracy_Table.tex")
)

message("Prediction Accuracy Summary:")
print(summary_table, row.names = FALSE)

message("Prediction accuracy table saved successfully.")

# =============================================================================
# 8. Visualization and Diagnostics
# =============================================================================

message("Generating validation plots...")

# -----------------------------------------------------------------------------
# 8.1 Distribution of mean differences
# -----------------------------------------------------------------------------

plot_mean_difference_histogram(
  results_df = RES_df,
  save_path = file.path(output_dir, "mean_difference_plot_clean.png")
)

# -----------------------------------------------------------------------------
# 8.2 Scatter plots: adjusted p-values and log2FC
# -----------------------------------------------------------------------------

scatter_plots <- plot_metric_scatter_panels(
  results_df = RES_df,
  save_path_combined = file.path(output_dir, "Scatter_plots_combined.png"),
  save_path_pvals = file.path(output_dir, "p1_scatter_pvals.png"),
  save_path_fc = file.path(output_dir, "p2_scatter_log2fc.png")
)

# -----------------------------------------------------------------------------
# 8.3 Mean-variance / CV relationship
# -----------------------------------------------------------------------------

plot_mean_variance_cv(
  results_df = RES_df,
  save_path = file.path(output_dir, "MeanVariance_CV.png")
)

# -----------------------------------------------------------------------------
# 8.4 Density plot for one representative statistical artifact gene
# -----------------------------------------------------------------------------

plot_df <- RES_df
plot_df$adj_p_t <- p.adjust(plot_df$t_p_val, method = "BH")
plot_df$adj_p_wilcox <- p.adjust(plot_df$Wilcox_p_val, method = "BH")

target_gene_row <- plot_df %>%
  dplyr::filter(
    p.adjust(Wilcox_p_val_rounded, method = "BH") < 0.05,
    abs(log2FC_counts) < 1
  ) %>%
  dplyr::slice(1)

if (nrow(target_gene_row) > 0) {
  target_gene <- target_gene_row$Names
  target_tp <- max(df$TP, na.rm = TRUE)
  
  counts_real <- as.numeric(df[df$TP == target_tp, target_gene])
  counts_pred <- as.numeric(df_pred[df_pred$Gene == target_gene, 2:(N_SIM + 1)])
  
  plot_real_predicted_density(
    gene = target_gene,
    real_counts = counts_real,
    pred_counts = counts_pred,
    save_path = file.path(output_dir, paste0("Density_", target_gene, ".png")),
    title_prefix = "Distribution Comparison"
  )
}

# -----------------------------------------------------------------------------
# 8.5 Failure-mode analysis
# -----------------------------------------------------------------------------

failing_genes <- plot_df %>%
  dplyr::filter(adj_p_t < 0.05 & abs(log2FC) >= 1) %>%
  dplyr::arrange(dplyr::desc(abs(log2FC)))

write.csv(
  failing_genes,
  file.path(output_dir, "S7_Ranked_Failure_Genes.csv"),
  row.names = FALSE
)

marker_status <- plot_df %>%
  dplyr::filter(Names %in% da_markers) %>%
  dplyr::select(Names, mean_real_count, adj_p_t, log2FC) %>%
  dplyr::mutate(Status = ifelse(adj_p_t >= 0.05 & abs(log2FC) < 1, "Success", "Fail"))

write.csv(
  marker_status,
  file.path(output_dir, "DA_Marker_Validation_Status.csv"),
  row.names = FALSE
)

print(marker_status)

# -----------------------------------------------------------------------------
# 8.6 Stratified accuracy by expression decile
# -----------------------------------------------------------------------------

stratified_accuracy <- plot_stratified_accuracy(
  results_df = plot_df,
  save_path = file.path(output_dir, "Stratified_Accuracy_Deciles.png"),
  n_bins = 10
)

write.csv(
  stratified_accuracy$stats,
  file.path(output_dir, "Stratified_Accuracy_Stats.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 8.7 Volcano plot and full failure list
# -----------------------------------------------------------------------------

volcano_result <- plot_volcano_failures(
  results_df = plot_df,
  save_path = file.path(output_dir, "Volcano_Fail_Analysis.png"),
  label_top_n = 10
)

write.csv(
  volcano_result$fail_list,
  file.path(output_dir, "Full_Fail_Genes_Analysis.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 8.8 Density plots for biological failure examples
# -----------------------------------------------------------------------------

fail_examples <- c("COL4A1", "FN1", "SLC1A3")
target_tp <- max(df$TP, na.rm = TRUE)

for (target_gene in fail_examples) {
  
  if (target_gene %in% colnames(df) && target_gene %in% df_pred$Gene) {
    counts_real <- as.numeric(df[df$TP == target_tp, target_gene])
    counts_pred <- as.numeric(df_pred[df_pred$Gene == target_gene, 2:(N_SIM + 1)])
    
    plot_real_predicted_density(
      gene = target_gene,
      real_counts = counts_real,
      pred_counts = counts_pred,
      save_path = file.path(output_dir, paste0("Supp_Fail_", target_gene, ".png")),
      title_prefix = "Non-linear Expression Shift"
    )
  }
}

message("Validation analysis complete. All plots saved to: ", output_dir)
