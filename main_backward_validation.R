#' @title Main Backward Validation Pipeline
#' @description Backward validation / backcasting pipeline for dynamic gene expression modelling.
#' @author Taras Lukashiv, Igor Malyk, Mathias Galati, Ahmed Hemedan and Venkata Satagopam
#' @date 2026-05-20

# =============================================================================
# 0. Clean Environment
# =============================================================================

rm(list = ls())
gc()

# =============================================================================
# 1. Load Modules
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
# 2. Initialization
# =============================================================================

set.seed(125)

output_dir <- "Results_Backward_Validation_New"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

N_SIM <- 100
n_opt <- 2
N_Steps <- 200
memory_decay <- 1
VALIDATION_MODE <- "backward"

# =============================================================================
# 3. Data Loading
# =============================================================================

message("Loading expression data...")

df <- process_cef_data("GSE76381_EmbryoMoleculeCounts.cef.txt")

if (!"Timepoint" %in% colnames(df)) {
  stop("Input data must contain a Timepoint column.")
}

df$TP <- suppressWarnings(as.numeric(substr(df$Timepoint, 6, nchar(df$Timepoint))))

if (any(is.na(df$TP))) {
  stop("Some TP values could not be extracted from Timepoint.")
}

all_tps <- sort(unique(na.omit(df$TP)))

if (length(all_tps) < 3) {
  stop("At least three timepoints are required for backward validation.")
}

message("Available timepoints: ", paste(all_tps, collapse = ", "))

message("Loading significant gene list...")

genes <- read.csv("all_with_sig_identif_new.csv", header = TRUE)

if (!"significant" %in% colnames(genes)) {
  stop("Gene list must contain a 'significant' column.")
}

genes <- genes[genes$significant == 1, ]

gene_names <- genes[, 1]
gene_names <- gene_names[gene_names %in% colnames(df)]

N_f <- length(gene_names)

if (N_f == 0) {
  stop("No significant genes from the gene list were found in the expression matrix.")
}

message("Number of significant genes used: ", N_f)

# =============================================================================
# 4. Backward Validation Setup
# =============================================================================

# Backcasting target:
# - train from the latest observed week backward to the second observed week;
# - predict the earliest observed week.
# Example for weeks 6–11: train on 11 -> 7 and predict 6.

test_tp  <- min(all_tps)                       # target, e.g. Week 6
start_tp <- max(all_tps)                       # latest observed timepoint, e.g. Week 11
end_tp   <- min(all_tps[all_tps > test_tp])    # next timepoint after target, e.g. Week 7

message("Backward validation setup:")
message("  start_tp     = ", start_tp)
message("  end_tp       = ", end_tp)
message("  prediction_tp = ", test_tp)

# =============================================================================
# 5. Initialize Results Storage
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
# 6. Main Backward Validation Loop
# =============================================================================

for (g in seq_len(N_f)) {
  
  gene <- gene_names[g]
  
  message("[Backward] Processing gene ", g, "/", N_f, " : ", gene)
  
  pred_res <- tryCatch(
    predict_gene_expression(
      df = df,
      gene = gene,
      start_tp = start_tp,
      end_tp = end_tp,
      prediction_tp = test_tp,
      n_opt = n_opt,
      N_SIM = N_SIM,
      N_Steps = N_Steps,
      memory_decay = memory_decay,
      direction = VALIDATION_MODE
    ),
    error = function(e) {
      warning("Prediction failed for gene ", gene, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(pred_res)) next
  
  sample_data_log <- log2(as.numeric(df[df$TP == test_tp, gene]) + 1)
  sample_simul_log <- pred_res$pred_log
  
  sample_data_counts <- 2^sample_data_log - 1
  sample_simul_counts <- pred_res$pred_counts
  
  metrics <- calculate_validation_metrics(
    sample_data_log = sample_data_log,
    sample_simul_log = sample_simul_log,
    sample_data_counts = sample_data_counts,
    sample_simul_counts = sample_simul_counts
  )
  
  RES_df[g, ] <- data.frame(
    Names = gene,
    t_p_val = metrics$t_p_val,
    Wilcox_p_val = metrics$Wilcox_p_val,
    W_log_dist = metrics$W_log_dist,
    W_orig_dist = metrics$W_orig_dist,
    Converg = 1,
    mean_real_count = metrics$mean_real_count,
    sd_real_count = metrics$sd_real_count,
    log2FC = metrics$log2FC,
    log2FC_counts = metrics$log2FC_counts,
    mean_pred_count = metrics$mean_pred_count,
    sd_pred_count = metrics$sd_pred_count,
    Wilcox_p_val_rounded = metrics$Wilcox_p_val_rounded,
    mean_diff = metrics$mean_diff,
    stringsAsFactors = FALSE
  )
  
  df_pred_log[g, 1] <- gene
  df_pred_log[g, 2:(N_SIM + 1)] <- sample_simul_log
  
  df_pred[g, 1] <- gene
  df_pred[g, 2:(N_SIM + 1)] <- round(sample_simul_counts, digits = 0)
}

RES_df <- RES_df[!is.na(RES_df$Names), ]
df_pred <- df_pred[!is.na(df_pred$Gene), ]
df_pred_log <- df_pred_log[!is.na(df_pred_log$Gene), ]

message("Successfully processed genes: ", nrow(RES_df))

if (nrow(RES_df) == 0) {
  stop("RES_df is empty. Check warnings from predict_gene_expression().")
}

# =============================================================================
# 7. Save Core Results
# =============================================================================

write.csv(
  df_pred,
  file.path(output_dir, paste0("Prediction", test_tp, "_backward_int.csv")),
  row.names = FALSE
)

write.csv(
  df_pred_log,
  file.path(output_dir, paste0("Prediction", test_tp, "_backward_log.csv")),
  row.names = FALSE
)

write.csv(
  RES_df,
  file.path(output_dir, paste0("Validation_results", test_tp, "_backward.csv")),
  row.names = FALSE
)

message("Backward validation core results saved to: ", output_dir)

# =============================================================================
# 8. Marker and Performance Analysis
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
# 9. Prediction Accuracy Summary Table
# =============================================================================

message("Generating prediction accuracy summary table...")

summary_table <- generate_prediction_accuracy_table(
  results_df = RES_df,
  save_path = file.path(output_dir, "Prediction_Accuracy_Summary.csv"),
  latex_path = file.path(output_dir, "Prediction_Accuracy_Table.tex"),
  caption = "Backward Validation Accuracy Metrics",
  label = "tab:backward-validation-accuracy"
)

message("Prediction Accuracy Summary:")
print(summary_table, row.names = FALSE)

# =============================================================================
# 10. Visualization and Diagnostics
# =============================================================================

message("Generating backward validation plots...")

# -----------------------------------------------------------------------------
# 10.1 Distribution of mean differences
# -----------------------------------------------------------------------------

plot_mean_difference_histogram(
  results_df = RES_df,
  save_path = file.path(output_dir, "mean_difference_plot_clean.png")
)

# -----------------------------------------------------------------------------
# 10.2 Scatter plots: adjusted p-values and log2FC
# -----------------------------------------------------------------------------

scatter_plots <- plot_metric_scatter_panels(
  results_df = RES_df,
  save_path_combined = file.path(output_dir, "Scatter_plots_combined.png"),
  save_path_pvals = file.path(output_dir, "p1_scatter_pvals.png"),
  save_path_fc = file.path(output_dir, "p2_scatter_log2fc.png")
)

# -----------------------------------------------------------------------------
# 10.3 Mean-variance / CV relationship
# -----------------------------------------------------------------------------

plot_mean_variance_cv(
  results_df = RES_df,
  save_path = file.path(output_dir, "MeanVariance_CV.png")
)

# -----------------------------------------------------------------------------
# 10.4 Density plot for one representative statistical artifact gene
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
  target_tp <- test_tp
  
  counts_real <- as.numeric(df[df$TP == target_tp, target_gene])
  counts_pred <- as.numeric(df_pred[df_pred$Gene == target_gene, 2:(N_SIM + 1)])
  
  plot_real_predicted_density(
    gene = target_gene,
    real_counts = counts_real,
    pred_counts = counts_pred,
    save_path = file.path(output_dir, paste0("Density_", target_gene, ".png")),
    title_prefix = "Backward Distribution Comparison"
  )
}

# -----------------------------------------------------------------------------
# 10.5 Failure-mode analysis
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
# 10.6 Stratified accuracy by expression decile
# -----------------------------------------------------------------------------

stratified_accuracy <- plot_stratified_accuracy(
  results_df = RES_df,
  save_path = file.path(output_dir, "Stratified_Accuracy_Deciles.png"),
  n_bins = 10
)

write.csv(
  stratified_accuracy$stats,
  file.path(output_dir, "Stratified_Accuracy_Stats.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 10.7 Volcano plot
# -----------------------------------------------------------------------------

volcano_res <- plot_volcano_failures(
  results_df = RES_df,
  save_path = file.path(output_dir, "Volcano_Fail_Analysis.png"),
  label_top_n = 10
)

write.csv(
  volcano_res$fail_list,
  file.path(output_dir, "Full_Fail_Genes_Analysis.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 10.8 Density plots for selected biological failure examples
# -----------------------------------------------------------------------------

fail_examples <- c("COL4A1", "FN1", "SLC1A3")
target_tp <- test_tp

for (target_gene in fail_examples) {
  
  if (target_gene %in% colnames(df) && target_gene %in% df_pred$Gene) {
    
    counts_real <- as.numeric(df[df$TP == target_tp, target_gene])
    counts_pred <- as.numeric(df_pred[df_pred$Gene == target_gene, 2:(N_SIM + 1)])
    
    plot_real_predicted_density(
      gene = target_gene,
      real_counts = counts_real,
      pred_counts = counts_pred,
      save_path = file.path(output_dir, paste0("Supp_Fail_", target_gene, ".png")),
      title_prefix = "Backward Failure Example"
    )
  }
}

message("Backward validation analysis complete. All outputs saved to: ", output_dir)
