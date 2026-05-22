#' @title Genome-scale Validation Extensions
#' @description Intermediate-horizon accuracy tables, extended density comparisons,
#' stratified deviation distributions, and per-gene DE counts.
#' @author Taras Lukashiv, Igor Malyk, Mathias Galati, Ahmed Hemedan and Venkata Satagopam
#' @date 2026-05-06

# =============================================================================
# 0. Load Modules
# =============================================================================

source("requirements.R")

library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)

source("simulation_functions.R")
source("parameter_estimation.R")
source("validation_functions.R")
source("process_cef_data.R")
source("prediction_functions.R")
source("extended_validation_functions.R")
source("plotting_functions.R")
source("trajectory_analysis.R")

# =============================================================================
# 1. Initialization
# =============================================================================

set.seed(125)

output_dir <- "Results_Genome_Scale_Extensions"
if (!dir.exists(output_dir)) dir.create(output_dir)

N_SIM <- 100
n_opt <- 2
N_Steps <- 200
memory_decay <- 1

# =============================================================================
# 2. Data Loading
# =============================================================================

df <- process_cef_data("GSE76381_EmbryoMoleculeCounts.cef.txt")
df$TP <- as.numeric(substr(df$Timepoint, 6, nchar(df$Timepoint)))

all_tps <- sort(unique(na.omit(df$TP)))

genes <- read.csv("all_with_sig_identif_new.csv", header = TRUE)
genes <- genes[genes$significant == 1, ]
gene_names <- genes[, 1]
gene_names <- gene_names[gene_names %in% colnames(df)]

N_f <- length(gene_names)

message("Number of significant genes used: ", N_f)
message("Available timepoints: ", paste(all_tps, collapse = ", "))

# =============================================================================
# 3. Helper: Run Validation for One Horizon
# =============================================================================

run_horizon_validation <- function(df,
                                   gene_names,
                                   start_tp,
                                   end_tp,
                                   prediction_tp,
                                   direction,
                                   output_subdir) {
  
  if (!dir.exists(output_subdir)) dir.create(output_subdir, recursive = TRUE)
  
  N_f <- length(gene_names)
  
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
  
  for (g in seq_along(gene_names)) {
    
    gene <- gene_names[g]
    message("[", direction, "] ", start_tp, " -> ", prediction_tp,
            " | Processing gene ", g, "/", N_f, " : ", gene)
    
    pred_res <- tryCatch(
      predict_gene_expression(
        df = df,
        gene = gene,
        start_tp = start_tp,
        end_tp = end_tp,
        prediction_tp = prediction_tp,
        n_opt = n_opt,
        N_SIM = N_SIM,
        N_Steps = N_Steps,
        memory_decay = memory_decay,
        direction = direction
      ),
      error = function(e) {
        warning("Prediction failed for gene ", gene, ": ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(pred_res)) next
    
    sample_data_log <- log2(as.numeric(df[df$TP == prediction_tp, gene]) + 1)
    sample_simul_log <- pred_res$pred_log
    
    sample_data_counts <- 2^sample_data_log - 1
    sample_simul_counts <- pred_res$pred_counts
    
    metrics <- calculate_validation_metrics(
      sample_data_log = sample_data_log,
      sample_simul_log = sample_simul_log,
      sample_data_counts = sample_data_counts,
      sample_simul_counts = sample_simul_counts
    )
    
    RES_df[g, 1] <- gene
    RES_df[g, 2] <- metrics$t_p_val
    RES_df[g, 3] <- metrics$Wilcox_p_val
    RES_df[g, 4] <- metrics$W_log_dist
    RES_df[g, 5] <- metrics$W_orig_dist
    RES_df[g, 6] <- 1
    RES_df[g, 7] <- metrics$mean_real_count
    RES_df[g, 8] <- metrics$sd_real_count
    RES_df[g, 9] <- metrics$log2FC
    RES_df[g, 10] <- metrics$log2FC_counts
    RES_df[g, 11] <- metrics$mean_pred_count
    RES_df[g, 12] <- metrics$sd_pred_count
    RES_df[g, 13] <- metrics$Wilcox_p_val_rounded
    RES_df[g, 14] <- metrics$mean_diff
    
    df_pred_log[g, 1] <- gene
    df_pred_log[g, 2:(N_SIM + 1)] <- sample_simul_log
    
    df_pred[g, 1] <- gene
    df_pred[g, 2:(N_SIM + 1)] <- round(sample_simul_counts, digits = 0)
  }
  
  RES_df <- RES_df[!is.na(RES_df$Names), ]
  df_pred <- df_pred[!is.na(df_pred$Gene), ]
  df_pred_log <- df_pred_log[!is.na(df_pred_log$Gene), ]
  
  write.csv(
    RES_df,
    file.path(output_subdir, paste0("Validation_results_", direction, "_", start_tp, "_to_", prediction_tp, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    df_pred,
    file.path(output_subdir, paste0("Prediction_counts_", direction, "_", start_tp, "_to_", prediction_tp, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    df_pred_log,
    file.path(output_subdir, paste0("Prediction_log_", direction, "_", start_tp, "_to_", prediction_tp, ".csv")),
    row.names = FALSE
  )
  
  accuracy_table <- generate_prediction_accuracy_table(
    results_df = RES_df,
    save_path = file.path(output_subdir, paste0("Accuracy_summary_", direction, "_", start_tp, "_to_", prediction_tp, ".csv")),
    latex_path = file.path(output_subdir, paste0("Accuracy_summary_", direction, "_", start_tp, "_to_", prediction_tp, ".tex")),
    caption = paste0("Prediction Accuracy Metrics (", direction, ": Week ", start_tp, " to Week ", prediction_tp, ")"),
    label = paste0("tab:accuracy_", direction, "_", start_tp, "_to_", prediction_tp)
  )
  
  return(list(
    results = RES_df,
    df_pred = df_pred,
    df_pred_log = df_pred_log,
    accuracy_table = accuracy_table
  ))
}

# =============================================================================
# 4. Intermediate-horizon Accuracy Tables
# =============================================================================

message("Running intermediate-horizon validation...")

forward_results <- list()
backward_results <- list()

forward_start <- min(all_tps)

# Forward validation requires at least two observed training timepoints.
# For weeks 6--11, this starts at target Week 8: train on 6--7, predict 8.
forward_targets <- all_tps[seq_along(all_tps) >= 3]

for (target_tp in forward_targets) {
  
  end_tp <- max(all_tps[all_tps < target_tp])
  
  horizon_label <- paste0("Week_", forward_start, "_to_", target_tp)
  
  forward_results[[horizon_label]] <- run_horizon_validation(
    df = df,
    gene_names = gene_names,
    start_tp = forward_start,
    end_tp = end_tp,
    prediction_tp = target_tp,
    direction = "forward",
    output_subdir = file.path(output_dir, "Forward_Horizons", horizon_label)
  )
}

backward_start <- max(all_tps)

# Backward validation also requires at least two observed training timepoints.
# For weeks 6--11, this starts at target Week 9: train on 11--10, predict 9.
backward_targets <- rev(all_tps[seq_along(all_tps) <= (length(all_tps) - 2)])

for (target_tp in backward_targets) {
  
  end_tp <- min(all_tps[all_tps > target_tp])
  
  horizon_label <- paste0("Week_", backward_start, "_to_", target_tp)
  
  backward_results[[horizon_label]] <- run_horizon_validation(
    df = df,
    gene_names = gene_names,
    start_tp = backward_start,
    end_tp = end_tp,
    prediction_tp = target_tp,
    direction = "backward",
    output_subdir = file.path(output_dir, "Backward_Horizons", horizon_label)
  )
}

forward_RES_list <- lapply(forward_results, function(x) x$results)
backward_RES_list <- lapply(backward_results, function(x) x$results)

forward_horizon_table <- generate_horizon_accuracy_table(
  results_list = forward_RES_list,
  horizon_labels = names(forward_RES_list),
  save_path = file.path(output_dir, "Forward_Intermediate_Horizon_Accuracy.csv")
)

backward_horizon_table <- generate_horizon_accuracy_table(
  results_list = backward_RES_list,
  horizon_labels = names(backward_RES_list),
  save_path = file.path(output_dir, "Backward_Intermediate_Horizon_Accuracy.csv")
)

print(forward_horizon_table)
print(backward_horizon_table)

# =============================================================================
# 5. Prediction Decay Analysis
# =============================================================================

forward_decay <- analyze_prediction_decay(
  horizon_results = forward_RES_list,
  horizon_labels = names(forward_RES_list),
  save_path = file.path(output_dir, "Forward_Prediction_Decay.csv")
)

backward_decay <- analyze_prediction_decay(
  horizon_results = backward_RES_list,
  horizon_labels = names(backward_RES_list),
  save_path = file.path(output_dir, "Backward_Prediction_Decay.csv")
)

plot_prediction_decay(
  forward_decay,
  save_path = file.path(output_dir, "Forward_Prediction_Decay_Success.png"),
  metric = "Success_Rate"
)

plot_prediction_decay(
  backward_decay,
  save_path = file.path(output_dir, "Backward_Prediction_Decay_Success.png"),
  metric = "Success_Rate"
)

plot_prediction_decay(
  forward_decay,
  save_path = file.path(output_dir, "Forward_Prediction_Decay_Fail.png"),
  metric = "Fail_Rate"
)

plot_prediction_decay(
  backward_decay,
  save_path = file.path(output_dir, "Backward_Prediction_Decay_Fail.png"),
  metric = "Fail_Rate"
)

# =============================================================================
# 6. Extended Density Comparisons
# =============================================================================

message("Generating extended density comparisons...")

density_dir <- file.path(output_dir, "Extended_Density_Comparisons")
if (!dir.exists(density_dir)) dir.create(density_dir, recursive = TRUE)

# Use longest forward horizon as representative
main_forward_name <- tail(names(forward_results), 1)
main_forward <- forward_results[[main_forward_name]]

rep_genes <- select_representative_genes(
  results_df = main_forward$results,
  n_success = 5,
  n_failure = 5,
  n_small_diff = 5
)

density_gene_sets <- list(
  success = rep_genes$success$Names,
  failure = rep_genes$failure$Names,
  small_diff_significant = rep_genes$small_diff_significant$Names
)

target_tp <- max(all_tps)

for (group_name in names(density_gene_sets)) {
  
  group_dir <- file.path(density_dir, group_name)
  if (!dir.exists(group_dir)) dir.create(group_dir, recursive = TRUE)
  
  for (gene in density_gene_sets[[group_name]]) {
    
    if (is.na(gene) || !gene %in% colnames(df)) next
    if (!gene %in% main_forward$df_pred$Gene) next
    
    real_counts <- as.numeric(df[df$TP == target_tp, gene])
    
    pred_counts <- as.numeric(
      main_forward$df_pred[
        main_forward$df_pred$Gene == gene,
        2:ncol(main_forward$df_pred)
      ]
    )
    
    plot_extended_density_comparison(
      gene_name = gene,
      real_counts = real_counts,
      predicted_counts = pred_counts,
      title_prefix = paste("Extended density comparison", group_name),
      save_path = file.path(group_dir, paste0("Density_", gene, ".png"))
    )
  }
}

# =============================================================================
# 7. Stratified Deviation Distributions
# =============================================================================

message("Running stratified deviation analysis...")

strat_dir <- file.path(output_dir, "Stratified_Deviation_Distributions")
if (!dir.exists(strat_dir)) dir.create(strat_dir, recursive = TRUE)

strat_mean <- analyze_stratified_deviations(
  results_df = main_forward$results,
  n_bins = 10,
  stratify_by = "mean_real_count",
  save_path = file.path(strat_dir, "Stratified_Deviations_By_Mean_Expression.csv")
)

p_strat_fc <- plot_stratified_deviation_distribution(
  stratified_object = strat_mean,
  y_var = "abs_log2FC",
  title = "Deviation Distribution Stratified by Mean Expression"
)

ggsave(
  filename = file.path(strat_dir, "Stratified_abs_log2FC_By_Mean_Expression.png"),
  plot = p_strat_fc,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

p_strat_counts_fc <- plot_stratified_deviation_distribution(
  stratified_object = strat_mean,
  y_var = "abs_log2FC_counts",
  title = "Count-scale Deviation Distribution Stratified by Mean Expression"
)

ggsave(
  filename = file.path(strat_dir, "Stratified_abs_log2FC_counts_By_Mean_Expression.png"),
  plot = p_strat_counts_fc,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

p_strat_mean_diff <- plot_stratified_deviation_distribution(
  stratified_object = strat_mean,
  y_var = "abs_mean_diff",
  title = "Mean Difference Distribution Stratified by Mean Expression"
)

ggsave(
  filename = file.path(strat_dir, "Stratified_abs_mean_diff_By_Mean_Expression.png"),
  plot = p_strat_mean_diff,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

# =============================================================================
# 8. Per-gene DE Counts
# =============================================================================

message("Calculating per-gene DE counts...")

de_dir <- file.path(output_dir, "Per_Gene_DE_Counts")
if (!dir.exists(de_dir)) dir.create(de_dir, recursive = TRUE)

de_ttest <- calculate_per_gene_DE_counts(
  results_df = main_forward$results,
  p_col = "t_p_val",
  fc_col = "log2FC",
  save_path = file.path(de_dir, "Per_Gene_DE_Counts_Ttest_Log.csv")
)

de_wilcox <- calculate_per_gene_DE_counts(
  results_df = main_forward$results,
  p_col = "Wilcox_p_val",
  fc_col = "log2FC_counts",
  save_path = file.path(de_dir, "Per_Gene_DE_Counts_Wilcoxon_Raw.csv")
)

de_wilcox_rounded <- calculate_per_gene_DE_counts(
  results_df = main_forward$results,
  p_col = "Wilcox_p_val_rounded",
  fc_col = "log2FC_counts",
  save_path = file.path(de_dir, "Per_Gene_DE_Counts_Wilcoxon_Rounded.csv")
)

write.csv(
  de_ttest$per_gene,
  file.path(de_dir, "Per_Gene_DE_Status_Ttest_Log.csv"),
  row.names = FALSE
)

write.csv(
  de_wilcox$per_gene,
  file.path(de_dir, "Per_Gene_DE_Status_Wilcoxon_Raw.csv"),
  row.names = FALSE
)

write.csv(
  de_wilcox_rounded$per_gene,
  file.path(de_dir, "Per_Gene_DE_Status_Wilcoxon_Rounded.csv"),
  row.names = FALSE
)

print(de_ttest$summary)
print(de_wilcox$summary)
print(de_wilcox_rounded$summary)

# =============================================================================
# 9. Genome-wide Calibration
# =============================================================================

message("Evaluating genome-wide calibration...")

calibration_dir <- file.path(output_dir, "Calibration")
if (!dir.exists(calibration_dir)) dir.create(calibration_dir, recursive = TRUE)

calibration_df <- evaluate_genomewide_calibration(
  df = df,
  df_pred = main_forward$df_pred,
  target_tp = target_tp,
  genes = gene_names,
  lower_prob = 0.025,
  upper_prob = 0.975,
  save_path = file.path(calibration_dir, "Genomewide_Calibration_95CI.csv")
)

plot_genomewide_calibration(
  calibration_df,
  save_path = file.path(calibration_dir, "Genomewide_Calibration_95CI.png")
)

# =============================================================================
# 10. Save Workspace Summary
# =============================================================================

saveRDS(
  list(
    forward_horizon_table = forward_horizon_table,
    backward_horizon_table = backward_horizon_table,
    forward_decay = forward_decay,
    backward_decay = backward_decay,
    representative_genes = rep_genes,
    stratified_mean_expression = strat_mean,
    de_ttest = de_ttest,
    de_wilcox = de_wilcox,
    de_wilcox_rounded = de_wilcox_rounded,
    calibration = calibration_df
  ),
  file.path(output_dir, "Genome_Scale_Extensions_Objects.rds")
)

message("Genome-scale validation extensions complete. Results saved to: ", output_dir)