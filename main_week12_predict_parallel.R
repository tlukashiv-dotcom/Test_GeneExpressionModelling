#' @title Genome-scale Week-12 Prediction Pipeline
#' @author Taras Lukashiv, Igor Malyk, Mathias Galati, Ahmed Hemedan and Venkata Satagopam
#' @date 2026-05-06

# =============================================================================
# 0. Load modules
# =============================================================================

source("requirements.R")

library(ggplot2)
library(dplyr)
library(foreach)
library(doParallel)
library(parallel)
library(mclust)

source("simulation_functions.R")
source("parameter_estimation.R")
source("validation_functions.R")
source("process_cef_data.R")
source("prediction_functions.R")
source("trajectory_analysis.R")

# =============================================================================
# 1. Initialization
# =============================================================================

set.seed(125)

output_dir <- "Results_Week12_Forecast"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

gene_dir <- file.path(output_dir, "Per_Gene_Results")
if (!dir.exists(gene_dir)) dir.create(gene_dir, recursive = TRUE)

section37_dir <- file.path(output_dir, "Section_3_7_Coherence")
if (!dir.exists(section37_dir)) dir.create(section37_dir, recursive = TRUE)

N_SIM <- 100
n_opt <- 2
N_Steps <- 200
memory_decay <- 1
N_CORES <- 10

# =============================================================================
# 2. Data loading
# =============================================================================

df <- process_cef_data("GSE76381_EmbryoMoleculeCounts.cef.txt")
df$TP <- as.numeric(substr(df$Timepoint, 6, nchar(df$Timepoint)))

all_tps <- sort(unique(na.omit(df$TP)))

start_tp <- min(all_tps)
end_tp <- max(all_tps)
prediction_tp <- end_tp + 1

message("Available timepoints: ", paste(all_tps, collapse = ", "))
message("Forecasting from Week ", start_tp, "–", end_tp, " to Week ", prediction_tp)

genes <- read.csv("all_with_sig_identif_new.csv", header = TRUE)
genes <- genes[genes$significant == 1, ]

gene_names <- genes[, 1]
gene_names <- gene_names[gene_names %in% colnames(df)]

N_f <- length(gene_names)

message("Number of genes used for Week-12 forecasting: ", N_f)

# =============================================================================
# 3. Helper functions
# =============================================================================

clean_gene_symbol <- function(x) {
  x <- sub("_loc.*$", "", x)
  x <- sub("_p[0-9]+$", "", x)
  x
}

summarise_prediction <- function(pred_counts) {
  
  q <- quantile(
    pred_counts,
    probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975),
    na.rm = TRUE
  )
  
  data.frame(
    mean_pred = mean(pred_counts, na.rm = TRUE),
    sd_pred = sd(pred_counts, na.rm = TRUE),
    q_0.025 = q[1],
    q_0.05 = q[2],
    q_0.25 = q[3],
    q_0.50 = q[4],
    q_0.75 = q[5],
    q_0.95 = q[6],
    q_0.975 = q[7],
    stringsAsFactors = FALSE
  )
}

save_week12_gene_outputs <- function(gene,
                                     pred_res,
                                     pred_counts,
                                     pred_log,
                                     df,
                                     gene_dir,
                                     N_SIM,
                                     n_opt = 2) {
  
  safe_gene <- gsub("[^A-Za-z0-9_\\-]", "_", gene)
  gene_subdir <- file.path(gene_dir, safe_gene)
  if (!dir.exists(gene_subdir)) dir.create(gene_subdir, recursive = TRUE)
  
  RES <- pred_res$simulation
  RES$X_untransf <- 2^RES$X - 1
  
  # ---------------------------------------------------------------------------
  # 3.1 Violin plot for observed weeks
  # ---------------------------------------------------------------------------
  
  df_violin <- data.frame(
    TP = as.factor(df$TP),
    gene_val = log2(as.numeric(df[[gene]]) + 1)
  )
  
  p_violin <- ggplot(df_violin, aes(x = TP, y = gene_val)) +
    geom_violin(trim = FALSE, fill = "#A4A4A4", color = "darkred", alpha = 0.5) +
    geom_point(
      position = position_jitter(width = 0.15, height = 0, seed = 1),
      alpha = 0.4,
      size = 1
    ) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5) +
    labs(
      title = paste("Violin Plot for Gene:", gene),
      y = paste0("log2(", gene, " + 1)"),
      x = "Time, week"
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold")
    )
  
  ggsave(
    filename = file.path(gene_subdir, paste0("Violin_plot_", safe_gene, ".png")),
    plot = p_violin,
    width = 7,
    height = 5,
    dpi = 300,
    bg = "white"
  )
  
  # ---------------------------------------------------------------------------
  # 3.2 Mixture component plots for observed weeks
  # ---------------------------------------------------------------------------
  
  observed_tps <- sort(unique(na.omit(df$TP)))
  
  for (tp in observed_tps) {
    
    sam <- log2(as.numeric(df[df$TP == tp, gene]) + 1)
    sam <- sam[is.finite(sam)]
    
    if (length(sam) < 5 || length(unique(sam)) < 2) next
    
    mix_c <- tryCatch(
      mclust::Mclust(sam, G = 1:n_opt),
      error = function(e) NULL
    )
    
    if (is.null(mix_c)) next
    
    # -------------------------------------------------------------------------
    # Classification uncertainty plot
    # -------------------------------------------------------------------------
    
    file_name_class <- file.path(
      gene_subdir,
      paste0("Classification_Refined_", tp, "_", safe_gene, ".png")
    )
    
    png(file_name_class, width = 1000, height = 600, res = 120)
    
    vals <- as.numeric(mix_c$data)
    uncertainty <- 1 - apply(mix_c$z, 1, max)
    
    if (mix_c$G <= 2) {
      class_cols <- c("#E41A1C77", "#377EB877")
      mean_cols <- c("#E41A1C", "#377EB8")
    } else {
      base_cols <- grDevices::rainbow(mix_c$G)
      class_cols <- paste0(base_cols, "77")
      mean_cols <- base_cols
    }
    
    point_cols <- class_cols[mix_c$classification]
    
    par(mar = c(5, 4, 4, 2) + 0.1)
    
    plot(
      vals,
      uncertainty,
      type = "n",
      xlab = paste0("Expression Level of ", gene, ", week ", tp),
      ylab = "Classification Uncertainty",
      main = "Refined Classification & Error Map",
      frame.plot = FALSE
    )
    
    grid(col = "grey90", lty = "solid")
    
    h <- hist(vals, plot = FALSE)
    
    if (max(h$density, na.rm = TRUE) > 0 && max(uncertainty, na.rm = TRUE) > 0) {
      rect(
        h$breaks[-length(h$breaks)],
        0,
        h$breaks[-1],
        h$density * (max(uncertainty, na.rm = TRUE) / max(h$density, na.rm = TRUE)),
        col = "#F0F0F0",
        border = "white"
      )
    }
    
    deterministic_offset <- (seq_along(vals) %% 10 - 5) / 200
    vals_shifted <- vals + deterministic_offset
    
    points(
      vals_shifted,
      uncertainty,
      pch = 19,
      col = point_cols,
      cex = 0.6 + uncertainty * 2
    )
    
    abline(
      v = mix_c$parameters$mean,
      col = mean_cols[seq_len(mix_c$G)],
      lty = 3,
      lwd = 1.5
    )
    
    legend(
      "topright",
      legend = c(paste0("Cluster ", seq_len(mix_c$G)), "Total Density Area"),
      col = c(mean_cols[seq_len(mix_c$G)], "#F0F0F0"),
      pch = c(rep(19, mix_c$G), 15),
      pt.cex = 1.5,
      bty = "n"
    )
    
    dev.off()
    
    file_name_density <- file.path(
      gene_subdir,
      paste0("Mixure_plot_stage_", tp, "_", safe_gene, ".png")
    )
    
    png(file_name_density, width = 1000, height = 700, res = 120)
    
    vals <- as.numeric(mix_c$data)
    x_range <- seq(
      min(vals) - diff(range(vals)) * 0.1,
      max(vals) + diff(range(vals)) * 0.1,
      length.out = 400
    )
    
    hist(
      vals,
      freq = FALSE,
      col = "grey85",
      border = "white",
      xlab = paste0("Sample ", gene, ", week ", tp),
      main = "Density with Components"
    )
    
    grid(col = "lightgray", lty = "dotted")
    
    cols <- if (mix_c$G <= 2) c("red", "blue") else rainbow(mix_c$G)
    
    for (k in seq_len(mix_c$G)) {
      
      mu <- mix_c$parameters$mean[k]
      vars <- mix_c$parameters$variance$sigmasq
      sigma <- sqrt(if (length(vars) > 1) vars[k] else vars)
      pro <- mix_c$parameters$pro[k]
      
      y_comp <- pro * dnorm(x_range, mean = mu, sd = sigma)
      
      lines(x_range, y_comp, col = cols[k], lwd = 2.5, lty = 2)
    }
    
    y_total <- tryCatch(
      mclust::dens(
        modelName = mix_c$modelName,
        data = x_range,
        parameters = mix_c$parameters
      ),
      error = function(e) NULL
    )
    
    if (!is.null(y_total)) {
      lines(x_range, y_total, col = "black", lwd = 2)
    }
    
    comp_labels <- paste0(
      "Comp ",
      seq_len(mix_c$G),
      " (",
      round(mix_c$parameters$pro * 100, 1),
      "%)"
    )
    
    legend(
      "topright",
      legend = c("Total Density", comp_labels),
      col = c("black", cols[seq_len(mix_c$G)]),
      lty = c(1, rep(2, mix_c$G)),
      lwd = 2,
      bty = "n",
      cex = 0.9
    )
    
    dev.off()
  }
  
  # ---------------------------------------------------------------------------
  # 3.3 Save transition matrix and mixture parameters
  # ---------------------------------------------------------------------------
  
  write.csv(
    pred_res$transition_matrix,
    file.path(gene_subdir, paste0("Transition_Matrix_", safe_gene, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    pred_res$parameters,
    file.path(gene_subdir, paste0("Mixture_Parameters_", safe_gene, ".csv")),
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 3.4 Save prediction vectors
  # ---------------------------------------------------------------------------
  
  write.csv(
    data.frame(Predicted_Week12_Counts = pred_counts),
    file.path(gene_subdir, paste0("Predicted_Week12_Counts_", safe_gene, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(Predicted_Week12_Log = pred_log),
    file.path(gene_subdir, paste0("Predicted_Week12_Log_", safe_gene, ".csv")),
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 3.5 Save full trajectories
  # ---------------------------------------------------------------------------
  
  trajectories_df <- data.frame(t(RES$X_untransf))
  trajectories_df <- cbind(Time = RES$time, trajectories_df)
  colnames(trajectories_df)[-1] <- paste0(
    "Trajectory_",
    seq_len(ncol(trajectories_df) - 1)
  )
  
  write.csv(
    trajectories_df,
    file.path(gene_subdir, paste0("Trajectories_Week12_", safe_gene, ".csv")),
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 3.6 Save CI metrics at Week 12
  # ---------------------------------------------------------------------------
  
  final_values <- RES$X_untransf[, ncol(RES$X_untransf)]
  
  Quant <- quantile(
    final_values,
    probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975),
    na.rm = TRUE
  )
  
  ci_df <- data.frame(
    Mean = mean(final_values, na.rm = TRUE),
    Sd = sd(final_values, na.rm = TRUE),
    Q_0.025 = Quant[1],
    Q_0.05 = Quant[2],
    Q_0.25 = Quant[3],
    Q_0.5 = Quant[4],
    Q_0.75 = Quant[5],
    Q_0.95 = Quant[6],
    Q_0.975 = Quant[7]
  )
  
  write.csv(
    ci_df,
    file.path(gene_subdir, paste0("Results_CI_metrics_Week12_", safe_gene, ".csv")),
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 3.7 Plot stochastic trajectories + confidence band
  # ---------------------------------------------------------------------------
  
  Average_Trajectory <- apply(RES$X_untransf, 2, mean, na.rm = TRUE)
  Lower_Bound <- apply(RES$X_untransf, 2, quantile, probs = 0.025, na.rm = TRUE)
  Upper_Bound <- apply(RES$X_untransf, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  png(
    file.path(gene_subdir, paste0("Simulations_plot_comb_Week12_", safe_gene, ".png")),
    width = 1000,
    height = 700,
    res = 120
  )
  
  y_max_val <- min(max(RES$X_untransf, na.rm = TRUE), 250)
  y_range <- c(min(RES$X_untransf, na.rm = TRUE), y_max_val)
  
  plot(
    RES$time,
    Average_Trajectory,
    type = "n",
    ylim = y_range,
    xlab = "Time, week",
    ylab = paste0("Expression of ", gene),
    main = paste("Simulated Week-12 Trajectories for", gene),
    frame.plot = FALSE,
    col.main = "black",
    font.main = 2
  )
  
  grid(col = "grey85", lty = "solid")
  
  for (i in seq_len(nrow(RES$X_untransf))) {
    lines(
      RES$time,
      RES$X_untransf[i, ],
      col = rgb(93, 173, 226, maxColorValue = 255, alpha = 45),
      lwd = 0.8
    )
  }
  
  polygon(
    c(RES$time, rev(RES$time)),
    c(Lower_Bound, rev(Upper_Bound)),
    col = rgb(243, 156, 18, maxColorValue = 255, alpha = 60),
    border = NA
  )
  
  lines(RES$time, Average_Trajectory, col = "black", lwd = 2.5)
  lines(RES$time, Lower_Bound, col = "#D35400", lty = 3, lwd = 1)
  lines(RES$time, Upper_Bound, col = "#D35400", lty = 3, lwd = 1)
  
  legend(
    "topright",
    legend = c("Mean Trajectory", "95% Confidence Interval", "Stochastic Simulations"),
    col = c("black", "#F39C12", "#5DADE2"),
    lwd = c(2.5, 8, 1),
    bty = "n",
    cex = 0.9,
    text.font = 2
  )
  
  dev.off()
  
  # ---------------------------------------------------------------------------
  # 3.8 Plot rounded average integer trajectory
  # ---------------------------------------------------------------------------
  
  Average_Trajectory_int <- round(Average_Trajectory, digits = 0)
  
  png(
    file.path(gene_subdir, paste0("Simulations_plot_average_int_Week12_", safe_gene, ".png")),
    width = 1000,
    height = 700,
    res = 120
  )
  
  y_min <- min(RES$X_untransf, na.rm = TRUE)
  y_max <- min(max(RES$X_untransf, na.rm = TRUE), 100)
  
  plot(
    RES$time,
    Average_Trajectory_int,
    type = "n",
    ylim = c(y_min, y_max),
    xlab = "Time, week",
    ylab = gene,
    main = paste("Average Integer Week-12 Trajectory:", gene),
    frame.plot = FALSE
  )
  
  grid(col = "grey90", lty = "solid")
  
  polygon(
    c(min(RES$time), RES$time, max(RES$time)),
    c(y_min, Average_Trajectory_int, y_min),
    col = rgb(0, 115, 194, maxColorValue = 255, alpha = 40),
    border = NA
  )
  
  lines(
    RES$time,
    Average_Trajectory_int,
    col = "#0029C9",
    lwd = 3
  )
  
  dev.off()
  
  return(gene_subdir)
}

# =============================================================================
# 4. Parallel setup
# =============================================================================

n_available <- parallel::detectCores()
n_workers <- min(N_CORES, n_available)

message("Starting parallel backend with ", n_workers, " workers.")

cl <- parallel::makeCluster(n_workers)
doParallel::registerDoParallel(cl)

clusterEvalQ(cl, {
  library(ggplot2)
  library(dplyr)
  library(mclust)
  library(mixdist)
  library(mixtools)
  library(transport)
  library(DescTools)
})

clusterExport(
  cl,
  varlist = c(
    "df",
    "gene_names",
    "start_tp",
    "end_tp",
    "prediction_tp",
    "n_opt",
    "N_SIM",
    "N_Steps",
    "memory_decay",
    "gene_dir",
    "clean_gene_symbol",
    "summarise_prediction",
    "save_week12_gene_outputs"
  ),
  envir = environment()
)

clusterEvalQ(cl, {
  source("simulation_functions.R")
  source("parameter_estimation.R")
  source("validation_functions.R")
  source("prediction_functions.R")
  source("trajectory_analysis.R")
})

# =============================================================================
# 5. Parallel Week-12 prediction
# =============================================================================

forecast_list <- foreach(
  g = seq_along(gene_names),
  .packages = c(
    "ggplot2",
    "dplyr",
    "mclust",
    "mixtools",
    "transport",
    "DescTools"
  ),
  .errorhandling = "pass"
) %dopar% {
  
  gene <- gene_names[g]
  
  tryCatch({
    
    pred_res <- predict_gene_expression(
      df = df,
      gene = gene,
      start_tp = start_tp,
      end_tp = end_tp,
      prediction_tp = prediction_tp,
      n_opt = n_opt,
      N_SIM = N_SIM,
      N_Steps = N_Steps,
      memory_decay = memory_decay,
      direction = "forward"
    )
    
    pred_log <- pred_res$pred_log
    pred_counts <- pred_res$pred_counts
    
    gene_output_dir <- save_week12_gene_outputs(
      gene = gene,
      pred_res = pred_res,
      pred_counts = pred_counts,
      pred_log = pred_log,
      df = df,
      gene_dir = gene_dir,
      N_SIM = N_SIM,
      n_opt = n_opt
    )
    
    week11_log <- log2(as.numeric(df[df$TP == end_tp, gene]) + 1)
    week11_counts <- 2^week11_log - 1
    
    ci_summary <- summarise_prediction(pred_counts)
    
    mean_week11_log <- mean(week11_log, na.rm = TRUE)
    mean_week12_log <- mean(pred_log, na.rm = TRUE)
    
    median_week11_counts <- median(week11_counts, na.rm = TRUE)
    median_week12_counts <- median(pred_counts, na.rm = TRUE)
    
    t_p <- tryCatch(
      t.test(week11_log, pred_log)$p.value,
      error = function(e) NA_real_
    )
    
    wilcox_p <- tryCatch(
      wilcox.test(week11_counts, pred_counts)$p.value,
      error = function(e) NA_real_
    )
    
    result_row <- data.frame(
      Gene = gene,
      GeneSymbol = clean_gene_symbol(gene),
      mean_week11_log = mean_week11_log,
      mean_week12_log = mean_week12_log,
      log2FC_log_mean = mean_week12_log - mean_week11_log,
      median_week11_count = median_week11_counts,
      median_week12_count = median_week12_counts,
      log2FC_count_median = log2((median_week12_counts + 1) / (median_week11_counts + 1)),
      t_p_val = t_p,
      Wilcox_p_val = wilcox_p,
      mean_pred_count = ci_summary$mean_pred,
      sd_pred_count = ci_summary$sd_pred,
      q_0.025 = ci_summary$q_0.025,
      q_0.05 = ci_summary$q_0.05,
      q_0.25 = ci_summary$q_0.25,
      q_0.50 = ci_summary$q_0.50,
      q_0.75 = ci_summary$q_0.75,
      q_0.95 = ci_summary$q_0.95,
      q_0.975 = ci_summary$q_0.975,
      gene_output_dir = gene_output_dir,
      stringsAsFactors = FALSE
    )
    
    list(
      gene = gene,
      success = TRUE,
      result_row = result_row,
      pred_log = pred_log,
      pred_counts = pred_counts,
      gene_output_dir = gene_output_dir
    )
    
  }, error = function(e) {
    
    list(
      gene = gene,
      success = FALSE,
      error = e$message
    )
  })
}

parallel::stopCluster(cl)

message("Parallel Week-12 forecasting complete.")

# =============================================================================
# 6. Collect results
# =============================================================================

successful <- forecast_list[
  sapply(forecast_list, function(x) isTRUE(x$success))
]

failed <- forecast_list[
  sapply(forecast_list, function(x) !isTRUE(x$success))
]

message("Successfully predicted genes: ", length(successful))
message("Failed genes: ", length(failed))

if (length(successful) == 0) {
  stop("No successful Week-12 predictions.")
}

forecast_summary <- bind_rows(
  lapply(successful, function(x) x$result_row)
)

forecast_summary$adj_p_t <- p.adjust(forecast_summary$t_p_val, method = "BH")
forecast_summary$adj_p_wilcox <- p.adjust(forecast_summary$Wilcox_p_val, method = "BH")

forecast_summary <- forecast_summary %>%
  mutate(
    DE_t_test = adj_p_t < 0.05 & abs(log2FC_log_mean) >= 1,
    DE_wilcox = adj_p_wilcox < 0.05 & abs(log2FC_count_median) >= 1,
    Direction_log = case_when(
      log2FC_log_mean > 0 ~ "Predicted_up",
      log2FC_log_mean < 0 ~ "Predicted_down",
      TRUE ~ "No_change"
    ),
    Direction_count = case_when(
      log2FC_count_median > 0 ~ "Predicted_up",
      log2FC_count_median < 0 ~ "Predicted_down",
      TRUE ~ "No_change"
    )
  )

write.csv(
  forecast_summary,
  file.path(output_dir, "Week12_Forecast_Summary.csv"),
  row.names = FALSE
)

if (length(failed) > 0) {
  failed_df <- data.frame(
    Gene = sapply(failed, function(x) x$gene),
    Error = sapply(failed, function(x) x$error),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    failed_df,
    file.path(output_dir, "Week12_Forecast_Failed_Genes.csv"),
    row.names = FALSE
  )
}

# =============================================================================
# 7. Save prediction matrices
# =============================================================================

df_pred <- data.frame(
  Gene = sapply(successful, function(x) x$gene),
  do.call(rbind, lapply(successful, function(x) round(x$pred_counts, digits = 0)))
)

df_pred_log <- data.frame(
  Gene = sapply(successful, function(x) x$gene),
  do.call(rbind, lapply(successful, function(x) x$pred_log))
)

colnames(df_pred)[-1] <- paste0("X_", seq_len(N_SIM))
colnames(df_pred_log)[-1] <- paste0("X_", seq_len(N_SIM))

write.csv(
  df_pred,
  file.path(output_dir, "Prediction12_int.csv"),
  row.names = FALSE
)

write.csv(
  df_pred_log,
  file.path(output_dir, "Prediction12_log.csv"),
  row.names = FALSE
)

# =============================================================================
# 8. Save pathfindR input tables for Section 3.7
# =============================================================================

pathfindr_input_ttest <- forecast_summary %>%
  transmute(
    Gene_symbol = GeneSymbol,
    logFC = log2FC_log_mean,
    adj.P.Val = adj_p_t,
    p_val = t_p_val
  ) %>%
  filter(
    !is.na(Gene_symbol),
    !is.na(logFC),
    !is.na(adj.P.Val),
    Gene_symbol != ""
  )

pathfindr_input_wilcox <- forecast_summary %>%
  transmute(
    Gene_symbol = GeneSymbol,
    logFC = log2FC_count_median,
    adj.P.Val = adj_p_wilcox,
    p_val = Wilcox_p_val
  ) %>%
  filter(
    !is.na(Gene_symbol),
    !is.na(logFC),
    !is.na(adj.P.Val),
    Gene_symbol != ""
  )

write.csv(
  pathfindr_input_ttest,
  file.path(section37_dir, "pathfindR_input_Week12_vs_Week11_Ttest.csv"),
  row.names = FALSE
)

write.csv(
  pathfindr_input_wilcox,
  file.path(section37_dir, "pathfindR_input_Week12_vs_Week11_Wilcoxon.csv"),
  row.names = FALSE
)

# Also save in root output directory
write.csv(
  pathfindr_input_ttest,
  file.path(output_dir, "pathfindR_input_Week12_vs_Week11_Ttest.csv"),
  row.names = FALSE
)

write.csv(
  pathfindr_input_wilcox,
  file.path(output_dir, "pathfindR_input_Week12_vs_Week11_Wilcoxon.csv"),
  row.names = FALSE
)

# =============================================================================
# 9. Save Section 3.7 coherence support tables
# =============================================================================

section37_ranked_ttest <- forecast_summary %>%
  arrange(adj_p_t, desc(abs(log2FC_log_mean))) %>%
  select(
    Gene,
    GeneSymbol,
    log2FC_log_mean,
    adj_p_t,
    mean_week11_log,
    mean_week12_log,
    mean_pred_count,
    q_0.025,
    q_0.975,
    Direction_log
  )

section37_ranked_wilcox <- forecast_summary %>%
  arrange(adj_p_wilcox, desc(abs(log2FC_count_median))) %>%
  select(
    Gene,
    GeneSymbol,
    log2FC_count_median,
    adj_p_wilcox,
    median_week11_count,
    median_week12_count,
    mean_pred_count,
    q_0.025,
    q_0.975,
    Direction_count
  )

write.csv(
  section37_ranked_ttest,
  file.path(section37_dir, "Week12_vs_Week11_Ranked_Genes_Ttest.csv"),
  row.names = FALSE
)

write.csv(
  section37_ranked_wilcox,
  file.path(section37_dir, "Week12_vs_Week11_Ranked_Genes_Wilcoxon.csv"),
  row.names = FALSE
)

section37_summary <- data.frame(
  Metric = c(
    "Predicted Week 12 genes",
    "T-test DE genes Week12 vs Week11",
    "Wilcoxon DE genes Week12 vs Week11",
    "Mean predicted Week12 count",
    "Median predicted Week12 count",
    "Mean log2FC log-space",
    "Median log2FC count-space"
  ),
  Value = c(
    nrow(forecast_summary),
    sum(forecast_summary$DE_t_test, na.rm = TRUE),
    sum(forecast_summary$DE_wilcox, na.rm = TRUE),
    round(mean(forecast_summary$mean_pred_count, na.rm = TRUE), 4),
    round(median(forecast_summary$median_week12_count, na.rm = TRUE), 4),
    round(mean(forecast_summary$log2FC_log_mean, na.rm = TRUE), 4),
    round(median(forecast_summary$log2FC_count_median, na.rm = TRUE), 4)
  )
)

write.csv(
  section37_summary,
  file.path(section37_dir, "Section_3_7_Week12_Coherence_Summary.csv"),
  row.names = FALSE
)

# =============================================================================
# 10. Save complete RDS object
# =============================================================================

saveRDS(
  list(
    forecast_summary = forecast_summary,
    prediction_counts = df_pred,
    prediction_log = df_pred_log,
    successful_genes = sapply(successful, function(x) x$gene),
    failed_genes = failed,
    section37_summary = section37_summary,
    settings = list(
      N_SIM = N_SIM,
      n_opt = n_opt,
      N_Steps = N_Steps,
      memory_decay = memory_decay,
      start_tp = start_tp,
      end_tp = end_tp,
      prediction_tp = prediction_tp,
      n_workers = n_workers
    )
  ),
  file.path(output_dir, "Week12_Forecast_Objects.rds")
)

# =============================================================================
# 11. Console summary
# =============================================================================

message("Week-12 forecast completed.")
message("Output directory: ", output_dir)
message("Per-gene outputs saved to: ", gene_dir)
message("Section 3.7 support files saved to: ", section37_dir)

message("T-test DE genes Week12 vs Week11: ", sum(forecast_summary$DE_t_test, na.rm = TRUE))
message("Wilcoxon DE genes Week12 vs Week11: ", sum(forecast_summary$DE_wilcox, na.rm = TRUE))

message("Mean predicted Week12 expression: ", round(mean(forecast_summary$mean_pred_count, na.rm = TRUE), 3))
message("Median predicted Week12 expression: ", round(median(forecast_summary$median_week12_count, na.rm = TRUE), 3))