#' @title Parameter Estimation for Gene Expression Mixtures
#' @description Functions to estimate Gaussian Mixture Model (GMM) parameters 
#' across timepoints and construct transition matrices for Semi-Markov processes.
#' @author Taras Lukashiv, Igor Malyk, Mathias Galati, Ahmed Hemedan and Venkata Satagopam
#' @date 2026-05-06

library(mclust)
library(transport)
library(modeest)

#' Estimate GMM Parameters for a Single Gene
#' 
#' @description Fits a Gaussian Mixture Model to expression data at each timepoint.
#' Handles zero-inflated data by providing fallback parameters.
#' 
#' @param gene_values Numeric vector of expression values (log-transformed).
#' @param timepoints Vector of timepoints corresponding to the expression values.
#' @param n_opt Integer. The fixed number of mixture components to return.
#' 
#' @return A data frame with columns: TP (Timepoint), pi (Weights), mu (Means), sigma (Variance).
#' @export
estimate_mixture_params <- function(gene_values, timepoints, n_opt = 2) {
  
  unique_tp <- sort(unique(as.numeric(as.character(timepoints))))
  param_df <- data.frame()
  
  for (tp in unique_tp) {
    sam <- as.numeric(gene_values[as.numeric(as.character(timepoints)) == tp])
    sam <- sam[is.finite(sam)]
    
    # 1. Handle cases with too many zeros (Zero-inflation check)
    if (quantile(sam, 0.95, na.rm = TRUE) == 0) {
      curr_param <- data.frame(
        TP = rep(tp, n_opt),
        pi = c(0.99, rep(0.01 / (n_opt - 1), n_opt - 1)),
        mu = c(0, rep(max(sam), n_opt - 1)),
        sigma = c(0.01, rep(0.05, n_opt - 1))
      )
    } else {
      # 2. Fit Gaussian Mixture Model using EM algorithm
      # We force G = n_opt to maintain consistent state dimensionality
      mix_c <- mclust::Mclust(sam, G = 1:n_opt, verbose = FALSE)
      
      if (is.null(mix_c)) {
        # Fallback if Mclust fails to converge
        new_var <- data.frame(
          pi = rep(1/n_opt, n_opt),
          mu = rep(mean(sam), n_opt),
          sigma = rep(var(sam), n_opt)
        )
      } else {
        new_var <- data.frame(
          pi = mix_c$parameters$pro,
          mu = mix_c$parameters$mean,
          sigma = mix_c$parameters$variance$sigmasq
        )
        
        # Ensure we always have exactly n_opt rows
        d1 <- nrow(new_var)
        if (d1 < n_opt) {
          new_var[d1:n_opt, ] <- new_var[d1, ]
          new_var$pi[d1:n_opt] <- new_var$pi[d1] / (n_opt + 1 - d1)
        }
      }
      
      # 3. Constrain Sigma to avoid extreme variances
      curr_param <- cbind(TP = tp, new_var)
      colnames(curr_param)[1] <- "TP"
      
      if (n_opt >= 2) {
        curr_param$sigma[2] <- min(curr_param$sigma[2], 10)
      }
    }
    param_df <- rbind(param_df, curr_param)
  }
  
  return(param_df)
}

#' Construct Transition Matrix via Wasserstein Distance
#' 
#' @description Builds a transition matrix between mixture components of 
#' subsequent timepoints using the inverse of the 1D Wasserstein distance.
#' 
#' @param param_df Data frame containing estimated GMM parameters.
#' @param weights_tr Numeric vector of weights for each time transition.
#' @param state_dim Integer. Number of mixture components (states).
#' 
#' @return A square matrix of transition probabilities.
#' @export
build_transition_matrix <- function(param_df, weights_tr, state_dim) {
  
  unique_tps <- unique(param_df$TP)
  n_transitions <- length(unique_tps) - 1
  nm <- paste0("stage", 1:state_dim)
  
  all_transitions <- list()
  
  for (i in 1:n_transitions) {
    s1 <- param_df[param_df$TP == unique_tps[i], c("mu", "sigma")]
    s2 <- param_df[param_df$TP == unique_tps[i+1], c("mu", "sigma")]
    
    temp_mat <- matrix(0, nrow = state_dim, ncol = state_dim)
    
    for (j in 1:state_dim) {
      for (k in 1:state_dim) {
        # Generate samples to estimate Wasserstein distance
        samp1 <- rnorm(100, mean = s1[j, 1], sd = sqrt(s1[j, 2]))
        samp2 <- rnorm(100, mean = s2[k, 1], sd = sqrt(s2[k, 2]))
        
        # Probability is inversely proportional to distance
        dist <- transport::wasserstein1d(samp1, samp2)
        temp_mat[j, k] <- 1 / dist 
      }
      # Normalize row
      temp_mat[j, ] <- temp_mat[j, ] / sum(temp_mat[j, ])
    }
    all_transitions[[i]] <- temp_mat
  }
  
  # Final matrix as a weighted average of transitions across time
  final_trans_prob <- matrix(0, nrow = state_dim, ncol = state_dim)
  
  for (i in 1:n_transitions) {
    final_trans_prob <- final_trans_prob + (all_transitions[[i]] * weights_tr[i])
  }
  
  # Re-normalize final matrix
  final_trans_prob <- final_trans_prob / rowSums(final_trans_prob)
  colnames(final_trans_prob) <- nm
  rownames(final_trans_prob) <- nm
  
  return(final_trans_prob)
}

#' Calculate Weighted Checked Parameters
#' 
#' @description Aggregates parameters across all timepoints using 
#' exponential weights to create a "global" state profile for simulation.
#' 
#' @param param_df Data frame from estimate_mixture_params.
#' @param weights Numeric vector of weights for timepoints.
#' @param n_opt Integer. Number of components.
#' 
#' @return A matrix [2 x n_opt] with Mu (row 1) and Sigma (row 2).
#' @export
get_weighted_check_params <- function(param_df, weights, n_opt) {
  par_check <- matrix(NA, ncol = n_opt, nrow = 2)
  
  for (i in 1:n_opt) {
    # Extract every i-th component (for each timepoint)
    indices <- seq(i, nrow(param_df), by = n_opt)
    par_check[1, i] <- weighted.mean(param_df$mu[indices], weights)
    par_check[2, i] <- weighted.mean(param_df$sigma[indices], weights)
  }
  
  return(par_check)
}