#' @title Semi-Markov and Brownian Bridge Simulation Toolkit

#' Simulate a Brownian Bridge Between Multiple Timepoints
#' @export
Simul_BB <- function(TIMES, LEVELS, N_simul, N_steps, PLOT = FALSE) {

  if (length(TIMES) != length(LEVELS)) {
    stop("TIMES and LEVELS must have the same length.")
  }

  N_int <- length(TIMES) - 1
  X_out <- NULL
  time_out <- NULL

  for (int in seq_len(N_int)) {
    t0 <- TIMES[int]
    t1 <- TIMES[int + 1]
    a <- LEVELS[int]
    b <- LEVELS[int + 1]

    time <- seq(t0, t1, length.out = N_steps + 1)
    h <- time[2] - time[1]

    X <- matrix(a, nrow = N_simul, ncol = N_steps + 1)

    for (i in seq_len(N_simul)) {
      for (j in 2:(N_steps + 1)) {
        X[i, j] <- X[i, j - 1] +
          (b - X[i, j - 1]) / (t1 - time[j - 1]) * h +
          rnorm(1, sd = sqrt(abs(h)))
      }
    }

    if (int == 1) {
      X_out <- X
      time_out <- time
    } else {
      X_out <- cbind(X_out, X[, 2:ncol(X), drop = FALSE])
      time_out <- c(time_out, time[2:length(time)])
    }
  }

  if (PLOT) {
    plot(time_out, X_out[1, ], type = "l", ylim = c(min(X_out), max(X_out)))
    if (nrow(X_out) > 1) {
      for (i in 2:nrow(X_out)) lines(time_out, X_out[i, ])
    }
  }

  list(time = time_out, X = X_out)
}

#' Simulate States for a Semi-Markov Process
#' @export
Simul_SM <- function(P, s0, N, Times = c()) {
  states <- c(s0)

  if (N < 2) return(states)

  for (i in 2:N) {
    states <- c(states, sample(seq_len(nrow(P)), 1, prob = P[states[i - 1], ]))
  }

  states
}

#' Main Stochastic Simulation Procedure
#'
#' Compatible with both old calls using N and new calls using N_SIM.
#' @export
SIMUL_PROC_rand_init <- function(P,
                                 Par,
                                 Times_sim,
                                 Times_forec,
                                 N = NULL,
                                 N_SIM = 100,
                                 N_Steps = 200,
                                 state_probs = NULL) {

  if (!is.null(N)) N_SIM <- N

  P <- as.matrix(P)
  n_states <- nrow(P)

  Par <- as.matrix(Par)
  if (nrow(Par) == 2 && ncol(Par) == n_states) {
    Par <- t(Par)
  }

  if (nrow(Par) != n_states || ncol(Par) < 2) {
    stop("Par must be either [state_dim x 2] or [2 x state_dim].")
  }

  if (is.null(state_probs)) {
    state_probs <- rep(1 / n_states, n_states)
  }

  state_probs <- as.numeric(state_probs[seq_len(n_states)])
  state_probs <- state_probs / sum(state_probs)

  Times_gen <- c(Times_sim, Times_forec)
  N_g <- length(Times_gen)

  X_gen <- NULL
  time <- NULL

  for (i in seq_len(N_SIM)) {
    S0 <- sample(seq_len(n_states), 1, prob = state_probs)
    seq_states_SM <- Simul_SM(P = P, s0 = S0, N = N_g)

    Av_norm <- Par[seq_states_SM, 1]
    Sd_norm <- Par[seq_states_SM, 2]
    seq_norm_dist <- rnorm(N_g, Av_norm, sqrt(Sd_norm))

    bb <- Simul_BB(
      TIMES = Times_gen,
      LEVELS = seq_norm_dist,
      N_simul = 1,
      N_steps = N_Steps,
      PLOT = FALSE
    )

    if (i == 1) {
      X_gen <- bb$X
      time <- bb$time
    } else {
      X_gen <- rbind(X_gen, bb$X)
    }
  }

  list(X = X_gen, time = time)
}

#' Normalize matrix rows to sum to one
#' @export
normalize_rows <- function(A) {
  A <- as.matrix(A)
  row_sums <- rowSums(A, na.rm = TRUE)
  row_sums[row_sums == 0 | !is.finite(row_sums)] <- 1
  A / row_sums
}
