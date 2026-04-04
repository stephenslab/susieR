# =============================================================================
# FIXED MIXTURE PRIOR
#
# Shared implementations for estimate_prior_method = "fixed_mixture".
# Evaluates Bayes factors on a pre-specified variance grid with given mixture
# weights, computes mixture posterior moments, and stores per-grid BF matrix.
#
# These functions are data-type-agnostic: they operate on betahat and shat2
# produced by the type-specific compute_ser_statistics().
# =============================================================================

#' Compute mixture log-Bayes factors and posterior inclusion probabilities
#'
#' For each grid point k and variant j, computes the Wakefield approximate
#' Bayes factor (ABF), then forms the mixture BF as a weighted sum over grid
#' points. Stores the full p x K log-BF matrix in model$lbf_grid[[l]] for
#' downstream use (e.g., mixsqp M-step in susieAnn).
#'
#' @param params Params object with prior_variance_grid (K-vector) and
#'   mixture_weights (K-vector summing to 1)
#' @param model Current model object with pi (prior weights)
#' @param ser_stats List with betahat (p-vector) and shat2 (p-vector)
#' @param l Effect index
#'
#' @return Updated model with alpha[l,], lbf[l], lbf_variable[l,], lbf_grid[[l]]
#'
#' @keywords internal
loglik_mixture_common <- function(params, model, ser_stats, l) {
  grid <- params$prior_variance_grid   # K-vector
  w    <- params$mixture_weights       # K-vector
  K    <- length(grid)

  betahat <- ser_stats$betahat         # p-vector
  shat2   <- pmax(ser_stats$shat2, .Machine$double.eps)  # p-vector
  p       <- length(betahat)

  # Compute p x K matrix of log-BFs (Wakefield ABF at each grid point)
  # lbf[j,k] = -0.5 * log(1 + V_k/shat2_j) + 0.5 * betahat_j^2 * V_k / (shat2_j * (V_k + shat2_j))
  lbf_grid <- matrix(0, nrow = p, ncol = K)
  for (k in seq_len(K)) {
    V_k <- grid[k]
    lbf_grid[, k] <- -0.5 * log(1 + V_k / shat2) +
      0.5 * betahat^2 * V_k / (shat2 * (V_k + shat2))
  }

  # Handle infinite shat2 (no information)
  inf_idx <- is.infinite(shat2)
  if (any(inf_idx)) {
    lbf_grid[inf_idx, ] <- 0
  }

  # Mixture log-BF per variant: log(sum_k w_k * exp(lbf_jk))
  # Use log-sum-exp for numerical stability
  log_w <- log(w + .Machine$double.eps)
  lbf_mix <- apply(lbf_grid, 1, function(row) {
    shifted <- row + log_w
    max_val <- max(shifted)
    max_val + log(sum(exp(shifted - max_val)))
  })

  # Store per-grid BF matrix for M-step (e.g., mixsqp in susieAnn)
  if (is.null(model$lbf_grid)) {
    model$lbf_grid <- vector("list", nrow(model$alpha))
  }
  model$lbf_grid[[l]] <- lbf_grid

  # Cache ser_stats for calculate_posterior_moments_mixture_common
  model$.ser_stats <- ser_stats

  # Compute posterior inclusion probabilities using existing machinery
  stable_res  <- lbf_stabilization(lbf_mix, model$pi, shat2)
  weights_res <- compute_posterior_weights(stable_res$lpo)

  model$alpha[l, ]        <- weights_res$alpha
  model$lbf[l]            <- weights_res$lbf_model
  model$lbf_variable[l, ] <- stable_res$lbf

  return(model)
}


#' Compute mixture posterior moments
#'
#' For each grid point k, computes the conjugate normal posterior moments
#' (mean, variance) given prior variance V_k. Forms the mixture posterior
#' using responsibility weights r_jk = w_k * BF_jk / sum_k' w_k' * BF_jk'.
#'
#' Uses betahat and shat2 from ser_stats (produced by the data-type-specific
#' compute_ser_statistics), so this function is data-type-agnostic.
#'
#' @param params Params object with prior_variance_grid and mixture_weights
#' @param model Model with lbf_grid[[l]] (p x K), alpha[l,] already computed,
#'   and ser_stats cached from loglik_mixture_common
#' @param l Effect index
#'
#' @return Updated model with mu[l,] and mu2[l,]
#'
#' @keywords internal
calculate_posterior_moments_mixture_common <- function(params, model, l) {
  grid <- params$prior_variance_grid
  w    <- params$mixture_weights
  K    <- length(grid)

  lbf_grid <- model$lbf_grid[[l]]      # p x K
  betahat  <- model$.ser_stats$betahat  # cached by loglik_mixture_common
  shat2    <- pmax(model$.ser_stats$shat2, .Machine$double.eps)
  p        <- length(betahat)

  # Responsibility weights: r_jk = w_k * BF_jk / sum_k' w_k' * BF_jk'
  # Work in log space for stability
  log_w <- log(w + .Machine$double.eps)
  log_r <- sweep(lbf_grid, 2, log_w, "+")  # p x K
  log_r_max <- apply(log_r, 1, max)
  r <- exp(log_r - log_r_max)              # p x K, shifted
  r <- r / rowSums(r)                       # normalize to responsibilities

  # Per-grid posterior moments
  post_mean  <- matrix(0, p, K)
  post_mean2 <- matrix(0, p, K)
  for (k in seq_len(K)) {
    V_k <- grid[k]
    pv_k <- V_k * shat2 / (V_k + shat2)    # posterior variance
    pm_k <- pv_k / shat2 * betahat          # posterior mean
    post_mean[, k]  <- pm_k
    post_mean2[, k] <- pv_k + pm_k^2        # E[beta^2]
  }

  # Mixture posterior: weighted average over grid points
  model$mu[l, ]  <- rowSums(r * post_mean)
  model$mu2[l, ] <- rowSums(r * post_mean2)

  # Clean up cached ser_stats
  model$.ser_stats <- NULL

  return(model)
}
