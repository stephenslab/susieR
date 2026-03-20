# =============================================================================
# DATA INITIALIZATION & CONFIGURATION
#
# Functions for data object setup, configuration, and preprocessing.
# These prepare data objects for model fitting and handle data-specific
# configurations like unmappable effects.
#
# Functions: configure_data, get_var_y
# =============================================================================

# Configure data
#' @keywords internal
configure_data.rss_lambda <- function(data, params) {
  return(configure_data.default(data, params))
}

# Get variance of y
#' @keywords internal
get_var_y.rss_lambda <- function(data, ...) {
  return(1)
}

# =============================================================================
# MODEL INITIALIZATION & SETUP
#
# Functions for initializing model objects and setting up initial states.
# These create model matrices, initialize fitted values, and prepare
# the SuSiE model for iterative fitting.
#
# Functions: initialize_susie_model, initialize_fitted, validate_prior, track_ibss_fit
# =============================================================================

# Initialize SuSiE model
#' @keywords internal
initialize_susie_model.rss_lambda <- function(data, params, var_y, ...) {
  # Base model
  model <- initialize_matrices(data, params, var_y)

  # Initialize SinvRj and RjSinvRj
  eigen_R <- get_eigen_R(data, model)
  D    <- eigen_R$values
  V    <- eigen_R$vectors
  tV   <- t(V)
  Dinv <- compute_Dinv(model, data)

  model$SinvRj   <- V %*% (Dinv * D * tV)
  model$RjSinvRj <- colSums(tV * (Dinv * D^2 * tV))

  return(model)
}

# Initialize fitted values
#' @keywords internal
initialize_fitted.rss_lambda <- function(data, mat_init) {
  return(list(Rz = as.vector(compute_Rv(data, colSums(mat_init$alpha * mat_init$mu)))))
}

# Validate prior variance
#' @keywords internal
validate_prior.rss_lambda <- function(data, params, model, ...) {
  return(validate_prior.default(data, params, model, ...))
}

# Track core parameters for tracking
#' @keywords internal
track_ibss_fit.rss_lambda <- function(data, params, model, tracking, iter, elbo, ...) {
  return(track_ibss_fit.default(data, params, model, tracking, iter, elbo, ...))
}

# =============================================================================
# SINGLE EFFECT REGRESSION & ELBO
#
# Core functions for single effect regression computation and ELBO calculation.
# These handle the mathematical core of SuSiE including residual computation, SER
# statistics, posterior moments, and log-likelihood calculations for the ELBO.
#
# Functions: compute_residuals, compute_ser_statistics, SER_posterior_e_loglik,
# calculate_posterior_moments, compute_kl, get_ER2, Eloglik, loglik, neg_loglik
# =============================================================================

# Compute residuals for single effect regression
#' @keywords internal
compute_residuals.rss_lambda <- function(data, params, model, l, ...) {
  # Remove lth effect from fitted values
  Rz_without_l <- model$Rz - compute_Rv(data, model$alpha[l, ] * model$mu[l, ], model)

  # Compute residuals
  r <- data$z - Rz_without_l

  # Store unified residuals in model
  model$residuals         <- r
  model$fitted_without_l  <- Rz_without_l
  model$residual_variance <- 1  # RSS lambda uses normalized residual variance

  # Dynamic stochastic LD variance inflation (z-score scale)
  if (!is.null(data$stochastic_ld_B)) {
    b_minus_l <- colSums(model$alpha * model$mu) - model$alpha[l, ] * model$mu[l, ]
    model$shat2_inflation <- compute_shat2_inflation_rss(
      data, model, Rz_without_l, b_minus_l)
  }

  return(model)
}

# Compute SER statistics
#' @keywords internal
compute_ser_statistics.rss_lambda <- function(data, params, model, l, ...) {
  signal  <- as.vector(crossprod(model$SinvRj, model$residuals))
  shat2   <- 1 / model$RjSinvRj
  betahat <- signal * shat2

  # Apply stochastic LD inflation to shat2
  if (!is.null(model$shat2_inflation))
    shat2 <- shat2 * model$shat2_inflation

  # Optimization parameters
  optim_init   <- log(max(c(betahat^2 - shat2, 1e-6), na.rm = TRUE))
  optim_bounds <- c(-30, 15)
  optim_scale  <- "log"

  return(list(
    betahat      = betahat,
    shat2        = shat2,
    optim_init   = optim_init,
    optim_bounds = optim_bounds,
    optim_scale  = optim_scale
  ))
}

# SER posterior expected log-likelihood
#' @keywords internal
SER_posterior_e_loglik.rss_lambda <- function(data, params, model, l) {
  Eb     <- model$alpha[l, ] * model$mu[l, ]
  Eb2    <- model$alpha[l, ] * model$mu2[l, ]
  eigen_R <- get_eigen_R(data, model)
  V      <- eigen_R$vectors
  Dinv   <- compute_Dinv(model, data)
  rR     <- compute_Rv(data, model$residuals, model)
  SinvEb <- V %*% (Dinv * crossprod(V, Eb))

  return(-0.5 * (-2 * sum(rR * SinvEb) + sum(model$RjSinvRj * Eb2)))
}

# Calculate posterior moments for single effect regression
#' @keywords internal
calculate_posterior_moments.rss_lambda <- function(data, params, model, V, l, ...) {
  shat2 <- 1 / model$RjSinvRj
  if (!is.null(model$shat2_inflation))
    shat2 <- shat2 * model$shat2_inflation

  post_var  <- V * shat2 / (V + shat2)
  signal    <- as.vector(crossprod(model$SinvRj, model$residuals))
  betahat   <- signal * (1 / model$RjSinvRj)
  post_mean <- post_var / shat2 * betahat
  post_mean2 <- post_var + post_mean^2

  # Store posterior moments in model
  model$mu[l, ] <- post_mean
  model$mu2[l, ] <- post_mean2

  return(model)
}

# Calculate KL divergence
#' @keywords internal
compute_kl.rss_lambda <- function(data, params, model, l) {
  model <- compute_kl.default(data, params, model, l)
  return(model)
}

# Expected squared residuals
#' @keywords internal
get_ER2.rss_lambda <- function(data, model) {
  # Eigen decomposition components
  eigen_R <- get_eigen_R(data, model)
  D     <- eigen_R$values
  V     <- eigen_R$vectors
  Dinv  <- compute_Dinv(model, data)

  # Cached quantities
  Vtz   <- get_Vtz(data, model)
  zbar  <- model$zbar
  postb2 <- model$diag_postb2

  # z^T S^{-1} z (use model z_null_norm2 if omega changed, else data)
  z_null_norm2 <- if (!is.null(model$z_null_norm2)) model$z_null_norm2 else data$z_null_norm2
  zSinvz <- sum((Dinv * Vtz) * Vtz) + z_null_norm2 / data$lambda

  # -2 zbar^T S^{-1} z
  tmp <- V %*% (Dinv * (D * Vtz))
  term2 <- -2 * sum(tmp * zbar)

  # zbar^T R S^{-1} R zbar
  Vtzbar <- crossprod(V, zbar)
  term3 <- sum((Vtzbar^2) * (Dinv * D^2))

  # RZ2 = sum((Z %*% RSinvR) * Z)
  VtZ <- model$Z %*% V
  term4 <- sum((VtZ^2) %*% (Dinv * D^2))

  # diag(RSinvR)^T postb2
  diag_RSinvR <- rowSums((V^2) * rep(Dinv * D^2, each = nrow(V)))
  term5 <- sum(diag_RSinvR * postb2)

  return(zSinvz + term2 + term3 - term4 + term5)
}

# Expected log-likelihood
#' @keywords internal
Eloglik.rss_lambda <- function(data, model) {
  D <- get_eigen_R(data, model)$values
  d <- model$sigma2 * D + data$lambda
  return(-(length(data$z) / 2) * log(2 * pi) - 0.5 *
           sum(log(d)) - 0.5 * get_ER2.rss_lambda(data, model))
}

# Log-likelihood for RSS
#' @keywords internal
loglik.rss_lambda <- function(data, params, model, V, ser_stats, l = NULL, ...) {
  # Wakefield ABF using betahat/shat2 from ser_stats (supports inflation)
  shat2 <- pmax(ser_stats$shat2, .Machine$double.eps)
  lbf   <- -0.5 * log(1 + V / shat2) +
    0.5 * ser_stats$betahat^2 * V / (shat2 * (V + shat2))

  # Stabilize logged Bayes Factor
  stable_res <- lbf_stabilization(lbf, model$pi, ser_stats$shat2)

  # Compute posterior weights
  weights_res <- compute_posterior_weights(stable_res$lpo)

  # Store in model if l is provided, otherwise return lbf_model for prior variance optimization
  if (!is.null(l)) {
    model$alpha[l, ] <- weights_res$alpha
    model$lbf[l] <- weights_res$lbf_model
    model$lbf_variable[l, ] <- stable_res$lbf
    return(model)
  } else {
    return(weights_res$lbf_model)
  }
}

#' @keywords internal
neg_loglik.rss_lambda <- function(data, params, model, V_param, ser_stats, ...) {
  # Convert parameter to V based on optimization scale (always log for RSS lambda)
  V <- exp(V_param)
  lbf_model <- loglik.rss_lambda(data, params, model, V, ser_stats)
  return(-lbf_model)
}

# =============================================================================
# MODEL UPDATES & FITTING
#
# Functions for iterative model updates and variance component estimation.
# These handle the dynamic aspects of model fitting including fitted value
# updates and variance component estimation.
#
# Functions: update_fitted_values, update_variance_components, update_derived_quantities
# =============================================================================

# Update fitted values
#' @keywords internal
update_fitted_values.rss_lambda <- function(data, params, model, l, ...) {
  model$Rz <- model$fitted_without_l + as.vector(compute_Rv(data, model$alpha[l, ] * model$mu[l, ], model))
  model    <- precompute_rss_lambda_terms(data, model)

  return(model)
}

# Update model variance (override to always run omega update for multi-panel)
#' @keywords internal
update_model_variance.rss_lambda <- function(data, params, model) {
  need_sigma2 <- isTRUE(params$estimate_residual_variance)
  need_omega  <- !is.null(data$K) && data$K > 1

  if (!need_sigma2 && !need_omega) return(model)

  variance_result <- update_variance_components(data, params, model)
  model <- modifyList(model, variance_result)
  if (need_sigma2) {
    model$sigma2 <- min(max(model$sigma2, params$residual_variance_lowerbound),
                        params$residual_variance_upperbound)
  }
  model <- update_derived_quantities(data, params, model)

  return(model)
}

# Update variance components
#' @keywords internal
#' @importFrom stats optimize
update_variance_components.rss_lambda <- function(data, params, model, ...) {
  result <- list()

  # Sigma2 estimation (only if requested)
  if (isTRUE(params$estimate_residual_variance)) {
    upper_bound <- 1 - data$lambda
    objective <- function(sigma2) {
      temp_model        <- model
      temp_model$sigma2 <- sigma2
      Eloglik.rss_lambda(data, temp_model)
    }
    est_sigma2 <- optimize(objective, interval = c(1e-4, upper_bound),
                           maximum = TRUE)$maximum
    if (objective(est_sigma2) < objective(upper_bound))
      est_sigma2 <- upper_bound
    result$sigma2 <- est_sigma2
  }

  # Multi-panel omega update via profile Eloglik (M-step of variational EM)
  # Maximizing Eloglik over omega on the simplex is guaranteed to not decrease
  # the ELBO, since KL doesn't depend on omega. The log-likelihood is concave
  # in omega, so this is a convex optimization problem with a unique global optimum.
  #
  # K=2: Brent's method (optimize) on omega_1 in [0,1]
  # K>2: Frank-Wolfe (conditional gradient) on simplex
  #
  # Uses eval_omega_eloglik (Rcpp) which combines eigendecomposition + Eloglik
  # evaluation in a single C++ call, avoiding costly R/C++ boundary crossing
  # for the p x p eigenvector matrix. Falls back to pure R if Rcpp unavailable.
  if (!is.null(data$K) && data$K > 1 && !is.null(data$panel_R)) {
    sigma2_cur <- if (!is.null(result$sigma2)) result$sigma2 else model$sigma2
    omega_cur  <- if (!is.null(model$omega)) model$omega else rep(1 / data$K, data$K)

    # Combined eigen(R(omega)) + Eloglik evaluator
    # Prefers Rcpp (avoids p x p matrix marshaling); falls back to pure R
    has_rcpp <- exists("eval_omega_eloglik", mode = "function")
    eval_omega <- function(omega_vec) {
      if (has_rcpp) {
        eval_omega_eloglik(data$panel_R, omega_vec, data$z, model$zbar,
                           model$diag_postb2, model$Z, sigma2_cur, data$lambda,
                           data$K, data$p)
      } else {
        eval_omega_eloglik_R(data$panel_R, omega_vec, data$z, model$zbar,
                              model$diag_postb2, model$Z, sigma2_cur, data$lambda,
                              data$K, data$p)
      }
    }

    if (data$K == 2) {
      opt <- optimize(function(w1) eval_omega(c(w1, 1 - w1)),
                      interval = c(0, 1), maximum = TRUE)
      result$omega <- c(opt$maximum, 1 - opt$maximum)
    } else {
      # Frank-Wolfe on simplex: each iteration picks best vertex, line-searches
      omega <- omega_cur
      cur_val <- eval_omega(omega)

      for (fw_iter in seq_len(10)) {
        vertex_vals <- vapply(seq_len(data$K), function(k) {
          e_k <- rep(0, data$K); e_k[k] <- 1; eval_omega(e_k)
        }, numeric(1))

        k_star <- which.max(vertex_vals)
        s <- rep(0, data$K); s[k_star] <- 1

        opt_gamma <- optimize(
          function(gamma) eval_omega((1 - gamma) * omega + gamma * s),
          interval = c(0, 1), maximum = TRUE)

        omega_new <- (1 - opt_gamma$maximum) * omega + opt_gamma$maximum * s
        new_val <- opt_gamma$objective

        if (new_val - cur_val < 1e-6) break
        omega <- omega_new
        cur_val <- new_val
      }
      result$omega <- omega
    }
  }

  result
}

# Update derived quantities
#' @keywords internal
update_derived_quantities.rss_lambda <- function(data, params, model) {
  # Multi-panel omega update: re-eigendecompose R(omega) from cached panel R matrices
  if (!is.null(data$K) && data$K > 1 && !is.null(model$omega)) {
    if (!is.null(data$panel_R)) {
      # Fast path: eigendecompose R(omega) = sum_k omega_k R_k via Rcpp
      model$eigen_R <- eigen_R_omega(data$panel_R, model$omega, data$K, data$p)
    } else {
      # Fallback: form X_meta and SVD
      X_meta <- form_X_meta(data$X_list, model$omega)
      model$X_meta  <- X_meta
      model$eigen_R <- eigen_from_X(X_meta, data$p)
    }
    model$Vtz          <- crossprod(model$eigen_R$vectors, data$z)
    model$z_null_norm2 <- max(sum(data$z^2) - sum(model$Vtz^2), 0)
    model$stochastic_ld_B <- 1 / sum(model$omega^2 / data$B_list)
  }

  # Recalculate Dinv with updated sigma2 (and potentially updated eigen_R)
  eigen_R <- get_eigen_R(data, model)
  Dinv <- compute_Dinv(model, data)
  V    <- eigen_R$vectors
  D    <- eigen_R$values
  tV   <- t(V)

  # Update SinvRj and RjSinvRj
  model$SinvRj   <- V %*% (Dinv * D * tV)
  model$RjSinvRj <- colSums(tV * (Dinv * (D^2) * tV))

  return(model)
}

# =============================================================================
# OUTPUT GENERATION & POST-PROCESSING
#
# Functions for generating final results and summary statistics.
# These process fitted models into interpretable outputs including
# credible sets, variable names, and fitted values.
#
# Functions: get_scale_factors, get_intercept, get_fitted, get_cs,
# get_variable_names, get_zscore
# =============================================================================

# Get scale factors
#' @keywords internal
get_scale_factors.rss_lambda <- function(data, params) {
  return(rep(1, data$p))
}

# Get intercept
#' @keywords internal
get_intercept.rss_lambda <- function(data, params, model, ...) {
  return(data$intercept_value)
}

# Get fitted values
#' @keywords internal
get_fitted.rss_lambda <- function(data, params, model, ...) {
  return(get_fitted.default(data, params, model, ...))
}

# Get credible sets
#' @keywords internal
get_cs.rss_lambda <- function(data, params, model, ...) {
  if (is.null(params$coverage) || is.null(params$min_abs_corr)) {
    return(NULL)
  }

  # Use current X_meta (multi-panel) or data$X (single-panel)
  X <- if (!is.null(model$X_meta)) model$X_meta else data$X
  if (!is.null(X)) {
    return(susie_get_cs(model,
                        X               = X,
                        coverage        = params$coverage,
                        min_abs_corr    = params$min_abs_corr,
                        n_purity        = params$n_purity))
  }

  return(susie_get_cs(model,
                      Xcorr           = safe_cov2cor(data$R),
                      check_symmetric = FALSE,
                      coverage        = params$coverage,
                      min_abs_corr    = params$min_abs_corr,
                      n_purity        = params$n_purity))
}

# Get variable names
#' @keywords internal
get_variable_names.rss_lambda <- function(data, model, ...) {
  return(assign_names(data, model, names(data$z)))
}

# Get univariate z-scores
#' @keywords internal
get_zscore.rss_lambda <- function(data, params, model, ...) {
  return(get_zscore.default(data, params, model))
}

# Clean up model object for RSS lambda data
#' @keywords internal
cleanup_model.rss_lambda <- function(data, params, model, ...) {
  # Remove common fields
  model <- cleanup_model.default(data, params, model, ...)

  # Remove RSS-lambda-specific temporary fields
  rss_fields <- c("SinvRj", "RjSinvRj", "Rz", "Z", "zbar", "diag_postb2",
                   "X_meta", "eigen_R", "Vtz", "omega", "stochastic_ld_B")

  for (field in rss_fields) {
    if (field %in% names(model)) {
      model[[field]] <- NULL
    }
  }

  return(model)
}
