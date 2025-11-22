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
  D    <- data$eigen_R$values
  V    <- data$eigen_R$vectors
  Dinv <- compute_Dinv(model, data)

  model$SinvRj   <- V %*% (Dinv * D * t(V))
  model$RjSinvRj <- colSums(t(V) * (Dinv * D^2 * t(V)))

  return(model)
}

# Initialize fitted values
#' @keywords internal
initialize_fitted.rss_lambda <- function(data, mat_init) {
  return(list(Rz = as.vector(data$R %*% colSums(mat_init$alpha * mat_init$mu))))
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
  Rz_without_l <- model$Rz - data$R %*% (model$alpha[l, ] * model$mu[l, ])

  # Compute residuals
  r <- data$z - Rz_without_l

  # Store unified residuals in model
  model$residuals         <- r
  model$fitted_without_l  <- Rz_without_l
  model$residual_variance <- 1  # RSS lambda uses normalized residual variance

  return(model)
}

# Compute SER statistics
#' @keywords internal
compute_ser_statistics.rss_lambda <- function(data, params, model, l, ...) {
  shat2 <- 1 / model$RjSinvRj

  # Optimization parameters
  init_vals    <- sapply(1:data$p, function(j) sum(model$SinvRj[, j] * model$residuals)^2) - (1 / model$RjSinvRj)
  optim_init   <- log(max(c(init_vals, 1e-6), na.rm = TRUE))
  optim_bounds <- c(-30, 15)
  optim_scale  <- "log"

  return(list(
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
  V      <- data$eigen_R$vectors
  Dinv   <- compute_Dinv(model, data)
  rR     <- data$R %*% model$residuals
  SinvEb <- V %*% (Dinv * crossprod(V, Eb))

  return(-0.5 * (-2 * sum(rR * SinvEb) + sum(model$RjSinvRj * Eb2)))
}

# Calculate posterior moments for single effect regression
#' @keywords internal
calculate_posterior_moments.rss_lambda <- function(data, params, model, V, loglik_res, ...) {
  post_var   <- (model$RjSinvRj + 1 / V)^(-1)
  post_mean  <- sapply(1:data$p, function(j) {
    post_var[j] * sum(model$SinvRj[, j] * model$residuals)
  })
  post_mean2 <- post_var + post_mean^2

  return(list(
    post_mean  = post_mean,
    post_mean2 = post_mean2,
    post_var   = post_var,
    rv         = 1
  ))
}

# Calculate KL divergence
#' @keywords internal
compute_kl.rss_lambda <- function(data, params, model, l) {
  return(compute_kl.default(data, params, model, l))
}

# Expected squared residuals
#' @keywords internal
get_ER2.rss_lambda <- function(data, model) {
  # Eigen decomposition components
  D     <- data$eigen_R$values
  V     <- data$eigen_R$vectors
  Dinv  <- compute_Dinv(model, data)

  # Cached quantities
  Vtz   <- data$Vtz
  zbar  <- model$zbar
  postb2 <- model$diag_postb2

  # z^T S^{-1} z
  zSinvz <- sum((Dinv * Vtz) * Vtz)

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
  D <- data$eigen_R$values
  d <- model$sigma2 * D + data$lambda
  return(-(length(data$z) / 2) * log(2 * pi) - 0.5 *
           sum(log(d)) - 0.5 * get_ER2.rss_lambda(data, model))
}

# Log-likelihood for RSS
#' @keywords internal
loglik.rss_lambda <- function(data, params, model, V, ser_stats, ...) {
  # Compute log Bayes factors
  lbf <- -0.5 * log(1 + V / ser_stats$shat2) +
    0.5 * (V / (1 + V / ser_stats$shat2)) * (crossprod(model$SinvRj, model$residuals)^2)

  # Stabilize logged Bayes Factor
  stable_res <- lbf_stabilization(lbf, model$pi, ser_stats$shat2)

  # Compute posterior weights
  weights_res <- compute_posterior_weights(stable_res$lpo)

  return(list(
    lbf       = stable_res$lbf,
    lbf_model = weights_res$lbf_model,
    alpha     = weights_res$alpha
  ))
}

#' @keywords internal
neg_loglik.rss_lambda <- function(data, params, model, V_param, ser_stats, ...) {
  # Convert parameter to V based on optimization scale (always log for RSS lambda)
  V   <- exp(V_param)
  res <- loglik.rss_lambda(data, params, model, V, ser_stats)
  return(-res$lbf_model)
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
update_fitted_values.rss_lambda <- function(data, params, model, l, res, ...) {
  model$Rz <- model$fitted_without_l + as.vector(data$R %*% (model$alpha[l, ] * model$mu[l, ]))
  model    <- precompute_rss_lambda_terms(data, model)

  return(model)
}

# Update variance components
#' @keywords internal
#' @importFrom stats optimize
update_variance_components.rss_lambda <- function(data, params, model, ...) {
  upper_bound <- 1 - data$lambda

  objective <- function(sigma2) {
    temp_model        <- model
    temp_model$sigma2 <- sigma2
    Eloglik.rss_lambda(data, temp_model)
  }

  est_sigma2 <- optimize(objective, interval = c(1e-4, upper_bound), maximum = TRUE)$maximum

  if (objective(est_sigma2) < objective(upper_bound)) {
    est_sigma2 <- upper_bound
  }

  list(sigma2 = est_sigma2)
}

# Update derived quantities
#' @keywords internal
update_derived_quantities.rss_lambda <- function(data, params, model) {
  # Recalculate Dinv with updated sigma2
  Dinv <- compute_Dinv(model, data)
  V    <- data$eigen_R$vectors
  D    <- data$eigen_R$values

  # Update SinvRj and RjSinvRj
  model$SinvRj   <- V %*% (Dinv * D * t(V))
  model$RjSinvRj <- colSums(t(V) * (Dinv * (D^2) * t(V)))

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

  return(susie_get_cs(model,
                      Xcorr           = muffled_cov2cor(data$R),
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
  rss_fields <- c("SinvRj", "RjSinvRj", "Rz", "Z", "zbar", "diag_postb2")

  for (field in rss_fields) {
    if (field %in% names(model)) {
      model[[field]] <- NULL
    }
  }

  return(model)
}
