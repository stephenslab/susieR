# =============================================================================
# DATA INITIALIZATION & CONFIGURATION
#
# Functions for data object setup, configuration, and preprocessing.
# These prepare data objects for model fitting and handle data-specific
# configurations like unmappable effects.
#
# Functions: configure_data, get_var_y
# =============================================================================

# Configure individual data for specified method
#' @keywords internal
configure_data.individual <- function(data, params) {
  if (params$unmappable_effects == "none") {
    return(configure_data.default(data, params))
  } else {
    warning_message("Individual-level data converted to sufficient statistics for unmappable effects methods")
    return(convert_individual_to_ss(data, params))
  }
}

# Get variance of y
#' @keywords internal
#' @importFrom stats var
get_var_y.individual <- function(data, ...) {
  return(var(drop(data$y)))
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
initialize_susie_model.individual <- function(data, params, var_y, ...) {

  # Base model
  model <- initialize_matrices(data, params, var_y)

  # Append predictor weights
  model$predictor_weights <- attr(data$X, "d")

  # Initialize Servin-Stephens parameters
  if (params$use_servin_stephens) {
    model$rv <- rep(1, params$L)
    model$marginal_loglik <- rep(as.numeric(NA), params$L)
  }

  return(model)
}

# Initialize fitted values
#' @keywords internal
initialize_fitted.individual <- function(data, mat_init) {
  return(list(Xr = compute_Xb(data$X, colSums(mat_init$alpha * mat_init$mu))))
}

# Validate prior variance
#' @keywords internal
validate_prior.individual <- function(data, params, model, ...) {
  return(validate_prior.default(data, params, model, ...))
}

# Track core parameters across iterations
#' @keywords internal
track_ibss_fit.individual <- function(data, params, model, tracking, iter, elbo, ...) {
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
compute_residuals.individual <- function(data, params, model, l, ...) {
  # Remove lth effect from fitted values
  Xr_without_l <- model$Xr - compute_Xb(data$X, model$alpha[l, ] * model$mu[l, ])

  # Compute residuals
  R   <- data$y - Xr_without_l
  XtR <- compute_Xty(data$X, R)

  # Store unified residuals in model
  model$residuals         <- XtR
  model$fitted_without_l  <- Xr_without_l
  model$raw_residuals     <- R
  model$residual_variance <- model$sigma2  # Standard residual variance

  return(model)
}

# Compute SER statistics
#' @keywords internal
compute_ser_statistics.individual <- function(data, params, model, l, ...) {
  betahat <- (1 / model$predictor_weights) * model$residuals
  shat2   <- model$residual_variance / model$predictor_weights

  # Optimization parameters
  optim_init   <- log(max(c(betahat^2 - shat2, 1), na.rm = TRUE))
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

# Posterior expected log-likelihood for single effect regression
#' @keywords internal
SER_posterior_e_loglik.individual <- function(data, params, model, l) {
  Eb  <- model$alpha[l, ] * model$mu[l, ]
  Eb2 <- model$alpha[l, ] * model$mu2[l, ]
  return(-0.5 * data$n * log(2 * pi * model$sigma2) -
           0.5 / model$sigma2 * (sum(model$raw_residuals * model$raw_residuals)
                                 - 2 * sum(model$raw_residuals * compute_Xb(data$X, Eb)) +
                                   sum(model$predictor_weights * Eb2)))
}

# Calculate posterior moments for single effect regression
#' @keywords internal
calculate_posterior_moments.individual <- function(data, params, model, V, l, ...) {
  if (params$use_servin_stephens) {
    if (V <= 0) {
      # Zero variance case
      post_mean  <- rep(0, data$p)
      post_mean2 <- rep(0, data$p)
      model$rv[l] <- 1
    } else {
      # Compute posterior moments for NIG prior
      moments <- compute_posterior_moments_NIG(data$n, model$predictor_weights,
                                               model$residuals, sum(model$raw_residuals^2),
                                               drop(cor(data$X, model$raw_residuals)),
                                               V, params$alpha0, params$beta0)

      post_mean  <- moments$post_mean
      post_mean2 <- moments$post_mean2

      # Compute weighted average of residual variance modes using PIPs
      model$rv[l] <- sum(model$alpha[l, ] * moments$rv)
    }
  } else {
    # Standard Gaussian posterior calculations
    post_var   <- (1 / V + model$predictor_weights / model$residual_variance)^(-1)
    post_mean  <- (1 / model$residual_variance) * post_var * model$residuals
    post_mean2 <- post_var + post_mean^2
  }

  # Store posterior moments in model
  model$mu[l, ] <- post_mean
  model$mu2[l, ] <- post_mean2

  return(model)
}

# Calculate KL divergence
#' @keywords internal
compute_kl.individual <- function(data, params, model, l) {
  if (params$use_servin_stephens) {
    # KL divergence for Servin-Stephens is only valid for L=1
    if (params$L == 1) {
      kl <- compute_kl_NIG(model$alpha[l, ], model$mu[l, ], model$mu2[l, ], model$pi, model$V[l],
                           a0 = 1, b0 = 1, a_post = data$n / 2 + 1, b_post = data$n * model$sigma2 + 1)
    } else {
      kl <- 0
    }
  } else {
    # Standard Gaussian KL divergence
    loglik_term <- model$lbf[l] + sum(dnorm(model$raw_residuals, 0, sqrt(model$sigma2), log = TRUE))
    kl <- -loglik_term + SER_posterior_e_loglik(data, params, model, l)
  }

  # Store in model and return
  model$KL[l] <- kl
  return(model)
}

# Expected squared residuals
#' @keywords internal
get_ER2.individual <- function(data, model) {
  Xr_L <- compute_MXt(model$alpha * model$mu, data$X)
  postb2 <- model$alpha * model$mu2
  return(sum((data$y - model$Xr)^2) - sum(Xr_L^2) + sum(model$predictor_weights * t(postb2)))
}

# Expected log-likelihood
#' @keywords internal
Eloglik.individual <- function(data, model) {
  return(-data$n / 2 * log(2 * pi * model$sigma2) -
           1 / (2 * model$sigma2) * get_ER2(data, model))
}

#' @importFrom Matrix colSums
#' @importFrom stats dnorm
#' @importFrom stats cor
#' @keywords internal
loglik.individual <- function(data, params, model, V, ser_stats, l = NULL, ...) {
  # Check if using Servin-Stephens prior
  if (params$use_servin_stephens) {
    # Compute log Bayes factors for NIG prior
    lbf <- compute_lbf_NIG(data$n, model$predictor_weights,
                           model$residuals, sum(model$raw_residuals^2),
                           drop(cor(data$X, model$raw_residuals)),
                           V, params$alpha0, params$beta0)
  } else {
    # Standard Gaussian prior log Bayes factors
    lbf <- dnorm(ser_stats$betahat, 0, sqrt(V + ser_stats$shat2), log = TRUE) -
      dnorm(ser_stats$betahat, 0, sqrt(ser_stats$shat2), log = TRUE)
  }

  # Stabilize logged Bayes Factor
  stable_res  <- lbf_stabilization(lbf, model$pi, ser_stats$shat2)

  # Compute posterior weights
  weights_res <- compute_posterior_weights(stable_res$lpo)

  # Store in model if l is provided, otherwise return lbf_model for prior variance optimization
  if (!is.null(l)) {
    # Store results in model
    model$alpha[l, ] <- weights_res$alpha
    model$lbf[l] <- weights_res$lbf_model
    model$lbf_variable[l, ] <- stable_res$lbf

    # Compute and store marginal log-likelihood for NIG prior
    model$marginal_loglik[l] <- compute_marginal_loglik(weights_res$lbf_model, data$n,
                                                         sum(model$raw_residuals^2),
                                                         params$alpha0, params$beta0,
                                                         params$use_servin_stephens)
    return(model)
  } else {
    return(weights_res$lbf_model)
  }
}

#' @keywords internal
neg_loglik.individual <- function(data, params, model, V_param, ser_stats, ...) {
  # Convert parameter to V based on optimization scale (always log for individual)
  V <- exp(V_param)
  lbf_model <- loglik.individual(data, params, model, V, ser_stats)
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
update_fitted_values.individual <- function(data, params, model, l, ...) {
  model$Xr <- model$fitted_without_l + compute_Xb(data$X, model$alpha[l, ] * model$mu[l, ])

  return(model)
}

# Update variance components for individual data
#' @keywords internal
update_variance_components.individual <- function(data, params, model, ...) {
  return(update_variance_components.default(data, params, model, ...))
}

# Update derived quantities for individual data
#' @keywords internal
update_derived_quantities.individual <- function(data, params, model) {
  return(update_derived_quantities.default(data, params, model))
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

# Get column scale factors
#' @keywords internal
get_scale_factors.individual <- function(data, params) {
  return(attr(data$X, "scaled:scale"))
}

# Get intercept
#' @keywords internal
get_intercept.individual <- function(data, params, model, ...) {
  if (params$intercept) {
    return(data$mean_y - sum(attr(data$X, "scaled:center") *
                               (colSums(model$alpha * model$mu) / attr(data$X, "scaled:scale"))))
  } else {
    return(0)
  }
}

# Get Fitted Values
#' @keywords internal
get_fitted.individual <- function(data, params, model, ...) {
  if (params$intercept) {
    fitted <- model$Xr + data$mean_y
  } else {
    fitted <- model$Xr
  }

  fitted <- drop(fitted)
  names(fitted) <- `if`(is.null(names(data$y)), rownames(data$X), names(data$y))

  return(fitted)
}

# Get Credible Sets
#' @keywords internal
get_cs.individual <- function(data, params, model, ...) {
  if (is.null(params$coverage) || is.null(params$min_abs_corr)) {
    return(NULL)
  }

  return(susie_get_cs(model,
                      X            = data$X,
                      coverage     = params$coverage,
                      min_abs_corr = params$min_abs_corr,
                      n_purity     = params$n_purity))
}


# Get Variable Names
#' @keywords internal
get_variable_names.individual <- function(data, model, ...) {
  return(assign_names(data, model, colnames(data$X)))
}

# Get univariate z-score
#' @keywords internal
get_zscore.individual <- function(data, params, model, ...) {
  if (isFALSE(params$compute_univariate_zscore)) {
    return(get_zscore.default(data, params, model))
  }

  X <- data$X

  if (!is.matrix(X)) {
    warning_message(
      "Calculation of univariate regression z-scores is not ",
      "implemented specifically for sparse or trend filtering ",
      "matrices, so this step may be slow if the matrix is large; ",
      "to skip this step set compute_univariate_zscore = FALSE"
    )
  }
  if (!is.null(model$null_weight) && model$null_weight != 0) {
    X <- X[, 1:(ncol(X) - 1)]
  }

  return(calc_z(X, data$y, center = params$intercept, scale = params$standardize))
}

# Clean up model object for individual data
#' @keywords internal
cleanup_model.individual <- function(data, params, model, ...) {
  # Remove common fields
  model <- cleanup_model.default(data, params, model, ...)

  # Remove individual-specific temporary fields
  model$raw_residuals <- NULL

  # Remove Servin-stephens specific temporary fields
  if (params$use_servin_stephens) {
    model$marginal_loglik <- NULL
    if (nrow(model$alpha) > 1) model$elbo <- NULL
  }

  return(model)
}
