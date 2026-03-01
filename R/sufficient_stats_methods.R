# =============================================================================
# DATA INITIALIZATION & CONFIGURATION
#
# Functions for data object setup, configuration, and preprocessing.
# These prepare data objects for model fitting and handle data-specific
# configurations like unmappable effects.
#
# Functions: configure_data, get_var_y
# =============================================================================

# Configure ss data for specified method
#' @keywords internal
configure_data.ss <- function(data, params) {
  if (params$unmappable_effects == "inf") {
    return(add_eigen_decomposition(data, params))
  } else {
    return(configure_data.default(data, params))
  }
}

# Get variance of y
#' @keywords internal
get_var_y.ss <- function(data, ...) {
  return(data$yty / (data$n - 1))
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
initialize_susie_model.ss <- function(data, params, var_y, ...) {

  # Base model
  model <- initialize_matrices(data, params, var_y)

  # Append predictor weights and initialize non-sparse quantities
  if (params$unmappable_effects == "inf") {
    # Initialize omega quantities for unmappable effects
    omega_res               <- compute_omega_quantities(data, tau2 = 0, sigma2 = var_y)
    model$omega_var         <- omega_res$omega_var
    model$predictor_weights <- omega_res$diagXtOmegaX
    model$XtOmegay          <- data$eigen_vectors %*% (data$VtXty / omega_res$omega_var)

    # Initialize unmappable variance component and coefficients
    model$tau2  <- 0
    model$theta <- rep(0, data$p)

  } else if (params$unmappable_effects == "ash") {
    pm <- if (!is.null(data$XtX)) data$XtX else data$X
    model$predictor_weights <- attr(pm, "d")
    model <- init_ash_fields(model, data$n, data$p, params$L, is_individual = FALSE)
  } else {
    pm <- if (!is.null(data$XtX)) data$XtX else data$X
    model$predictor_weights <- attr(pm, "d")

    # Initialize Servin-Stephens (NIG) parameters
    if (params$use_servin_stephens) {
      model$rv <- rep(1, params$L)
      model$marginal_loglik <- rep(as.numeric(NA), params$L)
    }
  }

  return(model)
}

# Initialize fitted values
#' @keywords internal
initialize_fitted.ss <- function(data, mat_init) {
  return(list(XtXr = compute_Rv(data, colSums(mat_init$alpha * mat_init$mu))))
}

# Validate Prior Variance
#' @keywords internal
validate_prior.ss <- function(data, params, model, ...) {
  if (isTRUE(params$check_prior)) {
    if (is.null(data$zm)) {
      bhat <- data$Xty / model$predictor_weights
      shat <- sqrt(model$sigma2 / model$predictor_weights)
      z <- bhat / shat
      data$zm <- max(abs(z[!is.nan(z)]))
    }
    if (any(model$V > 100 * (data$zm^2))) {
      stop(
        "Estimated prior variance is unreasonably large.\n",
        "This usually caused by mismatch between the summary statistics and the LD Matrix.\n",
        "Please check the input."
      )
    }
  }
  return(validate_prior.default(data, params, model, ...))
}

# Track core parameters across iterations
#' @keywords internal
track_ibss_fit.ss <- function(data, params, model, tracking, iter, elbo, ...) {
  if (params$unmappable_effects %in% c("inf", "ash")) {
    # Append non-sparse variance component to tracking
    tracking <- track_ibss_fit.default(data, params, model, tracking, iter, elbo, ...)
    if (isTRUE(params$track_fit)) {
      tracking[[iter]]$tau2 <- model$tau2
    }
    return(tracking)
  } else {
    # Use default for standard SS case
    return(track_ibss_fit.default(data, params, model, tracking, iter, elbo, ...))
  }
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
compute_residuals.ss <- function(data, params, model, l, ...) {
  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: Omega-weighted residuals
    b <- colSums(model$mu * model$alpha) - model$mu[l, ] * model$alpha[l, ]

    # Compute Residuals
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    XtOmegay <- data$eigen_vectors %*% (data$VtXty / omega_res$omega_var)
    XtOmegaXb <- data$eigen_vectors %*% ((t(data$eigen_vectors) %*% b) * data$eigen_values / omega_res$omega_var)
    XtOmegar <- XtOmegay - XtOmegaXb

    # Store residuals and parameters (unmappable case)
    model$residuals         <- XtOmegar
    model$predictor_weights <- omega_res$diagXtOmegaX  # Update for this iteration
    model$residual_variance <- 1                       # Already incorporated in Omega

    return(model)

  } else if (params$unmappable_effects == "ash") {
    # SuSiE-ash: explicit residualization
    XtXr_without_l <- model$XtXr - compute_Rv(data, model$alpha[l, ] * model$mu[l, ])

    # Subtract X'X*theta from residuals
    XtR <- data$Xty - model$XtX_theta - XtXr_without_l

    # Store residuals and parameters
    model$residuals         <- XtR
    model$fitted_without_l  <- XtXr_without_l
    model$residual_variance <- model$sigma2
    

    return(model)

  } else {
    # Remove lth effect from fitted values
    XtXr_without_l <- model$XtXr - compute_Rv(data, model$alpha[l, ] * model$mu[l, ])

    # Compute Residuals
    XtR <- data$Xty - XtXr_without_l

    # Store residuals and parameters (standard case)
    model$residuals         <- XtR
    model$fitted_without_l  <- XtXr_without_l # For fitted update
    model$residual_variance <- model$sigma2  # Standard residual variance

    # Compute r'r for Servin-Stephens (NIG) prior
    if (params$use_servin_stephens) {
      b_minus_l <- colSums(model$alpha * model$mu) - model$alpha[l, ] * model$mu[l, ]
      model$yy_residual <- as.numeric(
        data$yty - 2 * sum(b_minus_l * data$Xty) + sum(b_minus_l * XtXr_without_l))
      model$yy_residual <- max(model$yy_residual, .Machine$double.eps)
    }

    return(model)
  }
}

# Compute SER statistics
#' @keywords internal
compute_ser_statistics.ss <- function(data, params, model, l, ...) {
  betahat <- (1 / model$predictor_weights) * model$residuals
  shat2   <- model$residual_variance / model$predictor_weights

  # Inflate shat2 for stochastic LD correction (tau_j^2 / sigma_{j,0}^2)
  if (!is.null(data$shat2_inflation))
    shat2 <- shat2 * data$shat2_inflation

  # Optimization parameters
  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: optimize on linear scale
    optim_init   <- model$V[l]
    optim_bounds <- c(0, 1)
    optim_scale  <- "linear"
  } else {
    # Standard SuSiE and SuSiE-ash: optimize on log scale
    optim_init   <- log(max(c(betahat^2 - shat2, 1), na.rm = TRUE))
    optim_bounds <- c(-30, 15)
    optim_scale  <- "log"
  }

  return(list(
    betahat      = betahat,
    shat2        = shat2,
    optim_init   = optim_init,
    optim_bounds = optim_bounds,
    optim_scale  = optim_scale
  ))
}

# Posterior expected log-likelihood for a single effect regression
#' @keywords internal
SER_posterior_e_loglik.ss <- function(data, params, model, l) {
  Eb  <- model$alpha[l, ] * model$mu[l, ]
  Eb2 <- model$alpha[l, ] * model$mu2[l, ]

  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: Omega-weighted likelihood
    return(-0.5 * (-2 * sum(Eb * model$residuals) + sum(model$predictor_weights * as.vector(Eb2))))
  } else {
    # Standard SuSiE and SuSiE-ash
    return(-0.5 / model$residual_variance * (-2 * sum(Eb * model$residuals) + sum(model$predictor_weights * as.vector(Eb2))))
  }
}

# Calculate posterior moments for single effect regression
#' @keywords internal
calculate_posterior_moments.ss <- function(data, params, model, V, l, ...) {
  if (params$use_servin_stephens) {
    # Servin-Stephens (NIG) posterior moments
    if (V <= 0) {
      post_mean  <- rep(0, data$p)
      post_mean2 <- rep(0, data$p)
      model$rv[l] <- 1
    } else {
      nig_ss <- get_nig_sufficient_stats(data, model)
      moments <- compute_posterior_moments_NIG(data$n, model$predictor_weights,
                                               model$residuals, nig_ss$yy, nig_ss$sxy,
                                               V, params$alpha0, params$beta0, nig_ss$tau)
      post_mean  <- moments$post_mean
      post_mean2 <- moments$post_mean2
      model$rv[l] <- sum(model$alpha[l, ] * moments$rv)
    }
  } else {
    # Standard Gaussian posterior calculations
    shat2 <- model$residual_variance / model$predictor_weights
    if (!is.null(data$shat2_inflation))
      shat2 <- shat2 * data$shat2_inflation

    post_var   <- V * shat2 / (V + shat2)
    post_mean  <- V * (model$residuals / model$predictor_weights) / (V + shat2)
    post_mean2 <- post_var + post_mean^2
  }

  # Store posterior moments in model
  model$mu[l, ] <- post_mean
  model$mu2[l, ] <- post_mean2

  return(model)
}

# Calculate KL divergence
#' @keywords internal
compute_kl.ss <- function(data, params, model, l) {
  if (params$use_servin_stephens) {
    if (params$L == 1) {
      model$KL[l] <- compute_kl_NIG(model$alpha[l, ], model$mu[l, ], model$mu2[l, ],
                                     model$pi, model$V[l],
                                     a0 = 1, b0 = 1,
                                     a_post = data$n / 2 + 1,
                                     b_post = data$n * model$sigma2 + 1)
    } else {
      model$KL[l] <- 0
    }
  } else {
    model <- compute_kl.default(data, params, model, l)
  }
  return(model)
}

# Expected Squared Residuals
#' @keywords internal
get_ER2.ss <- function(data, model) {
  B       <- model$alpha * model$mu
  XB2     <- sum(compute_BR(data, B) * B)
  betabar <- colSums(B)
  postb2  <- model$alpha * model$mu2 # Posterior second moment.

  return(data$yty - 2 * sum(betabar * data$Xty) + sum(betabar * compute_Rv(data, betabar)) -
           XB2 + sum(model$predictor_weights * t(postb2)))
}

# Expected log-likelihood
#' @keywords internal
Eloglik.ss <- function(data, model) {
  # Standard log-likelihood computation
  return(-data$n / 2 * log(2 * pi * model$sigma2) -
           1 / (2 * model$sigma2) * get_ER2(data, model))
}

#' @importFrom Matrix colSums
#' @importFrom stats dnorm
#' @keywords internal
loglik.ss <- function(data, params, model, V, ser_stats, l = NULL, ...) {
  if (params$use_servin_stephens) {
    # Servin-Stephens (NIG) log Bayes factors
    nig_ss <- get_nig_sufficient_stats(data, model)
    lbf <- compute_lbf_NIG(data$n, model$predictor_weights,
                            model$residuals, nig_ss$yy, nig_ss$sxy,
                            V, params$alpha0, params$beta0, nig_ss$tau)
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
    model$alpha[l, ] <- weights_res$alpha
    model$lbf[l] <- weights_res$lbf_model
    model$lbf_variable[l, ] <- stable_res$lbf

    # Compute and store marginal log-likelihood for NIG prior
    if (params$use_servin_stephens) {
      model$marginal_loglik[l] <- compute_marginal_loglik(weights_res$lbf_model, data$n,
                                                           nig_ss$yy, params$alpha0, params$beta0,
                                                           TRUE)
    }
    return(model)
  } else {
    return(weights_res$lbf_model)
  }
}

#' @keywords internal
neg_loglik.ss <- function(data, params, model, V_param, ser_stats, ...) {
  # Convert parameter to V based on optimization scale
  V <- if (ser_stats$optim_scale == "log") exp(V_param) else V_param

  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: Omega-weighted objective with logSumExp trick
    # Apply stochastic LD inflation: effective pw = pw / inflation
    pw   <- model$predictor_weights
    infl <- if (!is.null(data$shat2_inflation)) data$shat2_inflation else 1
    return(-matrixStats::logSumExp(
      -0.5 * log(1 + V * pw / infl) +
        V * model$residuals^2 / (2 * infl * (1 + V * pw / infl)) +
        log(model$pi + sqrt(.Machine$double.eps))
    ))
  } else {
    # Standard SuSiE and SuSiE-ash: standard objective
    lbf_model <- loglik.ss(data, params, model, V, ser_stats)
    return(-lbf_model)
  }
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
update_fitted_values.ss <- function(data, params, model, l, ...) {
  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: include theta in fitted values
    model$XtXr <- as.vector(compute_Rv(data, colSums(model$alpha * model$mu) + model$theta))
  } else {
    # Standard SuSiE and SuSiE-ash: sparse component only
    model$XtXr <- model$fitted_without_l + as.vector(compute_Rv(data, model$alpha[l, ] * model$mu[l, ]))
  }
  return(model)
}

# Update variance components for ss data
#' @keywords internal
update_variance_components.ss <- function(data, params, model, ...) {
  if (params$unmappable_effects == "inf") {
    # Calculate omega
    L         <- nrow(model$alpha)
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    omega     <- matrix(rep(omega_res$diagXtOmegaX, L), nrow = L, ncol = data$p, byrow = TRUE) +
      matrix(rep(1 / model$V, data$p), nrow = L, ncol = data$p, byrow = FALSE)

    # Compute theta for infinitesimal effects.
    theta <- compute_theta_blup(data, model)

    # Sigma2 and tau2 update
    if (params$estimate_residual_method == "MLE") {
      mle_result <- mle_unmappable(data, params, model, omega)
      return(list(sigma2 = mle_result$sigma2,
                  tau2   = mle_result$tau2,
                  theta  = theta))
    } else {
      mom_result <- mom_unmappable(data, params, model, omega, model$tau2)
      return(list(sigma2 = mom_result$sigma2,
                  tau2   = mom_result$tau2,
                  theta  = theta))
    }
  } else if (params$unmappable_effects == "ash") {
    # SuSiE-ash: shared update dispatches to mr.ash.rss for SS data
    return(update_ash_variance_components(data, model, params))
  } else {
    # Use default method for standard SuSiE
    return(update_variance_components.default(data, params, model))
  }
}

# Update derived quantities for ss data
#' @keywords internal
update_derived_quantities.ss <- function(data, params, model) {
  if (params$unmappable_effects == "inf") {
    # Update omega quantities for next iteration
    omega_res               <- compute_omega_quantities(data, model$tau2, model$sigma2)
    model$omega_var         <- omega_res$omega_var
    model$predictor_weights <- omega_res$diagXtOmegaX
    model$XtOmegay          <- data$eigen_vectors %*% (data$VtXty / omega_res$omega_var)
    # Update fitted values to include theta
    b          <- colSums(model$alpha * model$mu)
    model$XtXr <- compute_Rv(data, b + model$theta)
    return(model)
  } else {
    return(update_derived_quantities.default(data, params, model))
  }
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
get_scale_factors.ss <- function(data, params) {
  pm <- if (!is.null(data$XtX)) data$XtX else data$X
  return(attr(pm, "scaled:scale"))
}

# Get intercept
#' @keywords internal
get_intercept.ss <- function(data, params, model, ...) {
  return(data$y_mean - sum(data$X_colmeans * (colSums(model$alpha * model$mu) / model$X_column_scale_factors)))
}

# Get Fitted Values
#' @keywords internal
get_fitted.ss <- function(data, params, model, ...) {
  return(get_fitted.default(data, params, model, ...))
}

# Get Credible Sets
#' @keywords internal
get_cs.ss <- function(data, params, model, ...) {
  if (is.null(params$coverage) || is.null(params$min_abs_corr)) {
    return(NULL)
  }

  if (!is.null(data$X)) {
    # Low-rank X path: data$X is B x p, columns are variables
    return(susie_get_cs(model,
                        X               = data$X,
                        coverage        = params$coverage,
                        min_abs_corr    = params$min_abs_corr,
                        n_purity        = params$n_purity))
  }

  if (any(!(diag(data$XtX) %in% c(0, 1)))) {
    Xcorr <- safe_cov2cor(data$XtX)
  } else {
    Xcorr <- data$XtX
  }

  return(susie_get_cs(model,
                      Xcorr           = Xcorr,
                      check_symmetric = FALSE,
                      coverage        = params$coverage,
                      min_abs_corr    = params$min_abs_corr,
                      n_purity        = params$n_purity))
}

# Get Variable Names
#' @keywords internal
get_variable_names.ss <- function(data, model, ...) {
  pm <- if (!is.null(data$XtX)) data$XtX else data$X
  return(assign_names(data, model, colnames(pm)))
}

# Get univariate z-score
#' @keywords internal
get_zscore.ss <- function(data, params, model, ...) {
  return(get_zscore.default(data, params, model))
}

# Clean up model object for sufficient statistics data
#' @keywords internal
cleanup_model.ss <- function(data, params, model, ...) {
  # Remove common fields
  model <- cleanup_model.default(data, params, model, ...)

  # FIXME: for non-standard fields please connect them to "runtime_xx" where xx is unmappable effect option
  
  # Remove SS-specific fields for unmappable effects
  if (!is.null(params$unmappable_effects) && params$unmappable_effects == "inf") {
    unmappable_fields <- c("omega_var", "XtOmegay")
    
    for (field in unmappable_fields) {
      if (field %in% names(model)) {
        model[[field]] <- NULL
      }
    }
  } else if (!is.null(params$unmappable_effects) && params$unmappable_effects == "ash") {
    model <- cleanup_ash_fields(model)
  }
  
  return(model)
}
