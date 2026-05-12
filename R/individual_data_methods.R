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
  if (params$unmappable_effects == "inf") {
    # Stay in the individual-data path: thin SVD of standardized X, no XtX.
    return(add_eigen_decomposition(data, params))
  }
  return(configure_data.default(data, params))
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

  # SuSiE-inf: predictor_weights are the diag of X'OmegaX (eigenspace),
  # not attr(X, "d").  Also cache omega_var and XtOmegay for the first
  # SER pass.
  if (params$unmappable_effects == "inf") {
    omega_res               <- compute_omega_quantities(data, tau2 = 0, sigma2 = var_y)
    model$omega_var         <- omega_res$omega_var
    model$predictor_weights <- omega_res$diagXtOmegaX
    model$XtOmegay          <- as.vector(data$eigen_vectors %*%
                                           (data$VtXty / omega_res$omega_var))

    model$tau2  <- 0
    model$theta <- rep(0, data$p)
    return(model)
  }

  # Append predictor weights
  model$predictor_weights <- attr(data$X, "d")

  # Initialize NIG parameters
  if (params$use_NIG) {
    model$rv <- rep(1, params$L)
    model$marginal_loglik <- rep(as.numeric(NA), params$L)
  }

  # Initialize ash (Mr.ASH) tracking fields
  if (params$unmappable_effects == "ash") {
    model <- init_ash_fields(model, data$n, data$p, params$L, is_individual = TRUE)
  } else if (params$unmappable_effects == "ash_filter_archived") {
    model <- init_ash_fields_filter_archived(model, data$n, data$p, params$L, is_individual = TRUE)
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
  if (params$unmappable_effects %in% c("ash", "ash_filter_archived")) {
    return(track_ibss_fit.default(data, params, model, tracking, iter, elbo, ...))
  }
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
  sw_l <- get_slot_weight(model, l)

  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: Omega-weighted residuals computed in eigenspace.
    # The full b_minus_l is built from (alpha, mu) directly here since
    # fitted_without_l (n-space) is not maintained on the inf path.
    # model$XtOmegay, model$omega_var, model$predictor_weights are cached
    # at iter boundaries.
    sw <- if (!is.null(model$slot_weights)) model$slot_weights else rep(1, nrow(model$alpha))
    b_minus_l <- colSums(sw * model$alpha * model$mu) - sw_l * model$alpha[l, ] * model$mu[l, ]

    # Compute V' b_minus_l once; reused for XtOmegaXb (Omega-weighted) and
    # the inflation-tail XtXr_without_l (un-weighted).  Saves one O(pr)
    # crossprod per SER step.
    Vtb_minus_l <- as.vector(crossprod(data$eigen_vectors, b_minus_l))

    XtOmegaXb <- as.vector(data$eigen_vectors %*%
                             (Vtb_minus_l * data$eigen_values / model$omega_var))
    model$residuals         <- model$XtOmegay - XtOmegaXb
    model$residual_variance <- 1   # Already incorporated in Omega

    # R inflation uses standard (non-Omega) quantities.
    XtXr_without_l <- as.vector(data$eigen_vectors %*%
                                  (data$eigen_values * Vtb_minus_l))
    r_vec <- data$Xty - XtXr_without_l
    infl_state <- compute_shat2_inflation(data, model, XtXr_without_l,
                                          b_minus_l, r_vec)
    model <- apply_inflation_state(model, infl_state)
    return(model)
  }

  # Remove lth effect from fitted values (scaled by slot weight)
  Xr_without_l <- model$Xr - sw_l * compute_Xb(data$X, model$alpha[l, ] * model$mu[l, ])

  # Compute residuals
  if (params$unmappable_effects %in% c("ash", "ash_filter_archived")) {
    # Subtract both sparse effects (without l) and ash theta
    R <- data$y - Xr_without_l - model$X_theta
  } else {
    R <- data$y - Xr_without_l
  }
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

  # Inflate shat2 for finite-reference R variance tracking (set by the
  # SuSiE-inf branch in compute_residuals.individual via shat2_inflation).
  if (!is.null(model$shat2_inflation))
    shat2 <- shat2 * model$shat2_inflation

  # Optimization parameters
  if (params$unmappable_effects == "inf") {
    optim_init   <- model$V[l]
    optim_bounds <- c(0, 1)
    optim_scale  <- "linear"
  } else {
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

# Posterior expected log-likelihood for single effect regression
#' @keywords internal
SER_posterior_e_loglik.individual <- function(data, params, model, l) {
  Eb  <- model$alpha[l, ] * model$mu[l, ]
  Eb2 <- model$alpha[l, ] * model$mu2[l, ]

  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: Omega-weighted likelihood (matches SER_posterior_e_loglik.ss).
    # raw_residuals is not maintained on the inf path; the marginal y-density
    # constant is absorbed by the standard compute_kl.default rather than
    # being canceled here.
    return(-0.5 * (-2 * sum(Eb * model$residuals) +
                     sum(model$predictor_weights * as.vector(Eb2))))
  }

  return(-0.5 * data$n * log(2 * pi * model$sigma2) -
           0.5 / model$sigma2 * (sum(model$raw_residuals * model$raw_residuals)
                                 - 2 * sum(model$raw_residuals * compute_Xb(data$X, Eb)) +
                                   sum(model$predictor_weights * Eb2)))
}

# Calculate posterior moments for single effect regression
#' @keywords internal
calculate_posterior_moments.individual <- function(data, params, model, V, l, ...) {
  if (params$use_NIG) {
    if (V <= 0) {
      # Zero variance case
      post_mean  <- rep(0, data$p)
      post_mean2 <- rep(0, data$p)
      model$rv[l] <- 1
    } else {
      # Compute posterior moments for NIG prior
      nig_ss <- get_nig_sufficient_stats(data, model)
      moments <- compute_posterior_moments_NIG(data$n, model$predictor_weights,
                                               model$residuals, nig_ss$yy, nig_ss$sxy,
                                               V, params$alpha0, params$beta0, nig_ss$tau)

      post_mean  <- moments$post_mean
      post_mean2 <- moments$post_mean2

      # Compute weighted average of residual variance modes using PIPs
      model$rv[l] <- sum(model$alpha[l, ] * moments$rv)
    }
    moments <- list(post_mean = post_mean, post_mean2 = post_mean2)
  } else {
    # Standard Gaussian posterior calculations
    shat2 <- model$residual_variance / model$predictor_weights
    if (!is.null(model$shat2_inflation))
      shat2 <- shat2 * model$shat2_inflation

    betahat <- model$residuals / model$predictor_weights
    moments <- gaussian_ser_moments(betahat, shat2, V)
  }

  # Store posterior moments in model
  store_ser_moments(model, l, moments)
}

# Calculate KL divergence
#' @keywords internal
compute_kl.individual <- function(data, params, model, l) {
  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: standard form KL = -lbf + SER_posterior_e_loglik (no
    # marginal-y constant to cancel since SER_posterior on the inf path
    # already omits it).
    return(compute_kl.default(data, params, model, l))
  }

  if (params$use_NIG) {
    # NIG KL only valid for L=1 (gIBSS for L>1 has no coherent ELBO; supp. line 503)
    if (params$L == 1) {
      ki <- nig_kl_inputs(data, params, model, l)
      kl <- compute_kl_NIG(model$alpha[l, ], model$mu[l, ], model$mu2[l, ],
                           model$pi, model$V[l],
                           a0 = params$alpha0 / 2, b0 = params$beta0 / 2,
                           a_post = ki$a_post, b_post = ki$b_post,
                           s_j_sq = ki$s_j_sq)
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
  # For ash, subtract theta contribution from residuals
  y_adj <- if (!is.null(model$X_theta)) data$y - model$X_theta else data$y
  # Slot-weight correction: E[||y - sum_l c_l X beta^(l)||^2] under Bern(chat_l)
  # = ||y - Xr||^2 + sum_l chat_l * E[||X b^(l)||^2] - chat_l^2 * ||X bbar_l||^2
  # When slot_weights is NULL (all weights = 1), reduces to the standard formula.
  sw <- if (!is.null(model$slot_weights)) model$slot_weights else rep(1, nrow(model$alpha))
  per_slot_Eb2 <- as.vector(postb2 %*% model$predictor_weights)  # L-vector
  per_slot_Xb2 <- rowSums(Xr_L^2)                                # L-vector
  return(sum((y_adj - model$Xr)^2) + sum(sw * per_slot_Eb2 - sw^2 * per_slot_Xb2))
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
  # Check if using NIG prior
  if (params$use_NIG) {
    # Compute log Bayes factors for NIG prior
    nig_ss <- get_nig_sufficient_stats(data, model)
    lbf <- compute_lbf_NIG(data$n, model$predictor_weights,
                           model$residuals, nig_ss$yy, nig_ss$sxy,
                           V, params$alpha0, params$beta0, nig_ss$tau)
  } else {
    # Standard Gaussian prior log Bayes factors
    lbf <- gaussian_ser_lbf(ser_stats$betahat, ser_stats$shat2, V)
  }

  ser_res <- apply_ser_lbf(model, lbf, ser_stats$shat2, l)

  # Store in model if l is provided, otherwise return lbf_model for prior variance optimization
  if (!is.null(l)) {
    model <- ser_res$model

    # Compute and store marginal log-likelihood for NIG prior
    if (params$use_NIG) {
      model$marginal_loglik[l] <- compute_marginal_loglik(ser_res$lbf_model, data$n,
                                                           nig_ss$yy, params$alpha0, params$beta0,
                                                           TRUE)
    }
    return(model)
  } else {
    return(ser_res$lbf_model)
  }
}

#' @keywords internal
neg_loglik.individual <- function(data, params, model, V_param, ser_stats, ...) {
  # Convert parameter to V based on optimization scale.  SuSiE-inf optimizes
  # on the linear scale; the rest (and the previous individual default) use
  # the log scale.
  V <- if (ser_stats$optim_scale == "log") exp(V_param) else V_param

  if (params$unmappable_effects == "inf") {
    pw   <- model$predictor_weights
    infl <- if (!is.null(model$shat2_inflation)) model$shat2_inflation else 1
    return(-matrixStats::logSumExp(
      -0.5 * log(1 + V * pw / infl) +
        V * model$residuals^2 / (2 * infl * (1 + V * pw / infl)) +
        log(model$pi + sqrt(.Machine$double.eps))
    ))
  }

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
  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: include theta in fitted values; recompute from scratch
    # because fitted_without_l is not maintained on the inf path.
    sw <- if (!is.null(model$slot_weights)) model$slot_weights else rep(1, nrow(model$alpha))
    b_total <- colSums(sw * model$alpha * model$mu) + model$theta
    model$Xr <- as.vector(compute_Xb(data$X, b_total))
    return(model)
  }

  sw_l <- get_slot_weight(model, l)
  model$Xr <- model$fitted_without_l + sw_l * compute_Xb(data$X, model$alpha[l, ] * model$mu[l, ])

  return(model)
}

# Update variance components for individual data
#' @keywords internal
update_variance_components.individual <- function(data, params, model, ...) {
  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: identical math to update_variance_components.ss; all eigenspace.
    # model$predictor_weights == diagXtOmegaX is cached at iter boundaries.
    L         <- nrow(model$alpha)
    omega     <- matrix(rep(model$predictor_weights, L), nrow = L, ncol = data$p, byrow = TRUE) +
      matrix(rep(1 / model$V, data$p), nrow = L, ncol = data$p, byrow = FALSE)

    theta <- compute_theta_blup(data, model)

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
  } else if (params$unmappable_effects == "ash_filter_archived") {
    # Original filter-based masking (archived for internal diagnostics)
    return(update_ash_variance_components_filter_archived(data, model, params))
  } else if (params$unmappable_effects == "ash") {
    # c_hat + 3 LD-interference heuristics
    return(update_ash_variance_components(data, model, params))
  }
  return(update_variance_components.default(data, params, model, ...))
}

# Update derived quantities for individual data
#' @keywords internal
update_derived_quantities.individual <- function(data, params, model) {
  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: refresh omega caches with new (tau2, sigma2) and update
    # the fitted vector to include theta.  Mirrors update_derived_quantities.ss.
    omega_res               <- compute_omega_quantities(data, model$tau2, model$sigma2)
    model$omega_var         <- omega_res$omega_var
    model$predictor_weights <- omega_res$diagXtOmegaX
    model$XtOmegay          <- as.vector(data$eigen_vectors %*%
                                           (data$VtXty / omega_res$omega_var))

    b <- colSums(model$alpha * model$mu)
    model$Xr <- as.vector(compute_Xb(data$X, b + model$theta))
    return(model)
  }
  if (params$unmappable_effects %in% c("ash", "ash_filter_archived")) {
    # For ash, recompute full Xr including sparse effects only
    # (theta is tracked separately via X_theta)
    # Use slot_weights (c_hat) if available to maintain consistency
    # with the c_hat-weighted Xr from ibss_fit + update_c_hat.
    sw <- if (!is.null(model$slot_weights)) model$slot_weights else rep(1, nrow(model$alpha))
    b <- colSums(sw * model$alpha * model$mu)
    model$Xr <- as.vector(compute_Xb(data$X, b))
    return(model)
  }
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

  # Include ash theta contribution
  if (!is.null(model$X_theta)) {
    fitted <- fitted + model$X_theta
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
      "to skip this step set compute_univariate_zscore = FALSE.",
      style = "hint"
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

  # Remove NIG specific temporary fields
  if (params$use_NIG) {
    model$marginal_loglik <- NULL
  }

  # Remove ash-specific runtime fields
  if (!is.null(params$unmappable_effects) && params$unmappable_effects == "ash") {
    model <- cleanup_ash_fields(model)
  } else if (!is.null(params$unmappable_effects) && params$unmappable_effects == "ash_filter_archived") {
    model <- cleanup_ash_fields_filter_archived(model)
  }

  return(model)
}
