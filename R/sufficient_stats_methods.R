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
  } else if (params$unmappable_effects == "ash_filter_archived") {
    pm <- if (!is.null(data$XtX)) data$XtX else data$X
    model$predictor_weights <- attr(pm, "d")
    model <- init_ash_fields_filter_archived(model, data$n, data$p, params$L, is_individual = FALSE)
  } else {
    pm <- if (!is.null(data$XtX)) data$XtX else data$X
    model$predictor_weights <- attr(pm, "d")

    # Initialize NIG parameters
    if (params$use_NIG) {
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
        "This usually caused by mismatch between the summary statistics and the R matrix.\n",
        "Please check the input."
      )
    }
  }
  return(validate_prior.default(data, params, model, ...))
}

# Track core parameters across iterations
#' @keywords internal
track_ibss_fit.ss <- function(data, params, model, tracking, iter, elbo, ...) {
  if (params$unmappable_effects %in% c("inf", "ash", "ash_filter_archived")) {
    return(track_ibss_fit.default(data, params, model, tracking, iter, elbo, ...))
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
  # Weighted sum of effects excluding l (slot_weights scale each effect's contribution)
  sw_l <- get_slot_weight(model, l)
  sw <- if (!is.null(model$slot_weights)) model$slot_weights else rep(1, nrow(model$alpha))
  b_minus_l <- colSums(sw * model$alpha * model$mu) - sw_l * model$alpha[l, ] * model$mu[l, ]

  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: Omega-weighted residuals.  model$XtOmegay, model$omega_var,
    # and model$predictor_weights (= diagXtOmegaX) are cached at iter
    # boundaries; no recompute here.
    XtOmegaXb <- as.vector(data$eigen_vectors %*%
                             ((crossprod(data$eigen_vectors, b_minus_l)) *
                                data$eigen_values / model$omega_var))

    model$residuals         <- model$XtOmegay - XtOmegaXb
    model$residual_variance <- 1   # Already incorporated in Omega

    # R inflation uses standard (non-Omega) quantities
    XtXr_without_l <- compute_Rv(data, b_minus_l)
    r <- data$Xty - XtXr_without_l
    infl_state <- compute_shat2_inflation(data, model, XtXr_without_l,
                                          b_minus_l, r)
    model <- apply_inflation_state(model, infl_state)
    return(model)
  }

  # Below are SuSiE, SuSiE-ASH and SuSiE-SS

  # Remove lth effect from fitted values (scaled by slot weight)
  XtXr_without_l <- model$XtXr - sw_l * compute_Rv(data, model$alpha[l, ] * model$mu[l, ])

  # Compute residuals (ash subtracts unmappable effect X'X*theta).
  is_ash <- params$unmappable_effects %in% c("ash", "ash_filter_archived")
  if (is_ash) {
    model$residuals <- data$Xty - model$XtX_theta - XtXr_without_l
  } else {
    model$residuals <- data$Xty - XtXr_without_l
  }

  model$fitted_without_l  <- XtXr_without_l
  model$residual_variance <- model$sigma2

  # NIG prior: compute residual sum of squares
  if (params$use_NIG) {
    model$yy_residual <- as.numeric(
      data$yty - 2 * sum(b_minus_l * data$Xty) + sum(b_minus_l * XtXr_without_l))
    model$yy_residual <- max(model$yy_residual, .Machine$double.eps)
  }

  # ASH path: residual subtracts theta (line 167), so the variance scale
  # s = eta^2 + v_g must also be built from b_minus_l + theta or the
  # data and variance model disagree on what has been removed.
  if (is_ash && !is.null(model$theta)) {
    XtX_theta <- if (!is.null(model$XtX_theta))
                   model$XtX_theta
                 else compute_Rv(data, model$theta)
    b_for_infl    <- b_minus_l + model$theta
    XtXr_for_infl <- XtXr_without_l + XtX_theta
  } else {
    b_for_infl    <- b_minus_l
    XtXr_for_infl <- XtXr_without_l
  }
  infl_state <- compute_shat2_inflation(data, model, XtXr_for_infl,
                                        b_for_infl, model$residuals)
  model <- apply_inflation_state(model, infl_state)

  return(model)
}

# compute_shat2_inflation moved to R/rss_mismatch.R.

# Compute SER statistics
#' @keywords internal
compute_ser_statistics.ss <- function(data, params, model, l, ...) {
  betahat <- (1 / model$predictor_weights) * model$residuals
  shat2   <- model$residual_variance / model$predictor_weights

  # Inflate shat2 for finite-reference R variance tracking (tau_j^2 / sigma^2)
  if (!is.null(model$shat2_inflation))
    shat2 <- shat2 * model$shat2_inflation

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
    return(gaussian_ser_posterior_e_loglik(
      model$alpha[l, ], model$mu[l, ], model$mu2[l, ],
      model$residuals / model$predictor_weights,
      model$residual_variance / model$predictor_weights))
  }
}

# Calculate posterior moments for single effect regression
#' @keywords internal
calculate_posterior_moments.ss <- function(data, params, model, V, l, ...) {
  if (params$use_NIG) {
    # NIG posterior moments
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
compute_kl.ss <- function(data, params, model, l) {
  if (params$use_NIG) {
    # NIG KL only valid for L=1 (gIBSS for L>1 has no coherent ELBO; supp. line 503)
    if (params$L == 1) {
      ki <- nig_kl_inputs(data, params, model, l)
      model$KL[l] <- compute_kl_NIG(model$alpha[l, ], model$mu[l, ], model$mu2[l, ],
                                     model$pi, model$V[l],
                                     a0 = params$alpha0 / 2, b0 = params$beta0 / 2,
                                     a_post = ki$a_post, b_post = ki$b_post,
                                     s_j_sq = ki$s_j_sq)
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
  postb2  <- model$alpha * model$mu2 # Posterior second moment.
  # Slot-weight correction: E[||y - sum_l c_l X beta^(l)||^2] under Bern(chat_l)
  # = y'y - 2 betabar_w' X'y + betabar_w' X'X betabar_w
  #   + sum_l chat_l * E[b^(l)' X'X b^(l)] - chat_l^2 * bbar_l' X'X bbar_l
  # When slot_weights is NULL (all weights = 1), reduces to the standard formula.
  sw <- if (!is.null(model$slot_weights)) model$slot_weights else rep(1, nrow(B))
  betabar <- colSums(sw * B)                                      # c_hat-weighted mean
  per_slot_XB2 <- rowSums(compute_BR(data, B) * B)                # bbar_l' R bbar_l
  per_slot_Eb2 <- as.vector(postb2 %*% model$predictor_weights)   # diag(X'X)' (alpha*mu2)_l

  return(data$yty - 2 * sum(betabar * data$Xty) + sum(betabar * compute_Rv(data, betabar)) -
           sum(sw^2 * per_slot_XB2) + sum(sw * per_slot_Eb2))
}

# Expected log-likelihood for the sufficient-stats path.  Without inflation,
# the standard regression log-likelihood under sigma2 (matches Eloglik.individual).
# With finite-R inflation, the SER posterior fits a betahat-scale augmented
# model; switch to the matching data-fit term.  Affects ELBO only; PIP/CS/
# sigma2 (which goes through est_residual_variance, not Eloglik) are unchanged.
#' @keywords internal
Eloglik.ss <- function(data, model) {
  if (!is.null(model$shat2_inflation))
    return(compute_augmented_eloglik_ss(data, model))
  -data$n / 2 * log(2 * pi * model$sigma2) -
    1 / (2 * model$sigma2) * get_ER2(data, model)
}

# Variational expectation of the augmented betahat-scale Gaussian log-
# likelihood under finite-R inflation.  Form derived in
# ld_mismatch_generativemodel.tex Sec. "Variational ELBO under the
# augmented variance".  The Var_q[(X'X beta)_j] correction requires
# (X'X)^2 element-wise; formed on each call.
#' @keywords internal
compute_augmented_eloglik_ss <- function(data, model) {
  pw     <- data$predictor_weights
  infl   <- model$shat2_inflation
  sigma2 <- model$sigma2
  p      <- length(pw)

  sw  <- if (!is.null(model$slot_weights)) model$slot_weights else rep(1, nrow(model$alpha))
  am  <- model$alpha * model$mu
  am2 <- model$alpha * model$mu2
  betabar  <- colSums(sw * am)

  res_mean <- data$Xty - compute_Rv(data, betabar)

  XtX <- if (!is.null(data$XtX)) data$XtX else crossprod(data$X)
  XtX_sq <- XtX * XtX
  F_mat <- am  %*% XtX
  G_mat <- am2 %*% XtX_sq
  var_corr <- as.vector(crossprod(G_mat - F_mat^2, sw^2))

  -p / 2 * log(2 * pi) -
    0.5 * sum(log(sigma2 * infl / pw)) -
    0.5 * sum((res_mean^2 + var_corr) / (pw * sigma2 * infl))
}

#' @importFrom Matrix colSums
#' @importFrom stats dnorm
#' @keywords internal
loglik.ss <- function(data, params, model, V, ser_stats, l = NULL, ...) {
  if (params$use_NIG) {
    # NIG log Bayes factors
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
    if (!params$use_NIG)
      model <- record_R_bf_attenuation(model, ser_stats, lbf, V, l)

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
neg_loglik.ss <- function(data, params, model, V_param, ser_stats, ...) {
  # Convert parameter to V based on optimization scale
  V <- if (ser_stats$optim_scale == "log") exp(V_param) else V_param

  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: Omega-weighted objective with logSumExp trick
    # Apply finite-reference R inflation: effective pw = pw / inflation
    pw   <- model$predictor_weights
    infl <- if (!is.null(model$shat2_inflation)) model$shat2_inflation else 1
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
  sw_l <- get_slot_weight(model, l)
  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: include theta in fitted values
    sw <- if (!is.null(model$slot_weights)) model$slot_weights else rep(1, nrow(model$alpha))
    model$XtXr <- as.vector(compute_Rv(data, colSums(sw * model$alpha * model$mu) + model$theta))
  } else {
    # Standard SuSiE and SuSiE-ash: sparse component only
    model$XtXr <- model$fitted_without_l + sw_l * as.vector(compute_Rv(data, model$alpha[l, ] * model$mu[l, ]))
  }
  return(model)
}

# Update variance components for ss data
#' @keywords internal
update_variance_components.ss <- function(data, params, model, ...) {
  if (params$unmappable_effects == "inf") {
    # Calculate omega.  model$predictor_weights == diagXtOmegaX is cached at
    # iter boundaries; no recompute needed.
    L         <- nrow(model$alpha)
    omega     <- matrix(rep(model$predictor_weights, L), nrow = L, ncol = data$p, byrow = TRUE) +
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
  } else if (params$unmappable_effects == "ash_filter_archived") {
    # Original filter-based masking (archived for internal diagnostics)
    return(update_ash_variance_components_filter_archived(data, model, params))
  } else if (params$unmappable_effects == "ash") {
    # c_hat + 3 LD-interference heuristics
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
  } else if (!is.null(params$unmappable_effects) && params$unmappable_effects == "ash_filter_archived") {
    model <- cleanup_ash_fields_filter_archived(model)
  }
  
  # Remove NIG specific temporary fields
  if (params$use_NIG) {
    model$marginal_loglik <- NULL
  }
  
  return(model)
}
