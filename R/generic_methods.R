# =============================================================================
# DATA INITIALIZATION & CONFIGURATION
#
# S3 generics dispatched on data objects setup, configuration, and preprocessing.
# These prepare data objects for model fitting and handle data-specific
# configurations like unmappable effects.
#
# Functions: configure_data, get_var_y
# =============================================================================

# Configure data object for specified method
#' @keywords internal
configure_data <- function(data, params) {
  UseMethod("configure_data")
}
#' @keywords internal
#' @export
#' @noRd
configure_data.default <- function(data, params) {
  return(data)
}

# Get variance of y
#' @keywords internal
get_var_y <- function(data, ...) {
  UseMethod("get_var_y")
}
#' @keywords internal
#' @export
#' @noRd
get_var_y.default <- function(data, ...) {
  stop("get_var_y: no method for class '", class(data)[1], "'")
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

# Initialize susie model object
#' @keywords internal
initialize_susie_model <- function(data, params, ...) {
  UseMethod("initialize_susie_model")
}
#' @keywords internal
#' @export
#' @noRd
initialize_susie_model.default <- function(data, params, ...) {
  stop("initialize_susie_model: no method for class '", class(data)[1], "'")
}

# Initialize fitted values
#' @keywords internal
initialize_fitted <- function(data, mat_init) {
  UseMethod("initialize_fitted")
}
#' @keywords internal
#' @export
#' @noRd
initialize_fitted.default <- function(data, mat_init, ...) {
  stop("initialize_fitted: no method for class '", class(data)[1], "'")
}

# Validate prior variance
#' @keywords internal
validate_prior <- function(data, params, model, ...) {
  UseMethod("validate_prior")
}
#' @keywords internal
#' @export
#' @noRd
validate_prior.default <- function(data, params, model, ...) {
  invisible(TRUE)
}

# Track core parameters of a susie fit across iterations
#' @keywords internal
track_ibss_fit <- function(data, params, model, tracking, iter, ...) {
  UseMethod("track_ibss_fit")
}
#' @keywords internal
#' @export
#' @noRd
track_ibss_fit.default <- function(data, params, model, tracking, iter, elbo, ...) {
  # Store iteration snapshot if tracking is enabled.
  if (isTRUE(params$track_fit)) {
    tracking[[iter]] <- make_track_snapshot(model, iter - 1L)
  }
  return(tracking)
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

#' Get the slot weight for effect l
#'
#' Returns the weight by which effect l's contribution to the fitted
#' values is scaled. When \code{model$slot_weights} is NULL (the default),
#' all effects have weight 1 and standard SuSiE behavior is recovered.
#'
#' Slot weights enable a natural mechanism for adaptively estimating the
#' number of effects: each slot l can have a weight in [0,1] reflecting
#' the posterior probability that the slot is active. With a suitable
#' prior on the number of active effects, this generalizes SuSiE's fixed
#' L to a data-driven estimate.
#'
#' @param model SuSiE model object.
#' @param l Effect index.
#'
#' @return Scalar weight (default 1).
#'
#' @keywords internal
get_slot_weight <- function(model, l) {
  if (is.null(model$slot_weights)) 1 else model$slot_weights[l]
}

# Compute residuals for single effect regression
#' @keywords internal
compute_residuals <- function(data, params, model, l, ...) {
  UseMethod("compute_residuals")
}
#' @keywords internal
#' @export
#' @noRd
compute_residuals.default <- function(data, params, model, l, ...) {
  stop("compute_residuals: no method for class '", class(data)[1], "'")
}

# Compute SER statistics (betahat, shat2)
#' @keywords internal
compute_ser_statistics <- function(data, params, model, l, ...) {
  UseMethod("compute_ser_statistics")
}
#' @keywords internal
#' @export
#' @noRd
compute_ser_statistics.default <- function(data, params, model, l, ...) {
  stop("compute_ser_statistics: no method for class '", class(data)[1], "'")
}

# Single effect regression posterior expected log-likelihood
#' @keywords internal
SER_posterior_e_loglik <- function(data, params, model, l) {
  UseMethod("SER_posterior_e_loglik")
}
#' @keywords internal
#' @export
#' @noRd
SER_posterior_e_loglik.default <- function(data, params, model, l) {
  stop("SER_posterior_e_loglik: no method for class '", class(data)[1], "'")
}

# Calculate posterior moments for single effect regression
#' @keywords internal
calculate_posterior_moments <- function(data, params, model, V, l, ...) {
  UseMethod("calculate_posterior_moments")
}
#' @keywords internal
#' @export
#' @noRd
calculate_posterior_moments.default <- function(data, params, model, V, l = NULL, ...) {
  stop("calculate_posterior_moments: no method for class '", class(data)[1], "'")
}

# Calculate KL divergence
#' @keywords internal
compute_kl <- function(data, params, model, l) {
  UseMethod("compute_kl")
}
#' @keywords internal
#' @export
#' @noRd
compute_kl.default <- function(data, params, model, l) {
  model$KL[l] <- -model$lbf[l] + SER_posterior_e_loglik(data, params, model, l)
  return(model)
}

# Expected squared residuals
#' @keywords internal
get_ER2 <- function(data, model) {
  UseMethod("get_ER2")
}
#' @keywords internal
#' @export
#' @noRd
get_ER2.default <- function(data, model) {
  stop("get_ER2: no method for class '", class(data)[1], "'")
}

# Expected log-likelihood
#' @keywords internal
Eloglik <- function(data, model) {
  UseMethod("Eloglik")
}
#' @keywords internal
#' @export
#' @noRd
Eloglik.default <- function(data, model) {
  stop("Eloglik: no method for class '", class(data)[1], "'")
}

# Variational E_q[log p(y|b, sigma^2)] under SuSiE-NIG. Non-S3 helper called
# from get_objective so we don't break the Eloglik(data, model) signature
# that downstream packages override (mvsusieR, mfsusieR).
# Decomposition: E_q[||y-Xb||^2 | sigma^2] = A + sigma^2 B,
#   B = sum_l sum_j alpha^(l)_j * r0^(l)_j * tau_j,
#   A = get_ER2 - E[sigma^2] * B.
# Eloglik = -n/2 log(2 pi) - n/2 (log b - digamma(a)) - 0.5 (A * a/b + B).
#' @keywords internal
nig_eloglik <- function(data, params, model) {
  n         <- data$n
  ERSS_marg <- get_ER2(data, model)
  a_post    <- (params$alpha0 + n) / 2
  b_post    <- (params$beta0 + ERSS_marg) / 2

  tau_v <- if (!is.null(model$shat2_inflation)) model$shat2_inflation else 1
  pw    <- model$predictor_weights
  B     <- 0
  for (l in seq_len(nrow(model$alpha))) {
    r0_l <- model$V[l] / (model$V[l] + tau_v / pw)
    B    <- B + sum(model$alpha[l, ] * r0_l * tau_v)
  }
  A <- ERSS_marg - (b_post / (a_post - 1)) * B
  -n / 2 * log(2 * pi) - n / 2 * (log(b_post) - digamma(a_post)) -
    0.5 * (A * a_post / b_post + B)
}

# Log-likelihood and posterior moments for fixed mixture prior
# (estimate_prior_method = "fixed_mixture"). Evaluates BFs on a
# pre-specified variance grid with given mixture weights.
#' @keywords internal
loglik_mixture <- function(data, params, model, ser_stats, l, ...) {
  UseMethod("loglik_mixture")
}
#' @keywords internal
#' @export
#' @noRd
loglik_mixture.default <- function(data, params, model, ser_stats, l, ...) {
  # Shared implementation for all data types.
  # compute_ser_statistics() (type-specific) has already produced betahat and shat2.
  model <- loglik_mixture_common(params, model, ser_stats, l)
  return(model)
}

#' @keywords internal
calculate_posterior_moments_mixture <- function(data, params, model, l, ...) {
  UseMethod("calculate_posterior_moments_mixture")
}
#' @keywords internal
#' @export
#' @noRd
calculate_posterior_moments_mixture.default <- function(data, params, model, l, ...) {
  # Shared implementation: mixture posterior from stored lbf_grid and ser_stats
  model <- calculate_posterior_moments_mixture_common(params, model, l)
  return(model)
}

# Log-likelihood for prior variance optimization
#' @keywords internal
loglik <- function(data, params, model, V, ser_stats, l = NULL, ...) {
  UseMethod("loglik")
}
#' @keywords internal
#' @export
#' @noRd
loglik.default <- function(data, params, model, V, ser_stats, l = NULL, ...) {
  stop("loglik: no method for class '", class(data)[1], "'")
}

# Negative log-likelihood for optimization (handles both log and linear scales)
#' @keywords internal
neg_loglik <- function(data, params, model, V_param, ser_stats, ...) {
  UseMethod("neg_loglik")
}
#' @keywords internal
#' @export
#' @noRd
neg_loglik.default <- function(data, params, model, V_param, ser_stats, ...) {
  stop("neg_loglik: no method for class '", class(data)[1], "'")
}

# EM update for prior variance
#' @keywords internal
em_update_prior_variance <- function(data, params, model, alpha, moments, V_init) {
  UseMethod("em_update_prior_variance")
}
#' @keywords internal
#' @export
#' @noRd
em_update_prior_variance.default <- function(data, params, model, alpha, moments, V_init) {
  if (!is.null(params$use_NIG) && params$use_NIG) {
    nig_ss <- get_nig_sufficient_stats(data, model)
    return(update_prior_variance_NIG_EM(data$n, model$predictor_weights,
                                         model$residuals, nig_ss$yy, nig_ss$sxy,
                                         alpha, V_init, params$alpha0, params$beta0,
                                         nig_ss$tau))
  }
  # Standard EM update
  sum(alpha * moments$post_mean2)
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
update_fitted_values <- function(data, params, model, l, ...) {
  UseMethod("update_fitted_values")
}
#' @keywords internal
#' @export
#' @noRd
update_fitted_values.default <- function(data, params, model, l, ...) {
  stop("update_fitted_values: no method for class '", class(data)[1], "'")
}

# Update variance components
#' @keywords internal
update_variance_components <- function(data, params, model, ...) {
  UseMethod("update_variance_components")
}
#' @keywords internal
#' @export
#' @noRd
update_variance_components.default <- function(data, params, model, ...) {
  if (isTRUE(params$use_NIG)) {
    # Posterior mean of IG((alpha0+n)/2, (beta0+ERSS)/2)
    sigma2 <- (params$beta0 + get_ER2(data, model)) /
              (params$alpha0 + data$n - 2)
  } else {
    sigma2 <- est_residual_variance(data, model)
  }
  return(list(sigma2 = sigma2))
}

# Update derived quantities after variance component changes
#' @keywords internal
update_derived_quantities <- function(data, params, model) {
  UseMethod("update_derived_quantities")
}
#' @keywords internal
#' @export
#' @noRd
update_derived_quantities.default <- function(data, params, model) {
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
# get_variable_names, get_zscore, cleanup_model
# =============================================================================

# Get column scale factors
#' @keywords internal
get_scale_factors <- function(data, params, ...) {
  UseMethod("get_scale_factors")
}
#' @keywords internal
#' @export
#' @noRd
get_scale_factors.default <- function(data, params, ...) {
  stop("get_scale_factors: no method for class '", class(data)[1], "'")
}

# Get intercept
#' @keywords internal
get_intercept <- function(data, params, model, ...) {
  UseMethod("get_intercept")
}
#' @keywords internal
#' @export
#' @noRd
get_intercept.default <- function(data, params, model, ...) {
  stop("get_intercept: no method for class '", class(data)[1], "'")
}

# Get fitted values
#' @keywords internal
get_fitted <- function(data, params, model, ...) {
  UseMethod("get_fitted")
}
#' @keywords internal
#' @export
#' @noRd
get_fitted.default <- function(data, params, model, ...) {
  return(NULL)
}

# Get credible sets
#' @keywords internal
get_cs <- function(data, params, model, ...) {
  UseMethod("get_cs")
}
#' @keywords internal
#' @export
#' @noRd
get_cs.default <- function(data, params, model, ...) {
  stop("get_cs: no method for class '", class(data)[1], "'")
}

# Get variable names
#' @keywords internal
get_variable_names <- function(data, model, ...) {
  UseMethod("get_variable_names")
}
#' @keywords internal
#' @export
#' @noRd
get_variable_names.default <- function(data, model, ...) {
  stop("get_variable_names: no method for class '", class(data)[1], "'")
}

# Get univariate z-scores
#' @keywords internal
get_zscore <- function(data, params, model, ...) {
  UseMethod("get_zscore")
}
#' @keywords internal
#' @export
#' @noRd
get_zscore.default <- function(data, params, model, ...) {
  return(NULL)
}

# Clean up model object by removing temporary computational fields
#' @keywords internal
cleanup_model <- function(data, params, model, ...) {
  UseMethod("cleanup_model")
}

#' Class-specific extra fields to strip in cleanup_model.default
#'
#' Default returns `character(0)`. Subclasses (e.g., mfsusieR's
#' `raw_residuals`, mvsusieR's `Y_imputed`/`llik_cache`) override
#' to add their per-class scratch fields. Result is unioned with
#' the standard temp_fields list inside `cleanup_model.default`.
#' @keywords internal
cleanup_extra_fields <- function(data) {
  UseMethod("cleanup_extra_fields")
}
#' @keywords internal
#' @export
#' @noRd
cleanup_extra_fields.default <- function(data) {
  character(0)
}

#' @keywords internal
#' @export
#' @noRd
cleanup_model.default <- function(data, params, model, ...) {
  # Remove temporary fields common to all data types
  temp_fields <- c("null_weight", "predictor_weights", "runtime",
                   "prev_elbo", "prev_alpha",
                   "residuals", "fitted_without_l", "residual_variance",
                   "shat2_inflation", "R_bf_attenuation",
                   "R_mismatch_ser_model",
                   cleanup_extra_fields(data))

  for (field in temp_fields) {
    if (field %in% names(model)) {
      model[[field]] <- NULL
    }
  }

  return(model)
}
