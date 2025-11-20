# =============================================================================
# DATA INITIALIZATION & CONFIGURATION
#
# Functions for data object setup, configuration, and preprocessing.
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
configure_data.default <- function(data, params) {
  return(data)
}

# Get variance of y
#' @keywords internal
get_var_y <- function(data, ...) {
  UseMethod("get_var_y")
}
#' @keywords internal
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
initialize_susie_model.default <- function(data, params, ...) {
  stop("initialize_susie_model: no method for class '", class(data)[1], "'")
}

# Initialize fitted values
#' @keywords internal
initialize_fitted <- function(data, mat_init) {
  UseMethod("initialize_fitted")
}
#' @keywords internal
initialize_fitted.default <- function(data, mat_init, ...) {
  stop("initialize_fitted: no method for class '", class(data)[1], "'")
}

# Validate prior variance
#' @keywords internal
validate_prior <- function(data, params, model, ...) {
  UseMethod("validate_prior")
}
#' @keywords internal
validate_prior.default <- function(data, params, model, ...) {
  invisible(TRUE)
}

# Track core parameters of a susie fit across iterations
#' @keywords internal
track_ibss_fit <- function(data, params, model, tracking, iter, ...) {
  UseMethod("track_ibss_fit")
}
#' @keywords internal
track_ibss_fit.default <- function(data, params, model, tracking, iter, elbo, ...) {
  # Initialize tracking structure on first iteration
  if (iter == 1) {
    tracking$convergence <- list(
      prev_elbo  = -Inf,
      prev_alpha = model$alpha
    )
  } else {
    # Update previous values for convergence checking
    tracking$convergence$prev_elbo  <- elbo[iter]
    tracking$convergence$prev_alpha <- model$alpha
  }

  # Store additional tracking information if requested
  if (isTRUE(params$track_fit)) {
    tracking[[iter]] <- list(
      alpha  = model$alpha,
      niter  = iter,
      V      = model$V,
      sigma2 = model$sigma2
    )
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

# Compute residuals for single effect regression
#' @keywords internal
compute_residuals <- function(data, params, model, l, ...) {
  UseMethod("compute_residuals")
}
#' @keywords internal
compute_residuals.default <- function(data, params, model, l, ...) {
  stop("compute_residuals: no method for class '", class(data)[1], "'")
}

# Compute SER statistics (betahat, shat2)
#' @keywords internal
compute_ser_statistics <- function(data, params, model, l, ...) {
  UseMethod("compute_ser_statistics")
}
#' @keywords internal
compute_ser_statistics.default <- function(data, params, model, l, ...) {
  stop("compute_ser_statistics: no method for class '", class(data)[1], "'")
}

# Single effect regression posterior expected log-likelihood
#' @keywords internal
SER_posterior_e_loglik <- function(data, params, model, l) {
  UseMethod("SER_posterior_e_loglik")
}
#' @keywords internal
SER_posterior_e_loglik.default <- function(data, params, model, l) {
  stop("SER_posterior_e_loglik: no method for class '", class(data)[1], "'")
}

# Calculate posterior moments for single effect regression
#' @keywords internal
calculate_posterior_moments <- function(data, params, model, V, ...) {
  UseMethod("calculate_posterior_moments")
}
#' @keywords internal
calculate_posterior_moments.default <- function(data, params, model, V, ...) {
  stop("calculate_posterior_moments: no method for class '", class(data)[1], "'")
}

# Calculate KL divergence
#' @keywords internal
compute_kl <- function(data, params, model, l) {
  UseMethod("compute_kl")
}
#' @keywords internal
compute_kl.default <- function(data, params, model, l) {
  return(-model$lbf[l] + SER_posterior_e_loglik(data, params, model, l))
}

# Expected squared residuals
#' @keywords internal
get_ER2 <- function(data, model) {
  UseMethod("get_ER2")
}
#' @keywords internal
get_ER2.default <- function(data, model) {
  stop("get_ER2: no method for class '", class(data)[1], "'")
}

# Expected log-likelihood
#' @keywords internal
Eloglik <- function(data, model) {
  UseMethod("Eloglik")
}
#' @keywords internal
Eloglik.default <- function(data, model) {
  stop("Eloglik: no method for class '", class(data)[1], "'")
}

# Log-likelihood for prior variance optimization
#' @keywords internal
loglik <- function(data, params, model, V, ser_stats, ...) {
  UseMethod("loglik")
}
#' @keywords internal
loglik.default <- function(data, params, model, V, ser_stats, ...) {
  stop("loglik: no method for class '", class(data)[1], "'")
}

# Negative log-likelihood for optimization (handles both log and linear scales)
#' @keywords internal
neg_loglik <- function(data, params, model, V_param, ser_stats, ...) {
  UseMethod("neg_loglik")
}
#' @keywords internal
neg_loglik.default <- function(data, params, model, V_param, ser_stats, ...) {
  stop("neg_loglik: no method for class '", class(data)[1], "'")
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
update_fitted_values <- function(data, params, model, l) {
  UseMethod("update_fitted_values")
}
#' @keywords internal
update_fitted_values.default <- function(data, params, model, l) {
  stop("update_fitted_values: no method for class '", class(data)[1], "'")
}

# Update variance components
#' @keywords internal
update_variance_components <- function(data, params, model, ...) {
  UseMethod("update_variance_components")
}
#' @keywords internal
update_variance_components.default <- function(data, params, model, ...) {
  # Standard MLE estimation (MLE and MoM are equivalent for standard SuSiE)
  sigma2 <- est_residual_variance(data, model)
  return(list(sigma2 = sigma2))
}

# Update derived quantities after variance component changes
#' @keywords internal
update_derived_quantities <- function(data, params, model) {
  UseMethod("update_derived_quantities")
}
#' @keywords internal
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
get_scale_factors.default <- function(data, params, ...) {
  stop("get_scale_factors: no method for class '", class(data)[1], "'")
}

# Get intercept
#' @keywords internal
get_intercept <- function(data, params, model, ...) {
  UseMethod("get_intercept")
}
#' @keywords internal
get_intercept.default <- function(data, params, model, ...) {
  stop("get_intercept: no method for class '", class(data)[1], "'")
}

# Get fitted values
#' @keywords internal
get_fitted <- function(data, params, model, ...) {
  UseMethod("get_fitted")
}
#' @keywords internal
get_fitted.default <- function(data, params, model, ...) {
  return(NULL)
}

# Get credible sets
#' @keywords internal
get_cs <- function(data, params, model, ...) {
  UseMethod("get_cs")
}
#' @keywords internal
get_cs.default <- function(data, params, model, ...) {
  stop("get_cs: no method for class '", class(data)[1], "'")
}

# Get variable names
#' @keywords internal
get_variable_names <- function(data, model, ...) {
  UseMethod("get_variable_names")
}
#' @keywords internal
get_variable_names.default <- function(data, model, ...) {
  stop("get_variable_names: no method for class '", class(data)[1], "'")
}

# Get univariate z-scores
#' @keywords internal
get_zscore <- function(data, params, model, ...) {
  UseMethod("get_zscore")
}
#' @keywords internal
get_zscore.default <- function(data, params, model, ...) {
  return(NULL)
}

# Clean up model object by removing temporary computational fields
#' @keywords internal
cleanup_model <- function(data, params, model, ...) {
  UseMethod("cleanup_model")
}
#' @keywords internal
cleanup_model.default <- function(data, params, model, ...) {
  # Remove temporary fields common to all data types
  temp_fields <- c("null_weight", "predictor_weights", "prev_elbo", "prev_alpha",
                   "residuals", "fitted_without_l", "residual_variance")

  for (field in temp_fields) {
    if (field %in% names(model)) {
      model[[field]] <- NULL
    }
  }

  return(model)
}
