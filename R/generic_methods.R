# =============================================================================
#' @section DATA INITIALIZATION & CONFIGURATION
#'
#' Functions for data object setup, configuration, and preprocessing.
#' These prepare data objects for model fitting and handle data-specific
#' configurations like unmappable effects.
#'
#' Functions: configure_data, get_var_y
# =============================================================================

# Configure data object for specified method
#' @keywords internal
configure_data <- function(data) {
  UseMethod("configure_data")
}
configure_data.default <- function(data) {
  return(data)
}

# Get variance of y
#' @keywords internal
get_var_y <- function(data, ...) {
  UseMethod("get_var_y")
}
get_var_y.default <- function(data, ...) {
  stop("get_var_y: no method for class '", class(data)[1], "'")
}

# =============================================================================
#' @section MODEL INITIALIZATION & SETUP
#'
#' Functions for initializing model objects and setting up initial states.
#' These create model matrices, initialize fitted values, and prepare
#' the SuSiE model for iterative fitting.
#'
#' Functions: initialize_susie_model, initialize_fitted, validate_prior, track_ibss_fit
# =============================================================================

# Initialize susie model object
#' @keywords internal
initialize_susie_model <- function(data, ...) {
  UseMethod("initialize_susie_model")
}
initialize_susie_model.default <- function(data, ...) {
  stop("initialize_susie_model: no method for class '", class(data)[1], "'")
}

# Initialize fitted values
#' @keywords internal
initialize_fitted <- function(data, alpha, mu) {
  UseMethod("initialize_fitted")
}
initialize_fitted.default <- function(data, alpha, mu, ...) {
  stop("initialize_fitted: no method for class '", class(data)[1], "'")
}

# Validate prior variance
#' @keywords internal
validate_prior <- function(data, model, check_prior, ...) {
  UseMethod("validate_prior")
}
validate_prior.default <- function(data, model, check_prior, ...) {
  invisible(TRUE)
}

# Track core parameters of a susie fit across iterations
#' @keywords internal
track_ibss_fit <- function(data, model, tracking, iter, track_fit, ...) {
  UseMethod("track_ibss_fit")
}
track_ibss_fit.default <- function(data, model, tracking, iter, track_fit, ...) {
  if (isTRUE(track_fit)) {
    tracking[[iter]] <- list(
      alpha = model$alpha,
      niter = iter,
      V = model$V,
      sigma2 = model$sigma2
    )
  }
  return(tracking)
}

# =============================================================================
#' @section SINGLE EFFECT REGRESSION & ELBO
#'
#' Core functions for single effect regression computation and ELBO calculation.
#' These handle the mathematical core of SuSiE including residual computation, SER
#' statistics, posterior moments, and log-likelihood calculations for the ELBO.
#'
#' Functions: compute_residuals, compute_ser_statistics, SER_posterior_e_loglik,
#' calculate_posterior_moments, compute_kl, get_ER2, Eloglik, loglik, neg_loglik
# =============================================================================

# Compute residuals for single effect regression
#' @keywords internal
compute_residuals <- function(data, model, l, ...) {
  UseMethod("compute_residuals")
}
compute_residuals.default <- function(data, model, l, ...) {
  stop("compute_residuals: no method for class '", class(data)[1], "'")
}

# Compute SER statistics (betahat, shat2)
#' @keywords internal
compute_ser_statistics <- function(data, model, residual_variance, ...) {
  UseMethod("compute_ser_statistics")
}
compute_ser_statistics.default <- function(data, model, residual_variance, ...) {
  stop("compute_ser_statistics: no method for class '", class(data)[1], "'")
}

# Single effect regression posterior expected log-likelihood
#' @keywords internal
SER_posterior_e_loglik <- function(data, model, Eb, Eb2) {
  UseMethod("SER_posterior_e_loglik")
}
SER_posterior_e_loglik.default <- function(data, model, Eb, Eb2) {
  stop("SER_posterior_e_loglik: no method for class '", class(data)[1], "'")
}

# Calculate posterior moments for single effect regression
#' @keywords internal
calculate_posterior_moments <- function(data, ...) {
  UseMethod("calculate_posterior_moments")
}
calculate_posterior_moments.default <- function(data, ...) {
  stop("calculate_posterior_moments: no method for class '", class(data)[1], "'")
}

# Calculate KL divergence
#' @keywords internal
compute_kl <- function(data, model, l) {
  UseMethod("compute_kl")
}
compute_kl.default <- function(data, model, l) {
  return(-model$lbf[l] + SER_posterior_e_loglik(data, model,
                                                model$alpha[l, ] * model$mu[l, ],
                                                model$alpha[l, ] * model$mu2[l, ]))
}

# Expected squared residuals
#' @keywords internal
get_ER2 <- function(data, model) {
  UseMethod("get_ER2")
}
get_ER2.default <- function(data, model) {
  stop("get_ER2: no method for class '", class(data)[1], "'")
}

# Expected log-likelihood
#' @keywords internal
Eloglik <- function(data, model) {
  UseMethod("Eloglik")
}
Eloglik.default <- function(data, model) {
  stop("Eloglik: no method for class '", class(data)[1], "'")
}

# Log-likelihood for prior variance optimization
#' @keywords internal
loglik <- function(data, ...) {
  UseMethod("loglik")
}
loglik.default <- function(data, ...) {
  stop("loglik: no method for class '", class(data)[1], "'")
}

# Negative log-likelihood for optimization (handles both log and linear scales)
#' @keywords internal
neg_loglik <- function(data, ...) {
  UseMethod("neg_loglik")
}
neg_loglik.default <- function(data, ...) {
  stop("neg_loglik: no method for class '", class(data)[1], "'")
}

# =============================================================================
#' @section MODEL UPDATES & FITTING
#'
#' Functions for iterative model updates and variance component estimation.
#' These handle the dynamic aspects of model fitting including fitted value
#' updates and variance component estimation.
#'
#' Functions: update_fitted_values, update_variance_components, update_derived_quantities
# =============================================================================

# Update fitted values
#' @keywords internal
update_fitted_values <- function(data, model, l) {
  UseMethod("update_fitted_values")
}
update_fitted_values.default <- function(data, model, l) {
  stop("update_fitted_values: no method for class '", class(data)[1], "'")
}

# Update variance components
#' @keywords internal
update_variance_components <- function(data, model, estimate_method = "MLE") {
  UseMethod("update_variance_components")
}
update_variance_components.default <- function(data, model, estimate_method = "MLE") {
  stop("update_variance_components: no method for class '", class(data)[1], "'")
}

# Update derived quantities after variance component changes
#' @keywords internal
update_derived_quantities <- function(data, model) {
  UseMethod("update_derived_quantities")
}
update_derived_quantities.default <- function(data, model) {
  return(model)
}

# =============================================================================
#' @section OUTPUT GENERATION & POST-PROCESSING
#'
#' Functions for generating final results and summary statistics.
#' These process fitted models into interpretable outputs including
#' credible sets, variable names, and fitted values.
#'
#' Functions: get_scale_factors, get_intercept, get_fitted, get_cs,
#' get_variable_names, get_zscore
# =============================================================================

# Get column scale factors
#' @keywords internal
get_scale_factors <- function(data, ...) {
  UseMethod("get_scale_factors")
}
get_scale_factors.default <- function(data, ...) {
  stop("get_scale_factors: no method for class '", class(data)[1], "'")
}

# Get intercept
#' @keywords internal
get_intercept <- function(data, model, ...) {
  UseMethod("get_intercept")
}
get_intercept.default <- function(data, model, ...) {
  stop("get_intercept: no method for class '", class(data)[1], "'")
}

# Get fitted values
#' @keywords internal
get_fitted <- function(data, model, ...) {
  UseMethod("get_fitted")
}
get_fitted.default <- function(data, model, ...) {
  stop("get_fitted: no method for class '", class(data)[1], "'")
}

# Get credible sets
#' @keywords internal
get_cs <- function(data, model, ...) {
  UseMethod("get_cs")
}
get_cs.default <- function(data, model, ...) {
  stop("get_cs: no method for class '", class(data)[1], "'")
}

# Get variable names
#' @keywords internal
get_variable_names <- function(data, model, ...) {
  UseMethod("get_variable_names")
}
get_variable_names.default <- function(data, model, ...) {
  stop("get_variable_names: no method for class '", class(data)[1], "'")
}

# Get univariate z-scores
#' @keywords internal
get_zscore <- function(data, model, ...) {
  UseMethod("get_zscore")
}
get_zscore.default <- function(data, model, ...) {
  return(NULL)
}
