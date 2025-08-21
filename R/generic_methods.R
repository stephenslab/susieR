# Generic methods for S3 dispatch

# Single effect regression posterior expected log-likelihood
SER_posterior_e_loglik <- function(data, model, R, Eb, Eb2) {
  UseMethod("SER_posterior_e_loglik")
}
SER_posterior_e_loglik.default <- function(data, model, R, Eb, Eb2) {
  stop("SER_posterior_e_loglik: no method for class '", class(data)[1], "'")
}

# Expected squared residuals
get_ER2 <- function(data, model) {
  UseMethod("get_ER2")
}
get_ER2.default <- function(data, model) {
  stop("get_ER2: no method for class '", class(data)[1], "'")
}

# Single effect update
single_effect_update <- function(data, model, l, ...) {
  UseMethod("single_effect_update")
}
single_effect_update.default <- function(data, model, l, ...) {
  stop("single_effect_update: no method for class '", class(data)[1], "'")
}

# Get column scale factors
get_scale_factors <- function(data, ...) {
  UseMethod("get_scale_factors")
}
get_scale_factors.default <- function(data, ...) {
  stop("get_scale_factors: no method for class '", class(data)[1], "'")
}

# Get intercept
get_intercept <- function(data, model, ...) {
  UseMethod("get_intercept")
}
get_intercept.default <- function(data, model, ...) {
  stop("get_intercept: no method for class '", class(data)[1], "'")
}

# Get fitted values
get_fitted <- function(data, model, ...) {
  UseMethod("get_fitted")
}
get_fitted.default <- function(data, model, ...) {
  stop("get_fitted: no method for class '", class(data)[1], "'")
}

# Get credible sets
get_cs <- function(data, model, ...) {
  UseMethod("get_cs")
}
get_cs.default <- function(data, model, ...) {
  stop("get_cs: no method for class '", class(data)[1], "'")
}

# Get variable names
get_variable_names <- function(data, model, ...) {
  UseMethod("get_variable_names")
}
get_variable_names.default <- function(data, model, ...) {
  stop("get_variable_names: no method for class '", class(data)[1], "'")
}

# Get univariate z-scores
get_zscore <- function(data, model, ...) {
  UseMethod("get_zscore")
}
get_zscore.default <- function(data, model, ...) {
  stop("get_zscore: no method for class '", class(data)[1], "'")
}

# Initialize fitted values
initialize_fitted <- function(data, alpha, mu) {
  UseMethod("initialize_fitted")
}
initialize_fitted.default <- function(data, alpha, mu, ...) {
  stop("initialize_fitted: no method for class '", class(data)[1], "'")
}

# Get variance of y
get_var_y <- function(data, ...) {
  UseMethod("get_var_y")
}
get_var_y.default <- function(data, ...) {
  stop("get_var_y: no method for class '", class(data)[1], "'")
}

# Validate prior variance
validate_prior <- function(data, model, check_prior, ...) {
  UseMethod("validate_prior")
}
validate_prior.default <- function(data, model, check_prior, ...) {
  stop("validate_prior: no method for class '", class(data)[1], "'")
}

# Extract core parameters of a susie fit across iterations
extract_core <- function(data, model, tracking, iter, track_fit, ...) {
  UseMethod("extract_core")
}
extract_core.default <- function(data, model, tracking, iter, track_fit, ...) {
  stop("extract_core: no method for class '", class(data)[1], "'")
}

# Initialize susie model
initialize_susie_model <- function(data, ...) {
  UseMethod("initialize_susie_model")
}
initialize_susie_model.default <- function(data, ...) {
  stop("initialize_susie_model: no method for class '", class(data)[1], "'")
}

# Configure data object for specified method
configure_data <- function(data) {
  UseMethod("configure_data")
}
configure_data.default <- function(data) {
  stop("configure_data: no method for class '", class(data)[1], "'")
}

# Add eigenvalue decomposition to data objects
add_eigen_decomposition <- function(data) {
  UseMethod("add_eigen_decomposition")
}
add_eigen_decomposition.default <- function(data) {
  stop("add_eigen_decomposition: no method for class '", class(data)[1], "'")
}

# Update variance components
update_variance_components <- function(data, model, estimate_method = "MLE") {
  UseMethod("update_variance_components")
}
update_variance_components.default <- function(data, model, estimate_method = "MLE") {
  stop("update_variance_components: no method for class '", class(data)[1], "'")
}

# Update derived quantities after variance component changes
update_derived_quantities <- function(data, model) {
  UseMethod("update_derived_quantities")
}
update_derived_quantities.default <- function(data, model) {
  stop("update_derived_quantities: no method for class '", class(data)[1], "'")
}

# Expected log-likelihood
Eloglik <- function(data, model) {
  UseMethod("Eloglik")
}
Eloglik.default <- function(data, model) {
  stop("Eloglik: no method for class '", class(data)[1], "'")
}

# Log-likelihood for prior variance optimization
loglik <- function(data, ...) {
  UseMethod("loglik")
}
loglik.default <- function(data, ...) {
  stop("loglik: no method for class '", class(data)[1], "'")
}

# Negative log-likelihood in log scale for optimization
neg_loglik_logscale <- function(data, ...) {
  UseMethod("neg_loglik_logscale")
}
neg_loglik_logscale.default <- function(data, ...) {
  stop("neg_loglik_logscale: no method for class '", class(data)[1], "'")
}
