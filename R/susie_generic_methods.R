# This file contains all of the generic methods that will be used for
# different data types.

# Posterior Expected Loglikelihood for a SER
SER_posterior_e_loglik <- function(data, model, R, Eb, Eb2)
  UseMethod("SER_posterior_e_loglik")
SER_posterior_e_loglik.default <- function(data, R, model, Eb, Eb2)
  stop("SER_posterior_e_loglik: no method for class '", class(data)[1], "'")

# Expected Squared Residuals
get_ER2 <- function(data, model) UseMethod("get_ER2")
get_ER2.default <- function(data, model)
  stop("get_ER2: no method for class '", class(data)[1], "'")

# Expected Log-likelihood
Eloglik <- function(data, model) UseMethod("Eloglik")
Eloglik.default <- function(data, model)
  stop("Eloglik: no method for class '", class(data)[1], "'")

# Get Objective
get_objective <- function(data, model) UseMethod("get_objective")
get_objective.default <- function(data, model)
  stop("get_objective: no method for class '", class(data)[1], "'")

# Estimate Residual Variance Objective
est_residual_variance <- function(data, model) UseMethod("est_residual_variance")
est_residual_variance.default <- function(data, model)
  stop("est_residual_variance: no method for class '", class(data)[1], "'")

# Single Effect Update
single_effect_update <- function(data, model, l, ...)
  UseMethod("single_effect_update")
single_effect_update.default <- function(data, ...)
  stop("single_effect_update: no method for class '", class(data)[1], "'")

# Get Column Scale factors
get_scale_factors <- function(data, ...)
  UseMethod("get_scale_factors")
get_scale_factors.default <- function(data, ...)
  stop("get_scale_factors: no method for class '", class(data)[1], "'")

# Get Intercept
get_intercept <- function(data, model, ...)
  UseMethod("get_intercept")
get_intercept.default <- function(data, ...)
  stop("get_intercept: no method for class '", class(data)[1], "'")

# Get Fitted
get_fitted <- function(data, model, ...)
  UseMethod("get_fitted")
get_fitted.default <- function(data, ...)
  stop("get_fitted: no method for class '", class(data)[1], "'")

# Get CS
get_cs <- function(data, model, ...)
  UseMethod("get_cs")
get_cs.default <- function(data, ...)
  stop("get_cs: no method for class '", class(data)[1], "'")

# Get PIP
get_pip <- function(data, model, ...)
  UseMethod("get_pip")
get_pip.default <- function(data, model, ...)
  stop("get_pip: no method for class '", class(data)[1], "'")

# Set Variable Names
get_variable_names <- function(data, model, ...)
  UseMethod("get_variable_names")
get_variable_names.default <- function(data, model, ...)
  stop("get_variable_names: no method for class '", class(data)[1], "'")

# Get univariate z-scores
get_zscore <- function(data, model, ...)
  UseMethod("get_zscore")
get_zscore.default <- function(data, model, ...)
  stop("get_zscore: no method for class '", class(data)[1], "'")
