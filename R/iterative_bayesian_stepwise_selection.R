# =============================================================================
# IBSS INITIALIZATION
#
# Initializes the SuSiE model object for Iterative Bayesian Stepwise Selection.
# Sets up model matrices, handles model_init, and prepares for IBSS.
# =============================================================================
#'
#' @param data Data object (individual, ss, or rss_lambda)
#' @param params Validated params object
#'
#' @return Initialized SuSiE model object with alpha, mu, mu2, V, sigma2, and fitted values
#'
#' @keywords internal
#' Initialize the IBSS model object
#'
#' Creates and initializes the model object for the IBSS algorithm,
#' including allocation of matrices for posterior quantities (alpha, mu, mu2),
#' prior variances, and fitted values. This is a building-block function
#' used by \code{\link{susie_workhorse}} and by downstream packages
#' that implement custom IBSS loops.
#'
#' @param data Data object (individual, ss, or rss_lambda).
#' @param params Validated params object.
#'
#' @return Initialized model object ready for the IBSS iteration loop.
#'
#' @importFrom utils modifyList
#' @export
#' @keywords internal
ibss_initialize <- function(data, params) {
  UseMethod("ibss_initialize")
}

#' @export
#' @keywords internal
ibss_initialize.default <- function(data, params) {

  # Set var(y)
  var_y <- get_var_y(data)

  # Adjust number of single effects if needed
  if (data$p < params$L) {
    params$L <- data$p
  }

  # Check & validate residual variance
  if (is.null(params$residual_variance)) {
    params$residual_variance <- var_y
  }
  # For multivariate models, residual_variance can be a matrix
  if (!is.matrix(params$residual_variance)) {
    if (!is.numeric(params$residual_variance)) {
      stop("Input residual variance sigma2 must be numeric.")
    }
    params$residual_variance <- as.numeric(params$residual_variance)
    if (length(params$residual_variance) != 1) {
      stop("Input residual variance sigma2 must be a scalar.")
    }
    if (params$residual_variance <= 0) {
      stop("Residual variance sigma2 must be positive (is your var(Y) zero?).")
    }
  }

  # Handle model initialization
  if (!is.null(params$model_init)) {
    # Validate the contents of model_init
    validate_init(data, params)

    # Prune effects with zero prior variance
    model_init_pruned <- prune_single_effects(params$model_init)

    # Adjust the number of effects
    adjustment <- adjust_L(params, model_init_pruned, var_y)
    params$L   <- adjustment$L

    # Create base model with all required fields
    mat_init <- initialize_susie_model(data, params, var_y)

    # Merge with adjusted model_init
    mat_init <- modifyList(mat_init, adjustment$model_init)

    # Reset iteration-specific values
    mat_init$KL  <- rep(as.numeric(NA), params$L)
    mat_init$lbf <- rep(as.numeric(NA), params$L)
  } else {
    # Create fresh model
    mat_init <- initialize_susie_model(data, params, var_y)
  }

  # Initialize fitted values and null index
  fitted     <- initialize_fitted(data, mat_init)
  null_index <- initialize_null_index(data, mat_init)

  # Preserve model class set by initialize_susie_model (e.g., "mvsusie")
  model_class <- class(mat_init)

  # Return assembled SuSiE object
  model <- c(mat_init,
             list(null_index = null_index),
             fitted)

  # Use the class from initialize_susie_model if it inherits from "susie",
  # otherwise default to "susie"
  if (inherits(mat_init, "susie")) {
    class(model) <- model_class
  } else {
    class(model) <- "susie"
  }
  model$converged <- FALSE

  return(model)
}

# =============================================================================
# IBSS FITTING
#
# Updates all L single effects in the SuSiE model for one IBSS iteration.
# Calls single_effect_update for each effect and validates prior variance estimates.
# =============================================================================
#'
#' @param data Data object (individual, ss, or rss_lambda)
#' @param params Validated params object
#' @param model Current SuSiE model object
#'
#' @return Updated SuSiE model object with new alpha, mu, mu2, V, lbf, KL, and
#' fitted values
#'
#' @keywords internal
#' @noRd
ibss_fit <- function(data, params, model) {

  # Repeat for each effect to update
  L <- nrow(model$alpha)
  if (L > 0) {
    for (l in seq_len(L)) {
      model <- single_effect_update(data, params, model, l)
    }
  }

  # Validate prior variance is reasonable
  validate_prior(data, params, model)

  return(model)
}

# =============================================================================
# IBSS FINALIZATION
#
# Finalizes the SuSiE model after convergence or maximum number of iterations
# reached. Computes credible sets, PIPs, intercept, fitted values, and z-scores.
# =============================================================================
#'
#' @param data Data object (individual, ss, or rss_lambda)
#' @param params Validated params object
#' @param model Final SuSiE model object from iterations
#' @param ... Additional finalization parameters
#'
#' @return Complete SuSiE model object with credible sets, PIPs, and summary statistics
#'
#' Finalize the IBSS model after convergence
#'
#' Computes credible sets, PIPs, z-scores, and cleans up temporary
#' fields from the model object. Building-block function for downstream
#' packages implementing custom IBSS loops.
#'
#' @param data Data object.
#' @param params Params object.
#' @param model Converged model object.
#' @param elbo ELBO values (optional).
#' @param iter Number of iterations completed.
#' @param tracking Tracking data (optional).
#'
#' @return Finalized model object with credible sets and PIPs.
#'
#' @export
#' @keywords internal
ibss_finalize <- function(data, params, model, elbo = NULL, iter = NA_integer_,
                          tracking = NULL) {

  # Append ELBO & iteration count to model output
  model$niter <- iter

  # Intercept & Fitted Values
  model$X_column_scale_factors <- get_scale_factors(data, params)
  model$intercept              <- get_intercept(data, params, model)
  model$fitted                 <- get_fitted(data, params, model)

  # Posterior Inclusion Probabilities, credible sets, z-scores
  model$sets <- get_cs(data, params, model)
  model$pip  <- susie_get_pip(model, prior_tol = params$prior_tol)
  model$z    <- get_zscore(data, params, model)

  # Tracking Across Iterations
  if (params$track_fit) model$trace <- tracking

  # Assign Variable Names
  model <- get_variable_names(data, model)

  # Stochastic LD diagnostics (from data → model, following sets/pip/z pattern)
  if (!is.null(data$stochastic_ld_diagnostics)) {
    model$stochastic_ld_diagnostics <- data$stochastic_ld_diagnostics
    # Store final-iteration per-variable penalty: v_j/sigma^2 = inflation - 1
    if (!is.null(model$shat2_inflation))
      model$stochastic_ld_diagnostics$per_variable_penalty <- model$shat2_inflation - 1
  }

  # Multi-panel omega weights
  if (!is.null(model$omega))
    model$omega_weights <- model$omega

  # Clean up temporary computational fields
  model <- cleanup_model(data, params, model)

  return(model)
}
