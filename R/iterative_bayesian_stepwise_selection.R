# =============================================================================
#' @section IBSS INITIALIZATION
#'
#' Initializes the SuSiE model object for Iterative Bayesian Stepwise Selection.
#' Sets up model matrices, handles model_init, and prepares for IBSS.
# =============================================================================
#'
#' @param data Data object (individual, ss, or rss_lambda)
#' @param params Validated params object
#'
#' @return Initialized SuSiE model object with alpha, mu, mu2, V, sigma2, and fitted values
#'
#' @keywords internal
#' @noRd
ibss_initialize <- function(data, params) {

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
  null_index <- initialize_null_index(data)

  # Return assembled SuSiE object
  model <- c(mat_init,
             list(null_index = null_index),
             fitted)

  class(model) <- "susie"

  return(model)
}

# =============================================================================
#' @section IBSS FITTING
#'
#' Updates all L single effects in the SuSiE model for one IBSS iteration.
#' Calls single_effect_update for each effect and validates prior variance estimates.
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
#' @section IBSS FINALIZATION
#'
#' Finalizes the SuSiE model after convergence or maximum number of iterations
#' reached. Computes credible sets, PIPs, intercept, fitted values, and z-scores.
# =============================================================================
#'
#' @param data Data object (individual, ss, or rss_lambda)
#' @param params Validated params object
#' @param model Final SuSiE model object from iterations
#' @param ... Additional finalization parameters
#'
#' @return Complete SuSiE model object with credible sets, PIPs, and summary statistics
#'
#' @keywords internal
#' @noRd
ibss_finalize <- function(data, params, model, elbo = NULL, iter = NA_integer_,
                          tracking = NULL) {

  # Append ELBO & iteration count to model output
  model$niter <- iter

  # Intercept & Fitted Values
  model$X_column_scale_factors <- get_scale_factors(data, params)
  model$intercept <- get_intercept(data, params, model, params$intercept)
  model$fitted <- get_fitted(data, params, model, params$intercept)

  # Tracking Across Iterations
  if (params$track_fit) {
    model$trace <- tracking
  }

  # Credible Sets
  model$sets <- get_cs(data, model, params$coverage, params$min_abs_corr, params$n_purity)

  # Posterior Inclusion Probabilities
  model$pip <- susie_get_pip(model, prune_by_cs = FALSE, prior_tol = params$prior_tol)

  # Set pi field from prior_weights
  if (is.null(model$pi)) {
    model$pi <- model$prior_weights
  }

  # Assign Variable Names
  model <- get_variable_names(data, model, data$null_weight)

  # Compute z-scores
  model$z <- get_zscore(data, model, params$compute_univariate_zscore,
                        params$intercept, params$standardize, data$null_weight)

  return(model)
}
