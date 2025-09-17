# Core IBSS functions: ibss_initialize(), ibss_fit(), and ibss_finalize()

ibss_initialize <- function(data, params) {

  # Define p and var_y
  p <- data$p
  var_y <- get_var_y(data)
  L <- params$L

  # Adjust number of single effects if needed
  if (p < params$L) {
    params$L <- p
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

  # Initialize prior weights if needed
  if (is.null(data$prior_weights)) {
    prior_weights <- rep(1 / p, p)
  } else {
    prior_weights <- data$prior_weights
  }

  # Handle model initialization
  if (!is.null(params$model_init)) {
    # Validate the contents of model_init
    validate_init(params$model_init, params$L, data$null_weight)

    # Prune effects with zero prior variance
    model_init_pruned <- susie_prune_single_effects(params$model_init)

    # Adjust the number of effects
    adjustment <- adjust_L(model_init_pruned, params$L, V = rep(params$scaled_prior_variance * var_y, params$L))
    params$L <- adjustment$L

    # Create base model with all required fields  
    mat_init <- initialize_susie_model(data, params, params$L, params$scaled_prior_variance, var_y,
                                       params$residual_variance, prior_weights)

    # Merge with adjusted model_init
    mat_init <- modifyList(mat_init, adjustment$model_init)

    # Reset iteration-specific values
    mat_init$KL <- rep(as.numeric(NA), params$L)
    mat_init$lbf <- rep(as.numeric(NA), params$L)
  } else {
    # Create fresh model
    mat_init <- initialize_susie_model(data, params, params$L, params$scaled_prior_variance, var_y,
                                       params$residual_variance, prior_weights)
  }

  # Initialize fitted values
  fitted <- initialize_fitted(data, params, mat_init$alpha, mat_init$mu)

  # Set null index
  null_index <- initialize_null_index(data$null_weight, p)

  # Return assembled SuSiE object
  model <- c(mat_init,
             list(null_index = null_index),
             fitted)

  class(model) <- "susie"

  return(model)
}

ibss_fit <- function(data, params, model) {

  estimate_prior_method <- match.arg(params$estimate_prior_method, c("optim", "EM", "simple"))

  # Set prior variance estimation if NA
  if (!params$estimate_prior_variance) {
    estimate_prior_method <- "none"
  }

  # Repeat for each effect to update
  L <- nrow(model$alpha)
  if (L > 0) {
    for (l in seq_len(L)) {
      model <- single_effect_update(data, params, model, l,
                                     optimize_V = estimate_prior_method,
                                     check_null_threshold = params$check_null_threshold)
    }
  }

  # Validate prior variance is reasonable
  validate_prior(data, params, model)

  return(model)
}

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
