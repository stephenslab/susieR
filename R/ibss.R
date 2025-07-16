# Core IBSS functions: ibss_initialize(), ibss_fit(), and ibss_finalize()

ibss_initialize <- function(data,
                            L                     = min(10, data$p),
                            scaled_prior_variance = 0.2,
                            residual_variance     = NULL,
                            prior_weights         = NULL,
                            null_weight           = NULL,
                            model_init            = NULL,
                            prior_tol             = 1e-9,
                            residual_variance_upperbound = Inf){
  # Define p and var_y
  p <- data$p
  var_y <- get_var_y(data)

  # Validate prior tolerance threshold
  if (!is.numeric(prior_tol) || length(prior_tol) != 1)
    stop("prior_tol must be a numeric scalar")
  if (prior_tol < 0)
    stop("prior_tol must be non-negative")
  if (prior_tol > 1)
    stop("prior_tol cannot exceed 1 (a single effect cannot account for more than 100% of outcome variance)")

  # Validate residual_variance_upperbound
  if (!is.numeric(residual_variance_upperbound) || length(residual_variance_upperbound) != 1)
    stop("residual_variance_upperbound must be a numeric scalar")
  if (residual_variance_upperbound <= 0)
    stop("residual_variance_upperbound must be positive")

  # Check prior variance
  if (!is.numeric(scaled_prior_variance) || scaled_prior_variance < 0)
    stop("Scaled prior variance should be positive number")

  # Check prior weights
  if (is.null(prior_weights))
    prior_weights <- rep(1 / p, p)

  # Check number of single effects
  if (p < L)
    L = p

  # Check & validate residual variance
  if (is.null(residual_variance))
    residual_variance <- var_y
  if (!is.numeric(residual_variance))
    stop("Input residual variance sigma2 must be numeric")
  residual_variance <- as.numeric(residual_variance)
  if (length(residual_variance) != 1)
    stop("Input residual variance sigma2 must be a scalar")
  if (residual_variance <= 0)
    stop("Residual variance sigma2 must be positive (is your var(Y) zero?)")

  # Check for pre-initialized model
  if (!missing(model_init) && !is.null(model_init)) {

    # Validate the contents of model_init
    validate_init(model_init, L, null_weight)

    # Prune effects with zero prior variance
    model_init <- susie_prune_single_effects(model_init)

    # Adjust the number of effects
    adjustment <- adjust_L(model_init, L, V = rep(scaled_prior_variance * var_y, L))

    # Return adjusted initialized model and number of effects
    model_init <- adjustment$model_init
    L          <- adjustment$L
  }

  # Build blank matrices
  mat_init <- initialize_matrices(data, L,
                                  scaled_prior_variance, var_y,
                                  residual_variance, prior_weights)

  # Handle model initialization using modifyList
  if (!missing(model_init) && !is.null(model_init)) {
    mat_init_list <- as.list(mat_init)
    model_init_list <- as.list(model_init)

    # Map field names between old and new conventions
    if (!is.null(model_init_list$sigma2)) {
      model_init_list$residual_variance <- model_init_list$sigma2
    }
    if (!is.null(model_init_list$pi)) {
      model_init_list$prior_weights <- model_init_list$pi
      model_init_list$pi <- NULL
    }

    mat_init_list <- modifyList(mat_init_list, model_init_list)
    mat_init <- mat_init_list

    # Reset KL and lbf to NA
    mat_init$KL <- rep(as.numeric(NA), L)
    mat_init$lbf <- rep(as.numeric(NA), L)
  }

  # Initialize fitted values
  fitted <- initialize_fitted(data, mat_init$alpha, mat_init$mu)

  # Set null index (for refine step)
  null_index <- initialize_null_index(null_weight, p)

  # Return assembled SuSiE object
  model <- c(
    mat_init,
    list(null_index = null_index),
    fitted
  )

  class(model) <- "susie"

  return(model)
}

ibss_fit = function(data, model,
                    estimate_prior_variance = TRUE,
                    estimate_prior_method   = c("optim","EM","simple"),
                    check_null_threshold    = 0,
                    check_prior             = FALSE){
  estimate_prior_method <- match.arg(estimate_prior_method)

  # Set prior variance estimation if NA
  if (!estimate_prior_variance)
    estimate_prior_method <- "none"

  # Repeat for each effect to update
  L <- nrow(model$alpha)
  if (L > 0)
    for (l in seq_len(L)){
      model <- single_effect_update(data, model, l,
                                  optimize_V = estimate_prior_method,
                                  check_null_threshold = check_null_threshold)
    }

  # Validate prior variance is reasonable
  validate_prior(data, model, check_prior)

  return(model)
}

ibss_finalize <- function(data,
                          model,
                          coverage               = 0.95,
                          min_abs_corr           = 0.50,
                          median_abs_corr        = NULL,
                          prior_tol              = 1e-9,
                          n_purity               = 100,
                          compute_univariate_zscore = FALSE,
                          intercept              = TRUE,
                          standardize            = TRUE,
                          elbo                   = NULL,
                          iter                   = NA_integer_,
                          null_weight            = NULL,
                          track_fit              = FALSE,
                          tracking               = NULL) {

  # Append ELBO & iteration count to model output
  model$niter <- iter

  # Intercept & Fitted Values
  model$X_column_scale_factors <- get_scale_factors(data)
  model$intercept              <- get_intercept(data, model, intercept)
  model$fitted                 <- get_fitted(data, model, intercept)

  # Tracking Across Iterations
  if (track_fit)
    model$trace = tracking

  # Credible Sets
  model$sets <- get_cs(data, model, coverage, min_abs_corr, n_purity)

  # Posterior Inclusion Probabilities
  model$pip <- get_pip(data, model, coverage, min_abs_corr, prior_tol)

  # Set pi field from prior_weights
  if (is.null(model$pi)) {
    model$pi <- model$prior_weights
  }

  # Assign Variable Names
  model <- get_variable_names(data, model, null_weight)

  # Compute z-scores
  model$z <- get_zscore(data, model, compute_univariate_zscore,
                        intercept, standardize, null_weight)

  return(model)
}
