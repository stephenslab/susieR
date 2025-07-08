
susie_engine <- function(data,
                         L, intercept = TRUE, standardize = TRUE,
                         scaled_prior_variance, residual_variance,
                         prior_weights, null_weight, model_init = NULL,
                         estimate_prior_variance,
                         estimate_prior_method = c("optim", "EM", "simple"),
                         check_null_threshold,
                         estimate_residual_variance,
                         residual_variance_lowerbound,
                         residual_variance_upperbound,
                         max_iter, tol, verbose, track_fit,
                         coverage, min_abs_corr,
                         prior_tol, n_purity, compute_univariate_zscore = FALSE,
                         check_prior = FALSE,
                         refine = FALSE) {

  # Validate method argument
  estimate_prior_method <- match.arg(estimate_prior_method)
  
  # Validate prior_tol
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

  # Initialize model object
  model <- ibss_initialize(data = data, L = L,
                           scaled_prior_variance = scaled_prior_variance,
                           residual_variance = residual_variance,
                           prior_weights = prior_weights,
                           null_weight = null_weight,
                           model_init = model_init)

  # Initialize tracking
  elbo <- rep(as.numeric(NA), max_iter + 1)
  elbo[1] <- -Inf
  tracking <- list()

  # Main IBSS iteration loop
  for (iter in seq_len(max_iter)) {
    # Track iteration progress
    tracking <- susie_extract_core(data, model, tracking, iter, track_fit)

    # Store previous model for convergence check
    model_prev <- model

    # Update all L effects
    model <- ibss_fit(data, model,
                      estimate_prior_variance = estimate_prior_variance,
                      estimate_prior_method = estimate_prior_method,
                      check_null_threshold = check_null_threshold,
                      check_prior = check_prior)

    # Calculate objective for tracking
    elbo[iter + 1] <- get_objective(data, model)

    # Handle convergence and variance updates
    result <- handle_convergence_and_variance(data, model, model_prev, elbo[iter], elbo[iter + 1],
                                               tol, estimate_residual_variance,
                                               residual_variance_lowerbound, residual_variance_upperbound)
    data <- result$data
    model <- result$model

    if (result$converged) {
      model$converged <- TRUE
      break
    }
  }

  # Check final convergence status
  if (is.null(model$converged)) {
    warning(paste("IBSS algorithm did not converge in", max_iter, "iterations!"))
    model$converged <- FALSE
  }

  model$elbo <- elbo[2:(iter + 1)]

  rerun <- function(prior_weights, model_init = NULL) {
    susie_engine(data                      = data,
                 L                         = L,
                 intercept                 = intercept,
                 standardize               = standardize,
                 scaled_prior_variance     = scaled_prior_variance,
                 residual_variance         = residual_variance,
                 prior_weights             = prior_weights,
                 null_weight               = null_weight,
                 model_init                = model_init,
                 estimate_prior_variance   = estimate_prior_variance,
                 estimate_prior_method     = estimate_prior_method,
                 check_null_threshold      = check_null_threshold,
                 estimate_residual_variance = estimate_residual_variance,
                 residual_variance_lowerbound = residual_variance_lowerbound,
                 residual_variance_upperbound = residual_variance_upperbound,
                 max_iter                  = max_iter,
                 tol                       = tol,
                 verbose                   = verbose,
                 track_fit                 = FALSE,
                 coverage                  = coverage,
                 min_abs_corr              = min_abs_corr,
                 prior_tol                 = prior_tol,
                 n_purity                  = n_purity,
                 compute_univariate_zscore = compute_univariate_zscore,
                 check_prior               = check_prior,
                 refine                    = FALSE)
  }

  if (refine) {
    model <- run_refine(data, model, rerun)
  }

  model <- ibss_finalize(data, model,
                         coverage = coverage,
                         min_abs_corr = min_abs_corr,
                         prior_tol = prior_tol,
                         n_purity = n_purity,
                         compute_univariate_zscore = compute_univariate_zscore,
                         intercept = intercept,
                         standardize = standardize,
                         elbo = elbo,
                         iter = iter,
                         null_weight = null_weight,
                         track_fit = track_fit,
                         tracking = tracking)

  return(model)
}
