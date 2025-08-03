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
                         check_prior = FALSE) {

  # Validate method argument
  estimate_prior_method <- match.arg(estimate_prior_method)


  # Initialize model object
  model <- ibss_initialize(data = data, L = L,
                           scaled_prior_variance = scaled_prior_variance,
                           residual_variance = residual_variance,
                           prior_weights = prior_weights,
                           null_weight = null_weight,
                           model_init = model_init,
                           prior_tol = prior_tol,
                           residual_variance_upperbound = residual_variance_upperbound)


  # Initialize tracking
  elbo <- rep(as.numeric(NA), max_iter + 1)
  elbo[1] <- -Inf
  tracking <- list()

  # Main IBSS iteration loop
  for (iter in seq_len(max_iter)) {
    # Track iteration progress
    tracking <- extract_core(data, model, tracking, iter, track_fit)

    # Store previous model for convergence check
    model_prev <- model

    # Update all L effects
    model <- ibss_fit(data, model,
                      estimate_prior_variance = estimate_prior_variance,
                      estimate_prior_method = estimate_prior_method,
                      check_null_threshold = check_null_threshold,
                      check_prior = check_prior)

    # Calculate objective for tracking
    elbo[iter + 1] <- get_objective(data, model, verbose = verbose)

    # Handle convergence and variance updates. I chose this way because susie-inf
    # updates variance parameters then checks for convergence vs original susie
    # which checks for convergence then updates residual variance.
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

  # Set ELBO from iterations
  model$elbo <- elbo[2:(iter + 1)]

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
