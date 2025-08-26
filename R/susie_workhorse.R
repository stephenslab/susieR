susie_workhorse <- function(data,
                            L, intercept = TRUE, standardize = TRUE,
                            scaled_prior_variance, residual_variance,
                            prior_weights, null_weight, model_init = NULL,
                            estimate_prior_variance,
                            estimate_prior_method = c("optim", "EM", "simple"),
                            check_null_threshold,
                            estimate_residual_variance,
                            estimate_residual_method = "MLE",
                            residual_variance_lowerbound,
                            residual_variance_upperbound,
                            max_iter, tol, verbose, track_fit,
                            coverage, min_abs_corr,
                            prior_tol, n_purity, compute_univariate_zscore = FALSE,
                            check_prior = FALSE, convergence_method = "elbo",
                            alpha0 = 0, beta0 = 0, refine = FALSE) {
  # Validate method argument
  estimate_prior_method <- match.arg(estimate_prior_method)

  # Initialize model object
  model <- ibss_initialize(
    data = data, L = L,
    scaled_prior_variance = scaled_prior_variance,
    residual_variance = residual_variance,
    prior_weights = prior_weights,
    null_weight = null_weight,
    model_init = model_init,
    prior_tol = prior_tol,
    residual_variance_upperbound = residual_variance_upperbound
  )

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
      check_prior = check_prior
    )

    # Calculate objective for tracking
    elbo[iter + 1] <- get_objective(data, model, verbose = verbose)

    # Check for convergence
    converged <- check_convergence(model_prev, model, elbo, tol, convergence_method, iter)

    if (converged) {
      model$converged <- TRUE
      break
    }

    # Update variance components if not converged and estimation is requested
    if (estimate_residual_variance) {
      result <- update_model_variance(
        data, model, residual_variance_lowerbound,
        residual_variance_upperbound, estimate_residual_method
      )
      data <- result$data
      model <- result$model
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
    tracking = tracking
  )

  # Run refinement if requested
  if (refine && !is.null(model$sets) && length(model$sets$cs) > 0) {
    model <- run_refine(
      model = model,
      data = data,
      L = L,
      intercept = intercept,
      standardize = standardize,
      scaled_prior_variance = scaled_prior_variance,
      residual_variance = residual_variance,
      prior_weights = prior_weights,
      null_weight = null_weight,
      model_init = model_init,
      estimate_prior_variance = estimate_prior_variance,
      estimate_prior_method = estimate_prior_method,
      check_null_threshold = check_null_threshold,
      estimate_residual_variance = estimate_residual_variance,
      estimate_residual_method = estimate_residual_method,
      residual_variance_lowerbound = residual_variance_lowerbound,
      residual_variance_upperbound = residual_variance_upperbound,
      max_iter = max_iter,
      tol = tol,
      verbose = verbose,
      track_fit = track_fit,
      coverage = coverage,
      min_abs_corr = min_abs_corr,
      prior_tol = prior_tol,
      n_purity = n_purity,
      compute_univariate_zscore = compute_univariate_zscore,
      check_prior = check_prior,
      convergence_method = convergence_method,
      alpha0 = alpha0,
      beta0 = beta0
    )
  }

  return(model)
}
