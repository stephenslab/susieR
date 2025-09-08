susie_workhorse <- function(data, L, intercept = TRUE, standardize = TRUE,
                            scaled_prior_variance, residual_variance,
                            prior_weights, null_weight, model_init = NULL,
                            estimate_prior_variance,
                            estimate_prior_method = c("optim", "EM", "simple"),
                            check_null_threshold, estimate_residual_variance,
                            estimate_residual_method = "MLE",
                            residual_variance_lowerbound,
                            residual_variance_upperbound,
                            max_iter, tol, verbose, track_fit,
                            coverage, min_abs_corr, prior_tol, n_purity,
                            compute_univariate_zscore = FALSE,
                            check_prior = FALSE, convergence_method = "elbo",
                            alpha0 = 0, beta0 = 0, refine = FALSE) {

  # Validate method argument
  estimate_prior_method <- match.arg(estimate_prior_method)

  # Initialize model object
  model <- ibss_initialize(
    data, L, scaled_prior_variance, residual_variance,
    prior_weights, null_weight, model_init,
    prior_tol, residual_variance_upperbound
  )

  # Initialize tracking
  elbo <- rep(as.numeric(NA), max_iter + 1)
  elbo[1] <- -Inf
  tracking <- list()

  # Main IBSS iteration loop
  for (iter in seq_len(max_iter)) {
    # Track iteration progress
    tracking <- track_ibss_fit(data, model, tracking, iter, track_fit)

    # Store previous model parameters for convergence check
    model$prev_elbo  <- elbo[iter]
    model$prev_alpha <- model$alpha

    # Update all L effects
    model <- ibss_fit(data, model, estimate_prior_variance,
                      estimate_prior_method, check_null_threshold,
                      check_prior)

    # Calculate objective for tracking
    elbo[iter + 1] <- get_objective(data, model, verbose = verbose)

    # Check for convergence
    converged <- check_convergence(model, elbo, tol, convergence_method, iter)

    if (converged) {
      model$converged <- TRUE
      break
    }

    # Update variance components if not converged and estimation is requested
    if (estimate_residual_variance) {
      model <- update_model_variance(data, model, residual_variance_lowerbound,
                                     residual_variance_upperbound, estimate_residual_method)
    }

  }

  # Check final convergence status
  if (is.null(model$converged)) {
    warning(paste("IBSS algorithm did not converge in", max_iter, "iterations!\n"))
    model$converged <- FALSE
  }

  # Set ELBO from iterations
  model$elbo <- elbo[2:(iter + 1)]

  model <- ibss_finalize(data, model, coverage, min_abs_corr,
                         median_abs_corr = NULL, prior_tol, n_purity,
                         compute_univariate_zscore, intercept, standardize,
                         elbo, iter, null_weight, track_fit, tracking)

  # Run refinement if requested
  if (refine && !is.null(model$sets) && length(model$sets$cs) > 0) {
    model <- run_refine(
      model, data, L, intercept, standardize,
      scaled_prior_variance, residual_variance,
      prior_weights, null_weight, model_init,
      estimate_prior_variance, estimate_prior_method,
      check_null_threshold, estimate_residual_variance,
      estimate_residual_method, residual_variance_lowerbound,
      residual_variance_upperbound, max_iter, tol,
      verbose, track_fit, coverage, min_abs_corr,
      prior_tol, n_purity, compute_univariate_zscore,
      check_prior, convergence_method, alpha0, beta0
    )
  }

  return(model)
}
