# The main workhorse function

susie_engine = function(data,
                        L, intercept = TRUE, standardize = TRUE,
                        scaled_prior_variance, residual_variance,
                        prior_weights, null_weight, model_init = NULL,
                        estimate_prior_variance,
                        estimate_prior_method = "optim",
                        check_null_threshold,
                        estimate_residual_variance,
                        residual_variance_lowerbound,
                        residual_variance_upperbound,
                        max_iter, tol, verbose, track_fit,
                        coverage, min_abs_corr,
                        prior_tol, n_purity, compute_univariate_zscore = FALSE,
                        check_prior = FALSE){

  # Prior Variance Method
  estimate_prior_method <- match.arg(estimate_prior_method)

  # Initialize Model Object
  model <- ibss_initialize(data                  = data,
                           L                     = L,
                           scaled_prior_variance = scaled_prior_variance,
                           residual_variance     = residual_variance,
                           prior_weights         = prior_weights,
                           null_weight           = null_weight,
                           model_init            = model_init)

  # Initialize ELBO & Tracking
  elbo <- rep(as.numeric(NA), max_iter + 1)
  elbo[1] <- -Inf
  tracking <- vector("list", max_iter)

  # IBSS Loop
  for(iter in seq_len(max_iter)){

    # Save alpha, prior variance, and residual variance if tracking
    # if (track_fit) tracking[[iter]] <- susie_slim(model)
    # TODO: Needs to adjust util function. It calls res$...

    # Update all L effects
    model <- ibss_fit(data, model,
                      estimate_prior_variance = estimate_prior_variance,
                      estimate_prior_method   = estimate_prior_method,
                      check_null_threshold    = check_null_threshold)

    # Validate prior variance is reasonable
    validate_prior(data, model, check_prior)

    # SS has a force iterate option. It appears this is used when IBSS is
    # still resolving zR discrepency. Since we no longer include that, seems
    # like we don't need.

    # Get Objective & Check Convergence
    elbo[iter + 1] <- get_objective(data, model)

    if (elbo[iter + 1] - elbo[iter] < tol) {
      model$converged <- TRUE
      break
    }

    # Update Residual Variance
    if (estimate_residual_variance) {
      model$sigma2 <- max(residual_variance_lowerbound,
                          est_residual_variance(data, model))
      model$sigma2 <- min(model$sigma2, residual_variance_upperbound)
    }
  }

  if (is.null(model$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!"))
    model$converged = FALSE
  }

  model <- ibss_finalize(data, model,
                         coverage                  = coverage,
                         min_abs_corr              = min_abs_corr,
                         median_abs_corr           = median_abs_corr,
                         prior_tol                 = prior_tol,
                         n_purity                  = n_purity,
                         compute_univariate_zscore = compute_univariate_zscore,
                         intercept                 = intercept,
                         standardize               = standardize,
                         elbo                      = elbo,
                         iter                      = iter,
                         null_weight               = null_weight,
                         track_fit                 = track_fit,
                         tracking                  = tracking)

  return(model)
}
