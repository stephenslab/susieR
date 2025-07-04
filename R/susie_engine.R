# The main workhorse function

susie_engine = function(data,
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
                        refine = FALSE){

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
  tracking <- list()

  # IBSS Loop
  for(iter in seq_len(max_iter)){

    # Save core values for each iteration

    # This could be absorbed into ibss_fit, but i dont believe it should. I think
    # it makes sense to absorb the various checks but when we are adding additional
    # features to the model output it should remain more transparent.
    tracking <- susie_extract_core(data, model, tracking, iter, track_fit)

    # Store previous model for convergence check
    model_prev <- model

    # Update all L effects
    model <- ibss_fit(data, model,
                      estimate_prior_variance = estimate_prior_variance,
                      estimate_prior_method   = estimate_prior_method,
                      check_null_threshold    = check_null_threshold,
                      check_prior             = check_prior)

    # Get Objective (for tracking)
    elbo[iter + 1] <- get_objective(data, model)

    # Handle convergence and variance updates based on method type
    if (update_variance_before_convergence(data)) {
      # Non-sparse: Update variance first, then check convergence
      if (estimate_residual_variance) {
        variance_result <- update_variance_components(data, model)
        model$sigma2 <- max(residual_variance_lowerbound, variance_result$sigma2)
        model$sigma2 <- min(model$sigma2, residual_variance_upperbound)
        
        # Update additional variance components if they exist
        if (!is.null(variance_result$tausq)) {
          model$tausq <- variance_result$tausq
        }
        
        # Update derived quantities after variance component changes
        data <- update_derived_quantities(data, model)
        
        # Transfer theta from data to model if computed (for non-sparse methods)
        if (!is.null(data$theta)) {
          model$theta <- data$theta
          
          # Update fitted values to include theta: XtXr = XtX %*% (b + theta)
          b <- colSums(model$alpha * model$mu)
          model$XtXr <- data$XtX %*% (b + model$theta)
        }
      }
      
      # Check convergence after variance update
      if (check_convergence(data, model_prev, model, elbo[iter], elbo[iter + 1], tol)) {
        model$converged <- TRUE
        break
      }
    } else {
      # Standard: Check convergence first, then update variance
      if (check_convergence(data, model_prev, model, elbo[iter], elbo[iter + 1], tol)) {
        model$converged <- TRUE
        break
      }
      
      # Update variance components after convergence check
      if (estimate_residual_variance) {
        variance_result <- update_variance_components(data, model)
        model$sigma2 <- max(residual_variance_lowerbound, variance_result$sigma2)
        model$sigma2 <- min(model$sigma2, residual_variance_upperbound)
        
        # Update additional variance components if they exist
        if (!is.null(variance_result$tausq)) {
          model$tausq <- variance_result$tausq
        }
        
        # Update derived quantities after variance component changes
        data <- update_derived_quantities(data, model)
        
        # Transfer theta from data to model if computed (for non-sparse methods)
        if (!is.null(data$theta)) {
          model$theta <- data$theta
          
          # Update fitted values to include theta: XtXr = XtX %*% (b + theta)
          b <- colSums(model$alpha * model$mu)
          model$XtXr <- data$XtX %*% (b + model$theta)
        }
      }
    }
  }

  # Update non-sparse
  #mr.ash()....

  if (is.null(model$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!"))
    model$converged = FALSE
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
