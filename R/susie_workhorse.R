susie_workhorse <- function(data, params) {

  # Validate method argument
  estimate_prior_method <- match.arg(params$estimate_prior_method, c("optim", "EM", "simple"))

  # Initialize model object
  model <- ibss_initialize(data, params)

  # Initialize tracking
  elbo <- rep(as.numeric(NA), params$max_iter + 1)
  elbo[1] <- -Inf
  tracking <- list()

  # Main IBSS iteration loop
  for (iter in seq_len(params$max_iter)) {
    cat("=== NEW SUSIE2.0: ITER", iter, "===\n")
    cat("Before ibss_fit: sigma2 =", model$sigma2, "\n")
    cat("Before ibss_fit: V =", paste(model$V, collapse=", "), "\n")
    if (exists("XtXr", model)) {
      cat("Before ibss_fit: XtXr[1:3] =", paste(model$XtXr[1:3], collapse=", "), "\n")
    }
    
    # Track iteration progress
    tracking <- track_ibss_fit(data, params, model, tracking, iter, params$track_fit)

    # Store previous model parameters for convergence check
    model$prev_elbo  <- elbo[iter]
    model$prev_alpha <- model$alpha

    # Update all L effects
    model <- ibss_fit(data, params, model)
    
    cat("After ibss_fit: V =", paste(model$V, collapse=", "), "\n")
    cat("After ibss_fit: alpha[1,1:3] =", paste(model$alpha[1,1:3], collapse=", "), "\n")
    cat("After ibss_fit: mu[1,1:3] =", paste(model$mu[1,1:3], collapse=", "), "\n")
    if (exists("XtXr", model)) {
      cat("After ibss_fit: XtXr[1:3] =", paste(model$XtXr[1:3], collapse=", "), "\n")
    }

    # Calculate objective for tracking
    elbo[iter + 1] <- get_objective(data, model, verbose = params$verbose)
    cat("ELBO =", elbo[iter + 1], "\n")

    # Check for convergence
    converged <- check_convergence(model, elbo, params$tol, params$convergence_method, iter)

    if (converged) {
      model$converged <- TRUE
      break
    }

    # Update variance components if not converged and estimation is requested
    if (params$estimate_residual_variance) {
      model <- update_model_variance(data, params, model, params$residual_variance_lowerbound,
                                     params$residual_variance_upperbound, params$estimate_residual_method)
      cat("After variance update: sigma2 =", model$sigma2, "\n")
    }

  }

  # Check final convergence status
  if (is.null(model$converged)) {
    warning_message(paste("IBSS algorithm did not converge in", params$max_iter, "iterations!"))
    model$converged <- FALSE
  }

  # Set ELBO from iterations
  model$elbo <- elbo[2:(iter + 1)]

  model <- ibss_finalize(data, params, model, elbo, iter, tracking)

  # Run refinement if requested
  if (params$refine && !is.null(model$sets) && length(model$sets$cs) > 0) {
    model <- run_refine(model, data, params)
  }

  return(model)
}
