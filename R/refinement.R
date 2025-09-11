# Function to run refinement iterations
#' @keywords internal
run_refine <- function(model, data, L, intercept, standardize,
                       scaled_prior_variance, residual_variance,
                       prior_weights, null_weight, model_init,
                       estimate_prior_variance, estimate_prior_method,
                       check_null_threshold, estimate_residual_variance,
                       estimate_residual_method, residual_variance_lowerbound,
                       residual_variance_upperbound, max_iter, tol,
                       verbose, track_fit, coverage, min_abs_corr,
                       prior_tol, n_purity, compute_univariate_zscore,
                       check_prior, convergence_method, alpha0, beta0) {

  # Warn if model_init was provided initially
  if (!is.null(model_init)) {
    warning_message("The given model_init is not used in refinement")
  }

  if (verbose) {
    message("Starting refinement process...")
    message("Initial ELBO:", susie_get_objective(model))
    message("Number of credible sets to refine:", length(model$sets$cs))
  }

  # Extract prior weights for refinement
  pw_s <- extract_prior_weights(model, data$null_weight)

  # Main refinement loop
  conti <- TRUE
  refine_iteration <- 0

  while (conti && !is.null(model$sets) && length(model$sets$cs) > 0) {
    refine_iteration <- refine_iteration + 1
    if (verbose) {
      message("\nRefinement iteration ", refine_iteration)
    }
    candidate_models <- list()

    for (cs_idx in seq_along(model$sets$cs)) {
      # Create modified prior weights (zero out CS variables)
      pw_cs <- pw_s
      pw_cs[model$sets$cs[[cs_idx]]] <- 0

      # Skip if all weights are zero
      if (all(pw_cs == 0)) {
        if (verbose) {
          message("  Skipping CS ", cs_idx, " - all prior weights would be zero")
        }
        break
      }

      if (verbose) {
        message("  Refining CS ", cs_idx, " with variables: ", paste(model$sets$cs[[cs_idx]], collapse = " "))
      }

      # Step 1: Create data with modified prior weights
      data_modified <- modify_prior_weights(data, pw_cs)

      # Step 2: Fit with modified weights, no initialization
      model_step1 <- susie_workhorse(
        data = data_modified,
        L = L,
        intercept = intercept,
        standardize = standardize,
        scaled_prior_variance = scaled_prior_variance,
        residual_variance = residual_variance,
        prior_weights = data_modified$prior_weights,
        null_weight = data_modified$null_weight,
        model_init = NULL,
        estimate_prior_variance = estimate_prior_variance,
        estimate_prior_method = estimate_prior_method,
        check_null_threshold = check_null_threshold,
        estimate_residual_variance = estimate_residual_variance,
        estimate_residual_method = estimate_residual_method,
        residual_variance_lowerbound = residual_variance_lowerbound,
        residual_variance_upperbound = residual_variance_upperbound,
        max_iter = max_iter,
        tol = tol,
        verbose = FALSE,
        track_fit = FALSE,
        coverage = coverage,
        min_abs_corr = min_abs_corr,
        prior_tol = prior_tol,
        n_purity = n_purity,
        compute_univariate_zscore = compute_univariate_zscore,
        check_prior = check_prior,
        convergence_method = convergence_method,
        alpha0 = alpha0,
        beta0 = beta0,
        refine = FALSE
      )

      # Step 3: Extract initialization from step1
      init_from_step1 <- list(
        alpha = model_step1$alpha,
        mu = model_step1$mu,
        mu2 = model_step1$mu2
      )
      class(init_from_step1) <- "susie"

      # Step 4: Create data with original prior weights
      data_original <- modify_prior_weights(data, pw_s)

      # Step 5: Fit with original weights using initialization
      model_step2 <- susie_workhorse(
        data = data_original,
        L = L,
        intercept = intercept,
        standardize = standardize,
        scaled_prior_variance = scaled_prior_variance,
        residual_variance = residual_variance,
        prior_weights = data_original$prior_weights,
        null_weight = data_original$null_weight,
        model_init = init_from_step1,
        estimate_prior_variance = estimate_prior_variance,
        estimate_prior_method = estimate_prior_method,
        check_null_threshold = check_null_threshold,
        estimate_residual_variance = estimate_residual_variance,
        estimate_residual_method = estimate_residual_method,
        residual_variance_lowerbound = residual_variance_lowerbound,
        residual_variance_upperbound = residual_variance_upperbound,
        max_iter = max_iter,
        tol = tol,
        verbose = FALSE,
        track_fit = FALSE,
        coverage = coverage,
        min_abs_corr = min_abs_corr,
        prior_tol = prior_tol,
        n_purity = n_purity,
        compute_univariate_zscore = compute_univariate_zscore,
        check_prior = check_prior,
        convergence_method = convergence_method,
        alpha0 = alpha0,
        beta0 = beta0,
        refine = FALSE
      )

      candidate_models <- c(candidate_models, list(model_step2))

      if (verbose) {
        message("    ELBO for refined CS ", cs_idx, ": ", susie_get_objective(model_step2))
      }
    }

    # Evaluate candidates
    if (length(candidate_models) == 0) {
      conti <- FALSE
      if (verbose) {
        message("No candidate models generated, stopping refinement")
      }
    } else {
      elbos <- sapply(candidate_models, susie_get_objective)
      current_elbo <- susie_get_objective(model)

      if (verbose) {
        message("\nCurrent ELBO: ", current_elbo)
        message("Candidate ELBOs: ", paste(elbos, collapse = " "))
        message("Best candidate ELBO: ", max(elbos))
      }

      if (max(elbos) - current_elbo > tol) {
        # Found improvement beyond tolerance
        best_idx <- which.max(elbos)
        model <- candidate_models[[best_idx]]
        if (verbose) {
          message("ELBO improved! Selected model from CS ", best_idx, " refinement")
          message("Improvement: ", max(elbos) - current_elbo)
        }
      } else {
        # No improvement beyond tolerance, stop
        if (verbose) {
          message("No improvement in ELBO beyond tolerance (", tol, "), stopping refinement")
        }
        conti <- FALSE
      }
    }
  }

  if (verbose) {
    message("\nRefinement complete!")
    message("Final ELBO: ", susie_get_objective(model))
    message("Total refinement iterations: ", refine_iteration)
  }

  return(model)
}

