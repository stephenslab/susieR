# Function to run refinement iterations
#' @keywords internal
run_refine <- function(model, data, params) {

  # Warn if model_init was provided initially
  if (!is.null(params$model_init)) {
    warning_message("The given model_init is not used in refinement")
  }

  if (params$verbose) {
    message("Starting refinement process...")
    message("Initial ELBO:", susie_get_objective(model))
    message("Number of credible sets to refine:", length(model$sets$cs))
  }

  # Extract prior weights for refinement
  pw_s <- extract_prior_weights(model)

  # Main refinement loop
  conti <- TRUE
  refine_iteration <- 0

  while (conti && !is.null(model$sets) && length(model$sets$cs) > 0) {
    refine_iteration <- refine_iteration + 1
    if (params$verbose) {
      message("\nRefinement iteration ", refine_iteration)
    }
    candidate_models <- list()

    for (cs_idx in seq_along(model$sets$cs)) {
      # Create modified prior weights (zero out CS variables)
      pw_cs <- pw_s
      pw_cs[model$sets$cs[[cs_idx]]] <- 0

      # Skip if all weights are zero
      if (all(pw_cs == 0)) {
        if (params$verbose) {
          message("  Skipping CS ", cs_idx, " - all prior weights would be zero")
        }
        break
      }

      if (params$verbose) {
        message("  Refining CS ", cs_idx, " with variables: ",
                paste(model$sets$cs[[cs_idx]], collapse = " "))
      }

      # Step 1: Fit with modified weights, no initialization
      params_step1                <- params
      params_step1$prior_weights  <- reconstruct_full_weights(pw_cs, model$null_weight)
      params_step1$null_weight    <- model$null_weight
      params_step1$model_init     <- NULL
      params_step1$verbose        <- FALSE
      params_step1$track_fit      <- FALSE
      params_step1$refine         <- FALSE

      model_step1 <- susie_workhorse(data, params_step1)

      # Step 3: Extract initialization from step1
      init_from_step1 <- list(
        alpha = model_step1$alpha,
        mu    = model_step1$mu,
        mu2   = model_step1$mu2
      )
      class(init_from_step1) <- "susie"

      # Step 2: Fit with original weights using initialization
      params_step2                <- params
      params_step2$prior_weights  <- reconstruct_full_weights(pw_s, model$null_weight)
      params_step2$null_weight    <- model$null_weight
      params_step2$model_init     <- init_from_step1
      params_step2$verbose        <- FALSE
      params_step2$track_fit      <- FALSE
      params_step2$refine         <- FALSE

      model_step2 <- susie_workhorse(data, params_step2)

      candidate_models <- c(candidate_models, list(model_step2))

      if (params$verbose) {
        message("    ELBO for refined CS ", cs_idx, ": ", susie_get_objective(model_step2))
      }
    }

    # Evaluate candidates
    if (length(candidate_models) == 0) {
      conti <- FALSE
      if (params$verbose) {
        message("No candidate models generated, stopping refinement")
      }
    } else {
      elbos <- sapply(candidate_models, susie_get_objective)
      current_elbo <- susie_get_objective(model)

      if (params$verbose) {
        message("\nCurrent ELBO: ", current_elbo)
        message("Candidate ELBOs: ", paste(elbos, collapse = " "))
        message("Best candidate ELBO: ", max(elbos))
      }

      if (max(elbos) - current_elbo > params$tol) {
        # Found improvement beyond tolerance
        best_idx <- which.max(elbos)
        model <- candidate_models[[best_idx]]
        if (params$verbose) {
          message("ELBO improved! Selected model from CS ", best_idx, " refinement")
          message("Improvement: ", max(elbos) - current_elbo)
        }
      } else {
        # No improvement beyond tolerance, stop
        if (params$verbose) {
          message("No improvement in ELBO beyond tolerance (", params$tol, "), stopping refinement")
        }
        conti <- FALSE
      }
    }
  }

  if (params$verbose) {
    message("\nRefinement complete!")
    message("Final ELBO: ", susie_get_objective(model))
    message("Total refinement iterations: ", refine_iteration)
  }

  return(model)
}

