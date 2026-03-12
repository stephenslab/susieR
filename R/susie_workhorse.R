#' SuSiE workhorse function
#'
#' Main orchestration function for the IBSS algorithm. Coordinates the complete
#' fitting pipeline from initialization through iteration to finalization.
#'
#' @param data Data object (individual, ss, or rss_lambda)
#' @param params Validated params object
#'
#' @return Complete fitted SuSiE model object with credible sets, PIPs, and
#'   other results.
#'
#' @export
#' @keywords internal
susie_workhorse <- function(data, params) {

  # Initialize model object
  model <- ibss_initialize(data, params)

  # Initialize ELBO & tracking
  elbo     <- rep(as.numeric(NA), params$max_iter + 1)
  elbo[1]  <- -Inf
  tracking <- list()

  # Initialize runtime state (convergence tracking, cleaned up at finalization)
  model$runtime <- list(
    prev_elbo     = -Inf,
    prev_alpha    = model$alpha,
    prev_pip_diff = NULL
  )

  # Main IBSS iteration loop
  for (iter in seq_len(params$max_iter)) {
    # Store iteration snapshot for track_fit
    tracking <- track_ibss_fit(data, params, model, tracking, iter, elbo)

    # Update all L effects
    model <- ibss_fit(data, params, model)

    # Calculate objective and check convergence
    elbo[iter + 1] <- get_objective(data, params, model)
    model <- check_convergence(data, params, model, elbo, iter)

    # Update convergence state for next iteration
    model$runtime$prev_elbo  <- elbo[iter + 1]
    model$runtime$prev_alpha <- model$alpha

    if (model$converged) {
      break
    }

    # Update variance components if not converged.
    # The method itself checks params to decide what to update,
    # allowing S3 overrides to update additional model parameters
    model <- update_model_variance(data, params, model)

  }

  # Check final convergence status
  if (!model$converged) {
    warning_message(paste("IBSS algorithm did not converge in", params$max_iter, "iterations!"))
  }

  # Set ELBO from iterations
  model$elbo <- elbo[2:(iter + 1)]

  # For Servin-Stephens (NIG prior), scale prior variance by residual
  # variance mode
  if (params$use_servin_stephens)
    model$V <- model$V * model$rv
    
  # Zero out effects with negligible prior variance
  model <- trim_null_effects(data, params, model)

  model <- ibss_finalize(data, params, model, elbo, iter, tracking)

  # Run refinement if requested
  if (params$refine && !is.null(model$sets) && length(model$sets$cs) > 0) {
    model <- run_refine(model, data, params)
  }

  return(model)
}
