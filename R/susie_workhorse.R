# =============================================================================
# SUSIE WORKHORSE
#
# Main orchestration function for the IBSS algorithm. Coordinates the complete
# fitting pipeline from initialization through iteration to finalization.
# =============================================================================
#'
#' @param data Data object (individual, ss, or rss_lambda)
#' @param params Validated params object
#'
#' @return Complete fitted SuSiE model object with credible sets, PIPs, and
#'   other results.
#'
#' @keywords internal
#' @noRd
susie_workhorse <- function(data, params) {

  # Initialize model object
  model <- ibss_initialize(data, params)

  # Initialize ELBO & tracking
  elbo     <- rep(as.numeric(NA), params$max_iter + 1)
  elbo[1]  <- -Inf
  tracking <- list()

  # Main IBSS iteration loop
  for (iter in seq_len(params$max_iter)) {
    # Track iteration progress and store previous values for convergence checking
    tracking <- track_ibss_fit(data, params, model, tracking, iter, elbo)

    # Update all L effects
    model <- ibss_fit(data, params, model)

    # Calculate objective and check convergence
    elbo[iter + 1] <- get_objective(data, params, model)
    converged      <- check_convergence(params, model, elbo, iter, tracking)

    if (converged) {
      model$converged <- TRUE
      break
    }

    # Update variance components if not converged and estimation is requested
    if (params$estimate_residual_variance) {
      model <- update_model_variance(data, params, model)
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
