#' Block coordinate ascent for iterative model refinement.
#'
#' Generic framework for post-convergence refinement of fitted models.
#' Repeatedly applies a block update (\code{step_fn}) to a fitted model,
#' monitoring ELBO for convergence and stability.  Both CS refinement
#' (\code{\link{run_refine}}) and residual variance estimation (mvsusieR)
#' are instances of block coordinate ascent over different parameter blocks.
#'
#' Convergence is reached when either:
#' \itemize{
#'   \item \code{step_fn} returns \code{converged = TRUE}
#'         (the block update signals no further improvement), or
#'   \item the relative ELBO change falls below \code{tol}
#'         (ELBO stabilized across block updates).
#' }
#'
#' A warning is issued if the ELBO decreases between iterations.
#'
#' @param model Fitted model (e.g., from \code{susie_workhorse} or
#'   \code{mvsusie_workhorse}).
#' @param data Data object passed to \code{step_fn}.
#' @param step_fn A function with signature
#'   \code{function(model, data, iter)} that performs one block coordinate
#'   update.  Must return a named list with elements:
#'   \describe{
#'     \item{model}{(required) The updated model.}
#'     \item{data}{(optional) Updated data object, e.g. after changing
#'       residual variance.  If \code{NULL}, the data is unchanged.}
#'     \item{converged}{(optional) Logical; if \code{TRUE}, stop
#'       iterating.}
#'     \item{log_msg}{(optional) Character string appended to verbose
#'       output.}
#'   }
#' @param max_iter Maximum number of block ascent iterations
#'   (default 100).
#' @param tol Convergence tolerance for relative ELBO change
#'   (default 1e-3).
#' @param verbose If \code{TRUE}, print progress each iteration
#'   (default \code{FALSE}).
#'
#' @return The refined model, with \code{model$converged} set to
#'   \code{TRUE} or \code{FALSE}.
#'
#' @export
block_coordinate_ascent <- function(model, data, step_fn,
                                     max_iter = 100, tol = 1e-3,
                                     verbose = FALSE) {
  prev_elbo  <- susie_get_objective(model)
  prev_model <- model
  converged  <- FALSE

  for (iter in seq_len(max_iter)) {
    result <- step_fn(model, data, iter)
    model <- result$model
    if (!is.null(result$data)) data <- result$data

    current_elbo <- susie_get_objective(model)
    elbo_change <- current_elbo - prev_elbo

    # If ELBO decreased, the step did not improve the objective.
    # Revert to the previous model and treat as converged.
    if (elbo_change < 0) {
      if (verbose)
        message(sprintf(
          "Block ascent iter %d: update did not improve ELBO (change=%.4g); ",
          iter, elbo_change),
          "rejecting update and stopping.")
      model <- prev_model
      converged <- TRUE
      break
    }

    if (verbose) {
      msg <- sprintf("Block ascent iter %d: ELBO=%.4f, change=%.4g",
                     iter, current_elbo, elbo_change)
      if (!is.null(result$log_msg))
        msg <- paste0(msg, ", ", result$log_msg)
      message(msg)
    }

    # Convergence: step_fn signals done, or ELBO stabilized
    if (isTRUE(result$converged)) {
      converged <- TRUE
      break
    }
    elbo_delta <- abs(elbo_change) / max(1, abs(current_elbo))
    if (elbo_delta < tol) {
      converged <- TRUE
      break
    }

    prev_elbo  <- current_elbo
    prev_model <- model
  }

  model$converged <- converged
  if (!converged)
    warning_message("Block coordinate ascent did not converge in ",
                    max_iter, " iterations")
  return(model)
}


# Credible set refinement via block coordinate ascent.
#
# For each credible set, perturbs prior weights (zeroing out the CS
# variables) and re-fits via a two-step procedure:
#   Step 1: fit with zeroed CS weights (explores alternative signals)
#   Step 2: re-fit with original weights, warm-started from Step 1
# The best candidate (highest ELBO) is accepted if it improves beyond
# tolerance.  This is repeated until no further improvement.
#
# @keywords internal
run_refine <- function(model, data, params) {

  if (!is.null(params$model_init))
    warning_message("The given model_init is not used in refinement")

  pw_s <- extract_prior_weights(model)

  # One block coordinate step: try refining each CS, pick best candidate.
  cs_refine_step <- function(model, data, iter) {
    if (is.null(model$sets) || length(model$sets$cs) == 0)
      return(list(model = model, converged = TRUE))

    candidates <- list()
    for (cs_idx in seq_along(model$sets$cs)) {
      # Zero out prior weights for variables in this CS
      pw_cs <- pw_s
      pw_cs[model$sets$cs[[cs_idx]]] <- 0
      if (all(pw_cs == 0)) break

      # Step 1: fit with zeroed CS weights (no initialization)
      p1 <- params
      p1$prior_weights <- reconstruct_full_weights(pw_cs, model$null_weight)
      p1$null_weight   <- model$null_weight
      p1$model_init    <- NULL
      p1$verbose       <- FALSE
      p1$track_fit     <- FALSE
      p1$refine        <- FALSE
      m1 <- susie_workhorse(data, p1)

      # Step 2: re-fit with original weights, warm-started from Step 1
      init <- list(alpha = m1$alpha, mu = m1$mu, mu2 = m1$mu2)
      class(init) <- "susie"
      p2 <- params
      p2$prior_weights <- reconstruct_full_weights(pw_s, model$null_weight)
      p2$null_weight   <- model$null_weight
      p2$model_init    <- init
      p2$verbose       <- FALSE
      p2$track_fit     <- FALSE
      p2$refine        <- FALSE
      candidates <- c(candidates, list(susie_workhorse(data, p2)))
    }

    if (length(candidates) == 0)
      return(list(model = model, converged = TRUE))

    elbos <- sapply(candidates, susie_get_objective)
    current_elbo <- susie_get_objective(model)

    if (max(elbos) - current_elbo > params$tol) {
      # Accept best candidate
      list(model = candidates[[which.max(elbos)]])
    } else {
      # No improvement beyond tolerance — converged
      list(model = model, converged = TRUE)
    }
  }

  block_coordinate_ascent(model, data, cs_refine_step,
                          max_iter = 100, tol = params$tol,
                          verbose = params$verbose)
}
