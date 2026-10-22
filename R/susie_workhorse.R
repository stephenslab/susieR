#' SuSiE workhorse function
#'
#' Main orchestration for the IBSS algorithm. When `params$L_greedy`
#' is non-NULL, runs a greedy outer loop that grows `L` in linear
#' steps of `params$L_greedy` until the fit has at least one empty
#' slot (`min(lbf) < params$greedy_lbf_cutoff`, default `0.1`) or `L` reaches
#' `params$L`. With `params$L_greedy = NULL` (default), runs a
#' single fixed-`L` IBSS, output bit-identical to prior susieR.
#'
#' @param data Data object (individual, ss, or rss_lambda).
#' @param params Validated params object.
#' @return Complete fitted SuSiE model.
#'
#' @export
#' @keywords internal
susie_workhorse <- function(data, params) {

  # Greedy-L outer loop. Saturation detected when any one slot's
  # lbf falls below greedy_lbf_cutoff (slot-invariant, single-round verdict).
  # Warm-start across rounds via params$model_init.
  if (!is.null(params$L_greedy)) {
    L_max   <- params$L
    L_step  <- params$L_greedy
    greedy_lbf_cutoff <- if (is.null(params$greedy_lbf_cutoff)) 0.1 else params$greedy_lbf_cutoff
    verbose <- isTRUE(params$verbose)
    history <- list()

    current_L <- min(L_step, L_max)
    fit       <- NULL
    round_n   <- 0L

    repeat {
      round_n <- round_n + 1L
      params_round          <- params
      params_round$L_greedy <- NULL                # avoid recursion
      params_round$L        <- current_L
      if (!is.null(fit)) params_round$model_init <- fit
      fit <- susie_workhorse(data, params_round)

      min_lbf <- min(fit$lbf, na.rm = TRUE)
      action  <- if (current_L >= L_max)    "L_max reached"
                 else if (min_lbf < greedy_lbf_cutoff) "saturated"
                 else                         "grow"
      history[[round_n]] <- list(L = current_L, min_lbf = min_lbf,
                                 action = action)
      if (action != "grow") break
      current_L <- min(current_L + L_step, L_max)
    }
    if (verbose) {
      message(sprintf("[L_greedy] %d round%s, greedy_lbf_cutoff=%.3f, final L=%d",
                      round_n, if (round_n == 1L) "" else "s",
                      greedy_lbf_cutoff, current_L))
      message(sprintf("%-6s %-5s %-10s %s",
                      "round", "L", "min(lbf)", "action"))
      for (i in seq_along(history)) {
        h <- history[[i]]
        message(sprintf("%-6d %-5d %-10.3f %s",
                        i, h$L, h$min_lbf, h$action))
      }
    }
    return(fit)
  }

  # Initialize model object
  model <- ibss_initialize(data, params)

  # Initialize ELBO & tracking
  elbo     <- rep(as.numeric(NA), params$max_iter + 1)
  elbo[1]  <- -Inf
  tracking <- list()

  # Initialize runtime state (convergence tracking, cleaned up at finalization)
  model$runtime <- list(
    prev_elbo  = -Inf,
    prev_alpha = model$alpha
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

  # SuSiE-ash: final Mr.ASH pass on the sparse fit the loop ended with.
  if (isTRUE(params$unmappable_effects %in% c("ash", "ash_filter_archived")))
    model <- run_final_ash_pass(data, params, model)

  # Check final convergence status. Emit both the cosmetic colored
  # message and a real R warning() so callers can `suppressWarnings()`
  # or `expect_warning()` on it (mirrors the Q_art artifact pattern in
  # rss_mismatch.R).
  if (!model$converged) {
    msg <- paste("IBSS algorithm did not converge in", params$max_iter, "iterations!")
    warning_message(msg)
    warning(msg, call. = FALSE)
  }

  # Set ELBO from iterations
  model$elbo <- elbo[2:(iter + 1)]

  # For NIG prior, scale prior variance by residual variance mode.
  if (isTRUE(params$use_NIG))
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
