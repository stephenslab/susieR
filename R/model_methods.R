# =============================================================================
# MODEL-LEVEL S3 METHODS

# S3 generics dispatched on model class (model field access, initialization, 
# convergence, ELBO)
# =============================================================================

#' Get prior variance for effect l
#' @keywords internal
get_prior_variance_l <- function(model, l) {
  UseMethod("get_prior_variance_l")
}
#' @keywords internal
get_prior_variance_l.default <- function(model, l) {
  model$V[l]
}

#' Set prior variance for effect l
#' @keywords internal
set_prior_variance_l <- function(model, l, V) {
  UseMethod("set_prior_variance_l")
}
#' @keywords internal
set_prior_variance_l.default <- function(model, l, V) {
  model$V[l] <- V
  model
}

#' Get posterior inclusion probabilities for effect l
#' @keywords internal
get_alpha_l <- function(model, l) {
  UseMethod("get_alpha_l")
}
#' @keywords internal
get_alpha_l.default <- function(model, l) {
  model$alpha[l, ]
}

#' Get posterior moments for effect l (for EM prior variance update)
#' @keywords internal
get_posterior_moments_l <- function(model, l) {
  UseMethod("get_posterior_moments_l")
}
#' @keywords internal
get_posterior_moments_l.default <- function(model, l) {
  list(post_mean = model$mu[l, ], post_mean2 = model$mu2[l, ])
}

#' Get PIP-weighted posterior mean for effect l (alpha * mu)
#' @keywords internal
get_posterior_mean_l <- function(model, l) {
  UseMethod("get_posterior_mean_l")
}
#' @keywords internal
get_posterior_mean_l.default <- function(model, l) {
  model$alpha[l, ] * model$mu[l, ]
}

#' Get sum of PIP-weighted posterior means across all effects
#' @keywords internal
get_posterior_mean_sum <- function(model) {
  UseMethod("get_posterior_mean_sum")
}
#' @keywords internal
get_posterior_mean_sum.default <- function(model) {
  colSums(model$alpha * model$mu)
}

# =============================================================================
# MODEL INITIALIZATION
#
# Initialize core model matrices and parameter storage.
# =============================================================================

#' @keywords internal
initialize_matrices <- function(data, params, var_y) {
  UseMethod("initialize_matrices")
}

#' @keywords internal
initialize_matrices.default <- function(data, params, var_y) {
  L <- params$L
  # NIG uses the scaled slab variance s0 = sigma0^2 / sigma^2 during fitting.
  V <- expand_scaled_prior_variance(
    params$scaled_prior_variance, ifelse(isTRUE(params$use_NIG), 1, var_y), L)
  mat_init <- list(
    alpha             = matrix(1 / data$p, L, data$p),
    mu                = matrix(0, L, data$p),
    mu2               = matrix(0, L, data$p),
    V                 = V,
    KL                = rep(as.numeric(NA), L),
    lbf               = rep(as.numeric(NA), L),
    lbf_variable      = matrix(as.numeric(NA), L, data$p),
    sigma2            = params$residual_variance,
    pi                = params$prior_weights,
    null_weight       = params$null_weight,
    predictor_weights = rep(as.numeric(NA), data$p)
  )

  return(mat_init)
}

# =============================================================================
# VARIANCE UPDATE
#
# Update residual variance (and possibly other variance components) after
# each IBSS iteration.
# =============================================================================

#' @keywords internal
#' @importFrom utils modifyList
update_model_variance <- function(data, params, model) {
  UseMethod("update_model_variance")
}

#' @keywords internal
update_model_variance.default <- function(data, params, model) {
  if (!isTRUE(params$estimate_residual_variance)) return(model)
  # Update variance components
  variance_result <- update_variance_components(data, params, model)
  model           <- modifyList(model, variance_result)

  # Apply bounds to residual variance
  model$sigma2    <- min(max(model$sigma2, params$residual_variance_lowerbound),
                         params$residual_variance_upperbound)

  # Update derived quantities after variance component changes
  model           <- update_derived_quantities(data, params, model)

  return(model)
}

# =============================================================================
# CONVERGENCE CHECKING
# =============================================================================

#' @keywords internal
check_convergence <- function(data, params, model, elbo, iter) {
  UseMethod("check_convergence")
}

#' Format the per-iter sigma2 cell for verbose output
#'
#' Default returns the scalar sigma2 in `%.4f`. Subclasses
#' (e.g., mfsusieR's list-of-vectors sigma2; mvsusieR's
#' matrix sigma2) override to a compact summary string of
#' fixed width.
#' @keywords internal
format_sigma2_summary <- function(model) {
  UseMethod("format_sigma2_summary")
}
#' @keywords internal
format_sigma2_summary.default <- function(model) {
  sprintf("%.4f", model$sigma2)
}

#' Append class-specific extra-diag columns to the verbose row
#'
#' Default returns an empty string. Subclasses override to inject
#' columns such as `max_pi_null`, `max_KL_l`, alpha-entropy
#' n_eff. Output is appended after the V column in the per-iter
#' tabular line.
#' @keywords internal
format_extra_diag <- function(model) {
  UseMethod("format_extra_diag")
}
#' @keywords internal
format_extra_diag.default <- function(model) {
  if (is.null(model$lambda_bias))
    return("")
  lambda_infl <- model$lambda_bias
  # Zero-masking of small finite values happens at source in
  # estimate_lambda_bias; we just sanitize non-finite for display.
  lambda_infl[!is.finite(lambda_infl)] <- 0
  if (length(lambda_infl) != 1)
    stop("lambda_bias must be a scalar on the SS path.")
  lb <- paste0("lambda_infl=", format(lambda_infl, digits = 2,
                                      scientific = TRUE))
  if (!is.null(model$B_corrected)) {
    B_corrected <- model$B_corrected
    if (length(unique(B_corrected[is.finite(B_corrected)])) == 1) {
      lb <- paste0(lb, " B_eff=", format(B_corrected[which(is.finite(B_corrected))[1]],
                                            digits = 4, scientific = FALSE))
    }
  }
  lb
}

#' @keywords internal
check_convergence.default <- function(data, params, model, elbo, iter) {
  verbose <- isTRUE(params$verbose)
  V_str <- format_V_summary(model$V)
  chat_str <- format_chat_summary(model)
  sigma2_str <- format_sigma2_summary(model)
  extra_str  <- format_extra_diag(model)

  # Tabular verbose format (ELBO-convergence path): columns are
  # iter, ELBO, delta, sigma2, mem, V (variable-width, last),
  # plus optional class-specific extras after V.
  verbose_row_fmt <- "%4d   %11.4f   %9s   %-9s   %-7s  %s%s%s"
  verbose_header  <- sprintf("%-4s   %11s   %9s   %-9s   %-7s  %s%s",
                             "iter", "ELBO", "delta", "sigma2", "mem", "V",
                             if (nzchar(extra_str)) "  extras" else "")

  # Skip convergence check on first iteration
  if (iter == 1) {
    model$converged <- FALSE
    if (verbose) {
      elbo_val <- elbo[iter + 1]
      if (!is.na(elbo_val) && is.finite(elbo_val)) {
        message(verbose_header)
        message(sprintf(verbose_row_fmt,
                        iter, elbo_val, "-", sigma2_str,
                        sprintf("%.2f GB", mem_used_gb()),
                        paste0(V_str, chat_str),
                        if (nzchar(extra_str)) paste0("  ", extra_str) else "",
                        ""))
      } else {
        message(sprintf("iter %3d: sigma2=%s, V=%s%s [mem: %.2f GB]",
                        iter, sigma2_str, V_str, chat_str, mem_used_gb()))
      }
    }
    return(model)
  }

  # Calculate difference in ELBO values
  ELBO_diff   <- elbo[iter + 1] - model$runtime$prev_elbo
  ELBO_failed <- is.na(ELBO_diff) || is.infinite(ELBO_diff)

  if (params$convergence_method == "pip" || ELBO_failed) {
    if (ELBO_failed && params$convergence_method == "elbo") {
      warning_message(paste0("Iteration ", iter, " produced an NA/infinite ELBO",
                             " value. Using pip-based convergence this iteration."))
    }
    # PIP/alpha convergence. pip_stall_window is reused as the maximum
    # short-cycle lag; it is no longer a "no improvement" stop.
    model <- check_alpha_pip_cycle_convergence(data, params, model)
    pip_diff <- model$runtime$pip_diff
    lambda_diff <- if (!is.null(model$runtime$lambda_bias_diff))
                     model$runtime$lambda_bias_diff else 0
    # Coordinate EB guard: fit_R_mismatch runs after the SER sweep, so a material
    # lambda update must be consumed by one more sweep before convergence.
    if (isTRUE(model$converged) && lambda_diff > params$tol) {
      model$converged <- FALSE
      model$convergence_reason <- paste0("lambda_infl_changed(",
                                         format(lambda_diff, digits = 3,
                                                scientific = TRUE), ")")
    }
    if (verbose) {
      conv_tag <- if (model$converged)
        paste0(" -- converged (", model$convergence_reason, ")")
      else
        ""
      lambda_str <- if (lambda_diff > 0)
        paste0(", max|d(lambda_infl)|=", format(lambda_diff, digits = 3,
                                                scientific = TRUE))
      else ""
      message(sprintf("iter %3d: max|d(alpha,PIP)|=%.2e%s, V=%s%s%s%s [mem: %.2f GB]",
                      iter, pip_diff, lambda_str, V_str, chat_str,
                      if (nzchar(extra_str)) paste0(", ", extra_str) else "",
                      conv_tag, mem_used_gb()))
    }

    if (model$converged && !is.null(params$unmappable_effects) &&
        params$unmappable_effects %in% c("ash", "ash_filter_archived")) {
      model <- run_final_ash_pass(data, params, model)
    }
    return(model)
  }

  # Converge when ELBO stabilizes: small non-negative change.
  # A large negative ELBO_diff means the objective dropped, not convergence.
  if (ELBO_diff < -params$tol) {
    warning_message(sprintf("ELBO decreased by %.2e at iteration %d",
                            -ELBO_diff, iter))
  }
  model$converged <- (ELBO_diff >= 0 && ELBO_diff < params$tol)
  lambda_diff <- if (!is.null(model$runtime$lambda_bias_diff))
                   model$runtime$lambda_bias_diff else 0
  # Coordinate EB guard: fit_R_mismatch runs after the SER sweep, so a material
  # lambda update must be consumed by one more sweep before declaring convergence.
  if (isTRUE(model$converged) && lambda_diff > params$tol)
    model$converged <- FALSE

  if (verbose)
    message(sprintf(verbose_row_fmt,
                    iter, elbo[iter + 1],
                    sprintf("%.2e", ELBO_diff),
                    sigma2_str,
                    sprintf("%.2f GB", mem_used_gb()),
                    paste0(V_str, chat_str),
                    if (nzchar(extra_str)) paste0("  ", extra_str) else "",
                    if (model$converged) "  converged" else ""))

  if (model$converged && !is.null(params$unmappable_effects) &&
      params$unmappable_effects %in% c("ash", "ash_filter_archived")) {
    model <- run_final_ash_pass(data, params, model)
  }
  return(model)
}

# =============================================================================
# OBJECTIVE FUNCTION (ELBO)
# =============================================================================

#' Compute the SuSiE ELBO (evidence lower bound)
#'
#' Building-block function used by downstream packages implementing
#' custom IBSS loops.
#'
#' @param data Data object.
#' @param params Params object.
#' @param model Model object.
#'
#' @return Scalar ELBO value.
#'
#' @export
#' @keywords internal
get_objective <- function(data, params, model) {
  UseMethod("get_objective")
}

#' @export
#' @keywords internal
get_objective.default <- function(data, params, model) {
  if (!is.null(params$unmappable_effects) && params$unmappable_effects == "inf") {
    # Compute omega
    L         <- nrow(model$alpha)
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    omega     <- matrix(0, L, data$p)

    for (l in seq_len(L)) {
      omega[l, ] <- omega_res$diagXtOmegaX + 1 / model$V[l]
    }

    # Compute total ELBO for infinitesimal effects model
    objective <- compute_elbo_inf(
      model$alpha, model$mu, omega, model$lbf,
      model$sigma2, model$tau2, data$n, data$p,
      data$eigen_vectors, data$eigen_values,
      data$VtXty, data$yty,
      eigen_vectors_sq = data$eigen_vectors_sq
    )
  } else if (params$use_NIG && nrow(model$alpha) == 1) {
    objective <- model$marginal_loglik[1]
  } else if (isTRUE(params$use_NIG)) {
    # NIG L>1: KL[l] is gated to 0 (gIBSS has no coherent ELBO); use the
    # proper variational expected log-likelihood.
    objective <- nig_eloglik(data, params, model)
  } else {
    # Standard ELBO computation. `na.rm = TRUE` so subclasses that
    # leave KL[l] = NA on null-effect rows (mfsusieR, mvsusieR) do
    # not need to override get_objective just to skip NAs.
    objective <- Eloglik(data, model) - sum(model$KL, na.rm = TRUE)
  }

  # Add slot prior ELBO terms when c_hat is active.
  # Without these, the ELBO is missing the prior and entropy contributions
  # from the slot activity model.
  if (!is.null(model$c_hat_state)) {
    objective <- objective + slot_prior_elbo(model)
  }

  if (is.infinite(objective)) {
    stop("get_objective() produced an infinite ELBO value")
  }
  return(objective)
}

# =============================================================================
# EFFECT TRIMMING
#
# Zero out effects with negligible prior variance after convergence.
# =============================================================================

#' @keywords internal
trim_null_effects <- function(data, params, model) {
  UseMethod("trim_null_effects")
}

#' @keywords internal
trim_null_effects.default <- function(data, params, model) {
  null_idx <- which(model$V < params$prior_tol)
  if (length(null_idx) == 0) return(model)

  model$V[null_idx] <- 0
  model$alpha[null_idx, ] <- rep(model$pi, each = length(null_idx))
  model$mu[null_idx, ] <- 0
  model$mu2[null_idx, ] <- 0
  model$lbf_variable[null_idx, ] <- 0
  model$lbf[null_idx] <- 0
  model$KL[null_idx] <- 0

  return(model)
}
