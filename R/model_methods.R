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
  mat_init <- list(
    alpha             = matrix(1 / data$p, L, data$p),
    mu                = matrix(0, L, data$p),
    mu2               = matrix(0, L, data$p),
    V                 = rep(params$scaled_prior_variance * var_y, L),
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
check_convergence <- function(data, params, model, elbo, iter, tracking) {
  UseMethod("check_convergence")
}

#' @keywords internal
check_convergence.default <- function(data, params, model, elbo, iter, tracking) {
  # Skip convergence check on first iteration
  if(iter == 1) {
    model$converged <- FALSE
    return(model)
  }

  # Calculate difference in ELBO values
  ELBO_diff   <- elbo[iter + 1] - tracking$convergence$prev_elbo
  ELBO_failed <- is.na(ELBO_diff) || is.infinite(ELBO_diff)

  if (params$convergence_method == "pip" || ELBO_failed) {
    # Fallback to PIP-based convergence if ELBO calculation fails
    if (ELBO_failed && params$convergence_method == "elbo") {
      warning_message(paste0("Iteration ", iter, " produced an NA/infinite ELBO
                             value. Using pip-based convergence this iteration."))
    }

    # For Servin-Stephens prior, require at least 3 iterations and average convergence
    # over 2 consecutive iterations for more stable convergence
    if (!is.null(params$use_servin_stephens) && params$use_servin_stephens) {
      if (iter <= 2) {
        model$converged <- FALSE
        return(model)  # Require at least 3 iterations
      }

      # Current iteration PIP difference
      current_diff <- max(abs(tracking$convergence$prev_alpha - model$alpha))

      # Average with previous iteration's difference if available
      if (!is.null(tracking$convergence$prev_pip_diff)) {
        avg_diff <- (current_diff + tracking$convergence$prev_pip_diff) / 2
      } else {
        avg_diff <- current_diff
      }

      # Store current diff for next iteration
      tracking$convergence$prev_pip_diff <- current_diff

      if (params$verbose) {
        message("max |change in PIP| (avg): ", format(avg_diff, digits = 6))
      }

      model$converged <- (avg_diff < params$tol)
      if (model$converged && !is.null(params$unmappable_effects) &&
          params$unmappable_effects == "ash") {
        model <- run_final_ash_pass(data, params, model)
      }
      return(model)
    } else {
      # Standard PIP convergence
      PIP_diff <- max(abs(tracking$convergence$prev_alpha - model$alpha))

      if (params$verbose) {
        message("max |change in PIP|: ", format(PIP_diff, digits = 6))
      }

      model$converged <- (PIP_diff < params$tol)
      if (model$converged && !is.null(params$unmappable_effects) &&
          params$unmappable_effects == "ash") {
        model <- run_final_ash_pass(data, params, model)
      }
      return(model)
    }
  }

  if (params$verbose) {
    message("ELBO: ", format(elbo[iter + 1], digits = 6))
  }

  model$converged <- (ELBO_diff < params$tol)
  if (model$converged && !is.null(params$unmappable_effects) &&
      params$unmappable_effects == "ash") {
    model <- run_final_ash_pass(data, params, model)
  }
  return(model)
}

# =============================================================================
# OBJECTIVE FUNCTION (ELBO)
# =============================================================================

#' @keywords internal
get_objective <- function(data, params, model) {
  UseMethod("get_objective")
}

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
      data$VtXty, data$yty
    )
  } else if (params$use_servin_stephens && nrow(model$alpha) == 1) {
    objective <- model$marginal_loglik[1]
  } else {
    # Standard ELBO computation
    objective <- Eloglik(data, model) - sum(model$KL)
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
  model$alpha[null_idx, ] <- 1 / ncol(model$alpha)
  model$mu[null_idx, ] <- 0
  model$mu2[null_idx, ] <- 0
  model$lbf_variable[null_idx, ] <- 0
  model$lbf[null_idx] <- 0
  model$KL[null_idx] <- 0

  return(model)
}
