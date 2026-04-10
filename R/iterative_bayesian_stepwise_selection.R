# =============================================================================
# IBSS INITIALIZATION
#
# Initializes the SuSiE model object for Iterative Bayesian Stepwise Selection.
# Sets up model matrices, handles model_init, and prepares for IBSS.
# =============================================================================
#' Initialize IBSS model
#'
#' Creates and initializes the model object for the IBSS algorithm.
#'
#' @param data Data object (individual, ss, or rss_lambda)
#' @param params Validated params object
#'
#' @return Initialized model object ready for the IBSS iteration loop.
#' @importFrom utils modifyList
#' @export
#' @keywords internal
ibss_initialize <- function(data, params) {
  UseMethod("ibss_initialize")
}

#' @rdname ibss_initialize
#' @export
#' @keywords internal
ibss_initialize.default <- function(data, params) {

  # Set var(y)
  var_y <- get_var_y(data)

  # Adjust number of single effects if needed
  if (data$p < params$L) {
    params$L <- data$p
  }

  # Check & validate residual variance
  if (is.null(params$residual_variance)) {
    params$residual_variance <- var_y
  }
  # For multivariate models, residual_variance can be a matrix
  if (!is.matrix(params$residual_variance)) {
    if (!is.numeric(params$residual_variance)) {
      stop("Input residual variance sigma2 must be numeric.")
    }
    params$residual_variance <- as.numeric(params$residual_variance)
    if (length(params$residual_variance) != 1) {
      stop("Input residual variance sigma2 must be a scalar.")
    }
    if (params$residual_variance <= 0) {
      stop("Residual variance sigma2 must be positive (is your var(Y) zero?).")
    }
  }

  # Handle model initialization
  if (!is.null(params$model_init)) {
    # Validate the contents of model_init
    validate_init(data, params)

    # Prune effects with zero prior variance
    model_init_pruned <- prune_single_effects(params$model_init)

    # Adjust the number of effects
    adjustment <- adjust_L(params, model_init_pruned, var_y)
    params$L   <- adjustment$L

    # Create base model with all required fields
    mat_init <- initialize_susie_model(data, params, var_y)

    # Merge with adjusted model_init
    mat_init <- modifyList(mat_init, adjustment$model_init)

    # Reset iteration-specific values
    mat_init$KL  <- rep(as.numeric(NA), params$L)
    mat_init$lbf <- rep(as.numeric(NA), params$L)
  } else {
    # Create fresh model
    mat_init <- initialize_susie_model(data, params, var_y)
  }

  # Initialize fitted values and null index
  fitted     <- initialize_fitted(data, mat_init)
  null_index <- initialize_null_index(data, mat_init)

  # Preserve model class set by initialize_susie_model (e.g., "mvsusie")
  model_class <- class(mat_init)

  # Return assembled SuSiE object
  model <- c(mat_init,
             list(null_index = null_index),
             fitted)

  # Use the class from initialize_susie_model if it inherits from "susie",
  # otherwise default to "susie"
  if (inherits(mat_init, "susie")) {
    class(model) <- model_class
  } else {
    class(model) <- "susie"
  }
  model$converged <- FALSE

  # Initialize Gamma-Poisson slot activity (c_hat) if a slot_prior
  # is specified. c_hat[l] = posterior probability that slot l is active.
  # Prior: mu ~ Gamma(nu, nu/C), c_l | mu ~ Bern(mu/L).
  sp <- params$slot_prior
  if (!is.null(sp)) {
    if (!is.slot_prior(sp))
      stop("slot_prior must be created by slot_prior_binomial() or ",
           "slot_prior_poisson(). Got class: ", paste(class(sp), collapse=", "))
    C_val <- sp$C
    L <- nrow(model$alpha)
    nu <- sp$nu

    # Warm start from c_hat_init or uniform at C/L
    if (!is.null(sp$c_hat_init) && length(sp$c_hat_init) == L) {
      c_hat <- sp$c_hat_init
      a_g <- nu + sum(c_hat)
    } else {
      c_hat <- rep(min(C_val / L, 1 - 1e-10), L)
      a_g <- nu + C_val
    }
    b_g <- nu / max(C_val, 1e-6) + 1  # constant within block

    prior_type <- if (inherits(sp, "slot_prior_binomial")) "binomial"
                  else "poisson"

    model$slot_weights <- c_hat
    model$c_hat_state <- list(
      C = C_val, nu = nu, a_g = a_g, b_g = b_g,
      update_schedule = sp$update_schedule,
      prior_type = prior_type
    )

    # Adaptive skip threshold
    if (sp$skip_threshold_multiplier > 0) {
      baseline_logodds <- digamma(a_g) - log(b_g) - log(L)
      c_hat_baseline <- 1 / (1 + exp(-baseline_logodds))
      model$c_hat_state$skip_threshold <-
        sp$skip_threshold_multiplier * c_hat_baseline
    } else {
      model$c_hat_state$skip_threshold <- 0
    }

    # Recompute fitted values with slot weights.
    # ibss_initialize builds Xr/XtXr/Rz assuming weight=1 for all slots;
    # with c_hat < 1 we need to correct.
    model <- recompute_fitted_weighted(data, model)
  }

  return(model)
}

# =============================================================================
# IBSS FITTING
#
# Updates all L single effects in the SuSiE model for one IBSS iteration.
# Calls single_effect_update for each effect and validates prior variance estimates.
# =============================================================================
#'
#' @param data Data object (individual, ss, or rss_lambda)
#' @param params Validated params object
#' @param model Current SuSiE model object
#'
#' @return Updated SuSiE model object with new alpha, mu, mu2, V, lbf, KL, and
#' fitted values
#'
#' @keywords internal
#' @noRd
ibss_fit <- function(data, params, model) {

  L <- nrow(model$alpha)
  use_c_hat <- !is.null(model$c_hat_state)

  if (L > 0) {
    for (l in seq_len(L)) {
      # Skip inactive slots when c_hat is below the adaptive threshold.
      # (Faithfully ported from susieAnn ibss_ann.R:137-139)
      if (use_c_hat &&
          model$slot_weights[l] < model$c_hat_state$skip_threshold) {
        next
      }

      # Standard SER step (susieAnn ibss_ann.R:145).
      # compute_residuals and update_fitted_values use
      # model$slot_weights[l] to scale effect l's contribution.
      model <- single_effect_update(data, params, model, l)

      # Gamma-Poisson slot activity update (susieAnn ibss_ann.R:148-163)
      if (use_c_hat) {
        model <- update_c_hat(data, model, l)
      }
    }
  }

  # Batch Gamma shape update: once per sweep (standard CAVI schedule).
  if (use_c_hat && model$c_hat_state$update_schedule == "batch") {
    model$c_hat_state$a_g <- model$c_hat_state$nu + sum(model$slot_weights)
  }

  # Validate prior variance is reasonable
  validate_prior(data, params, model)

  return(model)
}

# =============================================================================
# GAMMA-POISSON SLOT ACTIVITY (c_hat) HELPERS
#
# These functions implement the Gamma-Poisson slot activity model from
# susieAnn's Algorithm 1. Each slot l has a posterior activity probability
# c_hat[l] = sigmoid(E_q[log mu] - log L + lbf[l]), where mu is the
# per-block causal rate with Gamma(nu, nu/C) prior.
#
# Ported faithfully from susieAnn/R/ibss_ann.R:109-163.
# =============================================================================

#' Update c_hat for one slot after its SER step
#'
#' Computes the Gamma-Poisson posterior slot activity, adjusts the
#' fitted-values field for the weight change, and updates the Gamma
#' shape parameter eagerly.
#'
#' @param data Data object.
#' @param model Current SuSiE model with c_hat_state.
#' @param l Slot index.
#' @return Updated model.
#' @keywords internal
#' @noRd
update_c_hat <- function(data, model, l) {
  st <- model$c_hat_state
  old_c <- model$slot_weights[l]
  L <- nrow(model$alpha)

  # Slot activity update: c_hat = sigmoid(E[log mu] - log L + lbf)
  # Binomial type adds exact correction: - log(1 - E[mu]/L)
  lbf_l <- model$lbf[l]
  if (is.na(lbf_l) || !is.finite(lbf_l)) lbf_l <- 0
  log_odds <- digamma(st$a_g) - log(st$b_g) - log(L) + lbf_l
  if (st$prior_type == "binomial") {
    mu_over_L <- st$a_g / (st$b_g * L)
    if (mu_over_L < 1) log_odds <- log_odds - log(1 - mu_over_L)
  }
  log_odds <- max(min(log_odds, 20), -20)  # numerical guard
  new_c <- 1 / (1 + exp(-log_odds))
  model$slot_weights[l] <- new_c

  # Correct fitted field for the change in slot weight.
  # After update_fitted_values, the field contains effect l at old_c.
  # Adjust by (new_c - old_c) * Rv(b_bar_l).
  # (susieAnn ibss_ann.R:156-160)
  if (abs(new_c - old_c) > 1e-15) {
    b_bar_l <- model$alpha[l, ] * model$mu[l, ]
    model <- adjust_fitted_for_c_hat(data, model, b_bar_l, new_c - old_c)
  }

  # Gamma shape update: sequential (per-slot) or deferred to batch (per-sweep).
  # Sequential converges faster per iteration but batch has stricter
  # ELBO monotonicity guarantee (standard CAVI).
  if (st$update_schedule == "sequential") {
    model$c_hat_state$a_g <- st$nu + sum(model$slot_weights)
  }

  return(model)
}

#' Adjust fitted-values field for a c_hat weight change
#'
#' Adds delta_weight * Rv(b_bar_l) to the appropriate fitted-values field
#' (Xr, XtXr, or Rz depending on data backend).
#'
#' @param data Data object.
#' @param model Current model.
#' @param b_bar_l Posterior mean vector alpha[l,] * mu[l,].
#' @param delta_weight Change in slot weight (new_c - old_c).
#' @return Updated model.
#' @keywords internal
#' @noRd
adjust_fitted_for_c_hat <- function(data, model, b_bar_l, delta_weight) {
  fitted_field <- detect_fitted_field(model)
  if (fitted_field == "Xr") {
    # Individual data: Xr += delta_weight * X %*% b_bar_l
    model$Xr <- model$Xr +
      delta_weight * as.vector(compute_Xb(data$X, b_bar_l))
  } else {
    # SS or RSS: fitted += delta_weight * R %*% b_bar_l
    Rv_delta <- as.vector(compute_Rv(data, b_bar_l, model$X_meta))
    model[[fitted_field]] <- model[[fitted_field]] + delta_weight * Rv_delta
  }
  return(model)
}

#' Detect which fitted-values field the model uses
#'
#' Returns "Rz" (rss_lambda path), "XtXr" (sufficient stats path),
#' or "Xr" (individual data path).
#'
#' @param model SuSiE model object.
#' @return Character string: fitted field name.
#' @keywords internal
#' @noRd
detect_fitted_field <- function(model) {
  # (susieAnn ibss_ann.R:109-112)
  if ("Rz" %in% names(model)) "Rz"
  else if ("XtXr" %in% names(model)) "XtXr"
  else if ("Xr" %in% names(model)) "Xr"
  else stop("Cannot detect fitted-values field on model object.")
}

#' Recompute fitted values with slot weights
#'
#' After ibss_initialize builds Xr/XtXr/Rz assuming weight=1 for all
#' slots, this function recomputes the field as the slot-weighted sum
#' sum_l(c_hat[l] * alpha[l,] * mu[l,]).
#'
#' @param data Data object.
#' @param model Model with slot_weights set.
#' @return Updated model.
#' @keywords internal
#' @noRd
recompute_fitted_weighted <- function(data, model) {
  # (susieAnn ibss_ann.R:120-124)
  # Recompute the fitted-values field as the slot-weighted sum.
  # For individual data: Xr = X %*% b_weighted (length n)
  # For SS/RSS: XtXr/Rz = R %*% b_weighted (length p)
  fitted_field <- detect_fitted_field(model)
  L <- nrow(model$alpha)
  c_hat <- model$slot_weights
  b_weighted <- rep(0, ncol(model$alpha))
  for (ll in seq_len(L)) {
    b_weighted <- b_weighted + c_hat[ll] * model$alpha[ll, ] * model$mu[ll, ]
  }
  if (fitted_field == "Xr") {
    # Individual data: compute X %*% b_weighted (length n)
    model$Xr <- as.vector(compute_Xb(data$X, b_weighted))
  } else {
    # SS or RSS: compute R %*% b_weighted or X'X %*% b_weighted (length p)
    model[[fitted_field]] <- as.vector(compute_Rv(data, b_weighted, model$X_meta))
  }
  return(model)
}

# =============================================================================
# IBSS FINALIZATION
#
# Finalizes the SuSiE model after convergence or maximum number of iterations
# reached. Computes credible sets, PIPs, intercept, fitted values, and z-scores.
# =============================================================================
#' Finalize IBSS model
#'
#' Computes credible sets, PIPs, z-scores, and cleans up temporary
#' fields from the model object.
#'
#' @param data Data object (individual, ss, or rss_lambda)
#' @param params Validated params object
#' @param model Converged model object
#' @param elbo ELBO values (optional)
#' @param iter Number of iterations completed
#' @param tracking Tracking data (optional)
#'
#' @return Finalized model object with credible sets and PIPs.
#' @export
#' @keywords internal
ibss_finalize <- function(data, params, model, elbo = NULL, iter = NA_integer_,
                          tracking = NULL) {

  # Append ELBO & iteration count to model output
  model$niter <- iter

  # Intercept & Fitted Values
  model$X_column_scale_factors <- get_scale_factors(data, params)
  model$intercept              <- get_intercept(data, params, model)
  model$fitted                 <- get_fitted(data, params, model)

  # Posterior Inclusion Probabilities, credible sets, z-scores
  model$sets <- get_cs(data, params, model)
  model$pip  <- susie_get_pip(model, prior_tol = params$prior_tol)
  model$z    <- get_zscore(data, params, model)

  # Tracking Across Iterations
  if (params$track_fit) model$trace <- tracking

  # Assign Variable Names
  model <- get_variable_names(data, model)

  # Sketch diagnostics (from data -> model, following sets/pip/z pattern)
  if (!is.null(data$sketch_diagnostics)) {
    model$sketch_diagnostics <- data$sketch_diagnostics
    # Store final-iteration per-variable penalty: v_j/sigma^2 = inflation - 1
    if (!is.null(model$shat2_inflation))
      model$sketch_diagnostics$per_variable_penalty <- model$shat2_inflation - 1
  }

  # Multi-panel omega weights
  if (!is.null(model$omega))
    model$omega_weights <- model$omega

  # Store Gamma-Poisson c_hat results on output for user access
  # and for susieAnn to extract (a_g, b_g needed for genome-wide nu update).
  if (!is.null(model$c_hat_state)) {
    model$c_hat <- model$slot_weights
    model$C_hat <- sum(model$slot_weights)
    model$a_g   <- model$c_hat_state$a_g
    model$b_g   <- model$c_hat_state$b_g
    model$c_hat_state <- NULL  # cleanup internal state
  }

  # Clean up temporary computational fields
  model <- cleanup_model(data, params, model)

  return(model)
}
