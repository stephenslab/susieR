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

  # Initialize slot activity (c_hat) if specified
  sp <- params$slot_prior
  if (!is.null(sp)) {
    if (!is.slot_prior(sp))
      stop("slot_prior must be created by slot_prior_betabinom() or ",
           "slot_prior_poisson(). ",
           "Got class: ", paste(class(sp), collapse=", "))
    L <- nrow(model$alpha)

    if (inherits(sp, "slot_prior_betabinom")) {
      a_beta <- sp$a_beta
      b_beta <- sp$b_beta
      prior_mean <- a_beta / (a_beta + b_beta)

      if (!is.null(sp$c_hat_init) && length(sp$c_hat_init) == L) {
        c_hat <- sp$c_hat_init
      } else {
        c_hat <- rep(min(prior_mean, 1 - 1e-10), L)
      }

      model$slot_weights <- c_hat
      model$c_hat_state <- list(
        prior_type = "betabinom", a_beta = a_beta, b_beta = b_beta,
        update_schedule = sp$update_schedule,
        skip_threshold_multiplier = sp$skip_threshold_multiplier,
        skip_threshold = 0
      )
    } else {
      C_val <- sp$C
      nu <- sp$nu

      if (!is.null(sp$c_hat_init) && length(sp$c_hat_init) == L) {
        c_hat <- sp$c_hat_init
        a_g <- nu + sum(c_hat)
      } else {
        c_hat <- rep(min(C_val / L, 1 - 1e-10), L)
        a_g <- nu + C_val
      }
      b_g <- nu / max(C_val, 1e-6) + 1

      model$slot_weights <- c_hat
      model$c_hat_state <- list(
        prior_type = "poisson", C = C_val, nu = nu, a_g = a_g, b_g = b_g,
        update_schedule = sp$update_schedule,
        skip_threshold_multiplier = sp$skip_threshold_multiplier,
        skip_threshold = 0
      )
    }

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
      if (use_c_hat &&
          model$slot_weights[l] < model$c_hat_state$skip_threshold) {
        next
      }

      model <- single_effect_update(data, params, model, l)

      if (use_c_hat) {
        model <- update_c_hat(data, model, l)
      }
    }
  }

  # Gamma-Poisson batch shape update (once per sweep)
  if (use_c_hat && model$c_hat_state$update_schedule == "batch" &&
      model$c_hat_state$prior_type != "betabinom") {
    model$c_hat_state$a_g <- model$c_hat_state$nu + sum(model$slot_weights)
  }

  # Adaptive skip threshold: baseline c_hat with lbf=0, scaled by multiplier
  if (use_c_hat && model$c_hat_state$skip_threshold_multiplier > 0) {
    st <- model$c_hat_state
    L_val <- nrow(model$alpha)
    if (st$prior_type == "betabinom") {
      # Self-consistent baseline: one Newton step for k_{-l} = k_total - baseline
      k_total <- sum(model$slot_weights)
      approx <- log(st$a_beta + k_total) - log(st$b_beta + L_val - 1 - k_total)
      k_others <- k_total - 1 / (1 + exp(-approx))
      baseline_logodds <- log(st$a_beta + k_others) -
                          log(st$b_beta + L_val - 1 - k_others)
    } else {
      baseline_logodds <- digamma(st$a_g) - log(st$b_g) - log(L_val)
    }
    c_hat_baseline <- 1 / (1 + exp(-baseline_logodds))
    model$c_hat_state$skip_threshold <-
      st$skip_threshold_multiplier * c_hat_baseline
  }

  # Validate prior variance is reasonable
  validate_prior(data, params, model)

  return(model)
}

# =============================================================================
# SLOT ACTIVITY (c_hat) HELPERS
#
# c_hat[l] = posterior probability that slot l is active.
# Beta-Binomial: logit(c_l) = log(a + k_{-l}) - log(b + L-1 - k_{-l}) + lbf_l
# Gamma-Poisson: logit(c_l) = psi(a_g) - log(b_g) - log(L) + lbf_l
# =============================================================================

#' Update c_hat for slot l after its SER step.
#' @keywords internal
#' @noRd
update_c_hat <- function(data, model, l) {
  st <- model$c_hat_state
  old_c <- model$slot_weights[l]
  L <- nrow(model$alpha)

  lbf_l <- model$lbf[l]
  if (is.na(lbf_l) || !is.finite(lbf_l)) lbf_l <- 0

  if (st$prior_type == "betabinom") {
    k_others <- sum(model$slot_weights[-l])
    log_odds <- log(st$a_beta + k_others) -
                log(st$b_beta + L - 1 - k_others) + lbf_l
  } else {
    log_odds <- digamma(st$a_g) - log(st$b_g) - log(L) + lbf_l
  }
  log_odds <- max(min(log_odds, 20), -20)
  new_c <- 1 / (1 + exp(-log_odds))
  model$slot_weights[l] <- new_c

  if (abs(new_c - old_c) > 1e-15) {
    b_bar_l <- model$alpha[l, ] * model$mu[l, ]
    model <- adjust_fitted_for_c_hat(data, model, b_bar_l, new_c - old_c)
  }

  # Gamma shape update (sequential mode; Beta-Binomial uses k_others directly)
  if (st$prior_type != "betabinom" && st$update_schedule == "sequential") {
    model$c_hat_state$a_g <- st$nu + sum(model$slot_weights)
  }

  return(model)
}

#' Add delta_weight * R*b_bar_l to the fitted-values field (Xr/XtXr/Rz).
#' @keywords internal
#' @noRd
adjust_fitted_for_c_hat <- function(data, model, b_bar_l, delta_weight) {
  fitted_field <- detect_fitted_field(model)
  if (fitted_field == "Xr") {
    model$Xr <- model$Xr +
      delta_weight * as.vector(compute_Xb(data$X, b_bar_l))
  } else {
    model[[fitted_field]] <- model[[fitted_field]] +
      delta_weight * as.vector(compute_Rv(data, b_bar_l, model$X_meta))
  }
  return(model)
}

#' @keywords internal
#' @noRd
detect_fitted_field <- function(model) {
  if ("Rz" %in% names(model)) "Rz"
  else if ("XtXr" %in% names(model)) "XtXr"
  else if ("Xr" %in% names(model)) "Xr"
  else stop("Cannot detect fitted-values field on model object.")
}

#' Recompute fitted = sum_l c_hat[l] * R * bbar[l] from scratch.
#' @keywords internal
#' @noRd
recompute_fitted_weighted <- function(data, model) {
  fitted_field <- detect_fitted_field(model)
  L <- nrow(model$alpha)
  c_hat <- model$slot_weights
  b_weighted <- rep(0, ncol(model$alpha))
  for (ll in seq_len(L)) {
    b_weighted <- b_weighted + c_hat[ll] * model$alpha[ll, ] * model$mu[ll, ]
  }
  if (fitted_field == "Xr") {
    model$Xr <- as.vector(compute_Xb(data$X, b_weighted))
  } else {
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
    if (model$c_hat_state$prior_type == "betabinom") {
      model$a_beta <- model$c_hat_state$a_beta
      model$b_beta <- model$c_hat_state$b_beta
    } else {
      model$a_g <- model$c_hat_state$a_g
      model$b_g <- model$c_hat_state$b_g
    }
    model$c_hat_state <- NULL  # cleanup internal state
  }

  # Clean up temporary computational fields
  model <- cleanup_model(data, params, model)

  return(model)
}
