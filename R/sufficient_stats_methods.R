# =============================================================================
# DATA INITIALIZATION & CONFIGURATION
#
# Functions for data object setup, configuration, and preprocessing.
# These prepare data objects for model fitting and handle data-specific
# configurations like unmappable effects.
#
# Functions: configure_data, get_var_y
# =============================================================================

# Configure ss data for specified method
#' @keywords internal
configure_data.ss <- function(data, params) {
  if (params$unmappable_effects == "inf") {
    return(add_eigen_decomposition(data, params))
  } else {
    return(configure_data.default(data, params))
  }
}

# Get variance of y
#' @keywords internal
get_var_y.ss <- function(data, ...) {
  return(data$yty / (data$n - 1))
}

# =============================================================================
# MODEL INITIALIZATION & SETUP
#
# Functions for initializing model objects and setting up initial states.
# These create model matrices, initialize fitted values, and prepare
# the SuSiE model for iterative fitting.
#
# Functions: initialize_susie_model, initialize_fitted, validate_prior, track_ibss_fit
# =============================================================================

# Initialize SuSiE model
#' @keywords internal
initialize_susie_model.ss <- function(data, params, var_y, ...) {

  # Base model
  model <- initialize_matrices(data, params, var_y)

  # Append predictor weights and initialize non-sparse quantities
  if (params$unmappable_effects == "inf") {
    # Initialize omega quantities for unmappable effects
    omega_res               <- compute_omega_quantities(data, tau2 = 0, sigma2 = var_y)
    model$omega_var         <- omega_res$omega_var
    model$predictor_weights <- omega_res$diagXtOmegaX
    model$XtOmegay          <- data$eigen_vectors %*% (data$VtXty / omega_res$omega_var)

    # Initialize unmappable variance component and coefficients
    model$tau2  <- 0
    model$theta <- rep(0, data$p)

  } else if (params$unmappable_effects == "ash") {
    model$predictor_weights <- attr(data$XtX, "d")
    model$tau2              <- 0
    model$theta             <- rep(0, data$p)
    model$XtX_theta         <- rep(0, data$p)
    model$masked         <- rep(FALSE, data$p)  # Track masked variants for LD-aware exclusion
    model$ash_iter          <- 0                # Track Mr.ASH iterations
    model$ash_pi          <- NULL
    model$diffuse_iter_count <- rep(0, params$L)
    model$prev_sentinel <- rep(0, params$L)  # Track sentinel varables
    model$unmask_candidate_iters <- rep(0, data$p)  
    model$ever_unmasked <- rep(FALSE, data$p)
    model$force_exposed_iter <- rep(0, data$p)     # When position was force-exposed (0 = never)
    model$ever_diffuse <- rep(0, params$L)
    model$second_chance_used <- rep(FALSE, data$p) # Permanent protection after second chance
    model$prev_case <- rep(0, params$L)            # Track previous case assignment for oscillation detection
  } else {
    model$predictor_weights <- attr(data$XtX, "d")
  }

  return(model)
}

# Initialize fitted values
#' @keywords internal
initialize_fitted.ss <- function(data, mat_init) {
  return(list(XtXr = as.vector(data$XtX %*% colSums(mat_init$alpha * mat_init$mu))))
}

# Validate Prior Variance
#' @keywords internal
validate_prior.ss <- function(data, params, model, ...) {
  if (isTRUE(params$check_prior)) {
    if (is.null(data$zm)) {
      bhat <- data$Xty / model$predictor_weights
      shat <- sqrt(model$sigma2 / model$predictor_weights)
      z <- bhat / shat
      data$zm <- max(abs(z[!is.nan(z)]))
    }
    if (any(model$V > 100 * (data$zm^2))) {
      stop(
        "Estimated prior variance is unreasonably large.\n",
        "This usually caused by mismatch between the summary statistics and the LD Matrix.\n",
        "Please check the input."
      )
    }
  }
  return(validate_prior.default(data, params, model, ...))
}

# Track core parameters across iterations
#' @keywords internal
track_ibss_fit.ss <- function(data, params, model, tracking, iter, elbo, ...) {
  if (params$unmappable_effects %in% c("inf", "ash")) {
    # Append non-sparse variance component to tracking
    tracking <- track_ibss_fit.default(data, params, model, tracking, iter, elbo, ...)
    if (isTRUE(params$track_fit)) {
      tracking[[iter]]$tau2 <- model$tau2
    }
    return(tracking)
  } else {
    # Use default for standard SS case
    return(track_ibss_fit.default(data, params, model, tracking, iter, elbo, ...))
  }
}

# =============================================================================
# SINGLE EFFECT REGRESSION & ELBO
#
# Core functions for single effect regression computation and ELBO calculation.
# These handle the mathematical core of SuSiE including residual computation, SER
# statistics, posterior moments, and log-likelihood calculations for the ELBO.
#
# Functions: compute_residuals, compute_ser_statistics, SER_posterior_e_loglik,
# calculate_posterior_moments, compute_kl, get_ER2, Eloglik, loglik, neg_loglik
# =============================================================================

# Compute residuals for single effect regression
#' @keywords internal
compute_residuals.ss <- function(data, params, model, l, ...) {
  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: Omega-weighted residuals
    b <- colSums(model$mu * model$alpha) - model$mu[l, ] * model$alpha[l, ]

    # Compute Residuals
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    XtOmegay <- data$eigen_vectors %*% (data$VtXty / omega_res$omega_var)
    XtOmegaXb <- data$eigen_vectors %*% ((t(data$eigen_vectors) %*% b) * data$eigen_values / omega_res$omega_var)
    XtOmegar <- XtOmegay - XtOmegaXb

    # Store residuals and parameters (unmappable case)
    model$residuals         <- XtOmegar
    model$predictor_weights <- omega_res$diagXtOmegaX  # Update for this iteration
    model$residual_variance <- 1                       # Already incorporated in Omega

    return(model)

  } else if (params$unmappable_effects == "ash") {
    # SuSiE-ash: explicit residualization
    XtXr_without_l <- model$XtXr - data$XtX %*% (model$alpha[l, ] * model$mu[l, ])

    # Subtract X'X*theta from residuals
    XtR <- data$Xty - model$XtX_theta - XtXr_without_l

    # Store residuals and parameters
    model$residuals         <- XtR
    model$fitted_without_l  <- XtXr_without_l
    model$residual_variance <- model$sigma2
    

    return(model)

  } else {
    # Remove lth effect from fitted values
    XtXr_without_l <- model$XtXr - data$XtX %*% (model$alpha[l, ] * model$mu[l, ])

    # Compute Residuals
    XtR <- data$Xty - XtXr_without_l

    # Store residuals and parameters (standard case)
    model$residuals         <- XtR
    model$fitted_without_l  <- XtXr_without_l # For fitted update
    model$residual_variance <- model$sigma2  # Standard residual variance

    return(model)
  }
}

# Compute SER statistics
#' @keywords internal
compute_ser_statistics.ss <- function(data, params, model, l, ...) {
  betahat <- (1 / model$predictor_weights) * model$residuals
  shat2   <- model$residual_variance / model$predictor_weights

  # Optimization parameters
  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: optimize on linear scale
    optim_init   <- model$V[l]
    optim_bounds <- c(0, 1)
    optim_scale  <- "linear"
  } else {
    # Standard SuSiE and SuSiE-ash: optimize on log scale
    optim_init   <- log(max(c(betahat^2 - shat2, 1), na.rm = TRUE))
    optim_bounds <- c(-30, 15)
    optim_scale  <- "log"
  }

  return(list(
    betahat      = betahat,
    shat2        = shat2,
    optim_init   = optim_init,
    optim_bounds = optim_bounds,
    optim_scale  = optim_scale
  ))
}

# Posterior expected log-likelihood for a single effect regression
#' @keywords internal
SER_posterior_e_loglik.ss <- function(data, params, model, l) {
  Eb  <- model$alpha[l, ] * model$mu[l, ]
  Eb2 <- model$alpha[l, ] * model$mu2[l, ]

  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: Omega-weighted likelihood
    return(-0.5 * (-2 * sum(Eb * model$residuals) + sum(model$predictor_weights * as.vector(Eb2))))
  } else {
    # Standard SuSiE and SuSiE-ash
    return(-0.5 / model$residual_variance * (-2 * sum(Eb * model$residuals) + sum(model$predictor_weights * as.vector(Eb2))))
  }
}

# Calculate posterior moments for single effect regression
#' @keywords internal
calculate_posterior_moments.ss <- function(data, params, model, V, l, ...) {
  # Standard Gaussian posterior calculations
  post_var   <- (1 / V + model$predictor_weights / model$residual_variance)^(-1)
  post_mean  <- (1 / model$residual_variance) * post_var * model$residuals
  post_mean2 <- post_var + post_mean^2

  # Store posterior moments in model
  model$mu[l, ] <- post_mean
  model$mu2[l, ] <- post_mean2

  return(model)
}

# Calculate KL divergence
#' @keywords internal
compute_kl.ss <- function(data, params, model, l) {
  model <- compute_kl.default(data, params, model, l)
  return(model)
}

# Expected Squared Residuals
#' @keywords internal
get_ER2.ss <- function(data, model) {
  B       <- model$alpha * model$mu
  XB2     <- sum((B %*% data$XtX) * B)
  betabar <- colSums(B)
  postb2  <- model$alpha * model$mu2 # Posterior second moment.

  return(data$yty - 2 * sum(betabar * data$Xty) + sum(betabar * (data$XtX %*% betabar)) -
           XB2 + sum(model$predictor_weights * t(postb2)))
}

# Expected log-likelihood
#' @keywords internal
Eloglik.ss <- function(data, model) {
  # Standard log-likelihood computation
  return(-data$n / 2 * log(2 * pi * model$sigma2) -
           1 / (2 * model$sigma2) * get_ER2(data, model))
}

#' @importFrom Matrix colSums
#' @importFrom stats dnorm
#' @keywords internal
loglik.ss <- function(data, params, model, V, ser_stats, l = NULL, ...) {
  # log(bf) for each SNP
  lbf <- dnorm(ser_stats$betahat, 0, sqrt(V + ser_stats$shat2), log = TRUE) -
    dnorm(ser_stats$betahat, 0, sqrt(ser_stats$shat2), log = TRUE)

  # Stabilize logged Bayes Factor
  stable_res  <- lbf_stabilization(lbf, model$pi, ser_stats$shat2)

  # Compute posterior weights
  weights_res <- compute_posterior_weights(stable_res$lpo)

  # Store in model if l is provided, otherwise return lbf_model for prior variance optimization
  if (!is.null(l)) {
    model$alpha[l, ] <- weights_res$alpha
    model$lbf[l] <- weights_res$lbf_model
    model$lbf_variable[l, ] <- stable_res$lbf
    return(model)
  } else {
    return(weights_res$lbf_model)
  }
}

#' @keywords internal
neg_loglik.ss <- function(data, params, model, V_param, ser_stats, ...) {
  # Convert parameter to V based on optimization scale
  V <- if (ser_stats$optim_scale == "log") exp(V_param) else V_param

  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: Omega-weighted objective with logSumExp trick
    return(-matrixStats::logSumExp(
      -0.5 * log(1 + V * model$predictor_weights) +
        V * model$residuals^2 / (2 * (1 + V * model$predictor_weights)) +
        log(model$pi + sqrt(.Machine$double.eps))
    ))
  } else {
    # Standard SuSiE and SuSiE-ash: standard objective
    lbf_model <- loglik.ss(data, params, model, V, ser_stats)
    return(-lbf_model)
  }
}

# =============================================================================
# MODEL UPDATES & FITTING
#
# Functions for iterative model updates and variance component estimation.
# These handle the dynamic aspects of model fitting including fitted value
# updates and variance component estimation.
#
# Functions: update_fitted_values, update_variance_components, update_derived_quantities
# =============================================================================

# Update fitted values
#' @keywords internal
update_fitted_values.ss <- function(data, params, model, l, ...) {
  if (params$unmappable_effects == "inf") {
    # SuSiE-inf: include theta in fitted values
    model$XtXr <- as.vector(data$XtX %*% (colSums(model$alpha * model$mu) + model$theta))
  } else {
    # Standard SuSiE and SuSiE-ash: sparse component only
    model$XtXr <- model$fitted_without_l + as.vector(data$XtX %*% (model$alpha[l, ] * model$mu[l, ]))
  }
  return(model)
}

# Update variance components for ss data
#' @keywords internal
update_variance_components.ss <- function(data, params, model, ...) {
  if (params$unmappable_effects == "inf") {
    # Calculate omega
    L         <- nrow(model$alpha)
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    omega     <- matrix(rep(omega_res$diagXtOmegaX, L), nrow = L, ncol = data$p, byrow = TRUE) +
      matrix(rep(1 / model$V, data$p), nrow = L, ncol = data$p, byrow = FALSE)

    # Compute theta for infinitesimal effects.
    theta <- compute_theta_blup(data, model)

    # Sigma2 and tau2 update
    if (params$estimate_residual_method == "MLE") {
      mle_result <- mle_unmappable(data, params, model, omega)
      return(list(sigma2 = mle_result$sigma2,
                  tau2   = mle_result$tau2,
                  theta  = theta))
    } else {
      mom_result <- mom_unmappable(data, params, model, omega, model$tau2)
      return(list(sigma2 = mom_result$sigma2,
                  tau2   = mom_result$tau2,
                  theta  = theta))
    }
  } else if (params$unmappable_effects == "ash") {
    # =========================================================================
    # SuSiE-ash: Hybrid sparse + adaptive shrinkage model
    # 
    # KEY INSIGHT: Protect SuSiE's sparse effects from Mr.ASH absorption,
    # but let Mr.ASH absorb unmappable and unreliable signals.
    # 
    # Two types of diffusion (both indicate unreliable effects):
    #   1. WITHIN-EFFECT: Low purity - spread across variants not in tight LD
    #   2. CROSS-EFFECT: Sentinel collision - multiple effects compete for 
    #      same position (composite signal, not clean single causal)
    #
    # Classification into three cases:
    #   CASE 1 (diffuse): purity < 0.1 - protect neighborhood loosely
    #   CASE 2 (uncertain): low purity OR ever_diffuse - expose to Mr.ASH
    #   CASE 3 (confident): purity >= 0.5 AND never diffuse - subtract from residuals
    #
    # Cross-effect diffusion tracking:
    #   - Detect via sentinel collision (sentinels in tight LD across effects)
    #   - Mark effect as ever_diffuse (sticky, effect-level)
    #   - ever_diffuse effects get zero protection permanently
    #   - Real signals survive Mr.ASH competition; composites get absorbed
    #
    # Low purity (non-diffuse) effects:
    #   - Use wait-then-expose mechanism
    #   - Second chance allows recovery after Mr.ASH testing
    # =========================================================================

    # --- Protection thresholds ---
    # cPIP > 40% in reasonable LD range is plausible signal that we should mask from mr.ash
    neighborhood_pip_threshold <- if (!is.null(params$neighborhood_pip_threshold)) params$neighborhood_pip_threshold else 0.4
    # >10% chance of signal from a standalone variant is something we are potentially interested in 
    direct_pip_threshold <- if (!is.null(params$direct_pip_threshold)) params$direct_pip_threshold else 0.1
    # |R| = 0.5 (R^2 = 0.25) is reasonable LD range, following SuSiE default
    ld_threshold <- if (!is.null(params$ld_threshold)) params$ld_threshold else 0.5
    
    # --- Purity thresholds ---
    # Tentative CS coverage, default to 0.9 because it is tentative we don't have to require 0.95
    cs_threshold <- if (!is.null(params$working_cs_threshold)) params$working_cs_threshold else 0.9
    # >10% CS formation threshold - if cannot form 10% CS with reasonable purity, effect is emerging
    cs_formation_threshold <- if (!is.null(params$cs_formation_threshold)) params$cs_formation_threshold else 0.1
    # |R| = 0.5 (R^2 = 0.25) is reasonable purity for confident effects, following SuSiE default
    purity_threshold <- if (!is.null(params$purity_threshold)) params$purity_threshold else 0.5
    
    # --- LD thresholds for collision and exposure ---
    # |R| = 0.9 (R^2 = 0.81) for collision detection - effects sharing 81% variance are clearly competing
    collision_ld_threshold <- if (!is.null(params$collision_ld_threshold)) params$collision_ld_threshold else 0.9
    # |R| = 0.95 for tight LD neighborhood in exposure - variants nearly indistinguishable for fine-mapping
    tight_ld_threshold <- if (!is.null(params$tight_ld_threshold)) params$tight_ld_threshold else 0.95
    
    # --- Iteration counters for CASE 2 ---
    # Wait before expose to mr.ash
    # Default to 2 to expose early (not over-protecting it)
    diffuse_iter_count <- if (!is.null(params$diffuse_iter_count)) params$diffuse_iter_count else 2
    track_sentinel <- if (!is.null(params$track_sentinel)) params$track_sentinel else TRUE

    # --- Second chance mechanism ---
    # After force-exposing, wait N iterations then restore protection
    # Default to 3 to allow enough time for SuSiE to land on something else
    second_chance_wait <- if (!is.null(params$second_chance_wait)) params$second_chance_wait else 3

    # --- Unmasking stability ---
    # A masked position unmasks when SuSiE loses interest (PIP drops).
    # We wait 2 iterations to confirm this is stable, not oscillation.
    delayed_unmask_iter <- 2  # set this to 2 to prevent oscillation which is the minimum we can do
    
    L <- nrow(model$alpha)
    p <- ncol(model$alpha)
    model$ash_iter <- model$ash_iter + 1

    if (any(!(diag(data$XtX) %in% c(0, 1)))) {
      Xcorr <- safe_cov2cor(data$XtX)
    } else {
      Xcorr <- data$XtX
    }

    # =========================================================================
    # First pass: Compute sentinels and purity
    # =========================================================================
    sentinels <- apply(model$alpha, 1, which.max)
    effect_purity <- rep(NA, L)
    
    for (l in 1:L) {
      alpha_order <- order(model$alpha[l,], decreasing = TRUE)
      cumsum_alpha <- cumsum(model$alpha[l, alpha_order])
      cs_size <- sum(cumsum_alpha <= cs_threshold) + 1
      cs_indices <- alpha_order[1:min(cs_size, p)]
      effect_purity[l] <- get_purity(cs_indices, X = NULL, Xcorr = Xcorr, use_rfast = FALSE)[1]
    }

    # Detect current collision and update ever_diffuse
    current_collision <- rep(FALSE, L)
    for (l in 1:L) {
      # Skip effects with uniform/near-uniform alpha (no meaningful signal)
      if (max(model$alpha[l,]) - min(model$alpha[l,]) < 5e-5) next
      
      sentinel_l <- sentinels[l]
      
      # Only compare with other effects that have meaningful signal
      for (other_l in (1:L)[-l]) {
        if (max(model$alpha[other_l,]) - min(model$alpha[other_l,]) < 5e-5) next
        # Use collision_ld_threshold (0.9) - effects sharing 81% variance are competing
        if (abs(Xcorr[sentinel_l, sentinels[other_l]]) > collision_ld_threshold) {
          current_collision[l] <- TRUE
        }
      }
      # FIXME: this only counts if it is ever diffused but not how many collision partners --- which might be alternatively informative
      # To do that, we can add this line instead to right below "current_collision[l] <- TRUE"
      model$ever_diffuse[l] <- model$ever_diffuse[l] + current_collision[l]
    }

    # Initialize per-iteration outputs
    b_confident <- rep(0, p)
    alpha_protected <- matrix(0, nrow = L, ncol = p)
    force_unmask <- rep(FALSE, p)
    force_mask <- rep(FALSE, p)
 
    # =========================================================================
    # Second pass: Classify effects and determine protection
    # =========================================================================
    current_case <- rep(0, L)  # Track current case for oscillation detection

    for (l in 1:L) {
      purity <- effect_purity[l]
      sentinel <- sentinels[l]

      # Reset counter if sentinel changed
      if (track_sentinel && sentinel != model$prev_sentinel[l] && model$prev_sentinel[l] > 0) {
        if (abs(Xcorr[sentinel, model$prev_sentinel[l]]) < tight_ld_threshold) {
          model$diffuse_iter_count[l] <- 0
        }
      }

      is_ever_diffuse <- model$ever_diffuse[l] > 0

      # Can enter CASE 3? Need: good purity AND never cross-effect diffuse
      can_be_confident <- (purity >= purity_threshold) && !is_ever_diffuse

      if (purity < cs_formation_threshold) {
        # =================================================================
        # CASE 1: Diffuse within effect (purity < 0.1)
        # =================================================================
        current_case[l] <- 1
        model$diffuse_iter_count[l] <- 0
        moderate_ld_with_sentinel <- abs(Xcorr[sentinel,]) > ld_threshold
        meaningful_alpha <- model$alpha[l,] > 5/p
        to_protect <- moderate_ld_with_sentinel | meaningful_alpha
        alpha_protected[l, to_protect] <- model$alpha[l, to_protect]
        force_mask <- force_mask | moderate_ld_with_sentinel
      } else if (!can_be_confident) {
        # =================================================================
        # CASE 2: Uncertain (current_collision OR ever_diffuse OR low purity)
        # =================================================================
        current_case[l] <- 2
        if (current_collision[l]) {
          # Active collision: NO protection (let Mr.ASH decide)
          model$diffuse_iter_count[l] <- 0
        } else {
          # ever_diffuse (no current collision) OR low purity: wait then expose
          # This keeps testing if signal is real
          model$diffuse_iter_count[l] <- model$diffuse_iter_count[l] + 1
          if (model$diffuse_iter_count[l] >= diffuse_iter_count) {
            # Use tight_ld_threshold (0.95) for exposure - only expose nearly indistinguishable variants
            tight_ld_with_sentinel <- abs(Xcorr[sentinel,]) > tight_ld_threshold

            newly_exposed <- tight_ld_with_sentinel &
                            !model$second_chance_used &
                            (model$force_exposed_iter == 0)
            model$force_exposed_iter[newly_exposed] <- model$ash_iter

            if (any(newly_exposed)) {
              model$diffuse_iter_count[l] <- 0
            }

            alpha_protected[l,] <- model$alpha[l,]
            expose_positions <- tight_ld_with_sentinel & !model$second_chance_used
            alpha_protected[l, expose_positions] <- 0
            force_unmask <- force_unmask | expose_positions
          } else {
            # Still waiting - protect fully
            alpha_protected[l,] <- model$alpha[l,]
          }
        }
      } else {
        # =================================================================
        # CASE 3: Confident (good purity AND never cross-effect diffuse)
        # =================================================================
        current_case[l] <- 3
        model$diffuse_iter_count[l] <- 0
        b_confident <- b_confident + model$alpha[l,] * model$mu[l,]
        alpha_protected[l,] <- model$alpha[l,]
      }

      model$prev_sentinel[l] <- sentinel
    }

    # =========================================================================
    # Oscillation detection: C2 <-> C3 transitions indicate unstable effects
    # When detected, mark as ever_diffuse to prevent future C3 classification
    # This breaks the feedback loop that causes non-convergence
    # =========================================================================
    for (l in 1:L) {
      prev <- model$prev_case[l]
      curr <- current_case[l]

      # Skip inactive effects (case 0) and first iteration (prev_case = 0)
      if (curr == 0 || prev == 0) next

      # Detect C2 <-> C3 oscillation (transition in either direction)
      if ((prev == 2 && curr == 3) || (prev == 3 && curr == 2)) {
        # Mark effect as oscillating - this makes ever_diffuse > 0
        # which permanently prevents C3 classification
        model$ever_diffuse[l] <- model$ever_diffuse[l] + 1

        # If we just classified as C3 but detected oscillation,
        # we need to undo the C3 treatment for this iteration
        if (curr == 3) {
          # Remove from b_confident
          b_confident <- b_confident - model$alpha[l,] * model$mu[l,]
          # Keep alpha_protected as is (will be treated as C2 effectively)
        }
      }
    }

    # Update prev_case for next iteration
    model$prev_case <- current_case

    # Compute residuals (with confident effects removed) and run Mr.ASH
    pip_protected <- susie_get_pip(alpha_protected)
    residuals <- data$y - data$X %*% b_confident

    # =========================================================================
    # Masking logic
    # Mask Mr.ASH theta at positions where SuSiE has signal (or nearby).
    # This prevents double-counting between SuSiE and Mr.ASH.
    # =========================================================================
    LD_adj <- abs(Xcorr) > ld_threshold
    neighborhood_pip <- as.vector(LD_adj %*% pip_protected)
    want_masked <- (neighborhood_pip > neighborhood_pip_threshold) | 
                   (pip_protected > direct_pip_threshold) |
                   force_mask

    # Track iterations wanting unmask (for stable unmask after N iterations)
    dont_want_mask <- !want_masked
    model$unmask_candidate_iters[model$masked & dont_want_mask] <- 
      model$unmask_candidate_iters[model$masked & dont_want_mask] + 1
    model$unmask_candidate_iters[want_masked | !model$masked] <- 0

    # Unmask if stable for N iterations (one-chance rule), or force unmask from CASE 2
    ready_to_unmask <- (model$masked & 
                       (model$unmask_candidate_iters >= delayed_unmask_iter) &
                       !model$ever_unmasked) |
                       (model$masked & force_unmask)

    model$ever_unmasked[ready_to_unmask] <- TRUE
    masked <- (model$masked | want_masked) & !ready_to_unmask & !model$ever_unmasked

    # =========================================================================
    # Second chance: restore protection after wait period
    # 
    # After CASE 2 force-exposes positions, Mr.ASH has N iterations to absorb
    # any noise. If positions are re-masked after this window, SuSiE can
    # potentially recover the signal if it was legitimate.
    # =========================================================================
    waited_long_enough <- (model$force_exposed_iter > 0) & 
                          (model$ash_iter - model$force_exposed_iter) >= second_chance_wait
    should_restore <- waited_long_enough & !model$second_chance_used
    
    if (any(should_restore)) {
      model$second_chance_used[should_restore] <- TRUE
      model$force_exposed_iter[should_restore] <- 0
      model$ever_unmasked[should_restore] <- FALSE
      masked[should_restore] <- TRUE
    }

    # Fit Mr.ASH: initialize from previous fit
    convtol <- if (model$ash_iter < 2) 1e-3 else 1e-4
    mrash_output <- mr.ash(
      X             = data$X,
      y             = residuals,
      intercept     = FALSE,
      standardize   = FALSE,
      sigma2        = model$sigma2,
      update.sigma2 = params$estimate_residual_variance,
      beta.init     = model$theta,
      pi            = model$ash_pi,
      tol = list(convtol = convtol, epstol = 1e-12),
      verbose = params$verbose,
      max.iter      = 1000
    )

    theta_new  <- mrash_output$beta
    sigma2_new <- mrash_output$sigma2
    tau2_new   <- sum(mrash_output$data$sa2 * mrash_output$pi) * mrash_output$sigma2
    
    # Zero out theta for masked variants
    theta_new[masked] <- 0
    if (FALSE) {
      diagnose_susie_ash_iter(data, model, Xcorr, mrash_output,
                      effect_purity, sentinels, current_collision,
                      alpha_protected, pip_protected, neighborhood_pip, 
                      masked, sigma2_new, tau2_new,
                      want_masked, ready_to_unmask, force_unmask, force_mask,
                      purity_threshold, cs_formation_threshold,
                      collision_ld_threshold, tight_ld_threshold,
                      ld_threshold, diffuse_iter_count)
    }

    XtX_theta_new <- as.vector(data$XtX %*% theta_new)

    return(list(
      sigma2              = sigma2_new,
      tau2                = tau2_new,
      theta               = theta_new,
      XtX_theta           = XtX_theta_new,
      ash_pi              = mrash_output$pi,
      sa2                 = mrash_output$data$sa2,
      ash_iter            = model$ash_iter,
      diffuse_iter_count = model$diffuse_iter_count,
      prev_sentinel       = model$prev_sentinel,
      masked              = masked,
      unmask_candidate_iters = model$unmask_candidate_iters,
      ever_unmasked       = model$ever_unmasked,
      force_exposed_iter  = model$force_exposed_iter,
      second_chance_used  = model$second_chance_used,
      ever_diffuse        = model$ever_diffuse,
      prev_case           = model$prev_case
      ))
  } else {
    # Use default method for standard SuSiE
    return(update_variance_components.default(data, params, model))
  }
}

# Update derived quantities for ss data
#' @keywords internal
update_derived_quantities.ss <- function(data, params, model) {
  if (params$unmappable_effects == "inf") {
    # Update omega quantities for next iteration
    omega_res               <- compute_omega_quantities(data, model$tau2, model$sigma2)
    model$omega_var         <- omega_res$omega_var
    model$predictor_weights <- omega_res$diagXtOmegaX
    model$XtOmegay          <- data$eigen_vectors %*% (data$VtXty / omega_res$omega_var)
    # Update fitted values to include theta
    b          <- colSums(model$alpha * model$mu)
    model$XtXr <- data$XtX %*% (b + model$theta)
    return(model)
  } else {
    return(update_derived_quantities.default(data, params, model))
  }
}

# =============================================================================
# OUTPUT GENERATION & POST-PROCESSING
#
# Functions for generating final results and summary statistics.
# These process fitted models into interpretable outputs including
# credible sets, variable names, and fitted values.
#
# Functions: get_scale_factors, get_intercept, get_fitted, get_cs,
# get_variable_names, get_zscore
# =============================================================================

# Get column scale factors
#' @keywords internal
get_scale_factors.ss <- function(data, params) {
  return(attr(data$XtX, "scaled:scale"))
}

# Get intercept
#' @keywords internal
get_intercept.ss <- function(data, params, model, ...) {
  return(data$y_mean - sum(data$X_colmeans * (colSums(model$alpha * model$mu) / model$X_column_scale_factors)))
}

# Get Fitted Values
#' @keywords internal
get_fitted.ss <- function(data, params, model, ...) {
  return(get_fitted.default(data, params, model, ...))
}

# Get Credible Sets
#' @keywords internal
get_cs.ss <- function(data, params, model, ...) {
  if (is.null(params$coverage) || is.null(params$min_abs_corr)) {
    return(NULL)
  }

  if (any(!(diag(data$XtX) %in% c(0, 1)))) {
    Xcorr <- safe_cov2cor(data$XtX)
  } else {
    Xcorr <- data$XtX
  }

  return(susie_get_cs(model,
                      Xcorr           = Xcorr,
                      check_symmetric = FALSE,
                      coverage        = params$coverage,
                      min_abs_corr    = params$min_abs_corr,
                      n_purity        = params$n_purity))
}

# Get Variable Names
#' @keywords internal
get_variable_names.ss <- function(data, model, ...) {
  return(assign_names(data, model, colnames(data$XtX)))
}

# Get univariate z-score
#' @keywords internal
get_zscore.ss <- function(data, params, model, ...) {
  return(get_zscore.default(data, params, model))
}

# Clean up model object for sufficient statistics data
#' @keywords internal
cleanup_model.ss <- function(data, params, model, ...) {
  # Remove common fields
  model <- cleanup_model.default(data, params, model, ...)

  # FIXME: for non-standard fields please connect them to "runtime_xx" where xx is unmappable effect option
  
  # Remove SS-specific fields for unmappable effects
  if (!is.null(params$unmappable_effects) && params$unmappable_effects == "inf") {
    unmappable_fields <- c("omega_var", "XtOmegay")
    
    for (field in unmappable_fields) {
      if (field %in% names(model)) {
        model[[field]] <- NULL
      }
    }
  } else if (!is.null(params$unmappable_effects) && params$unmappable_effects == "ash") {
    ash_fields <- c("XtX_theta", "masked", "ash_iter")
    
    for (field in ash_fields) {
      if (field %in% names(model)) {
        model[[field]] <- NULL
      }
    }
  }
  
  return(model)
}
