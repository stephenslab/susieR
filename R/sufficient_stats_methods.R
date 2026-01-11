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

    # SuSiE-ash: explicit residualization approach
    model$predictor_weights <- attr(data$XtX, "d")
    model$tau2              <- 0
    model$theta             <- rep(0, data$p)
    model$sigma2_tilde      <- var_y
    model$XtX_theta         <- rep(0, data$p)
    model$contested         <- rep(FALSE, data$p)  # Track contested variants for LD-aware exclusion
    model$ash_iter          <- 0                   # Track Mr.ASH iterations for purity tracking
    model$min_purity        <- rep(1, params$L)    # Track minimum purity observed per effect

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
    model$residual_variance <- model$sigma2_tilde

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
#' @importFrom mr.ash.alpha mr.ash
#' @keywords internal
update_variance_components.ss <- function(data, params, model, ...) {
  if (params$unmappable_effects == "inf") {
    # Calculate omega
    L         <- nrow(model$alpha)
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    omega     <- matrix(rep(omega_res$diagXtOmegaX, L), nrow = L, ncol = data$p,byrow = TRUE) +
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
    # LD-aware exclusion with purity-based selective protection
    # High-purity effects: subtract from residuals, fully protect from Mr.ASH  
    # Low-purity effects: keep in residuals, expose sentinel's tight LD block to Mr.ASH
    pip_threshold         <- if (!is.null(params$pip_threshold)) params$pip_threshold else 0.5
    ld_threshold          <- if (!is.null(params$ld_threshold)) params$ld_threshold else 0.5
    purity_threshold      <- if (!is.null(params$purity_threshold)) params$purity_threshold else 0.5
    sentinel_ld_threshold <- if (!is.null(params$sentinel_ld_threshold)) params$sentinel_ld_threshold else 0.9
    # n_purity_iterations: number of iterations to track minimum purity per effect
    # If > 0, effects that show low purity during these iterations stay exposed forever
    # If <= 0, use current purity only (no memory)
    n_purity_iterations   <- if (!is.null(params$n_purity_iterations)) params$n_purity_iterations else 1

    n <- data$n
    L <- nrow(model$alpha)
    p <- ncol(model$alpha)

    # Update iteration counter
    model$ash_iter <- model$ash_iter + 1

    # Get correlation matrix for purity computation
    if (any(!(diag(data$XtX) %in% c(0, 1)))) {
      Xcorr <- muffled_cov2cor(data$XtX)
    } else {
      Xcorr <- data$XtX
    }

    # Compute purity for each effect and build:
    # 1. b_confident: effects to subtract from residuals (high purity only)
    # 2. pip_protected: PIPs that contribute to protection
    b_confident <- rep(0, p)
    pip_protected <- rep(0, p)
    effect_purity <- rep(NA, L)
    
    for (l in 1:L) {
      # Get CS for this effect (variants with non-negligible alpha)
      cs_threshold <- 0.9 # params$coverage  # coverage threshold, a bit more lenient for this purpose, not default 95%
      alpha_order <- order(model$alpha[l,], decreasing = TRUE)
      cumsum_alpha <- cumsum(model$alpha[l, alpha_order])
      cs_size <- sum(cumsum_alpha <= cs_threshold) + 1
      cs_indices <- alpha_order[1:min(cs_size, p)]
      
      # Compute purity for this effect
      purity <- get_purity(cs_indices, X = NULL, Xcorr = Xcorr, use_rfast = FALSE)[1]  # min purity
      effect_purity[l] <- purity
      
      # Track minimum purity during warmup iterations
      if (n_purity_iterations > 0 && model$ash_iter <= n_purity_iterations) {
        model$min_purity[l] <- min(model$min_purity[l], purity)
      }
      
      # Determine purity to use for protection decision
      purity_for_decision <- if (n_purity_iterations > 0) model$min_purity[l] else purity
      
      if (purity_for_decision >= purity_threshold) {
        # High purity: subtract from residuals AND protect entire effect
        b_confident <- b_confident + model$alpha[l,] * model$mu[l,]
        pip_protected <- pip_protected + model$alpha[l,]
      } else {
        # Low purity: DON'T subtract, protect everything EXCEPT sentinel + tight LD neighbors
        sentinel <- which.max(model$alpha[l,])
        tight_ld_with_sentinel <- abs(Xcorr[sentinel,]) > sentinel_ld_threshold
        
        # Add alpha to protection, but zero out sentinel's tight LD block
        alpha_protected <- model$alpha[l,]
        alpha_protected[tight_ld_with_sentinel] <- 0
        pip_protected <- pip_protected + alpha_protected
      }
    }

    # Compute residuals using only high-purity effects
    residuals <- data$y - data$X %*% b_confident

    # Build ash grid
    sa2 <- create_ash_grid(data)

    mrash_output <- mr.ash(
      X             = data$X,
      y             = residuals,
      sa2           = sa2,
      intercept     = FALSE,
      standardize   = FALSE,
      sigma2        = model$sigma2,
      update.sigma2 = params$estimate_residual_variance,
      max.iter      = 1000
    )

    theta_new  <- mrash_output$beta
    sigma2_new <- mrash_output$sigma2
    tau2_new   <- sum(mrash_output$data$sa2 * mrash_output$pi) * mrash_output$sigma2

    # Compute neighborhood PIP using protected PIPs only
    LD_adj <- abs(data$XtX) > n * ld_threshold
    neighborhood_pip <- as.vector(LD_adj %*% pip_protected)

    # Update contested: once contested, always contested (monotonic to prevent cycling)
    new_contested <- neighborhood_pip > pip_threshold
    contested <- model$contested | new_contested

    # Store theta stats before exclusion for diagnostics
    theta_sum2_before <- sum(theta_new^2)
    theta_max_before <- max(abs(theta_new))

    # Zero out theta for contested variants
    theta_new[contested] <- 0

    if (FALSE) {
      # Use min_purity for counting if n_purity_iterations > 0
      purity_used <- if (n_purity_iterations > 0) model$min_purity else effect_purity
      n_high_purity <- sum(purity_used >= purity_threshold, na.rm = TRUE)
      n_low_purity <- L - n_high_purity
      b_full <- colSums(model$alpha * model$mu)
      w <- colSums(data$X^2)
      grid_scale <- n / median(w)
      pi_null <- mrash_output$pi[1]
      pi_nonnull <- 1 - pi_null
      
      cat(sprintf("SuSiE-ash iter %d:\n", model$ash_iter))
      cat(sprintf("  Purity: %d high (>=%.2f), %d low | current: %s\n", 
                  n_high_purity, purity_threshold, n_low_purity,
                  paste(round(effect_purity, 2), collapse=", ")))
      if (n_purity_iterations > 0) {
        cat(sprintf("  Min purity (iter 1-%d): %s\n", n_purity_iterations,
                    paste(round(model$min_purity, 2), collapse=", ")))
      }
      cat(sprintf("  Effect lbf: %s\n", paste(round(model$lbf, 1), collapse=", ")))
      cat(sprintf("  Residuals: var(y-Xb_conf)=%.4f | b_full²=%.2e, b_conf²=%.2e\n",
                  var(as.vector(residuals)), sum(b_full^2), sum(b_confident^2)))
      cat(sprintf("  Grid: scale=%.3f (n/med(w)) | sa2 range=[%.2e, %.2e]\n",
                  grid_scale, min(sa2), max(sa2)))
      cat(sprintf("  Mr.ASH: sigma2=%.4f, tau2=%.2e | pi_null=%.1f%%, pi_nonnull=%.1f%%\n",
                  sigma2_new, tau2_new, pi_null*100, pi_nonnull*100))
      cat(sprintf("  Theta before: sum²=%.2e, max|θ|=%.4f | after exclusion: sum²=%.2e\n",
                  theta_sum2_before, theta_max_before, sum(theta_new^2)))
      cat(sprintf("  LD-exclusion: %d newly contested, %d total contested\n",
                  sum(new_contested), sum(contested)))
    }

    # Compute derived quantities for explicit residualization
    sigma2_tilde_new <- sigma2_new + tau2_new
    XtX_theta_new    <- as.vector(data$XtX %*% theta_new)

    return(list(
      sigma2       = sigma2_new,
      tau2         = tau2_new,
      theta        = theta_new,
      sigma2_tilde = sigma2_tilde_new,
      XtX_theta    = XtX_theta_new,
      contested    = contested,
      ash_iter     = model$ash_iter,
      min_purity   = model$min_purity,
      ash_pi       = mrash_output$pi,
      sa2          = mrash_output$data$sa2
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
    Xcorr <- muffled_cov2cor(data$XtX)
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
  
  # Remove SS-specific fields for unmappable effects
  if (!is.null(params$unmappable_effects) && params$unmappable_effects == "inf") {
    unmappable_fields <- c("omega_var", "XtOmegay")
    
    for (field in unmappable_fields) {
      if (field %in% names(model)) {
        model[[field]] <- NULL
      }
    }
  } else if (!is.null(params$unmappable_effects) && params$unmappable_effects == "ash" && params$verbose) {
    ash_fields <- c("sigma2_tilde", "XtX_theta", "contested", "ash_iter", "min_purity")
    
    for (field in ash_fields) {
      if (field %in% names(model)) {
        model[[field]] <- NULL
      }
    }
  }
  
  return(model)
}
