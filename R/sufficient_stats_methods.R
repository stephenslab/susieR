# Sufficient statistics data backend methods

# Initialize fitted values
initialize_fitted.ss <- function(data, alpha, mu) {
  return(list(XtXr = compute_Xb(data$XtX, colSums(alpha * mu))))
}

# Initialize SuSiE model
initialize_susie_model.ss <- function(data, L, scaled_prior_variance, var_y,
                                      residual_variance, prior_weights, ...) {

  # Base model
  model <- initialize_matrices(data, L, scaled_prior_variance, var_y,
                               residual_variance, prior_weights)

  # Append predictor weights and initialize non-sparse quantities
  if (data$unmappable_effects %in% c("inf", "ash")) {

    # Initialize omega quantities for unmappable effects
    omega_res               <- compute_omega_quantities(data, tau2 = 0, sigma2 = 1)
    model$omega_var         <- omega_res$omega_var
    model$predictor_weights <- omega_res$diagXtOmegaX
    model$XtOmegay          <- data$eigen_vectors %*% (data$VtXty / omega_res$omega_var)

    # Initialize unmappable variance component and coefficients
    model$tau2  <- 0
    model$theta <- rep(0, data$p)
  } else {
    model$predictor_weights <- attr(data$XtX, "d")
  }

  return(model)
}

# Get variance of y
get_var_y.ss <- function(data, ...) {
  return(data$yty / (data$n - 1))
}

# Configure ss data for specified method
configure_data.ss <- function(data) {
  if (data$unmappable_effects == "none") {
    return(configure_data.default(data))
  } else {
    return(add_eigen_decomposition(data))
  }
}

# Track core parameters across iterations
track_ibss_fit.ss <- function(data, model, tracking, iter, track_fit, ...) {
  if (data$unmappable_effects %in% c("inf", "ash")) {
    # Append non-sparse variance component to tracking
    tracking <- track_ibss_fit.default(data, model, tracking, iter, track_fit, ...)
    if (isTRUE(track_fit)) {
      tracking[[iter]]$tau2 <- model$tau2
    }
    return(tracking)
  } else {
    # Use default for standard SS case
    return(track_ibss_fit.default(data, model, tracking, iter, track_fit, ...))
  }
}

# Validate Prior Variance
validate_prior.ss <- function(data, model, check_prior, ...) {
  if (isTRUE(check_prior)) {
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
  return(validate_prior.default(data, model, check_prior, ...))
}

# Expected Squared Residuals
get_ER2.ss <- function(data, model) {
  B       <- model$alpha * model$mu
  XB2     <- sum((B %*% data$XtX) * B)
  betabar <- colSums(B)
  postb2  <- model$alpha * model$mu2 # Posterior second moment.

  return(data$yty - 2 * sum(betabar * data$Xty) + sum(betabar * (data$XtX %*% betabar)) -
    XB2 + sum(model$predictor_weights * t(postb2)))
}

# Posterior expected log-likelihood for a single effect regression
SER_posterior_e_loglik.ss <- function(data, model, Eb, Eb2) {
  if (data$unmappable_effects == "none") {
    # Standard SuSiE
    return(-0.5 / model$sigma2 * (-2 * sum(Eb * model$residuals) + sum(model$predictor_weights * as.vector(Eb2))))
  } else {
    # Omega-weighted likelihood
    return(-0.5 * (-2 * sum(Eb * model$residuals) + sum(model$predictor_weights * as.vector(Eb2))))
  }
}

# Compute residuals for single effect regression
compute_residuals.ss <- function(data, model, l, ...) {
  if (!is.null(data$unmappable_effects) && data$unmappable_effects != "none") {
    # Remove lth effect from fitted values
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
compute_ser_statistics.ss <- function(data, model, residual_variance, l, ...) {
  betahat <- (1 / model$predictor_weights) * model$residuals
  shat2 <- residual_variance / model$predictor_weights

  # Optimization parameters
  if (data$unmappable_effects == "none") {
    # Standard SuSiE: optimize on log scale
    optim_init <- log(max(c(betahat^2 - shat2, 1), na.rm = TRUE))
    optim_bounds <- c(-30, 15)
    optim_scale <- "log"
  } else {
    # Unmappable effects: optimize on linear scale
    optim_init <- model$V[l]
    optim_bounds <- c(0, 1)
    optim_scale <- "linear"
  }

  return(list(
    betahat = betahat,
    shat2 = shat2,
    optim_init = optim_init,
    optim_bounds = optim_bounds,
    optim_scale = optim_scale
  ))
}

# Calculate KL divergence
compute_kl.ss <- function(data, model, l) {
  return(compute_kl.default(data, model, l))
}

# Update fitted values
update_fitted_values.ss <- function(data, model, l) {
  if (data$unmappable_effects != "none") {
    model$XtXr <- compute_Xb(data$XtX, colSums(model$alpha * model$mu) + model$theta)
  } else {
    model$XtXr <- model$fitted_without_l + compute_Xb(data$XtX, model$alpha[l, ] * model$mu[l, ])
  }
  return(model)
}

# Get column scale factors
get_scale_factors.ss <- function(data) {
  return(attr(data$XtX, "scaled:scale"))
}

# Get intercept
get_intercept.ss <- function(data, model, ...) {
  return(sum(data$X_colmeans * (colSums(model$alpha * model$mu) / model$X_column_scale_factors)))
}

# Get Fitted Values
get_fitted.ss <- function(data, model, ...) {
  return(NULL)
}

# Get Credible Sets
get_cs.ss <- function(data, model, coverage, min_abs_corr, n_purity) {
  if (is.null(coverage) || is.null(min_abs_corr)) {
    return(NULL)
  }

  if (any(!(diag(data$XtX) %in% c(0, 1)))) {
    Xcorr <- muffled_cov2cor(data$XtX)
  } else {
    Xcorr <- data$XtX
  }

  return(susie_get_cs(model,
                      coverage = coverage,
                      Xcorr = Xcorr,
                      min_abs_corr = min_abs_corr,
                      check_symmetric = FALSE,
                      n_purity = n_purity))
}

# Get Variable Names
get_variable_names.ss <- function(data, model, null_weight) {
  return(assign_names(model, colnames(data$XtX), null_weight, data$p))
}

# Get univariate z-score
get_zscore.ss <- function(data, model, ...) {
  return(get_zscore.default(data, model))
}

# Update variance components for ss data
update_variance_components.ss <- function(data, model, estimate_method = "MLE") {
  if (data$unmappable_effects == "inf") {
    # Calculate omega
    L         <- nrow(model$alpha)
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    omega     <- matrix(rep(omega_res$diagXtOmegaX, L), nrow = L, ncol = data$p,byrow = TRUE) +
                 matrix(rep(1 / model$V, data$p), nrow = L, ncol = data$p, byrow = FALSE)

    # Compute theta for infinitesimal effects.
    theta <- compute_theta_blup(data, model)

    # Sigma2 and tau2 update
    if (estimate_method == "MLE") {
      mle_result <- mle_unmappable(model$alpha, model$mu, omega, model$sigma2, model$tau2, data$n,
        data$eigen_vectors, data$eigen_values, data$VtXty, data$yty,
        est_sigma2 = TRUE, est_tau2 = TRUE,
        verbose = FALSE
      )
      return(list(sigma2 = mle_result$sigma2,
                  tau2   = mle_result$tau2,
                  theta  = theta))
    } else {
      mom_result <- mom_unmappable(model$alpha, model$mu, omega, model$sigma2, model$tau2, data$n,
        data$eigen_vectors, data$eigen_values, data$VtXty, data$Xty, data$yty,
        est_sigma2 = TRUE, est_tau2 = TRUE, verbose = FALSE
      )
      return(list(sigma2 = mom_result$sigma2,
                  tau2   = mom_result$tau2,
                  theta  = theta))
    }
  } else if (data$unmappable_effects == "ash") {
    # Compute omega from current iteration
    L         <- nrow(model$alpha)
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    omega     <- matrix(rep(omega_res$diagXtOmegaX, L), nrow = L, ncol = data$p, byrow = TRUE) +
                 matrix(rep(1 / model$V, data$p), nrow = L, ncol = data$p, byrow = FALSE)

    # Update the sparse effect variance
    sparse_var <- mean(colSums(model$alpha * model$V))

    # Update sigma2 and tau2 using sparse effect variance and omega from current iteration
    mom_result <- mom_unmappable(model$alpha, model$mu, omega,
                                 sigma2 = model$sigma2, tau2 = sparse_var, data$n,
                                 data$eigen_vectors, data$eigen_values, data$VtXty,
                                 data$Xty, data$yty,
                                 est_sigma2 = TRUE, est_tau2 = TRUE, verbose = FALSE)

    # Compute diagXtOmegaX and XtOmega for mr.ash using sparse effect variance and MoM residual variance
    omega_res <- compute_omega_quantities(data, sparse_var, mom_result$sigma2)
    XtOmega <- data$eigen_vectors %*% sweep(data$VtXt, 1, 1/omega_res$omega_var, `*`)

    # Compute variance grid using tau2 estimates from MoM
    est_sa2 <- 100 * mom_result$tau2 * (seq(0, 1, length.out = 10))^2

    # Call mr.ash directly with pre-computed quantities
    mrash_output <- mr.ash.alpha.mccreight::mr.ash(
      X             = data$X,
      y             = data$y,
      sa2           = est_sa2,
      intercept     = FALSE,
      standardize   = FALSE,
      sigma2        = mom_result$sigma2,
      update.sigma2 = FALSE,
      diagXtOmegaX  = omega_res$diagXtOmegaX,
      XtOmega       = XtOmega,
      V             = data$eigen_vectors,
      tausq         = sparse_var,
      sum_Dsq       = sum(data$eigen_values),
      Dsq           = data$eigen_values,
      VtXt          = data$VtXt
    )

    return(list(
      sigma2 = mrash_output$sigma2,
      tau2   = sum(est_sa2 * mrash_output$pi),
      theta  = mrash_output$beta,
      ash_pi = mrash_output$pi
    ))
  } else {
    # For standard SuSiE MLE and MoM are equivalent
    sigma2 <- est_residual_variance(data, model)
    return(list(sigma2 = sigma2))
  }
}

# Update derived quantities for ss data
update_derived_quantities.ss <- function(data, model) {
  if (data$unmappable_effects %in% c("inf", "ash")) {
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
    return(update_derived_quantities.default(data, model))
  }
}

# Expected log-likelihood
Eloglik.ss <- function(data, model) {
  # Standard log-likelihood computation
  return(-data$n / 2 * log(2 * pi * model$sigma2) -
    1 / (2 * model$sigma2) * get_ER2(data, model))
}

#' @importFrom Matrix colSums
#' @importFrom stats dnorm
loglik.ss <- function(data, model, V, ser_stats, ...) {
  # log(bf) for each SNP
  lbf <- dnorm(ser_stats$betahat, 0, sqrt(V + ser_stats$shat2), log = TRUE) -
    dnorm(ser_stats$betahat, 0, sqrt(ser_stats$shat2), log = TRUE)

  # Stabilize logged Bayes Factor
  stable_res  <- lbf_stabilization(lbf, model$pi, ser_stats$shat2)

  # Compute posterior weights
  weights_res <- compute_posterior_weights(stable_res$lpo)

  # Compute gradient
  gradient    <- compute_lbf_gradient(weights_res$alpha, ser_stats$betahat, ser_stats$shat2, V)

  return(list(
    lbf       = stable_res$lbf,
    lbf_model = weights_res$lbf_model,
    alpha     = weights_res$alpha,
    gradient  = gradient
  ))
}

neg_loglik.ss <- function(data, model, V_param, ser_stats, ...) {
  # Convert parameter to V based on optimization scale
  V <- if (ser_stats$optim_scale == "log") exp(V_param) else V_param

  if (data$unmappable_effects == "none") {
    # Standard objective
    res <- loglik.ss(data, model, V, ser_stats)
    return(-res$lbf_model)
  } else {
    # Unmappable objective with logSumExp trick
    return(-matrixStats::logSumExp(
      -0.5 * log(1 + V * model$predictor_weights) +
      V * model$residuals^2 / (2 * (1 + V * model$predictor_weights)) +
        log(model$pi + sqrt(.Machine$double.eps))
    ))
  }
}

# Calculate posterior moments for single effect regression
calculate_posterior_moments.ss <- function(data, model, V,
                                           residual_variance, ...) {
  # Standard Gaussian posterior calculations
  post_var   <- (1 / V + model$predictor_weights / residual_variance)^(-1)
  post_mean  <- (1 / residual_variance) * post_var * model$residuals
  post_mean2 <- post_var + post_mean^2

  return(list(
    post_mean  = post_mean,
    post_mean2 = post_mean2,
    post_var   = post_var
  ))
}
