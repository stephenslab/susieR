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

# FIXME: we don't have "configure_data" function here? if the only place we need for configre_data.individal then it is simply converting individual level data to suff stats which should not happen as a generic method. It should be done outside of the methods and we caution that it is temporary (due to lack of implementation of individual level updates for the unmappable effects)

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
      d <- attr(data$XtX, "d")
      bhat <- data$Xty / d
      shat <- sqrt(model$sigma2 / d)
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
  invisible(TRUE)
}

# Expected Squared Residuals
get_ER2.ss <- function(data, model) {
  B <- model$alpha * model$mu
  XB2 <- sum((B %*% data$XtX) * B)
  betabar <- colSums(B)
  d <- attr(data$XtX, "d")
  postb2 <- model$alpha * model$mu2 # Posterior second moment.
  return(data$yty - 2 * sum(betabar * data$Xty) + sum(betabar * (data$XtX %*% betabar)) -
    XB2 + sum(d * t(postb2)))
}

# Posterior expected log-likelihood for a single effect regression
SER_posterior_e_loglik.ss <- function(data, model, XtR, Eb, Eb2) {
  if (data$unmappable_effects == "none") {
    # Standard SuSiE
    return(-0.5 / model$sigma2 * (-2 * sum(Eb * XtR) + sum(model$predictor_weights * as.vector(Eb2))))
  } else {
    # Omega-weighted likelihood
    return(-0.5 * (-2 * sum(Eb * XtR) + sum(model$predictor_weights * as.vector(Eb2))))
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

    # Store unified residuals in model (unmappable case)
    model$residuals <- XtOmegar             # For SER & KL
    model$omega_res <- omega_res            # For predictor_weights update
    model$unmappable <- TRUE                # Flag for special handling
    
    return(model)
  } else {
    # Remove lth effect from fitted values
    XtXr_without_l <- model$XtXr - data$XtX %*% (model$alpha[l, ] * model$mu[l, ])

    # Compute Residuals
    XtR <- data$Xty - XtXr_without_l

    # Store unified residuals in model (standard case)
    model$residuals <- XtR                  # For SER & KL
    model$fitted_without_l <- XtXr_without_l # For fitted update
    model$unmappable <- FALSE               # Flag for special handling

    return(model)
  }
}

# Compute SER statistics
compute_ser_statistics.ss <- function(data, model, residual_variance, l, ...) {
  betahat <- (1 / model$predictor_weights) * model$residuals
  shat2 <- residual_variance / model$predictor_weights

  # Compute optimization parameters based on unmappable effects
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

# Single Effect Update
single_effect_update.ss <- function(
    data, model, l,
    optimize_V, check_null_threshold) {

  # Compute residuals and store in model
  model <- compute_residuals(data, model, l)

  # Update model predictor_weights and residual variance based on unmappable status
  if (model$unmappable) {
    model$predictor_weights <- model$omega_res$diagXtOmegaX
    residual_variance <- 1  # Already incorporated in Omega
  } else {
    residual_variance <- model$sigma2
  }

  # Run single effect regression
  res <- single_effect_regression(
    data                 = data,
    model                = model,
    l                    = l,
    residual_variance    = residual_variance,
    optimize_V           = optimize_V,
    check_null_threshold = check_null_threshold
  )

  # Store results
  model$alpha[l, ]        <- res$alpha
  model$mu[l, ]           <- res$mu
  model$mu2[l, ]          <- res$mu2
  model$V[l]              <- res$V
  model$lbf[l]            <- res$lbf_model
  model$lbf_variable[l, ] <- res$lbf

  # Compute KL divergence
  model$KL[l] <- -res$lbf_model +
    SER_posterior_e_loglik(data, model, model$residuals,
                           Eb  = res$alpha * res$mu,
                           Eb2 = res$alpha * res$mu2)

  # Update fitted values
  if (model$unmappable) {
    b <- colSums(model$alpha * model$mu)
    model$XtXr <- compute_Xb(data$XtX, b + model$theta)
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
  variable_names <- colnames(data$XtX)
  return(assign_names(model, variable_names, null_weight, data$p))
}

# Get univariate z-score
# FIXME: again for what's returned as NULL do we still have to define them or we can put in the default behavior with generic method?
get_zscore.ss <- function(data, model, ...) {
  return(NULL)
}

# Configure ss data for specified method
configure_data.ss <- function(data) {
  if (data$unmappable_effects == "none") {
    return(data)
  } else {
    return(add_eigen_decomposition(data))
  }
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
    return(model)
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
loglik.ss <- function(data, model, V, ser_stats, prior_weights, ...) {
  # log(bf) for each SNP
  lbf <- dnorm(ser_stats$betahat, 0, sqrt(V + ser_stats$shat2), log = TRUE) -
    dnorm(ser_stats$betahat, 0, sqrt(ser_stats$shat2), log = TRUE)

  # Stabilize logged Bayes Factor
  stable_res <- lbf_stabilization(lbf, prior_weights, ser_stats$shat2)

  # Compute posterior weights
  weights_res <- compute_posterior_weights(stable_res$lpo)

  # Compute gradient
  gradient <- compute_lbf_gradient(weights_res$alpha, ser_stats$betahat, ser_stats$shat2, V)

  return(list(
    lbf = stable_res$lbf,
    lbf_model = weights_res$lbf_model,
    alpha = weights_res$alpha,
    gradient = gradient
  ))
}

neg_loglik.ss <- function(data, model, V_param, ser_stats, prior_weights, ...) {
  # Convert parameter to V based on optimization scale
  V <- if (ser_stats$optim_scale == "log") exp(V_param) else V_param

  if (data$unmappable_effects == "none") {
    # Standard objective
    res <- loglik.ss(data, model, V, ser_stats, prior_weights)
    return(-res$lbf_model)
  } else {
    # Unmappable objective with logSumExp trick
    return(-matrixStats::logSumExp(
      -0.5 * log(1 + V * model$predictor_weights) +
      V * model$residuals^2 / (2 * (1 + V * model$predictor_weights)) +
        log(prior_weights + sqrt(.Machine$double.eps))
    ))
  }
}

# Calculate posterior moments for single effect regression
calculate_posterior_moments.ss <- function(data, model, V,
                                           residual_variance, ...) {
  # Standard Gaussian posterior calculations
  post_var <- (1 / V + model$predictor_weights / residual_variance)^(-1)
  post_mean <- (1 / residual_variance) * post_var * model$residuals
  post_mean2 <- post_var + post_mean^2

  return(list(
    post_mean = post_mean,
    post_mean2 = post_mean2,
    post_var = post_var
  ))
}
