# Sufficient statistics data backend methods

# Initialize fitted values
initialize_fitted.ss <- function(data, alpha, mu) {
  if (data$unmappable_effects %in% c("inf", "ash")) {
    b <- colSums(alpha * mu)
    XtXr <- data$XtX %*% b
    return(list(XtXr = XtXr))
  } else {
    return(list(XtXr = data$XtX %*% colSums(alpha * mu)))
  }
}

# Initialize susie model
initialize_susie_model.ss <- function(data, L, scaled_prior_variance, var_y,
                                      residual_variance, prior_weights, ...) {
  if (data$unmappable_effects %in% c("inf", "ash")) {
    return(initialize_matrices(data$p, L, scaled_prior_variance, var_y,
      residual_variance, prior_weights,
      include_unmappable = TRUE
    ))
  } else {
    return(initialize_matrices(
      data$p, L, scaled_prior_variance, var_y,
      residual_variance, prior_weights
    ))
  }
}

# Get variance of y
get_var_y.ss <- function(data, ...) {
  return(data$yty / (data$n - 1))
}

# Extract core parameters across iterations
extract_core.ss <- function(data, model, tracking, iter, track_fit, ...) {
  if (isTRUE(track_fit)) {
    if (data$unmappable_effects %in% c("inf", "ash")) {
      tracking_item <- list(
        alpha = model$alpha,
        niter = iter,
        V = model$V,
        sigma2 = model$sigma2,
        tau2 = model$tau2
      )
      # Add ash-specific tracking
      if (data$unmappable_effects == "ash" && !is.null(model$ash_pi)) {
        tracking_item$ash_pi <- model$ash_pi
        tracking_item$theta <- model$theta
      }
      tracking[[iter]] <- tracking_item
    } else {
      tracking[[iter]] <- list(
        alpha = model$alpha,
        niter = iter,
        V = model$V,
        sigma2 = model$sigma2
      )
    }
  }
  return(tracking)
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

# Posterior expected log-likelihood for a single effect regression
SER_posterior_e_loglik.ss <- function(data, model, XtR, Eb, Eb2) {
  return(-0.5 / model$sigma2 * (-2 * sum(Eb * XtR) + sum(attr(data$XtX, "d") * as.vector(Eb2))))
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

# Single Effect Update
single_effect_update.ss <- function(
    data, model, l,
    optimize_V, check_null_threshold) {
  if (data$unmappable_effects %in% c("inf", "ash")) {
    # Infinitesimal / ash single effect update

    # Remove lth effect
    b <- colSums(model$mu * model$alpha) - model$mu[l, ] * model$alpha[l, ]

    # Compute Residuals
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    XtOmegay <- data$eigen_vectors %*% (data$VtXty / omega_res$omega_var)
    XtOmegaXb <- data$eigen_vectors %*% ((t(data$eigen_vectors) %*% b) * data$eigen_values / omega_res$omega_var)
    XtOmegar <- XtOmegay - XtOmegaXb

    res <- single_effect_regression(
      data                 = data,
      Xty                  = XtOmegar,
      dXtX                 = omega_res$diagXtOmegaX,
      V                    = model$V[l],
      residual_variance    = 1, # Already incorporated in Omega
      prior_weights        = model$pi,
      optimize_V           = optimize_V,
      check_null_threshold = check_null_threshold,
      unmappable_effects   = TRUE
    )

    model$alpha[l, ] <- res$alpha
    model$mu[l, ] <- res$mu
    model$mu2[l, ] <- res$mu2
    model$V[l] <- res$V
    model$lbf[l] <- res$lbf_model
    model$lbf_variable[l, ] <- res$lbf

    # TODO: KL and mu2 for infinitesimal and ash models is not properly implemented at the moment.
    model$KL[l] <- -res$lbf_model + SER_posterior_e_loglik(data, model, XtOmegar,
      Eb  = model$alpha[l, ] * model$mu[l, ],
      Eb2 = model$alpha[l, ] * model$mu2[l, ]
    )

    b <- colSums(model$alpha * model$mu)
    model$XtXr <- data$XtX %*% (b + model$theta)

  } else {
    # Ordinary SuSiE single effect update

    # Remove lth effect
    model$XtXr <- model$XtXr - data$XtX %*% (model$alpha[l, ] * model$mu[l, ])

    # Compute Residuals
    XtR <- data$Xty - model$XtXr

    res <- single_effect_regression(
      data                 = data,
      Xty                  = XtR,
      dXtX                 = attr(data$XtX, "d"),
      V                    = model$V[l],
      residual_variance    = model$sigma2,
      prior_weights        = model$pi,
      optimize_V           = optimize_V,
      check_null_threshold = check_null_threshold,
      unmappable_effects   = FALSE
    )

    res$KL <- -res$lbf_model +
      SER_posterior_e_loglik(data, model, XtR,
        Eb  = res$alpha * res$mu,
        Eb2 = res$alpha * res$mu2
      )

    # Update alpha and mu for adding effect back
    model$alpha[l, ] <- res$alpha
    model$mu[l, ] <- res$mu
    model$mu2[l, ] <- res$mu2
    model$V[l] <- res$V
    model$lbf[l] <- res$lbf_model
    model$lbf_variable[l, ] <- res$lbf
    model$KL[l] <- res$KL

    model$XtXr <- model$XtXr + data$XtX %*% (model$alpha[l, ] * model$mu[l, ])
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
    L <- nrow(model$alpha)
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    omega <- matrix(rep(omega_res$diagXtOmegaX, L), nrow = L, ncol = data$p, byrow = TRUE) +
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
                  tau2 = mle_result$tau2,
                  theta = theta))
    } else {
      mom_result <- mom_unmappable(model$alpha, model$mu, omega, model$sigma2, model$tau2, data$n,
        data$eigen_vectors, data$eigen_values, data$VtXty, data$Xty, data$yty,
        est_sigma2 = TRUE, est_tau2 = TRUE, verbose = FALSE
      )
      return(list(sigma2 = mom_result$sigma2,
                  tau2 = mom_result$tau2,
                  theta = theta))
    }
  } else if (data$unmappable_effects == "ash") {
    # Compute omega from current iteration
    L <- nrow(model$alpha)
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    omega <- matrix(rep(omega_res$diagXtOmegaX, L), nrow = L, ncol = data$p, byrow = TRUE) +
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
    # TODO: Use Rcpp for this XtOmega computation
    omega_res <- compute_omega_quantities(data, sparse_var, mom_result$sigma2)
    XtOmega <- data$eigen_vectors %*% sweep(data$VtXt, 1, 1/omega_res$omega_var, `*`)

    # Compute variance grid using tau2 estimates from MoM
    est_sa2 <- 100 * mom_result$tau2 * (seq(0, 1, length.out = 10))^2

    # Call mr.ash directly with pre-computed quantities
    # TODO: We can fix the exact why to call the correction version of mr.ash,
    # but for now I did a direct call to my customized package.
    mrash_output <- mr.ash.alpha.mccreight::mr.ash(
      X = data$X,
      y = data$y,
      sa2 = est_sa2,
      intercept = FALSE,
      standardize = FALSE,
      sigma2 = mom_result$sigma2,
      update.sigma2 = FALSE,
      diagXtOmegaX = omega_res$diagXtOmegaX,
      XtOmega = XtOmega,
      V = data$eigen_vectors,
      tausq = sparse_var,
      sum_Dsq = sum(data$eigen_values),
      Dsq = data$eigen_values,
      VtXt = data$VtXt
    )

    # Extract results from mr.ash
    ash_result <- list(
      sigma2 = mrash_output$sigma2,
      tau2 = sum(est_sa2 * mrash_output$pi),
      theta = mrash_output$beta,
      pi = mrash_output$pi,
      est_sa2 = est_sa2
    )

    # Store theta, mixture weights, and grid in model object
    model$theta <- ash_result$theta
    model$ash_pi <- ash_result$pi
    model$ash_grid <- ash_result$est_sa2

    return(list(
      sigma2 = ash_result$sigma2,
      tau2 = ash_result$tau2,
      theta = ash_result$theta,
      ash_pi = ash_result$pi
    ))
  } else {
    # For standard SuSiE w/ ss data, MLE and MoM are equivalent
    sigma2 <- est_residual_variance(data, model)
    return(list(sigma2 = sigma2, tau2 = NULL))
  }
}

# Update derived quantities for ss data
update_derived_quantities.ss <- function(data, model) {
  if (data$unmappable_effects %in% c("inf", "ash")) {

    # Update var, diagXtOmegaX, and XtOmegay for next iteration.
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    data$var <- omega_res$omega_var
    data$diagXtOmegaX <- omega_res$diagXtOmegaX
    data$XtOmegay <- data$eigen_vectors %*% (data$VtXty / omega_res$omega_var)

    return(data)
  } else {
    return(data)
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
loglik.ss <- function(data, V, betahat, shat2, prior_weights) {
  # log(bf) for each SNP
  lbf <- dnorm(betahat, 0, sqrt(V + shat2), log = TRUE) -
    dnorm(betahat, 0, sqrt(shat2), log = TRUE)
  lpo <- lbf + log(prior_weights + sqrt(.Machine$double.eps))

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] <- 0
  lpo[is.infinite(shat2)] <- 0

  maxlpo <- max(lpo)
  w_weighted <- exp(lpo - maxlpo)
  weighted_sum_w <- sum(w_weighted)
  alpha <- w_weighted / weighted_sum_w

  # Compute gradient
  T2 <- betahat^2 / shat2
  grad_components <- 0.5 * (1 / (V + shat2)) * ((shat2 / (V + shat2)) * T2 - 1)
  grad_components[is.nan(grad_components)] <- 0
  gradient <- sum(alpha * grad_components)

  return(list(
    lbf_model = log(weighted_sum_w) + maxlpo,
    lbf = lbf,
    alpha = alpha,
    gradient = gradient
  ))
}

neg_loglik_logscale.ss <- function(data, lV, betahat, shat2, prior_weights) {
  res <- loglik.ss(data, exp(lV), betahat, shat2, prior_weights)
  return(-res$lbf_model)
}
