# Sufficient statistics data backend methods

# Initialize fitted values
initialize_fitted.ss <- function(data, alpha, mu) {
  if (data$unmappable_effects == "inf") {
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
  if (data$unmappable_effects != "none") {
    return(initialize_matrices(data$p, L, scaled_prior_variance, var_y,
                                  residual_variance, prior_weights, include_unmappable = TRUE))
  } else {
    return(initialize_matrices(data$p, L, scaled_prior_variance, var_y,
                                  residual_variance, prior_weights))
  }
}

# Get variance of y
get_var_y.ss <- function(data, ...) {
  return(data$yty / (data$n - 1))
}

# Extract core parameters across iterations
extract_core.ss <- function(data, model, tracking, iter, track_fit, ...) {
  if (isTRUE(track_fit)) {
    if (data$unmappable_effects != "none") {
      tracking[[iter]] <- list(alpha = model$alpha,
                               niter = iter,
                               V = model$V,
                               sigma2 = model$sigma2,
                               tau2 = model$tau2)
    } else {
      tracking[[iter]] <- list(alpha = model$alpha,
                               niter = iter,
                               V = model$V,
                               sigma2 = model$sigma2)
    }
  }
  return(tracking)
}

# Validate Prior Variance
validate_prior.ss <- function(data, model, check_prior, ...) {
  if (isTRUE(check_prior)) {
    if (is.null(data$zm)) {
      d    <- attr(data$XtX, "d")
      bhat <- data$Xty / d
      shat <- sqrt(model$sigma2 / d)
      z    <- bhat / shat
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
SER_posterior_e_loglik.ss <- function (data, model, XtR, Eb, Eb2)
  return(-0.5 / model$sigma2 * (-2 * sum(Eb * XtR) + sum(attr(data$XtX, "d") * as.vector(Eb2))))

# Expected Squared Residuals
get_ER2.ss <- function (data, model) {
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

  if (data$unmappable_effects == "inf") {
    # Unmappable effects implementation
    alpha <- model$alpha
    mu <- model$mu

    V <- data$eigen_vectors
    Dsq <- data$eigen_values
    VtXty <- data$VtXty
    Xty <- data$Xty
    n <- data$n
    p <- data$p
    L <- nrow(alpha)

    sigma2 <- model$sigma2
    tau2 <- if (is.null(model$tau2)) 0 else model$tau2

    if (!is.null(data$diagXtOmegaX) && !is.null(data$XtOmegay)) {
      diagXtOmegaX <- data$diagXtOmegaX
      XtOmegay <- data$XtOmegay
    } else {
      var <- tau2 * Dsq + sigma2
      diagXtOmegaX <- rowSums(sweep(V^2, 2, (Dsq / var), `*`))
      XtOmegay <- V %*% (VtXty / var)
    }

    # Remove lth effect
    b <- colSums(mu * alpha) - mu[l, ] * alpha[l, ]

    # Compute Residuals
    XtOmegaXb <- V %*% ((t(V) %*% b) * Dsq / (tau2 * Dsq + sigma2))
    XtOmegar <- XtOmegay - XtOmegaXb

    res <- single_effect_regression(
      data                 = data,
      Xty                  = XtOmegar,
      dXtX                 = diagXtOmegaX,
      V                    = model$V[l],
      residual_variance    = 1,  # Already incorporated in Omega
      prior_weights        = model$pi,
      optimize_V           = optimize_V,
      check_null_threshold = check_null_threshold,
      unmappable_effects   = TRUE
    )

    model$alpha[l,]         <- res$alpha
    model$mu[l, ]           <- res$mu
    model$mu2[l, ]          <- res$mu2
    model$V[l]              <- res$V
    model$lbf[l]            <- res$lbf_model
    model$lbf_variable[l, ] <- res$lbf

    # TODO: KL and mu2 for infinitesimal model is not properly implemented at the moment.
    model$KL[l] <- -res$lbf_model + SER_posterior_e_loglik(data, model, XtOmegar,
                                                            Eb  = model$alpha[l, ] * model$mu[l, ],
                                                            Eb2 = model$alpha[l, ] * model$mu2[l, ])

    b <- colSums(model$alpha * model$mu)
    if (!is.null(model$theta)) {
      model$XtXr <- data$XtX %*% (b + model$theta)
    } else {
      model$XtXr <- data$XtX %*% b
    }

  } else {
    # Standard approach
    # Remove lth effect
    model$XtXr <- model$XtXr - data$XtX %*% (model$alpha[l, ] * model$mu[l, ])

    # Compute Residuals
    XtR <- data$Xty - model$XtXr
    d <- attr(data$XtX, "d")

    res <- single_effect_regression(
      data                 = data,
      Xty                  = XtR,
      dXtX                 = d,
      V                    = model$V[l],
      residual_variance    = model$sigma2,
      prior_weights        = model$pi,
      optimize_V           = optimize_V,
      check_null_threshold = check_null_threshold,
      unmappable_effects   = FALSE)

    res$KL <- -res$lbf_model +
      SER_posterior_e_loglik(data, model, XtR,
                             Eb  = res$alpha * res$mu,
                             Eb2 = res$alpha * res$mu2)

    # Update alpha and mu for adding effect back
    model$alpha[l,]         <- res$alpha
    model$mu[l, ]           <- res$mu
    model$mu2[l, ]          <- res$mu2
    model$V[l]              <- res$V
    model$lbf[l]            <- res$lbf_model
    model$lbf_variable[l, ] <- res$lbf
    model$KL[l]             <- res$KL

    model$XtXr <- model$XtXr + data$XtX %*% (model$alpha[l, ] * model$mu[l, ])
  }

  return(model)
}

# Get column scale factors
get_scale_factors.ss <- function(data){
  return(attr(data$XtX,"scaled:scale"))
}

# Get intercept
get_intercept.ss <- function(data, model, ...){
  return(sum(data$X_colmeans * (colSums(model$alpha * model$mu)/model$X_column_scale_factors)))
}

# Get Fitted Values
get_fitted.ss <- function(data, model, ...){
  return(NULL)
}

# Get Credible Sets
get_cs.ss <- function(data, model, coverage, min_abs_corr, n_purity){

  if (is.null(coverage) || is.null(min_abs_corr)) return(NULL)

  if (data$unmappable_effects == "inf" && !is.null(data$eigen_vectors)) {
    # For unmappable effects, construct correlation matrix from eigen decomposition
    LD <- (data$eigen_vectors %*% diag(data$eigen_values)) %*% t(data$eigen_vectors) / data$n
    if (any(!(diag(LD) %in% c(0,1)))) {
      Xcorr <- muffled_cov2cor(LD)
    } else {
      Xcorr <- LD
    }
    return(susie_get_cs(model, coverage = coverage,
                        Xcorr = Xcorr,
                        min_abs_corr = min_abs_corr,
                        check_symmetric = FALSE,
                        n_purity = n_purity))
  } else {
    # Standard approach using XtX
    if(any(!(diag(data$XtX) %in% c(0,1)))){
      return(susie_get_cs(model,coverage = coverage,
                          Xcorr = muffled_cov2cor(data$XtX),
                          min_abs_corr = min_abs_corr,
                          check_symmetric = FALSE,
                          n_purity = n_purity))
    }else{
      return(susie_get_cs(model,coverage = coverage,
                          Xcorr = data$XtX,
                          min_abs_corr = min_abs_corr,
                          check_symmetric = FALSE,
                          n_purity = n_purity))
    }
  }
}


# Get Variable Names
get_variable_names.ss <- function(data, model, null_weight){
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
    return(data)  # No changes needed
  } else {
    # Add eigen decomposition for unmappable effects
    data <- add_eigen_decomposition(data)
    return(data)
  }
}

# Add eigen decomposition to ss objects (for unmappable effects methods)
add_eigen_decomposition.ss <- function(data) {
  # Compute eigen decomposition of correlation matrix
  eigen_decomp <- compute_eigen_decomposition(data$XtX, data$n)

  # Add eigen components to data object
  data$eigen_vectors <- eigen_decomp$V
  data$eigen_values  <- eigen_decomp$Dsq
  data$VtXty         <- t(eigen_decomp$V) %*% data$Xty  # Compute VtXty

  # Initialize derived quantities for unmappable effects methods
  # These will be updated when variance components change
  sigmasq <- 1  # Default initial value
  tausq <- 0    # Default initial value
  var <- tausq * data$eigen_values + sigmasq
  data$var <- var
  data$diagXtOmegaX <- rowSums(sweep(data$eigen_vectors^2, 2, (data$eigen_values / var), `*`))
  data$XtOmegay <- data$eigen_vectors %*% (data$VtXty / var)

  return(data)
}

# Update variance components for ss data
update_variance_components.ss <- function(data, model) {
  if (data$unmappable_effects == "inf") {
    # Method of Moments for unmappable effects
    alpha <- model$alpha
    mu <- model$mu
    p <- data$p
    L <- nrow(alpha)

    sigma2 <- model$sigma2
    tau2 <- if (is.null(model$tau2)) 0 else model$tau2

    var <- tau2 * data$eigen_values + sigma2
    diagXtOmegaX <- rowSums(sweep(data$eigen_vectors^2, 2, (data$eigen_values / var), `*`))
    omega <- matrix(rep(diagXtOmegaX, L), nrow = L, ncol = p, byrow = TRUE) +
             matrix(rep(1 / model$V, p), nrow = L, ncol = p, byrow = FALSE)

    mom_result <- MoM(alpha, mu, omega, sigma2, tau2, data$n,
                      data$eigen_vectors, data$eigen_values, data$VtXty,
                      data$Xty, data$yty,
                      est_sigma2 = TRUE, est_tau2 = TRUE, verbose = FALSE)

    return(list(sigma2 = mom_result$sigma2, tau2 = mom_result$tau2))
  } else {
    # Standard approach
    sigma2 <- est_residual_variance(data, model)
    return(list(sigma2 = sigma2, tausq = NULL))
  }
}

# Update derived quantities for ss data
update_derived_quantities.ss <- function(data, model) {
  if (data$unmappable_effects == "inf") {
    sigma2 <- model$sigma2
    tau2 <- if (is.null(model$tau2)) 0 else model$tau2

    var <- tau2 * data$eigen_values + sigma2
    data$var <- var
    data$diagXtOmegaX <- rowSums(sweep(data$eigen_vectors^2, 2, (data$eigen_values / var), `*`))
    data$XtOmegay <- data$eigen_vectors %*% (data$VtXty / var)

    data$theta <- compute_theta_blup(data, model)

    return(data)
  } else {
    return(data)  # No changes needed for standard ss data
  }
}

# Check convergence for ss data
check_convergence.ss <- function(data, model_prev, model_current, elbo_prev, elbo_current, tol) {
  if (data$unmappable_effects == "inf") {
    # PIP-based convergence for unmappable effects
    PIP_diff <- max(abs(model_prev$alpha - model_current$alpha))
    return(PIP_diff < tol)
  } else {
    # Standard ELBO-based convergence (uses pre-computed ELBO values)
    return(elbo_current - elbo_prev < tol)
  }
}

# Update variance before convergence check for ss data
update_variance_before_convergence.ss <- function(data) {
  if (data$unmappable_effects == "inf") {
    return(TRUE)  # Update variance before convergence check
  } else {
    # Standard behavior: Check convergence first, then update variance
    return(FALSE)
  }
}

# Handle convergence and variance updates for ss data
handle_convergence_and_variance.ss <- function(data, model, model_prev, elbo_prev, elbo_current,
                                                tol, estimate_residual_variance,
                                                residual_variance_lowerbound, residual_variance_upperbound) {
  if (data$unmappable_effects == "inf") {
    # Unmappable effects: Update variance first, then check convergence
    if (estimate_residual_variance) {
      result <- update_model_variance(data, model, residual_variance_lowerbound, residual_variance_upperbound)
      data <- result$data
      model <- result$model
    }

    converged <- check_convergence(data, model_prev, model, elbo_prev, elbo_current, tol)
  } else {
    # Standard: Check convergence first, then update variance
    converged <- check_convergence(data, model_prev, model, elbo_prev, elbo_current, tol)

    if (!converged && estimate_residual_variance) {
      result <- update_model_variance(data, model, residual_variance_lowerbound, residual_variance_upperbound)
      data <- result$data
      model <- result$model
    }
  }

  return(list(data = data, model = model, converged = converged))
}

# Expected log-likelihood
Eloglik.ss <- function(data, model) {
  return(-data$n / 2 * log(2 * pi * model$sigma2) -
         1 / (2 * model$sigma2) * get_ER2(data, model))
}
