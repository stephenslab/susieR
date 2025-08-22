# Individual-level data backend methods

# Initialize fitted values
initialize_fitted.individual <- function(data, alpha, mu) {
  return(list(Xr = compute_Xb(data$X, colSums(alpha * mu))))
}

# Initialize susie model
initialize_susie_model.individual <- function(data, L, scaled_prior_variance, var_y,
                                              residual_variance, prior_weights, ...) {
  return(initialize_matrices(
    data$p, L, scaled_prior_variance, var_y,
    residual_variance, prior_weights
  ))
}

# Get variance of y
get_var_y.individual <- function(data, ...) {
    return(var(drop(data$y)))
}

# Configure individual data for specified method
configure_data.individual <- function(data) {
  if (data$unmappable_effects == "none") {
    return(data)
  } else {
    warning("Individual-level data converted to sufficient statistics for unmappable effects methods\n")
    return(convert_individual_to_ss_unmappable(data))
  }
}

# Extract core parameters across iterations
extract_core.individual <- function(data, model, tracking, iter, track_fit, ...) {
  if (isTRUE(track_fit)) {
    tracking[[iter]] <- list(
      alpha = model$alpha,
      niter = iter,
      V = model$V,
      sigma2 = model$sigma2
    )
  }
  return(tracking)
}

# Validate prior variance
validate_prior.individual <- function(data, model, check_prior, ...) {
  invisible(TRUE)
}

# Posterior expected log-likelihood for single effect regression
SER_posterior_e_loglik.individual <- function(data, model, R, Eb, Eb2) {
  return(-0.5 * data$n * log(2 * pi * model$sigma2) -
    0.5 / model$sigma2 * (sum(R * R) - 2 * sum(R * compute_Xb(data$X, Eb)) +
      sum(attr(data$X, "d") * Eb2)))
}

# Expected squared residuals
get_ER2.individual <- function(data, model) {
  Xr_L <- compute_MXt(model$alpha * model$mu, data$X)
  postb2 <- model$alpha * model$mu2
  return(sum((data$y - model$Xr)^2) - sum(Xr_L^2) + sum(attr(data$X, "d") * t(postb2)))
}

# Compute residuals for single effect regression
compute_residuals.individual <- function(data, model, l, ...) {
  # Remove lth effect from fitted values
  Xr_without_l <- model$Xr - compute_Xb(data$X, model$alpha[l, ] * model$mu[l, ])

  # Compute residuals
  R <- data$y - Xr_without_l
  XtR <- compute_Xty(data$X, R)

  return(list(
    R = R,
    XtR = XtR,
    Xr_without_l = Xr_without_l
  ))
}

# Compute SER statistics
compute_ser_statistics.individual <- function(data, model, residuals, dXtX, residual_variance, ...) {
  betahat <- (1 / dXtX) * residuals
  shat2 <- residual_variance / dXtX
  
  # Compute initial value for optimization
  optim_init <- log(max(c(betahat^2 - shat2, 1), na.rm = TRUE))
  
  return(list(
    betahat = betahat,
    shat2 = shat2,
    optim_init = optim_init
  ))
}

# Single Effect Update
single_effect_update.individual <- function(
    data, model, l,
    optimize_V, check_null_threshold) {

  # Compute residuals
  residuals <- compute_residuals(data, model, l)

  # Update Xr (removing lth effect)
  model$Xr <- residuals$Xr_without_l

  # Append residual to data object (needed for Servin-Stephens)
  data$R <- residuals$R

  res <- single_effect_regression(
    data                 = data,
    model                = model,
    l                    = l,
    residuals            = residuals$XtR,
    dXtX                 = attr(data$X, "d"),
    optimize_V           = optimize_V,
    check_null_threshold = check_null_threshold,
    unmappable_effects   = FALSE
  )

  # log-likelihood term using current residual vector (not available in ss)
  res$loglik <- res$lbf_model +
    sum(dnorm(residuals$R, 0, sqrt(model$sigma2), log = TRUE))

  res$KL <- -res$loglik +
    SER_posterior_e_loglik(data, model, residuals$R,
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

  # Update Xr (add lth effect back)
  model$Xr <- model$Xr + compute_Xb(data$X, model$alpha[l, ] * model$mu[l, ])

  return(model)
}

# Get column scale factors
get_scale_factors.individual <- function(data) {
  return(attr(data$X, "scaled:scale"))
}

# Get intercept
get_intercept.individual <- function(data, model, intercept) {
  if (intercept) {
    return(data$mean_y - sum(attr(data$X, "scaled:center") *
      (colSums(model$alpha * model$mu) / attr(data$X, "scaled:scale"))))
  } else {
    return(0)
  }
}

# Get Fitted Values
get_fitted.individual <- function(data, model, intercept) {
  if (intercept) {
    fitted <- model$Xr + data$mean_y
  } else {
    fitted <- model$Xr
  }

  fitted <- drop(fitted)
  names(fitted) <- `if`(is.null(names(data$y)), rownames(data$X), names(data$y))

  return(fitted)
}

# Get Credible Sets
get_cs.individual <- function(data, model, coverage, min_abs_corr, n_purity) {
  if (is.null(coverage) || is.null(min_abs_corr)) {
    return(NULL)
  }

  return(susie_get_cs(model,
                      coverage = coverage,
                      X = data$X,
                      min_abs_corr = min_abs_corr,
                      n_purity = n_purity))
}


# Get Variable Names
get_variable_names.individual <- function(data, model, null_weight) {
  variable_names <- colnames(data$X)
  return(assign_names(model, variable_names, null_weight, data$p))
}

# Get univariate z-score
get_zscore.individual <- function(data, model, compute_univariate_zscore,
                                  intercept, standardize, null_weight) {
  if (isFALSE(compute_univariate_zscore)) {
    return(NULL)
  }

  X <- data$X

  if (!is.matrix(X)) {
    warning(
      "Calculation of univariate regression z-scores is not ",
      "implemented specifically for sparse or trend filtering ",
      "matrices, so this step may be slow if the matrix is large; ",
      "to skip this step set compute_univariate_zscore = FALSE"
    )
  }
  if (!is.null(null_weight) && null_weight != 0) {
    X <- X[, 1:(ncol(X) - 1)]
  }

  return(calc_z(X, data$y, center = intercept, scale = standardize))
}

# Update variance components for individual data
update_variance_components.individual <- function(data, model, estimate_method = "MLE") {
  # For standard SuSiE w/ individual data, MLE and MoM are equivalent
  sigma2 <- est_residual_variance(data, model)
  return(list(sigma2 = sigma2, tausq = NULL))
}

# Update derived quantities for individual data
update_derived_quantities.individual <- function(data, model) {
  return(data) # No changes needed for individual data
}

# Expected log-likelihood
Eloglik.individual <- function(data, model) {
  return(-data$n / 2 * log(2 * pi * model$sigma2) -
    1 / (2 * model$sigma2) * get_ER2(data, model))
}

#' @importFrom Matrix colSums
#' @importFrom stats dnorm
loglik.individual <- function(data, model = NULL, V, residuals, ser_stats, prior_weights, ...) {
  # Check if using Servin-Stephens prior
  if (data$use_servin_stephens) {
    # Calculate Servin-Stephens logged Bayes factors
    lbf <- do.call(c, lapply(1:data$p, function(j){
      compute_lbf_servin_stephens(x = data$X[,j],
                                  y = data$R,
                                  s0 = sqrt(V),
                                  alpha0 = data$alpha0,
                                  beta0 = data$beta0)}))

  } else {
    # Standard Gaussian prior log Bayes factors
    lbf <- dnorm(ser_stats$betahat, 0, sqrt(V + ser_stats$shat2), log = TRUE) -
      dnorm(ser_stats$betahat, 0, sqrt(ser_stats$shat2), log = TRUE)
  }

  # Stabilize logged Bayes Factor
  stable_res <- lbf_stabilization(lbf, prior_weights, ser_stats$shat2)

  # Compute posterior weights
  weights_res <- compute_posterior_weights(stable_res$lpo)

  # Compute gradient
  gradient <- compute_lbf_gradient(weights_res$alpha, ser_stats$betahat, ser_stats$shat2, V, data$use_servin_stephens)

  return(list(
    lbf = stable_res$lbf,
    lbf_model = weights_res$lbf_model,
    alpha = weights_res$alpha,
    gradient = gradient
  ))
}

neg_loglik_logscale.individual <- function(data, model = NULL, lV, residuals, ser_stats, prior_weights, ...) {
  res <- loglik.individual(data, model, exp(lV), residuals, ser_stats, prior_weights)
  return(-res$lbf_model)
}

# Calculate posterior moments for single effect regression
calculate_posterior_moments.individual <- function(data, model = NULL, V, residuals,
                                                   dXtX, residual_variance, ...) {
  # Initialize beta_1 as NULL (only used for Servin-Stephens)
  beta_1 <- NULL

  if (data$use_servin_stephens) {
    if (V <= 0) {
      # Zero variance case
      post_mean  <- rep(0, data$p)
      post_mean2 <- rep(0, data$p)
      post_var   <- rep(0, data$p)
      beta_1     <- rep(0, data$p)
    } else {
      # Calculate Servin Stephens Posterior Mean
      post_mean <- sapply(1:data$p, function(j){
        posterior_mean_servin_stephens(dXtX[j], residuals[j], V)
      })

      # Calculate Servin Stephens Posterior Variance
      var_result <- lapply(1:data$p, function(j){
        posterior_var_servin_stephens(dXtX[j], residuals[j],
                                      crossprod(data$R),
                                      data$n, V)
      })

      post_var <- sapply(var_result, function(x) x$post_var)
      beta_1 <- sapply(var_result, function(x) x$beta1)
      post_mean2 <- post_mean^2 + post_var
    }
  } else {
    # Standard Gaussian posterior calculations
    post_var <- (1 / V + dXtX / residual_variance)^(-1)
    post_mean <- (1 / residual_variance) * post_var * residuals
    post_mean2 <- post_var + post_mean^2
  }

  return(list(
    post_mean = post_mean,
    post_mean2 = post_mean2,
    post_var = post_var,
    beta_1 = beta_1
  ))
}
