# S3 method implementations for RSS Lambda data (class: rss_lambda)

# Initialize fitted values
initialize_fitted.rss_lambda <- function(data, alpha, mu) {
  return(list(Rz = compute_Xb(data$R, colSums(alpha * mu))))
}

# Initialize susie model
initialize_susie_model.rss_lambda <- function(data, L, scaled_prior_variance, var_y,
                                              residual_variance, prior_weights, ...) {
  model <- initialize_matrices(data$p, L, scaled_prior_variance, var_y,
    residual_variance, prior_weights,
    include_unmappable = FALSE
  )

  model$lambda <- data$lambda
  model$intercept <- data$intercept_value

  # Initialize SinvRj and RjSinvRj based on initial sigma2
  eigenS_values <- model$sigma2 * data$eigen_R$values + data$lambda
  Dinv <- 1 / eigenS_values
  Dinv[is.infinite(Dinv)] <- 0

  model$SinvRj <- data$eigen_R$vectors %*% (Dinv * data$eigen_R$values * t(data$eigen_R$vectors))

  tmp <- t(data$eigen_R$vectors)
  model$RjSinvRj <- colSums(tmp * (Dinv * (data$eigen_R$values^2) * tmp))

  return(model)
}

# Get variance of y
get_var_y.rss_lambda <- function(data, ...) {
  return(1)
}

# Configure data
configure_data.rss_lambda <- function(data) {
  return(data)
}

# Extract core parameters for tracking
extract_core.rss_lambda <- function(data, model, tracking, iter, track_fit, ...) {
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
validate_prior.rss_lambda <- function(data, model, check_prior, ...) {
  return(model)
}

# Expected squared residuals
get_ER2.rss_lambda <- function(data, model, sigma2 = NULL) {
  if (is.null(sigma2)) {
    sigma2 <- model$sigma2
  }

  d <- sigma2 * data$eigen_R$values + data$lambda
  Dinv <- 1 / d
  Dinv[is.infinite(Dinv)] <- 0

  SinvR <- data$eigen_R$vectors %*%
    ((Dinv * data$eigen_R$values) * t(data$eigen_R$vectors))
  Utz <- crossprod(data$eigen_R$vectors, data$z)
  zSinvz <- sum(Utz * (Dinv * Utz))

  Z <- model$alpha * model$mu
  RSinvR <- data$R %*% SinvR
  RZ2 <- sum((Z %*% RSinvR) * Z)

  zbar <- colSums(Z)
  postb2 <- model$alpha * model$mu2

  return(zSinvz - 2 * sum((SinvR %*% data$z) * zbar) +
    sum(zbar * (RSinvR %*% zbar)) -
    RZ2 + sum(diag(RSinvR) * t(postb2)))
}

# SER posterior expected log-likelihood
SER_posterior_e_loglik.rss_lambda <- function(data, model, r, Eb, Eb2, dXtX) {
  rR <- data$R %*% r

  d <- model$sigma2 * data$eigen_R$values + data$lambda
  Dinv <- 1 / d
  Dinv[is.infinite(Dinv)] <- 0
  SinvEb <- data$eigen_R$vectors %*% (Dinv * crossprod(data$eigen_R$vectors, Eb))

  return(-0.5 * (-2 * sum(rR * SinvEb) + sum(model$RjSinvRj * Eb2)))
}

# Expected log-likelihood
Eloglik.rss_lambda <- function(data, model) {
  d <- model$sigma2 * data$eigen_R$values + data$lambda
  result <- -(length(data$z) / 2) * log(2 * pi) - 0.5 * sum(log(d)) - 0.5 * get_ER2.rss_lambda(data, model)
  return(result)
}

# Log-likelihood for RSS
loglik.rss_lambda <- function(data, model, V, residuals, ser_stats, prior_weights, ...) {
  # Compute log Bayes factors
  lbf <- sapply(1:data$p, function(j) {
    -0.5 * log(1 + (V / ser_stats$shat2[j])) +
      0.5 * (V / (1 + (V / ser_stats$shat2[j]))) * sum(model$SinvRj[, j] * residuals)^2
  })

  # Stabilize logged Bayes Factor
  stable_res <- lbf_stabilization(lbf, prior_weights, ser_stats$shat2)

  # Compute posterior weights
  weights_res <- compute_posterior_weights(stable_res$lpo)

  return(list(
    lbf = stable_res$lbf,
    lbf_model = weights_res$lbf_model,
    alpha = weights_res$alpha,
    gradient = NA # Gradient not implemented for RSS with lambda > 0
  ))
}

neg_loglik.rss_lambda <- function(data, model, V_param, residuals, ser_stats, prior_weights, ...) {
  # Convert parameter to V based on optimization scale (always log for RSS lambda)
  V <- exp(V_param)
  res <- loglik.rss_lambda(data, model, V, residuals, ser_stats, prior_weights)
  return(-res$lbf_model)
}

# Compute residuals for single effect regression
compute_residuals.rss_lambda <- function(data, model, l, ...) {
  # Remove lth effect from fitted values
  Rz_without_l <- model$Rz - data$R %*% (model$alpha[l, ] * model$mu[l, ])

  # Compute residuals
  r <- data$z - Rz_without_l

  return(list(
    z_residual = r,
    Rz_without_l = Rz_without_l
  ))
}

# Compute SER statistics
compute_ser_statistics.rss_lambda <- function(data, model, residuals, dXtX, residual_variance, l, ...) {
  # For RSS-lambda, we only need shat2
  shat2 <- 1 / model$RjSinvRj
  
  # Compute initial value for optimization: max((z^T Sigma^{-1} R_j)^2 - 1/RjSinvRj, 1e-6)
  init_vals <- sapply(1:data$p, function(j) sum(model$SinvRj[, j] * residuals)^2) - (1 / model$RjSinvRj)
  optim_init <- log(max(c(init_vals, 1e-6), na.rm = TRUE))
  optim_bounds <- c(-30, 15)
  optim_scale <- "log"
  
  return(list(
    betahat = NULL,  # Not used for RSS-lambda
    shat2 = shat2,
    dXtX = dXtX,  # Store for consistency
    optim_init = optim_init,
    optim_bounds = optim_bounds,
    optim_scale = optim_scale
  ))
}

# Single effect update
single_effect_update.rss_lambda <- function(data, model, l,
                                            optimize_V, check_null_threshold) {

  # Compute residuals
  residuals <- compute_residuals(data, model, l)

  # Update Rz (removing lth effect)
  model$Rz <- residuals$Rz_without_l

  res <- single_effect_regression(
    data = data,
    model = model,
    l = l,
    residuals = residuals$z_residual,
    dXtX = NULL,  # Not used for RSS-lambda
    optimize_V = optimize_V,
    check_null_threshold = check_null_threshold
  )

  model$mu[l, ] <- res$mu
  model$alpha[l, ] <- res$alpha
  model$mu2[l, ] <- res$mu2
  model$V[l] <- res$V
  model$lbf[l] <- res$lbf_model
  model$lbf_variable[l, ] <- res$lbf

  model$KL[l] <- -res$lbf_model +
    SER_posterior_e_loglik(
      data, model, residuals$z_residual,
      res$alpha * res$mu, res$alpha * res$mu2,
      dXtX = diag(data$R)  # For RSS, diagonal of R matrix
    )

  # Add lth effect back
  model$Rz <- model$Rz + compute_Xb(data$R, model$alpha[l, ] * model$mu[l, ])

  return(model)
}

# Get scale factors
get_scale_factors.rss_lambda <- function(data) {
  return(rep(1, data$p))
}

# Get intercept
get_intercept.rss_lambda <- function(data, model, ...) {
  return(model$intercept)
}

# Get fitted values
get_fitted.rss_lambda <- function(data, model, ...) {
  return(model$Rz)
}

# Get credible sets
# FIXME: data$R has been used previously as residuals but here it is R the LD matrix. Even though these are different data types I still think for clarity we should, throughout the code base, use "R" for correlatin/LD matrix. You can use `yr` for residual or something else, in the individual level data objects.
get_cs.rss_lambda <- function(data, model, coverage, min_abs_corr, n_purity) {
  if (is.null(coverage) || is.null(min_abs_corr)) {
    return(NULL)
  }

  return(susie_get_cs(model,
                      coverage = coverage,
                      Xcorr = muffled_cov2cor(data$R),
                      min_abs_corr = min_abs_corr,
                      check_symmetric = FALSE,
                      n_purity = n_purity))
}

# Get variable names
get_variable_names.rss_lambda <- function(data, model, null_weight) {
  variable_names <- names(data$z)
  return(assign_names(model, variable_names, null_weight, data$p))
}

# Get univariate z-scores
get_zscore.rss_lambda <- function(data, model, ...) {
  return(data$z)
}

# Update variance components
update_variance_components.rss_lambda <- function(data, model, estimate_method = "MLE") {
  # For lambda != 0, use optimization to estimate sigma2
  upper_bound <- 1 - data$lambda

  Eloglik <- function(sigma2) {
    model_temp <- model
    model_temp$sigma2 <- sigma2

    eigenS_values <- sigma2 * data$eigen_R$values + data$lambda
    Dinv <- 1 / eigenS_values
    Dinv[is.infinite(Dinv)] <- 0

    model_temp$SinvRj <- data$eigen_R$vectors %*% (Dinv * data$eigen_R$values * t(data$eigen_R$vectors))
    tmp <- t(data$eigen_R$vectors)
    model_temp$RjSinvRj <- colSums(tmp * (Dinv * (data$eigen_R$values^2) * tmp))

    log_det_Sigma <- sum(log(eigenS_values[eigenS_values > 0]))
    ER2_term <- get_ER2.rss_lambda(data, model_temp, sigma2 = sigma2)

    return(-0.5 * log_det_Sigma - 0.5 * ER2_term)
  }

  opt_result <- optimize(Eloglik, interval = c(1e-4, upper_bound), maximum = TRUE)
  est_sigma2 <- opt_result$maximum

  if (Eloglik(est_sigma2) < Eloglik(upper_bound)) {
    est_sigma2 <- upper_bound
  }

  # Return just the sigma2 value - RSS lambda doesn't use tau2
  return(list(sigma2 = est_sigma2, tau2 = NULL))
}

# Update derived quantities
update_derived_quantities.rss_lambda <- function(data, model) {
  # Calculate updated SinvRj and RjSinvRj based on new sigma2
  eigenS_values <- model$sigma2 * data$eigen_R$values + data$lambda
  Dinv <- 1 / eigenS_values
  Dinv[is.infinite(Dinv)] <- 0

  data$SinvRj_temp <- data$eigen_R$vectors %*% (Dinv * data$eigen_R$values * t(data$eigen_R$vectors))

  tmp <- t(data$eigen_R$vectors)
  data$RjSinvRj_temp <- colSums(tmp * (Dinv * (data$eigen_R$values^2) * tmp))

  return(data)
}


# Calculate posterior moments for single effect regression
calculate_posterior_moments.rss_lambda <- function(data, model, V, residuals, ...) {
  # RSS-lambda specific posterior calculations
  post_var <- (model$RjSinvRj + 1 / V)^(-1)
  post_mean <- sapply(1:data$p, function(j) {
    post_var[j] * sum(model$SinvRj[, j] * residuals)
  })
  post_mean2 <- post_var + post_mean^2

  return(list(
    post_mean = post_mean,
    post_mean2 = post_mean2,
    post_var = post_var,
    beta_1 = NULL  # Not used for RSS-lambda
    # FIXME: when it is not used we should not have to include it as a NULL
    # Also beta_1 as a notation is a bit weird. What is it exactly from Servin_Stephens and can we clarify on why we need it?
  ))
}
