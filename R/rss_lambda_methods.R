# S3 method implementations for RSS Lambda data (class: rss_lambda)

# Initialize fitted values
initialize_fitted.rss_lambda <- function(data, alpha, mu) {
  Rz <- as.vector(data$R %*% colSums(alpha * mu))
  return(list(Rz = Rz))
}

# Initialize susie model
initialize_susie_model.rss_lambda <- function(data, L, scaled_prior_variance, var_y,
                                              residual_variance, prior_weights, ...) {
  model <- initialize_matrices(data$p, L, scaled_prior_variance, var_y,
                               residual_variance, prior_weights, include_unmappable = FALSE)

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

# Add eigen decomposition
add_eigen_decomposition.rss_lambda <- function(data) {
  return(data)
}

# Extract core parameters for tracking
extract_core.rss_lambda <- function(data, model, tracking, iter, track_fit, ...) {
  if (isTRUE(track_fit)) {
    tracking[[iter]] <- list(alpha = model$alpha,
                             niter = iter,
                             V = model$V,
                             sigma2 = model$sigma2)
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
  Dinv <- 1/d
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

  return(zSinvz - 2*sum((SinvR %*% data$z) * zbar) +
         sum(zbar * (RSinvR %*% zbar)) -
         RZ2 + sum(diag(RSinvR) * t(postb2)))
}

# SER posterior expected log-likelihood
SER_posterior_e_loglik.rss_lambda <- function(data, model, r, Eb, Eb2) {
  rR <- data$R %*% r

  d <- model$sigma2 * data$eigen_R$values + data$lambda
  Dinv <- 1/d
  Dinv[is.infinite(Dinv)] <- 0
  SinvEb <- data$eigen_R$vectors %*% (Dinv * crossprod(data$eigen_R$vectors, Eb))

  return(-0.5 * (-2 * sum(rR * SinvEb) + sum(model$RjSinvRj * Eb2)))
}

# Single effect update
single_effect_update.rss_lambda <- function(data, model, l,
                                            optimize_V, check_null_threshold) {

  model$Rz <- model$Rz - data$R %*% (model$alpha[l,] * model$mu[l,])

  r <- data$z - model$Rz

  res <- single_effect_regression(
    data = data,
    z = r,
    SinvRj = model$SinvRj,
    RjSinvRj = model$RjSinvRj,
    V = model$V[l],
    prior_weights = model$pi,
    optimize_V = optimize_V,
    check_null_threshold = check_null_threshold
  )

  model$mu[l,] <- res$mu
  model$alpha[l,] <- res$alpha
  model$mu2[l,] <- res$mu2
  model$V[l] <- res$V
  model$lbf[l] <- res$lbf_model
  model$lbf_variable[l,] <- res$lbf

  model$KL[l] <- -res$lbf_model +
    SER_posterior_e_loglik.rss_lambda(data, model, r,
                                      res$alpha * res$mu, res$alpha * res$mu2)

  model$Rz <- model$Rz + data$R %*% (model$alpha[l,] * model$mu[l,])

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
get_cs.rss_lambda <- function(data, model, coverage, min_abs_corr, n_purity) {
  if (is.null(coverage) || is.null(min_abs_corr)) return(NULL)

  return(susie_get_cs(model, coverage = coverage,
                      Xcorr = data$R,
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

# Check convergence
check_convergence.rss_lambda <- function(data, model_prev, model_current,
                                         elbo_prev, elbo_current, tol, convergence_method) {
  if (convergence_method == "pip") {
    # PIP-based convergence
    PIP_diff <- max(abs(model_prev$alpha - model_current$alpha))
    return(PIP_diff < tol)
  } else {
    # Standard ELBO-based convergence
    return((elbo_current - elbo_prev) < tol)
  }
}


# Expected log-likelihood
Eloglik.rss_lambda <- function(data, model) {
  d <- model$sigma2 * data$eigen_R$values + data$lambda
  result <- -(length(data$z)/2)*log(2*pi) - 0.5*sum(log(d)) - 0.5*get_ER2.rss_lambda(data, model)
  return(result)
}

# Log-likelihood for RSS
loglik.rss_lambda <- function(data, V, z, SinvRj, RjSinvRj, shat2, prior_weights) {
  p <- length(z)
  
  # Compute log Bayes factors
  lbf <- sapply(1:p, function(j)
    -0.5 * log(1 + (V/shat2[j])) +
     0.5 * (V/(1 + (V/shat2[j]))) * sum(SinvRj[,j] * z)^2
  )
  
  # Add log prior weights
  lpo <- lbf + log(prior_weights + sqrt(.Machine$double.eps))
  
  # Deal with special case of infinite shat2 (e.g., happens if X does not vary)
  lbf[is.infinite(shat2)] <- 0
  lpo[is.infinite(shat2)] <- 0
  
  # Compute log-sum-exp of weighted lbf
  maxlpo <- max(lpo)
  w_weighted <- exp(lpo - maxlpo)
  weighted_sum_w <- sum(w_weighted)
  alpha <- w_weighted / weighted_sum_w
  
  return(list(
    lbf_model = log(weighted_sum_w) + maxlpo,
    lbf = lbf,
    alpha = alpha,
    gradient = NA  # Gradient not computed for RSS
  ))
}

neg_loglik_logscale.rss_lambda <- function(data, lV, z, SinvRj, RjSinvRj, shat2, prior_weights) {
  res <- loglik.rss_lambda(data, exp(lV), z, SinvRj, RjSinvRj, shat2, prior_weights)
  return(-res$lbf_model)
}