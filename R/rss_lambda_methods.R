# S3 method implementations for RSS Lambda data (class: rss_lambda)

# Initialize fitted values
initialize_fitted.rss_lambda <- function(data, alpha, mu) {
  return(list(Rz = as.vector(data$R %*% colSums(alpha * mu))))
}

# Initialize SuSiE model
initialize_susie_model.rss_lambda <- function(data, L, scaled_prior_variance, var_y,
                                              residual_variance, prior_weights, ...) {
  # Base model
  model <- initialize_matrices(data, L, scaled_prior_variance, var_y,
                               residual_variance, prior_weights)

  # Initialize SinvRj and RjSinvRj
  D    <- data$eigen_R$values
  V    <- data$eigen_R$vectors
  Dinv <- compute_Dinv(model, data)

  model$SinvRj   <- V %*% (Dinv * D * t(V))
  model$RjSinvRj <- colSums(t(V) * (Dinv * D^2 * t(V)))

  return(model)
}

# Get variance of y
get_var_y.rss_lambda <- function(data, ...) {
  return(1)
}

# Configure data
configure_data.rss_lambda <- function(data) {
  return(configure_data.default(data))
}

# Track core parameters for tracking
track_ibss_fit.rss_lambda <- function(data, model, tracking, iter, track_fit, ...) {
  return(track_ibss_fit.default(data, model, tracking, iter, track_fit, ...))
}

# Validate prior variance
validate_prior.rss_lambda <- function(data, model, check_prior, ...) {
  return(validate_prior.default(data, model, check_prior, ...))
}

get_ER2.rss_lambda <- function(data, model) {
  # Eigen decomposition components
  D     <- data$eigen_R$values
  V     <- data$eigen_R$vectors
  Dinv  <- compute_Dinv(model, data)

  # Cached quantities
  Vtz   <- data$Vtz
  zbar  <- model$zbar
  postb2 <- model$diag_postb2

  # z^T S^{-1} z
  zSinvz <- sum((Dinv * Vtz) * Vtz)

  # -2 zbar^T S^{-1} z
  tmp <- V %*% (Dinv * (D * Vtz))
  term2 <- -2 * sum(tmp * zbar)

  # zbar^T R S^{-1} R zbar
  Vtzbar <- crossprod(V, zbar)
  term3 <- sum((Vtzbar^2) * (Dinv * D^2))

  # RZ2 = sum((Z %*% RSinvR) * Z)
  VtZ <- model$Z %*% V
  term4 <- sum((VtZ^2) %*% (Dinv * D^2))

  # diag(RSinvR)^T postb2
  diag_RSinvR <- rowSums((V^2) * rep(Dinv * D^2, each = nrow(V)))
  term5 <- sum(diag_RSinvR * postb2)

  return(zSinvz + term2 + term3 - term4 + term5)
}

# SER posterior expected log-likelihood
SER_posterior_e_loglik.rss_lambda <- function(data, model, Eb, Eb2) {
  V      <- data$eigen_R$vectors
  Dinv   <- compute_Dinv(model, data)
  rR     <- data$R %*% model$residuals
  SinvEb <- V %*% (Dinv * crossprod(V, Eb))

  return(-0.5 * (-2 * sum(rR * SinvEb) + sum(model$RjSinvRj * Eb2)))
}

# Expected log-likelihood
Eloglik.rss_lambda <- function(data, model) {
  D <- data$eigen_R$values
  d <- model$sigma2 * D + data$lambda
  return(-(length(data$z) / 2) * log(2 * pi) - 0.5 *
           sum(log(d)) - 0.5 * get_ER2.rss_lambda(data, model))
}

# Log-likelihood for RSS
loglik.rss_lambda <- function(data, model, V, ser_stats, ...) {
  # Compute log Bayes factors
  lbf <- -0.5 * log(1 + V / ser_stats$shat2) +
         0.5 * (V / (1 + V / ser_stats$shat2)) * (crossprod(model$SinvRj, model$residuals)^2)

  # Stabilize logged Bayes Factor
  stable_res <- lbf_stabilization(lbf, model$pi, ser_stats$shat2)

  # Compute posterior weights
  weights_res <- compute_posterior_weights(stable_res$lpo)

  return(list(
    lbf       = stable_res$lbf,
    lbf_model = weights_res$lbf_model,
    alpha     = weights_res$alpha
  ))
}

neg_loglik.rss_lambda <- function(data, model, V_param, ser_stats, ...) {
  # Convert parameter to V based on optimization scale (always log for RSS lambda)
  V   <- exp(V_param)
  res <- loglik.rss_lambda(data, model, V, ser_stats)
  return(-res$lbf_model)
}

# Compute residuals for single effect regression
compute_residuals.rss_lambda <- function(data, model, l, ...) {
  # Remove lth effect from fitted values
  Rz_without_l <- model$Rz - data$R %*% (model$alpha[l, ] * model$mu[l, ])

  # Compute residuals
  r <- data$z - Rz_without_l

  # Store unified residuals in model
  model$residuals         <- r
  model$fitted_without_l  <- Rz_without_l
  model$residual_variance <- 1  # RSS lambda uses normalized residual variance

  return(model)
}

# Compute SER statistics
compute_ser_statistics.rss_lambda <- function(data, model, residual_variance, l, ...) {
  shat2 <- 1 / model$RjSinvRj

  # Optimization parameters
  init_vals <- sapply(1:data$p, function(j) sum(model$SinvRj[, j] * model$residuals)^2) - (1 / model$RjSinvRj)
  optim_init <- log(max(c(init_vals, 1e-6), na.rm = TRUE))
  optim_bounds <- c(-30, 15)
  optim_scale <- "log"

  return(list(
    shat2 = shat2,
    optim_init = optim_init,
    optim_bounds = optim_bounds,
    optim_scale = optim_scale
  ))
}

# Calculate KL divergence
compute_kl.rss_lambda <- function(data, model, l) {
  return(compute_kl.default(data, model, l))
}

# Update fitted values
update_fitted_values.rss_lambda <- function(data, model, l) {
  model$Rz <- model$fitted_without_l + as.vector(data$R %*% (model$alpha[l, ] * model$mu[l, ]))
  model    <- precompute_rss_lambda_terms(data, model)

  return(model)
}

# Get scale factors
get_scale_factors.rss_lambda <- function(data) {
  return(rep(1, data$p))
}

# Get intercept
get_intercept.rss_lambda <- function(data, model, ...) {
  return(data$intercept_value)
}

# Get fitted values
get_fitted.rss_lambda <- function(data, model, ...) {
  return(model$Rz)
}

# Get credible sets
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
  return(get_zscore.default(data, model))
}

update_variance_components.rss_lambda <- function(data, model, estimate_method = "MLE") {
  upper_bound <- 1 - data$lambda

  objective <- function(sigma2) {
    temp_model        <- model
    temp_model$sigma2 <- sigma2
    Eloglik.rss_lambda(data, temp_model)
  }

  est_sigma2 <- optimize(objective, interval = c(1e-4, upper_bound), maximum = TRUE)$maximum

  if (objective(est_sigma2) < objective(upper_bound)) {
    est_sigma2 <- upper_bound
  }

  list(sigma2 = est_sigma2)
}

# Update derived quantities
update_derived_quantities.rss_lambda <- function(data, model) {
  # Recalculate Dinv with updated sigma2
  Dinv <- compute_Dinv(model, data)
  V <- data$eigen_R$vectors
  D <- data$eigen_R$values

  # Update SinvRj and RjSinvRj
  model$SinvRj <- V %*% (Dinv * D * t(V))
  model$RjSinvRj <- colSums(t(V) * (Dinv * (D^2) * t(V)))

  return(model)
}


# Calculate posterior moments for single effect regression
calculate_posterior_moments.rss_lambda <- function(data, model, V, ...) {
  post_var   <- (model$RjSinvRj + 1 / V)^(-1)
  post_mean  <- sapply(1:data$p, function(j) {
    post_var[j] * sum(model$SinvRj[, j] * model$residuals)
  })
  post_mean2 <- post_var + post_mean^2

  return(list(
    post_mean  = post_mean,
    post_mean2 = post_mean2,
    post_var   = post_var
  ))
}
