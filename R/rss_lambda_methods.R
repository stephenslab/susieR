# S3 method implementations for RSS Lambda data (class: rss_lambda)
# This file contains methods for handling RSS with correlated errors (Sigma = sigma2*R + lambda*I)

# Initialize fitted values
initialize_fitted.rss_lambda <- function(data, alpha, mu) {
  # Rz = R %*% colSums(alpha * mu)
  Rz <- as.vector(data$R %*% colSums(alpha * mu))
  return(list(Rz = Rz))
}

# Initialize susie model
initialize_susie_model.rss_lambda <- function(data, L, scaled_prior_variance, var_y,
                                              residual_variance, prior_weights, ...) {
  # Use standard initialization
  model <- initialize_matrices(data$p, L, scaled_prior_variance, var_y,
                               residual_variance, prior_weights, include_unmappable = FALSE)
  
  # Add RSS-lambda specific fields
  model$lambda <- data$lambda
  model$intercept <- data$intercept_value
  
  # Initialize placeholders for Sigma-related matrices
  model$SinvRj <- NULL
  model$RjSinvRj <- NULL
  
  # Compute initial Sigma-related matrices
  model <- update_variance_components.rss_lambda(data, model)
  
  return(model)
}

# Get variance of y
get_var_y.rss_lambda <- function(data, ...) {
  # For RSS, assume standardized scale
  return(1)
}

# Configure data
configure_data.rss_lambda <- function(data) {
  # RSS lambda data is already configured
  return(data)
}

# Add eigen decomposition (ghost function - not used for rss_lambda)
add_eigen_decomposition.rss_lambda <- function(data) {
  # Not used for RSS lambda - eigen decomposition already done in constructor
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
  # No special validation needed for RSS lambda
  return(model)
}

# Expected squared residuals
get_ER2.rss_lambda <- function(data, model, sigma2 = NULL) {
  # Following the original get_ER2_rss implementation from elbo_rss.R
  # Allow sigma2 to be specified for residual variance estimation
  if (is.null(sigma2)) {
    sigma2 <- model$sigma2
  }
  
  d <- sigma2 * data$eigen_R$values + data$lambda
  Dinv <- 1/d
  Dinv[is.infinite(Dinv)] <- 0
  
  # Compute key quantities using eigen decomposition
  # SinvR = U * diag(Dinv * eigenvalues) * U^T
  SinvR <- data$eigen_R$vectors %*%
          ((Dinv * data$eigen_R$values) * t(data$eigen_R$vectors))
  Utz <- crossprod(data$eigen_R$vectors, data$z)
  zSinvz <- sum(Utz * (Dinv * Utz))
  
  # Aggregate effects
  Z <- model$alpha * model$mu
  if (data$lambda == 0) {
    RSinvR <- data$R / sigma2
  } else {
    RSinvR <- data$R %*% SinvR
  }
  RZ2 <- sum((Z %*% RSinvR) * Z)
  
  zbar <- colSums(Z)
  postb2 <- model$alpha * model$mu2  # Posterior second moment
  
  return(zSinvz - 2*sum((SinvR %*% data$z) * zbar) +
         sum(zbar * (RSinvR %*% zbar)) -
         RZ2 + sum(diag(RSinvR) * t(postb2)))
}

# SER posterior expected log-likelihood
SER_posterior_e_loglik.rss_lambda <- function(data, model, r, Eb, Eb2) {
  # RSS-specific computation from original elbo_rss.R
  # E[log p(r | b)] propto -0.5 * E[(r - R*b)' * Sigma^(-1) * (r - R*b)]
  # Only need the parts that depend on b
  # r is the residuals, Eb is alpha*mu, Eb2 is alpha*mu2
  
  # Compute R %*% r
  rR <- data$R %*% r
  
  # Compute Sigma^(-1) %*% Eb using eigen decomposition
  d <- model$sigma2 * data$eigen_R$values + data$lambda
  Dinv <- 1/d
  Dinv[is.infinite(Dinv)] <- 0
  SinvEb <- data$eigen_R$vectors %*% (Dinv * crossprod(data$eigen_R$vectors, Eb))
  
  return(-0.5 * (-2 * sum(rR * SinvEb) + sum(model$RjSinvRj * Eb2)))
}

# Single effect update
single_effect_update.rss_lambda <- function(data, model, l, 
                                            optimize_V, check_null_threshold) {
  
  # Remove lth effect from fitted values
  model$Rz <- model$Rz - data$R %*% (model$alpha[l,] * model$mu[l,])
  
  # Compute residuals
  r <- data$z - model$Rz
  
  # Perform single effect regression using precomputed matrices
  res <- single_effect_regression_rss(
    z = r,
    SinvRj = model$SinvRj,
    RjSinvRj = model$RjSinvRj,
    V = model$V[l],
    prior_weights = model$pi,
    optimize_V = optimize_V,
    check_null_threshold = check_null_threshold
  )
  
  # Update model with results
  model$mu[l,] <- res$mu
  model$alpha[l,] <- res$alpha
  model$mu2[l,] <- res$mu2
  model$V[l] <- res$V
  model$lbf[l] <- res$lbf_model
  model$lbf_variable[l,] <- res$lbf
  
  # Compute KL divergence  
  # KL[q(b_l) || p(b_l)] = -ELBO_l = -log(BF_l) + E_q[log q(b_l)/p(b_l | r_l, sigma^2)]
  # where r_l is the residual excluding effect l
  model$KL[l] <- -res$lbf_model + 
    SER_posterior_e_loglik.rss_lambda(data, model, r,
                                      res$alpha * res$mu, res$alpha * res$mu2)
  
  # Update fitted values
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
  
  # Use correlation matrix for purity calculations
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
update_variance_components.rss_lambda <- function(data, model) {
  # This method updates sigma2 and recomputes all Sigma-dependent matrices
  # It's called both during initialization and during iterations
  
  # Compute Sigma = sigma2*R + lambda*I
  eigenS_values <- model$sigma2 * data$eigen_R$values + data$lambda
  Dinv <- 1 / eigenS_values
  Dinv[is.infinite(Dinv)] <- 0
  
  # Compute Sigma^(-1) R_j = U (sigma2 D + lambda)^(-1) D U^T e_j
  model$SinvRj <- data$eigen_R$vectors %*% (Dinv * data$eigen_R$values * t(data$eigen_R$vectors))
  
  # Compute diagonal of R * Sigma^(-1) * R
  if (data$lambda == 0) {
    model$RjSinvRj <- diag(data$R) / model$sigma2
  } else {
    tmp <- t(data$eigen_R$vectors)
    model$RjSinvRj <- colSums(tmp * (Dinv * (data$eigen_R$values^2) * tmp))
  }
  
  return(model)
}

# Update derived quantities (ghost function - not used)
update_derived_quantities.rss_lambda <- function(data, model) {
  # Not used for RSS lambda
  return(model)
}

# Check convergence
check_convergence.rss_lambda <- function(data, model_prev, model_current, 
                                         elbo_prev, elbo_current, tol) {
  # Standard ELBO-based convergence
  return((elbo_current - elbo_prev) < tol)
}

# Update variance before convergence check
update_variance_before_convergence.rss_lambda <- function(data) {
  # For RSS lambda, we update variance before convergence check
  return(TRUE)
}

# Handle convergence and variance updates
handle_convergence_and_variance.rss_lambda <- function(data, model, model_prev, 
                                                       elbo_prev, elbo_current, tol,
                                                       estimate_residual_variance,
                                                       residual_variance_lowerbound,
                                                       residual_variance_upperbound) {
  
  converged <- FALSE
  
  # Check convergence first
  if ((elbo_current - elbo_prev) < tol) {
    converged <- TRUE
  }
  
  # Update residual variance if requested
  if (estimate_residual_variance && !converged) {
    if (data$lambda == 0) {
      # Closed form solution when lambda = 0
      num_nonzero_eigenvals <- sum(data$eigen_R$values != 0)
      # Pass sigma2 = 1 as in the original implementation
      est_sigma2 <- (1/num_nonzero_eigenvals) * get_ER2.rss_lambda(data, model, sigma2 = 1)
      
      if (est_sigma2 < 0) {
        stop("Estimating residual variance failed: the estimated value is negative")
      }
      
      # Enforce upper bound of 1 for RSS
      if (est_sigma2 > 1) {
        est_sigma2 <- 1
      }
      
      # Enforce bounds
      est_sigma2 <- max(residual_variance_lowerbound, 
                        min(residual_variance_upperbound, est_sigma2))
    } else {
      # Use optimization when lambda > 0
      # The constraint is sigma2 + lambda <= 1
      upper_bound <- min(residual_variance_upperbound, 1 - data$lambda)
      
      # Define log-likelihood function for optimization
      Eloglik <- function(sigma2) {
        # Update model sigma2 temporarily
        model_temp <- model
        model_temp$sigma2 <- sigma2
        model_temp <- update_variance_components.rss_lambda(data, model_temp)
        
        # Compute expected log-likelihood
        # E[log p(z | b, sigma2)] = -n/2 * log(2*pi) - 1/2 * log|Sigma| - 1/2 * E[(z-Rb)'Sigma^(-1)(z-Rb)]
        
        # log|Sigma| = sum(log(eigenvalues of Sigma))
        eigenS_values <- sigma2 * data$eigen_R$values + data$lambda
        log_det_Sigma <- sum(log(eigenS_values[eigenS_values > 0]))
        
        # Expected squared residuals term
        ER2_term <- get_ER2.rss_lambda(data, model_temp, sigma2 = sigma2)
        
        return(-0.5 * log_det_Sigma - 0.5 * ER2_term)
      }
      
      # Optimize
      opt_result <- optimize(Eloglik, interval = c(residual_variance_lowerbound, upper_bound),
                             maximum = TRUE)
      est_sigma2 <- opt_result$maximum
      
      # Double-check boundary
      if (Eloglik(est_sigma2) < Eloglik(upper_bound)) {
        est_sigma2 <- upper_bound
      }
    }
    
    model$sigma2 <- est_sigma2
    model <- update_variance_components.rss_lambda(data, model)
  }
  
  return(list(data = data, model = model, converged = converged))
}

# Expected log-likelihood
Eloglik.rss_lambda <- function(data, model) {
  d <- model$sigma2 * data$eigen_R$values + data$lambda
  if (data$lambda == 0) {
    result <- -(sum(d != 0)/2) * log(2*pi*model$sigma2) - 0.5*get_ER2.rss_lambda(data, model)
  } else {
    result <- -(length(data$z)/2)*log(2*pi) - 0.5*sum(log(d)) - 0.5*get_ER2.rss_lambda(data, model)
  }
  return(result)
}