### Non-sparse backend methods ###
# Implements infinitesimal (ss_inf) and adaptive shrinkage (ss_ash) methods

# Compute eigenvalue decomposition for non-sparse methods
compute_eigen_decomposition <- function(XtX, n) {
  LD <- XtX / n
  eig <- eigen(LD, symmetric = TRUE)
  idx <- order(eig$values, decreasing = TRUE)

  list(
    V = eig$vectors[, idx],
    Dsq = pmax(eig$values[idx] * n, 0),
    VtXty = NULL
  )
}

# Method of Moments variance estimation for non-sparse methods
MoM <- function(alpha, mu, omega, sigma2, tau2, n, V, Dsq, VtXty, Xty, yty,
                est_sigma2, est_tau2, verbose) {
  # Subroutine to estimate sigma^2, tau^2 using MoM
  L <- nrow(mu)  # Note: mu is L×p format
  p <- ncol(mu)

  ### Compute A. corresponds to the matrix in equation (37) of the supplement:
  ### where Tr(X'X) = sum(Dsq) and Tr(X'X)^2 = sum(Dsq^2)
  A <- matrix(0, nrow = 2, ncol = 2)
  A[1, 1] <- n
  A[1, 2] <- sum(Dsq)
  A[2, 1] <- A[1, 2]
  A[2, 2] <- sum(Dsq^2)

  # Compute diag(V'MV)
  b <- colSums(mu * alpha)  # equation 48
  Vtb <- t(V) %*% b
  diagVtMV <- Vtb^2  # portion of equation 51 + 52
  tmpD <- rep(0, p)

  for (l in seq_len(L)) {
    bl <- mu[l, ] * alpha[l, ]
    Vtbl <- t(V) %*% bl
    diagVtMV <- diagVtMV - Vtbl^2
    tmpD <- tmpD + alpha[l, ] * (mu[l, ]^2 + 1 / omega[l, ])
  }

  diagVtMV <- diagVtMV + rowSums(sweep(t(V)^2, 2, tmpD, `*`))

  # Compute x
  x <- rep(0, 2)
  x[1] <- yty - 2 * sum(b * Xty) + sum(Dsq * diagVtMV)  # equation 51
  x[2] <- sum(Xty^2) - 2 * sum(Vtb * VtXty * Dsq) + sum(Dsq^2 * diagVtMV)  # equation 52

  # Solves system of equations from equation 37
  if (est_tau2) {
    sol <- solve(A, x)
    if (sol[1] > 0 && sol[2] > 0) {
      sigma2 <- sol[1]
      tau2 <- sol[2]
    } else {
      sigma2 <- x[1] / n
      tau2 <- 0
    }
    if (verbose) {
      cat(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigma2, tau2))
    }
  } else if (est_sigma2) {
    sigma2 <- (x[1] - A[1, 2] * tau2) / n
    if (verbose) {
      cat(sprintf("Update sigma^2 to %f\n", sigma2))
    }
  }
  return(list(sigma2 = sigma2, tau2 = tau2))
}

# Single Effect Update for ss_inf (MATHEMATICAL EQUIVALENT to existing implementation)
single_effect_update.ss_inf <- function(data, model, l,
                                       optimize_V, check_null_threshold) {
  
  # Use L×p format directly (same as susieR 2.0)
  alpha <- model$alpha   # L×p matrix (posterior inclusion probabilities)
  mu <- model$mu         # L×p matrix
  
  # Extract data components (use derived quantities if available)
  V <- data$eigen_vectors    # p×p matrix
  Dsq <- data$eigen_values   # p-vector
  VtXty <- data$VtXty        # p-vector
  Xty <- data$Xty            # p-vector
  n <- data$n
  p <- data$p
  L <- nrow(alpha)
  
  # Current variance components
  sigma2 <- model$sigma2
  tau2 <- if (is.null(model$tau2)) 0 else model$tau2
  
  # Use precomputed derived quantities if available, otherwise compute
  if (!is.null(data$diagXtOmegaX) && !is.null(data$XtOmegay)) {
    diagXtOmegaX <- data$diagXtOmegaX
    XtOmegay <- data$XtOmegay
  } else {
    # Compute variance and precision matrices
    var <- tau2 * Dsq + sigma2
    diagXtOmegaX <- rowSums(sweep(V^2, 2, (Dsq / var), `*`))
    XtOmegay <- V %*% (VtXty / var)
  }
  
  # Single Effect Regression for effect l (EXACT algorithm from susie_inf.R)
  b <- colSums(mu * alpha) - mu[l, ] * alpha[l, ]
  XtOmegaXb <- V %*% ((t(V) %*% b) * Dsq / (tau2 * Dsq + sigma2))
  XtOmegar <- XtOmegay - XtOmegaXb
  
  # Get current effect size variance
  V_l <- model$V[l]
  
  # Update Prior Variance V[l] using optimization (susie-inf only supports "optim")
  # This corresponds to est_ssq = TRUE in the original susie_inf implementation
  if (!is.null(optimize_V) && optimize_V == "optim") {
    # Prior probabilities (uniform for now)
    logpi0 <- rep(log(1 / p), p)
    
    # Objective function for prior variance optimization
    f <- function(x) {
      -matrixStats::logSumExp(-0.5 * log(1 + x * diagXtOmegaX) +
                                x * XtOmegar^2 / (2 * (1 + x * diagXtOmegaX)) +
                                logpi0)
    }
    
    # Optimize with bounds (equivalent to V_range in original)
    res <- optim(par = V_l,
                 fn = f,
                 method = "Brent",
                 lower = 0,      # V_range[1]
                 upper = 1)      # V_range[2]
    
    if (!is.null(res$par) && res$convergence == 0) {
      V_l <- res$par
    }
  }
  
  # Update omega, mu, and alpha for effect l
  omega_l <- diagXtOmegaX + 1 / V_l
  mu[l, ] <- XtOmegar / omega_l
  
  # Compute log Bayes factors
  lbf_variable_l <- XtOmegar^2 / (2 * omega_l) - 0.5 * log(omega_l * V_l)
  
  # Prior probabilities (uniform for now)
  logpi0 <- rep(log(1 / p), p)
  log_alpha <- lbf_variable_l + logpi0
  lbf_l <- matrixStats::logSumExp(log_alpha)
  alpha[l, ] <- exp(log_alpha - lbf_l)
  
  # Compute second moment (mu2) properly 
  mu2 <- mu[l, ]^2 + 1 / omega_l  # Second moment = E[B^2] = Var[B] + E[B]^2
  
  # Update model components (already in L×p format)
  model$alpha[l, ]         <- alpha[l, ]
  model$mu[l, ]            <- mu[l, ]
  model$mu2[l, ]           <- mu2  # Proper second moment
  model$V[l]               <- V_l  # Updated prior variance
  model$lbf[l]            <- lbf_l
  model$lbf_variable[l, ] <- lbf_variable_l
  
  # Calculate KL divergence (same as standard ss)
  model$KL[l] <- -lbf_l + SER_posterior_e_loglik(data, model, XtOmegar,
                                                 Eb  = model$alpha[l, ] * model$mu[l, ],
                                                 Eb2 = model$alpha[l, ] * model$mu2[l, ])
  
  # Update fitted values to include current theta (if available)
  # For infinitesimal model: XtXr = XtX %*% (b + theta)
  b <- colSums(model$alpha * model$mu)
  if (!is.null(model$theta)) {
    model$XtXr <- data$XtX %*% (b + model$theta)
  } else {
    model$XtXr <- data$XtX %*% b  # Fallback to main effects only
  }
  
  return(model)
}

# Single Effect Update for ss_ash (same implementation as ss_inf for now)
single_effect_update.ss_ash <- single_effect_update.ss_inf

# Initialize matrices for ss_inf (includes tau2 and theta)
initialize_matrices.ss_inf <- function(data, L, scaled_prior_variance, var_y,
                                       residual_variance, prior_weights, ...) {
  return(create_matrix_initialization(data$p, L, scaled_prior_variance, var_y, 
                                      residual_variance, prior_weights, include_non_sparse = TRUE))
}

# Initialize matrices for ss_ash (same as ss_inf)
initialize_matrices.ss_ash <- initialize_matrices.ss_inf

# Initialize fitted values for ss_inf
initialize_fitted.ss_inf <- function(data, alpha, mu) {
  b <- colSums(alpha * mu)
  XtXr <- data$XtX %*% b
  return(list(XtXr = XtXr))
}

# Initialize fitted values for ss_ash (same as ss_inf)
initialize_fitted.ss_ash <- initialize_fitted.ss_inf

# Update variance components for ss_inf (Method of Moments)
update_variance_components.ss_inf <- function(data, model) {
  # Use L×p format directly
  alpha <- model$alpha   # L×p matrix (posterior inclusion probabilities)
  mu <- model$mu         # L×p matrix
  
  # Initialize omega matrix (needed for MoM)
  p <- data$p
  L <- nrow(alpha)
  
  # Current variance components
  sigma2 <- model$sigma2
  tau2 <- if (is.null(model$tau2)) 0 else model$tau2
  
  # Compute omega matrix
  var <- tau2 * data$eigen_values + sigma2
  diagXtOmegaX <- rowSums(sweep(data$eigen_vectors^2, 2, (data$eigen_values / var), `*`))
  omega <- matrix(rep(diagXtOmegaX, L), nrow = L, ncol = p, byrow = TRUE) + 
           matrix(rep(1 / model$V, p), nrow = L, ncol = p, byrow = FALSE)
  
  # Method of Moments variance estimation
  mom_result <- MoM(alpha, mu, omega, sigma2, tau2, data$n, 
                    data$eigen_vectors, data$eigen_values, data$VtXty, 
                    data$Xty, data$yty, 
                    est_sigma2 = TRUE, est_tau2 = TRUE, verbose = FALSE)
  
  return(list(sigma2 = mom_result$sigma2, tau2 = mom_result$tau2))
}

# Update variance components for ss_ash (same as ss_inf for now)
update_variance_components.ss_ash <- update_variance_components.ss_inf

# Update derived quantities for ss_inf (recompute var, diagXtOmegaX, XtOmegay, theta)
update_derived_quantities.ss_inf <- function(data, model) {
  # Update variance-dependent quantities after variance component changes
  sigma2 <- model$sigma2
  tau2 <- if (is.null(model$tau2)) 0 else model$tau2
  
  # Recompute derived quantities
  var <- tau2 * data$eigen_values + sigma2
  data$var <- var
  data$diagXtOmegaX <- rowSums(sweep(data$eigen_vectors^2, 2, (data$eigen_values / var), `*`))
  data$XtOmegay <- data$eigen_vectors %*% (data$VtXty / var)
  
  # Compute theta (random effects) using BLUP and store in data for transfer to model
  data$theta <- compute_theta_blup(data, model)
  
  return(data)
}

# Update derived quantities for ss_ash (same as ss_inf)
update_derived_quantities.ss_ash <- update_derived_quantities.ss_inf

# Check convergence for ss_inf (uses PIP differences, not ELBO)
check_convergence.ss_inf <- function(data, model_prev, model_current, elbo_prev, elbo_current, tol) {
  # Non-sparse convergence based on maximum change in PIP (alpha) values
  # This matches the original susie_inf convergence criterion
  PIP_diff <- max(abs(model_prev$alpha - model_current$alpha))
  return(PIP_diff < tol)
}

# Check convergence for ss_ash (same as ss_inf)
check_convergence.ss_ash <- check_convergence.ss_inf

# Update variance before convergence check for ss_inf
update_variance_before_convergence.ss_inf <- function(data) {
  # Non-sparse behavior: Update variance first, then check convergence
  return(TRUE)
}

# Update variance before convergence check for ss_ash (same as ss_inf)
update_variance_before_convergence.ss_ash <- update_variance_before_convergence.ss_inf

# Handle convergence and variance updates for ss_inf (non-sparse behavior)
handle_convergence_and_variance.ss_inf <- function(data, model, model_prev, elbo_prev, elbo_current, 
                                                    tol, estimate_residual_variance, 
                                                    residual_variance_lowerbound, residual_variance_upperbound) {
  # Non-sparse: Update variance first, then check convergence
  if (estimate_residual_variance) {
    result <- update_model_variance(data, model, residual_variance_lowerbound, residual_variance_upperbound)
    data <- result$data
    model <- result$model
  }
  
  converged <- check_convergence(data, model_prev, model, elbo_prev, elbo_current, tol)
  
  return(list(data = data, model = model, converged = converged))
}

# Handle convergence and variance updates for ss_ash (same as ss_inf)
handle_convergence_and_variance.ss_ash <- handle_convergence_and_variance.ss_inf

# Compute theta (random effects) using BLUP for non-sparse methods
compute_theta_blup <- function(data, model) {
  # Use L×p format directly (same as susieR 2.0)
  alpha <- model$alpha   # L×p matrix (posterior inclusion probabilities)
  mu <- model$mu         # L×p matrix
  
  # Current variance components
  sigma2 <- model$sigma2
  tau2 <- if (is.null(model$tau2)) 0 else model$tau2
  
  # Compute posterior means of main effects: b = colSums(mu * alpha)
  b <- colSums(mu * alpha)
  
  # Compute XtOmegaXb using eigen decomposition
  var <- tau2 * data$eigen_values + sigma2
  XtOmegaXb <- data$eigen_vectors %*% ((t(data$eigen_vectors) %*% b) * data$eigen_values / var)
  
  # Compute residual: XtOmegar = XtOmegay - XtOmegaXb
  XtOmegar <- data$XtOmegay - XtOmegaXb
  
  # Compute theta using BLUP: theta = tau2 * XtOmegar
  theta <- tau2 * XtOmegar
  
  return(theta)
}

# Credible Sets for non-sparse methods (uses standard susie_get_cs format)
get_cs.ss_inf <- function(data, model, coverage, min_abs_corr, n_purity) {
  
  if (is.null(coverage) || is.null(min_abs_corr)) return(NULL)
  
  # For non-sparse methods, we need to construct correlation matrix from eigen decomposition
  if (!is.null(data$XtX)) {
    # Use XtX if available
    if (any(!(diag(data$XtX) %in% c(0,1)))) {
      Xcorr <- muffled_cov2cor(data$XtX)
    } else {
      Xcorr <- data$XtX
    }
  } else {
    # Reconstruct correlation matrix from eigen decomposition
    LD <- (data$eigen_vectors %*% diag(data$eigen_values)) %*% t(data$eigen_vectors) / data$n
    if (any(!(diag(LD) %in% c(0,1)))) {
      Xcorr <- muffled_cov2cor(LD)
    } else {
      Xcorr <- LD
    }
  }
  
  # Use standard susie_get_cs function to ensure consistent output format
  return(susie_get_cs(model, coverage = coverage,
                      Xcorr = Xcorr,
                      min_abs_corr = min_abs_corr,
                      check_symmetric = FALSE,
                      n_purity = n_purity))
}

# Credible Sets for ss_ash (same as ss_inf)
get_cs.ss_ash <- get_cs.ss_inf

# Get marginal PIP for non-sparse methods (uses standard susie_get_pip format)
get_pip.ss_inf <- function(data, model, coverage, min_abs_corr, prior_tol) {
  
  if (is.null(coverage) || is.null(min_abs_corr)) return(NULL)
  
  # Use standard susie_get_pip function to ensure consistent output format
  return(susie_get_pip(model, prune_by_cs = FALSE, prior_tol = prior_tol))
}

# Get marginal PIP for ss_ash (same as ss_inf)
get_pip.ss_ash <- get_pip.ss_inf
