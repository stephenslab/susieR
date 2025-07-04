### Non-sparse specific backend methods ###

# TODO: determine whether these first three functions should be in utils or here.

# Matrix conversion utilities for compatibility between susieR (L×p) and non-sparse code (p×L) formats
convert_to_nonsparse_format <- function(alpha, mu, mu2) {
  # Convert from susieR 2.0 format (L×p) to non-sparse format (p×L)
  list(
    PIP = t(alpha),  # p×L matrix
    mu = t(mu),      # p×L matrix
    mu2 = t(mu2)     # p×L matrix
  )
}

convert_from_nonsparse_format <- function(PIP, mu, mu2) {
  # Convert from non-sparse format (p×L) to susieR format (L×p)
  list(
    alpha = t(PIP),  # L×p matrix
    mu = t(mu),      # L×p matrix
    mu2 = t(mu2)     # L×p matrix
  )
}

# Helper function for eigen decomposition
compute_eigen_decomposition <- function(XtX, n) {
  # Compute LD matrix
  LD <- XtX / n

  # Eigen decomposition
  eig <- eigen(LD, symmetric = TRUE)

  # Order eigenvalues
  idx <- order(eig$values, decreasing = TRUE)

  # Return components in decreasing eigenvalue order
  list(
    V    = eig$vectors[, idx],                    # p×p eigenvectors
    Dsq  = pmax(eig$values[idx] * n, 0),         # p-vector eigenvalues (descending)
    VtXty = NULL  # initialize value to NULL
  )
}

# Method of Moments variance estimation
# TODO: determine whether this should be in utils
MoM <- function(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, Xty, yty,
                est_sigmasq, est_tausq, verbose) {
  # Subroutine to estimate sigma^2, tau^2 using MoM
  p <- nrow(mu)  # Note: mu is p×L in non-sparse format
  L <- ncol(mu)

  ### Compute A. corresponds to the matrix in equation (37) of the supplement:
  ### where Tr(X'X) = sum(Dsq) and Tr(X'X)^2 = sum(Dsq^2)
  A <- matrix(0, nrow = 2, ncol = 2)
  A[1, 1] <- n
  A[1, 2] <- sum(Dsq)
  A[2, 1] <- A[1, 2]
  A[2, 2] <- sum(Dsq^2)

  # Compute diag(V'MV)
  b <- rowSums(mu * PIP)  # equation 48
  Vtb <- t(V) %*% b
  diagVtMV <- Vtb^2  # portion of equation 51 + 52
  tmpD <- rep(0, p)

  for (l in seq_len(L)) {
    bl <- mu[, l] * PIP[, l]
    Vtbl <- t(V) %*% bl
    diagVtMV <- diagVtMV - Vtbl^2
    tmpD <- tmpD + PIP[, l] * (mu[, l]^2 + 1 / omega[, l])
  }

  diagVtMV <- diagVtMV + rowSums(sweep(t(V)^2, 2, tmpD, `*`))

  # Compute x
  x <- rep(0, 2)
  x[1] <- yty - 2 * sum(b * Xty) + sum(Dsq * diagVtMV)  # equation 51
  x[2] <- sum(Xty^2) - 2 * sum(Vtb * VtXty * Dsq) + sum(Dsq^2 * diagVtMV)  # equation 52

  # Solves system of equations from equation 37
  if (est_tausq) {
    sol <- solve(A, x)
    if (sol[1] > 0 && sol[2] > 0) {
      sigmasq <- sol[1]
      tausq <- sol[2]
    } else {
      sigmasq <- x[1] / n
      tausq <- 0
    }
    if (verbose) {
      cat(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigmasq, tausq))
    }
  } else if (est_sigmasq) {
    sigmasq <- (x[1] - A[1, 2] * tausq) / n
    if (verbose) {
      cat(sprintf("Update sigma^2 to %f\n", sigmasq))
    }
  }
  return(list(sigmasq = sigmasq, tausq = tausq))
}

# Single Effect Update for ss_inf (MATHEMATICAL EQUIVALENT to existing implementation)
single_effect_update.ss_inf <- function(data, model, l,
                                       optimize_V, check_null_threshold) {
  
  # Convert susieR 2.0 format (L×p) to non-sparse format (p×L) for calculations
  nonsparse_format <- convert_to_nonsparse_format(model$alpha, model$mu, model$mu2)
  PIP <- nonsparse_format$PIP   # p×L matrix
  mu <- nonsparse_format$mu     # p×L matrix
  
  # Extract data components (use derived quantities if available)
  V <- data$eigen_vectors    # p×p matrix
  Dsq <- data$eigen_values   # p-vector
  VtXty <- data$VtXty        # p-vector
  Xty <- data$Xty            # p-vector
  n <- data$n
  p <- data$p
  L <- ncol(PIP)
  
  # Current variance components
  sigmasq <- model$sigma2
  tausq <- if (is.null(model$tausq)) 0 else model$tausq
  
  # Use precomputed derived quantities if available, otherwise compute
  if (!is.null(data$diagXtOmegaX) && !is.null(data$XtOmegay)) {
    diagXtOmegaX <- data$diagXtOmegaX
    XtOmegay <- data$XtOmegay
  } else {
    # Compute variance and precision matrices
    var <- tausq * Dsq + sigmasq
    diagXtOmegaX <- rowSums(sweep(V^2, 2, (Dsq / var), `*`))
    XtOmegay <- V %*% (VtXty / var)
  }
  
  # Single Effect Regression for effect l (EXACT algorithm from susie_inf.R)
  b <- rowSums(mu * PIP) - mu[, l] * PIP[, l]
  XtOmegaXb <- V %*% ((t(V) %*% b) * Dsq / (tausq * Dsq + sigmasq))
  XtOmegar <- XtOmegay - XtOmegaXb
  
  # Get current effect size variance
  ssq_l <- model$V[l]
  
  # Update Prior Variance ssq[l] using optimization (susie-inf only supports "optim")
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
    
    # Optimize with bounds (equivalent to ssq_range in original)
    res <- optim(par = ssq_l,
                 fn = f,
                 method = "Brent",
                 lower = 0,      # ssq_range[1]
                 upper = 1)      # ssq_range[2]
    
    if (!is.null(res$par) && res$convergence == 0) {
      ssq_l <- res$par
    }
  }
  
  # Update omega, mu, and PIP for effect l
  omega_l <- diagXtOmegaX + 1 / ssq_l
  mu[, l] <- XtOmegar / omega_l
  
  # Compute log Bayes factors
  lbf_variable_l <- XtOmegar^2 / (2 * omega_l) - 0.5 * log(omega_l * ssq_l)
  
  # Prior probabilities (uniform for now)
  logpi0 <- rep(log(1 / p), p)
  logPIP <- lbf_variable_l + logpi0
  lbf_l <- matrixStats::logSumExp(logPIP)
  PIP[, l] <- exp(logPIP - lbf_l)
  
  # Compute second moment (mu2) properly 
  mu2 <- matrix(0, nrow = p, ncol = L)
  mu2[, l] <- mu[, l]^2 + 1 / omega_l  # Second moment = E[B^2] = Var[B] + E[B]^2
  
  # Convert back to susieR 2.0 format (L×p)
  susie_format <- convert_from_nonsparse_format(PIP, mu, mu2)
  
  # Update model components
  model$alpha[l, ]         <- susie_format$alpha[l, ]
  model$mu[l, ]            <- susie_format$mu[l, ]
  model$mu2[l, ]           <- susie_format$mu2[l, ]  # Proper second moment
  model$V[l]               <- ssq_l  # Updated prior variance
  model$lbf[l]            <- lbf_l
  model$lbf_variable[l, ] <- lbf_variable_l
  
  # Calculate KL divergence (same as standard ss)
  model$KL[l] <- -lbf_l + SER_posterior_e_loglik(data, model, XtOmegar,
                                                 Eb  = susie_format$alpha[l, ] * susie_format$mu[l, ],
                                                 Eb2 = susie_format$alpha[l, ] * susie_format$mu[l, ]^2)
  
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

# Initialize matrices for ss_inf (includes tausq and theta)
initialize_matrices.ss_inf <- function(data, L, scaled_prior_variance, var_y,
                                       residual_variance, prior_weights, ...){
  p <- data$p

  mat_init <- list(
    alpha        = matrix(1 / p, L, p),
    mu           = matrix(0,     L, p),
    mu2          = matrix(0,     L, p),
    V            = rep(scaled_prior_variance * var_y, L),
    KL           = rep(as.numeric(NA), L),
    lbf          = rep(as.numeric(NA), L),
    lbf_variable = matrix(as.numeric(NA), L, p),
    sigma2       = residual_variance,
    tausq        = 0,  # Initialize infinitesimal variance to 0
    theta        = rep(0, p),  # Initialize random effects to 0
    pi           = prior_weights
  )

  return(mat_init)
}

# Initialize matrices for ss_ash (same as ss_inf)
initialize_matrices.ss_ash <- initialize_matrices.ss_inf

# Initialize fitted values for ss_inf (includes random effects theta)
initialize_fitted.ss_inf <- function(data, alpha, mu) {
  # Main effects: b = colSums(alpha * mu)
  b <- colSums(alpha * mu)
  
  # For infinitesimal model: fitted = X(b + theta)
  # Since we're working with sufficient statistics, we compute XtXr = XtX %*% (b + theta)
  # Initially theta = 0, so this reduces to the standard calculation
  XtXr <- data$XtX %*% b
  
  return(list(XtXr = XtXr))
}

# Initialize fitted values for ss_ash (same as ss_inf)
initialize_fitted.ss_ash <- initialize_fitted.ss_inf

# Update variance components for ss_inf (Method of Moments)
update_variance_components.ss_inf <- function(data, model) {
  # Convert to non-sparse format for MoM calculation
  nonsparse_format <- convert_to_nonsparse_format(model$alpha, model$mu, model$mu2)
  PIP <- nonsparse_format$PIP   # p×L matrix
  mu <- nonsparse_format$mu     # p×L matrix
  
  # Initialize omega matrix (needed for MoM)
  p <- data$p
  L <- ncol(PIP)
  
  # Current variance components
  sigmasq <- model$sigma2
  tausq <- if (is.null(model$tausq)) 0 else model$tausq
  
  # Compute omega matrix
  var <- tausq * data$eigen_values + sigmasq
  diagXtOmegaX <- rowSums(sweep(data$eigen_vectors^2, 2, (data$eigen_values / var), `*`))
  omega <- matrix(rep(diagXtOmegaX, L), nrow = p, ncol = L) + 
           matrix(rep(1 / model$V, each = p), nrow = p, ncol = L)
  
  # Method of Moments variance estimation
  mom_result <- MoM(PIP, mu, omega, sigmasq, tausq, data$n, 
                    data$eigen_vectors, data$eigen_values, data$VtXty, 
                    data$Xty, data$yty, 
                    est_sigmasq = TRUE, est_tausq = TRUE, verbose = FALSE)
  
  return(list(sigma2 = mom_result$sigmasq, tausq = mom_result$tausq))
}

# Update variance components for ss_ash (same as ss_inf for now)
update_variance_components.ss_ash <- update_variance_components.ss_inf

# Update derived quantities for ss_inf (recompute var, diagXtOmegaX, XtOmegay, theta)
update_derived_quantities.ss_inf <- function(data, model) {
  # Update variance-dependent quantities after variance component changes
  sigmasq <- model$sigma2
  tausq <- if (is.null(model$tausq)) 0 else model$tausq
  
  # Recompute derived quantities
  var <- tausq * data$eigen_values + sigmasq
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

# Compute theta (random effects) using BLUP for non-sparse methods
compute_theta_blup <- function(data, model) {
  # Convert susieR 2.0 format (L×p) to non-sparse format (p×L) for calculations
  nonsparse_format <- convert_to_nonsparse_format(model$alpha, model$mu, model$mu2)
  PIP <- nonsparse_format$PIP   # p×L matrix
  mu <- nonsparse_format$mu     # p×L matrix
  
  # Current variance components
  sigmasq <- model$sigma2
  tausq <- if (is.null(model$tausq)) 0 else model$tausq
  
  # Compute posterior means of main effects: b = rowSums(mu * PIP)
  b <- rowSums(mu * PIP)
  
  # Compute XtOmegaXb using eigen decomposition
  var <- tausq * data$eigen_values + sigmasq
  XtOmegaXb <- data$eigen_vectors %*% ((t(data$eigen_vectors) %*% b) * data$eigen_values / var)
  
  # Compute residual: XtOmegar = XtOmegay - XtOmegaXb
  XtOmegar <- data$XtOmegay - XtOmegaXb
  
  # Compute theta using BLUP: theta = tausq * XtOmegar
  theta <- tausq * XtOmegar
  
  return(theta)
}
