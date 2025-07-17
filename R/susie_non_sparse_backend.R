# Unmappable effects backend methods: infinitesimal (ss_inf) and adaptive shrinkage (ss_ash)

# Optimize prior variance for unmappable effects methods
optimize_prior_variance_unmappable <- function(V_init, XtOmegar, diagXtOmegaX, prior_weights,
                                               bounds = c(0, 1)) {
  p <- length(XtOmegar)

  logpi0 <- if (!is.null(prior_weights)) {
    log(prior_weights + sqrt(.Machine$double.eps))
  } else {
    rep(log(1/p), p)
  }

  objective <- function(V) {
    -matrixStats::logSumExp(-0.5 * log(1 + V * diagXtOmegaX) +
                              V * XtOmegar^2 / (2 * (1 + V * diagXtOmegaX)) +
                              logpi0)
  }

  res <- optim(par = V_init,
               fn = objective,
               method = "Brent",
               lower = bounds[1],
               upper = bounds[2])

  if (!is.null(res$par) && res$convergence == 0) {
    return(res$par)
  } else {
    return(V_init)
  }
}

# Single Effect Update for ss_inf
single_effect_update.ss_inf <- function(data, model, l,
                                       optimize_V, check_null_threshold) {

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

  b <- colSums(mu * alpha) - mu[l, ] * alpha[l, ]
  XtOmegaXb <- V %*% ((t(V) %*% b) * Dsq / (tau2 * Dsq + sigma2))
  XtOmegar <- XtOmegay - XtOmegaXb

  V_l <- model$V[l]
  if (!is.null(optimize_V) && optimize_V == "optim") {
    V_l <- optimize_prior_variance_unmappable(V_l, XtOmegar, diagXtOmegaX, model$pi)
  }

  res <- single_effect_regression(
    Xty                  = XtOmegar,
    dXtX                 = diagXtOmegaX,
    V                    = V_l,
    residual_variance    = 1, # I hard coded residual variance to 1 as we
                              # account for this in the precision matrix, Omega.
    prior_weights        = model$pi,
    optimize_V           = "none",
    check_null_threshold = check_null_threshold
  )

  model$alpha[l, ]         <- res$alpha
  model$mu[l, ]            <- res$mu
  model$mu2[l, ]           <- res$mu2
  model$V[l]               <- res$V
  model$lbf[l]             <- res$lbf_model
  model$lbf_variable[l, ]  <- res$lbf

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

  return(model)
}

# Initialize susie model for ss_inf (includes tau2 and theta)
initialize_susie_model.ss_inf <- function(data, L, scaled_prior_variance, var_y,
                                       residual_variance, prior_weights, ...) {
  return(initialize_matrices(data$p, L, scaled_prior_variance, var_y,
                                residual_variance, prior_weights, include_unmappable = TRUE))
}
# Initialize fitted values for ss_inf
initialize_fitted.ss_inf <- function(data, alpha, mu) {
  b <- colSums(alpha * mu)
  XtXr <- data$XtX %*% b
  return(list(XtXr = XtXr))
}
# Extract core parameters across iterations for ss_inf (includes tau2)
extract_core.ss_inf <- function(data, model, tracking, iter, track_fit, ...) {
  if (isTRUE(track_fit)) {
    tracking[[iter]] <- list(alpha = model$alpha,
                             niter = iter,
                             V = model$V,
                             sigma2 = model$sigma2,
                             tau2 = model$tau2)
  }
  return(tracking)
}
# Update variance components for ss_inf (Method of Moments)
update_variance_components.ss_inf <- function(data, model) {
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
}
# Update derived quantities for ss_inf
update_derived_quantities.ss_inf <- function(data, model) {
  sigma2 <- model$sigma2
  tau2 <- if (is.null(model$tau2)) 0 else model$tau2

  var <- tau2 * data$eigen_values + sigma2
  data$var <- var
  data$diagXtOmegaX <- rowSums(sweep(data$eigen_vectors^2, 2, (data$eigen_values / var), `*`))
  data$XtOmegay <- data$eigen_vectors %*% (data$VtXty / var)

  data$theta <- compute_theta_blup(data, model)

  return(data)
}
# Check convergence for ss_inf (uses PIP differences, not ELBO)
check_convergence.ss_inf <- function(data, model_prev, model_current, elbo_prev, elbo_current, tol) {
  PIP_diff <- max(abs(model_prev$alpha - model_current$alpha))
  return(PIP_diff < tol)
}
# Update variance before convergence check for ss_inf
update_variance_before_convergence.ss_inf <- function(data) {
  return(TRUE)
}
# Handle convergence and variance updates for ss_inf
handle_convergence_and_variance.ss_inf <- function(data, model, model_prev, elbo_prev, elbo_current,
                                                    tol, estimate_residual_variance,
                                                    residual_variance_lowerbound, residual_variance_upperbound) {
  if (estimate_residual_variance) {
    result <- update_model_variance(data, model, residual_variance_lowerbound, residual_variance_upperbound)
    data <- result$data
    model <- result$model
  }

  converged <- check_convergence(data, model_prev, model, elbo_prev, elbo_current, tol)

  return(list(data = data, model = model, converged = converged))
}

# Credible Sets for unmappable effects methods
get_cs.ss_inf <- function(data, model, coverage, min_abs_corr, n_purity) {

  if (is.null(coverage) || is.null(min_abs_corr)) return(NULL)

  if (!is.null(data$XtX)) {
    if (any(!(diag(data$XtX) %in% c(0,1)))) {
      Xcorr <- muffled_cov2cor(data$XtX)
    } else {
      Xcorr <- data$XtX
    }
  } else {
    LD <- (data$eigen_vectors %*% diag(data$eigen_values)) %*% t(data$eigen_vectors) / data$n
    if (any(!(diag(LD) %in% c(0,1)))) {
      Xcorr <- muffled_cov2cor(LD)
    } else {
      Xcorr <- LD
    }
  }

  return(susie_get_cs(model, coverage = coverage,
                      Xcorr = Xcorr,
                      min_abs_corr = min_abs_corr,
                      check_symmetric = FALSE,
                      n_purity = n_purity))
}