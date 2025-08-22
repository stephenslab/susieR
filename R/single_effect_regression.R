#' Single Effect Regression
#'
#' Performs single effect regression (SER) for a single effect in the SuSiE model.
#' Uses S3 method dispatch based on data class for computation.
#'
#' @param data Data object with class determining computation method (individual, ss, or rss_lambda)
#' @param model Current SuSiE model containing alpha, mu, mu2, V, sigma2, and other parameters
#' @param l Integer index of the effect being updated (1 to L)
#' @param residuals Unified residuals input: XtR for individual/ss data, z_residual for RSS
#' @param dXtX A p-vector of diagonal elements of X'X (from data attributes)
#' @param residual_variance The residual variance (sigma^2)
#' @param prior_weights A p-vector of prior weights for each variable
#' @param optimize_V Method for optimizing prior variance: "none", "optim", "uniroot", "EM", or "simple"
#' @param check_null_threshold Threshold for setting V to zero for numerical stability
#' @param unmappable_effects Whether to use unmappable effects optimization (for infinitesimal/ash model)
#'
#' @return A list containing:
#' \item{alpha}{Posterior inclusion probabilities (p-vector)}
#' \item{mu}{Posterior means (p-vector)}
#' \item{mu2}{Posterior second moments (p-vector)}
#' \item{lbf}{Log Bayes factors for each variable (p-vector)}
#' \item{lbf_model}{Model log Bayes factor (scalar)}
#' \item{V}{Optimized prior variance (scalar)}
#'
#' @keywords internal
#' @noRd
single_effect_regression <-
  function(data, model, l,
           residuals,
           dXtX,
           residual_variance = NULL,
           prior_weights = NULL,
           optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
           check_null_threshold = 0,
           unmappable_effects = FALSE) {

    # Match prior variance optimization argument
    optimize_V <- match.arg(optimize_V)

    # Store Prior Variance Value for the lth Effect
    V <- model$V[l]

    # Extract weights
    if (!is.null(model$pi)) {
      prior_weights <- model$pi
    } else {
      prior_weights <- rep(1 / data$p, data$p)
    }

    # Set residual_variance if not provided
    if (is.null(residual_variance)) {
      residual_variance <- model$sigma2
    }

    # Compute SER statistics
    ser_stats <- compute_ser_statistics(data, model, residuals, dXtX, residual_variance)

    # Optimize Prior Variance of lth effect
    if (optimize_V != "EM" && optimize_V != "none") {
      # Specific prior variance optimization function for unmappable effects
      if (unmappable_effects && optimize_V == "optim") {
        V <- optimize_prior_variance_unmappable(V, residuals, dXtX, prior_weights)
      } else {
        V <- optimize_prior_variance(optimize_V, data, model, ser_stats, residuals,
          prior_weights, alpha = NULL, post_mean2 = NULL, V_init = V,
          check_null_threshold = check_null_threshold)
      }
    }

    # Use loglik generic to compute logged Bayes factors and alpha
    loglik_res <- loglik(data, model, V, residuals, ser_stats, prior_weights)
    lbf <- loglik_res$lbf
    alpha <- loglik_res$alpha
    lbf_model <- loglik_res$lbf_model

    # Compute posterior moments
    moments <- calculate_posterior_moments(data, model = model, V = V,
                                          residuals = residuals,
                                          dXtX = dXtX,
                                          residual_variance = residual_variance)
    post_mean <- moments$post_mean
    post_mean2 <- moments$post_mean2
    beta_1 <- moments$beta_1

    # Expectation-maximization prior variance update using posterior moments.
    if (optimize_V == "EM") {
      V <- optimize_prior_variance(optimize_V, data, model, ser_stats, residuals,
        prior_weights, alpha, post_mean2, V_init = NULL,
        check_null_threshold = check_null_threshold,
        use_servin_stephens = data$use_servin_stephens,
        beta_1 = moments$beta_1, n = data$n)
    }

    return(list(
      alpha = alpha,
      mu = post_mean,
      mu2 = post_mean2,
      lbf = lbf,
      lbf_model = lbf_model,
      V = V
    ))
  }

# Optimization functions
optimize_prior_variance <- function(optimize_V, data, model, ser_stats, residuals,
                                    prior_weights, alpha = NULL, post_mean2 = NULL,
                                    V_init = NULL, check_null_threshold = 0,
                                    use_servin_stephens = FALSE, beta_1 = NULL, n = NULL) {
  V <- V_init
  if (optimize_V != "simple") {
    if (optimize_V == "optim") {
      # Use unified interface for optimization
      lV <- optim(
        par = ser_stats$optim_init,
        fn = neg_loglik_logscale,
        data = data, model = model,
        residuals = residuals, ser_stats = ser_stats,
        prior_weights = prior_weights,
        method = "Brent", lower = -30, upper = 15
      )$par

      # Check if new estimate improves likelihood
      if (neg_loglik_logscale(data, model, lV, residuals, ser_stats, prior_weights) >
          neg_loglik_logscale(data, model, log(V), residuals, ser_stats, prior_weights)) {
        lV <- log(V)
      }
      V <- exp(lV)
    } else if (optimize_V == "uniroot") {
      V <- est_V_uniroot(data, model, residuals, ser_stats, prior_weights)
    } else if (optimize_V == "EM") {
      if (use_servin_stephens) {
        # Servin-Stephens EM update
        V <- sqrt(sum(alpha * (ser_stats$betahat^2 + beta_1/(n - 2) + ser_stats$shat2)))
      } else {
        # Standard EM update
        V <- sum(alpha * post_mean2)
      }
    } else {
      stop("Invalid option for optimize_V method")
    }
  }

  # TODO: Need to formally test the effect of this on a larger scale.

  # Set V exactly 0 if that beats the numerical value by
  # check_null_threshold in loglik. check_null_threshold = 0.1 is
  # exp(0.1) = 1.1 on likelihood scale; it means that for parsimony
  # reasons we set estimate of V to zero, if its numerical estimate is
  # only "negligibly" different from zero. We use a likelihood ratio
  # of exp(check_null_threshold) to define "negligible" in this
  # context. This is fairly modest condition compared to, say, a
  # formal LRT with p-value 0.05. But the idea is to be lenient to
  # non-zeros estimates unless they are indeed small enough to be
  # neglible. See more intuition at
  # https://stephens999.github.io/fiveMinuteStats/LR_and_BF.html
  if (loglik(data, model, 0, residuals, ser_stats, prior_weights)$lbf_model +
    check_null_threshold >= loglik(data, model, V, residuals, ser_stats, prior_weights)$lbf_model) {
    V <- 0
  }

  return(V)
}


optimize_prior_variance_unmappable <- function(V_init, XtOmegar, diagXtOmegaX, prior_weights,
                                               bounds = c(0, 1)) {
  p <- length(XtOmegar)

  logpi0 <- if (!is.null(prior_weights)) {
    log(prior_weights + sqrt(.Machine$double.eps))
  } else {
    rep(log(1 / p), p)
  }

  objective <- function(V) {
    -matrixStats::logSumExp(-0.5 * log(1 + V * diagXtOmegaX) +
      V * XtOmegar^2 / (2 * (1 + V * diagXtOmegaX)) +
      logpi0)
  }

  res <- optim(
    par = V_init,
    fn = objective,
    method = "Brent",
    lower = bounds[1],
    upper = bounds[2]
  )

  if (!is.null(res$par) && res$convergence == 0) {
    return(res$par)
  } else {
    return(V_init)
  }
}

# Estimate prior variance using uniroot
est_V_uniroot <- function(data, model, residuals, ser_stats, prior_weights) {
  # Define loglikelihood and gradient as function of lV:=log(V)
  # to improve numerical optimization
  neg_loglik_grad_logscale <- function(lV, data, model, residuals, ser_stats, prior_weights) {
    -exp(lV) * loglik(data, model, exp(lV), residuals, ser_stats, prior_weights)$gradient
  }

  V.u <- uniroot(neg_loglik_grad_logscale, c(-10, 10),
    extendInt = "upX",
    data = data, model = model, residuals = residuals,
    ser_stats = ser_stats, prior_weights = prior_weights
  )
  return(exp(V.u$root))
}

# Vector of gradients of logBF_j for each j, with respect to prior
# variance V.
lbf.grad <- function(V, shat2, T2) {
  l <- 0.5 * (1 / (V + shat2)) * ((shat2 / (V + shat2)) * T2 - 1)
  l[is.nan(l)] <- 0
  return(l)
}
