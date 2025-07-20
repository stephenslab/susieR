#' Single Effect Regression
#'
#' Performs single effect regression (SER) on sufficient statistics.
#' This is an internal function that fits a single effect in the SuSiE model.
#'
#' @param Xty A p-vector of X'y values (sufficient statistics)
#' @param dXtX A p-vector of diagonal elements of X'X (sufficient statistics)
#' @param V Prior variance for the single effect
#' @param residual_variance The residual variance (sigma^2)
#' @param prior_weights A p-vector of prior weights for each variable (default NULL for uniform)
#' @param optimize_V Method for optimizing prior variance: "none", "optim", "uniroot", "EM", or "simple"
#' @param check_null_threshold Threshold for setting V to zero for numerical stability
#'
#' @return A list containing:
#' \item{alpha}{Posterior inclusion probabilities}
#' \item{mu}{Posterior means}
#' \item{mu2}{Posterior second moments}
#' \item{lbf}{Log Bayes factors}
#' \item{V}{Optimized prior variance}
#' \item{lbf_model}{Model log Bayes factor}
#'
#' @keywords internal
#' @noRd
single_effect_regression <-
  function (Xty,
            dXtX,
            V,
            residual_variance = 1,
            prior_weights = NULL,
            optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
            check_null_threshold = 0,
            unmappable_effects = FALSE) {

    optimize_V <- match.arg(optimize_V)
    betahat <- (1 / dXtX) * Xty
    shat2 <- residual_variance / dXtX

    # Check prior weights
    if (is.null(prior_weights))
      prior_weights <- rep(1 / length(dXtX), length(dXtX))

    # Optimize Prior Variance of lth effect
    if (optimize_V != "EM" && optimize_V != "none") {
      if (unmappable_effects && optimize_V == "optim") {
        V <- optimize_prior_variance_unmappable(V, Xty, dXtX, prior_weights)
      } else {
        V <- optimize_prior_variance(optimize_V, betahat, shat2, prior_weights,
                                    alpha = NULL, post_mean2 = NULL, V_init = V,
                                    check_null_threshold = check_null_threshold)
      }
    }

    # log(bf) for each SNP
    lbf <- dnorm(betahat, 0, sqrt(V + shat2), log = TRUE) -
      dnorm(betahat, 0, sqrt(shat2), log = TRUE)
    lpo <- lbf + log(prior_weights + sqrt(.Machine$double.eps))

    # Deal with special case of infinite shat2 (e.g., happens if X does
    # not vary).
    lbf[is.infinite(shat2)] <- 0
    lpo[is.infinite(shat2)] <- 0
    maxlpo <- max(lpo)

    # w is proportional to BF, but subtract max for numerical stability.
    w_weighted <- exp(lpo - maxlpo)

    # Posterior prob for each SNP.
    weighted_sum_w <- sum(w_weighted)
    alpha <- w_weighted / weighted_sum_w
    post_var <- (1 / V + dXtX / residual_variance)^(-1) # Posterior variance.
    post_mean <- (1 / residual_variance) * post_var * Xty
    post_mean2 <- post_var + post_mean^2 # Second moment.

    # BF for single effect model.
    lbf_model <- maxlpo + log(weighted_sum_w)

    if(optimize_V == "EM")
      V <- optimize_prior_variance(optimize_V, betahat, shat2, prior_weights,
                                  alpha, post_mean2,
                                  check_null_threshold = check_null_threshold)

    return(list(alpha = alpha,
                mu = post_mean,
                mu2 = post_mean2,
                lbf = lbf,
                lbf_model = lbf_model,
                V = V))

  }

# Estimate prior variance.
est_V_uniroot <- function (betahat, shat2, prior_weights) {
  V.u <- uniroot(negloglik.grad.logscale, c(-10, 10), extendInt = "upX",
                betahat = betahat, shat2 = shat2, prior_weights = prior_weights)
  return(exp(V.u$root))
}

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

optimize_prior_variance <- function (optimize_V, betahat, shat2, prior_weights,
                                    alpha = NULL, post_mean2 = NULL,
                                    V_init = NULL, check_null_threshold = 0) {
  V <- V_init
  if (optimize_V != "simple") {
    if(optimize_V == "optim") {
      lV <- optim(par = log(max(c(betahat^2 - shat2, 1), na.rm = TRUE)),
                 fn = neg.loglik.logscale, betahat = betahat, shat2 = shat2,
                 prior_weights = prior_weights, method = "Brent", lower = -30,
                 upper = 15)$par
      ## if the estimated one is worse than current one, don't change it.
      if(neg.loglik.logscale(lV, betahat = betahat, shat2 = shat2, prior_weights = prior_weights) >
         neg.loglik.logscale(log(V), betahat = betahat,
                             shat2 = shat2, prior_weights = prior_weights)){
        lV <- log(V)
      }
      V <- exp(lV)
    } else if (optimize_V == "uniroot")
      V <- est_V_uniroot(betahat, shat2, prior_weights)
    else if (optimize_V == "EM")
      V <- sum(alpha * post_mean2)
    else
      stop("Invalid option for optimize_V method")
  }

  # Set V exactly 0 if that beats the numerical value by
  # check_null_threshold in loglik. check_null_threshold = 0.1 is
  # exp(0.1) = 1.1 on likelihood scale; it means that for parsimony
  # reasons we set estiate of V to zero, if its numerical estimate is
  # only "negligibly" different from zero. We use a likelihood ratio
  # of exp(check_null_threshold) to define "negligible" in this
  # context. This is fairly modest condition compared to, say, a
  # formal LRT with p-value 0.05. But the idea is to be lenient to
  # non-zeros estimates unless they are indeed small enough to be
  # neglible. See more intuition at
  # https://stephens999.github.io/fiveMinuteStats/LR_and_BF.html
  if (loglik(0, betahat, shat2, prior_weights) +
      check_null_threshold >= loglik(V, betahat, shat2, prior_weights))
    V <- 0

  return(V)
}

# In these functions, s2 represents residual_variance, and shat2 is an
# estimate of it.

# The log likelihood function for SER model (based on summary data
# betahat, shat2) as a function of prior variance V.
#
#' @importFrom Matrix colSums
#' @importFrom stats dnorm
loglik <- function (V, betahat, shat2, prior_weights) {

  #log(bf) for each SNP
  lbf <- dnorm(betahat, 0, sqrt(V + shat2), log = TRUE) -
    dnorm(betahat, 0, sqrt(shat2), log = TRUE)
  lpo <- lbf + log(prior_weights + sqrt(.Machine$double.eps))

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] <- 0
  lpo[is.infinite(shat2)] <- 0

  maxlpo <- max(lpo)
  w_weighted <- exp(lpo - maxlpo)
  weighted_sum_w <- sum(w_weighted)
  return(log(weighted_sum_w) + maxlpo)
}

neg.loglik.logscale <- function(lV, betahat, shat2, prior_weights)
  -loglik(exp(lV), betahat, shat2, prior_weights)

#' @importFrom Matrix colSums
#' @importFrom stats dnorm
loglik.grad <- function(V, betahat, shat2, prior_weights) {

  # log(bf) for each SNP.
  lbf <- dnorm(betahat, 0, sqrt(V + shat2), log = TRUE) -
    dnorm(betahat, 0, sqrt(shat2), log = TRUE)
  lpo <- lbf + log(prior_weights + sqrt(.Machine$double.eps))

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] <- 0
  lpo[is.infinite(shat2)] <- 0

  maxlpo <- max(lpo)
  w_weighted <- exp(lpo - maxlpo)
  weighted_sum_w <- sum(w_weighted)
  alpha <- w_weighted / weighted_sum_w
  return(sum(alpha * lbf.grad(V, shat2, betahat^2 / shat2)))
}

# Define loglikelihood and gradient as function of lV:=log(V)
# to improve numerical optimization
negloglik.grad.logscale <- function (lV, betahat, shat2, prior_weights)
  -exp(lV) * loglik.grad(exp(lV), betahat, shat2, prior_weights)

# Vector of gradients of logBF_j for each j, with respect to prior
# variance V.
lbf.grad <- function (V, shat2, T2) {
  l <- 0.5 * (1 / (V + shat2)) * ((shat2 / (V + shat2)) * T2 - 1)
  l[is.nan(l)] <- 0
  return(l)
}

lbf <- function (V, shat2, T2) {
  l <- 0.5 * log(shat2 / (V + shat2)) + 0.5 * T2 * (V / (V + shat2))
  l[is.nan(l)] <- 0
  return(l)
}
