#' Single Effect Regression
#'
#' Performs single effect regression (SER) on sufficient statistics.
#' This is an internal function that fits a single effect in the SuSiE model.
#'
#' @param data Data object for class checking (determines standard vs RSS computation)
#' @param Xty A p-vector of X'y values (sufficient statistics)
#' @param dXtX A p-vector of diagonal elements of X'X (sufficient statistics)
#' @param z A p-vector of z-scores (RSS-specific)
#' @param SinvRj A p×p matrix where column j is Σ^(-1)R_j (RSS-specific)
#' @param RjSinvRj A p-vector where element j is R_j'Σ^(-1)R_j (RSS-specific)
#' @param V Prior variance for the single effect
#' @param residual_variance The residual variance (sigma^2)
#' @param prior_weights A p-vector of prior weights for each variable (default NULL for uniform)
#' @param optimize_V Method for optimizing prior variance: "none", "optim", "uniroot", "EM", or "simple"
#' @param check_null_threshold Threshold for setting V to zero for numerical stability
#' @param unmappable_effects Whether to use unmappable effects optimization
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
  function(data = NULL,
           Xty = NULL,
           dXtX = NULL,
           z = NULL,
           SinvRj = NULL,
           RjSinvRj = NULL,
           V,
           residual_variance = 1,
           prior_weights = NULL,
           optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
           check_null_threshold = 0,
           unmappable_effects = FALSE) {

    # Match prior variance optimization argument
    optimize_V <- match.arg(optimize_V)

    # Check if input data is RSS with non-zero lambda
    if (inherits(data, "rss_lambda")) {

      # RSS-specific computation path
      p <- length(z)
      shat2 <- 1 / RjSinvRj

      if (is.null(prior_weights)) {
        prior_weights <- rep(1 / p, p)
      }

      if (optimize_V != "EM" && optimize_V != "none") {
        V <- optimize_prior_variance_rss(optimize_V, data, z, SinvRj, RjSinvRj, shat2,
          prior_weights,
          alpha = NULL, post_mean2 = NULL,
          V_init = V, check_null_threshold = check_null_threshold
        )
      }

      # Use loglik to compute log Bayes factors and alpha
      loglik_res <- loglik(data, V, z, SinvRj, RjSinvRj, shat2, prior_weights)
      lbf <- loglik_res$lbf
      alpha <- loglik_res$alpha
      lbf_model <- loglik_res$lbf_model

      # Compute posterior mean and variance
      post_var <- (RjSinvRj + 1 / V)^(-1)
      post_mean <- sapply(1:p, function(j) post_var[j] * sum(SinvRj[, j] * z))
      post_mean2 <- post_var + post_mean^2

      if (optimize_V == "EM") {
        V <- sum(alpha * post_mean2)
      }
    } else {
      # Computation path for individual, ss, and rss (with lambda = 0)
      betahat <- (1 / dXtX) * Xty
      shat2 <- residual_variance / dXtX
      beta_1 <- NULL

      # Check prior weights
      if (is.null(prior_weights)) {
        prior_weights <- rep(1 / length(dXtX), length(dXtX))
      }

      # Optimize Prior Variance of lth effect
      if (optimize_V != "EM" && optimize_V != "none") {

        # Specific prior variance optimization function for unmappable effects
        if (unmappable_effects && optimize_V == "optim") {
          V <- optimize_prior_variance_unmappable(V, Xty, dXtX, prior_weights)
        } else {
          # TODO: look into consolidating prior variance into singular function.
          V <- optimize_prior_variance(optimize_V, data, betahat, shat2, prior_weights,
            alpha = NULL, post_mean2 = NULL, V_init = V,
            check_null_threshold = check_null_threshold)
        }

      }

      # Use loglik generic to compute logged Bayes factors and alpha
      loglik_res <- loglik(data, V, betahat, shat2, prior_weights)
      lbf <- loglik_res$lbf
      alpha <- loglik_res$alpha
      lbf_model <- loglik_res$lbf_model

      # Compute posterior moments
      # TODO: We could convert posterior moment calculation to generics to reduce complexity here.
      # Alternatively, we could make a misc helper function in the utils file specifically for this
      # servin stephens. Something like if(servin){calclulate_servin_stephens_moments}else{original code}


      if (data$use_servin_stephens){
        if(V <= 0){
          post_mean  = rep(0, data$p)
          post_mean2 = rep(0, data$p)
          post_var   = rep(0, data$p)
          beta_1     = rep(0, data$p)
        } else {
          # Calculate Servin Stephens Posterior Mean
          post_mean = sapply(1:data$p, function(j){
            posterior_mean_servin_stephens(dXtX[j],
                                           Xty[j],
                                           V)})

          # Calculate Servin Stephens Posterior Variance
          var_result = lapply(1:data$p, function(j){
            posterior_var_servin_stephens(dXtX[j],
                                          Xty[j],
                                          crossprod(data$R),
                                          data$n,
                                          V)})

          post_var = sapply(var_result, function(x) x$post_var)
          beta_1 = sapply(var_result, function(x) x$beta1)
          post_mean2 = post_mean^2 + post_var

        }
      } else {
        # Standard Gaussian posterior calculations
        post_var <- (1 / V + dXtX / residual_variance)^(-1) # Posterior variance.
        post_mean <- (1 / residual_variance) * post_var * Xty
        post_mean2 <- post_var + post_mean^2 # Second moment.
      }

      # Expectation-maximization prior variance update using posterior moments.
      if (optimize_V == "EM") {
        V <- optimize_prior_variance(optimize_V, data, betahat, shat2, prior_weights,
          alpha, post_mean2,
          check_null_threshold = check_null_threshold,
          use_servin_stephens = data$use_servin_stephens,
          beta_1 = beta_1,
          n = data$n
        )
      }
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

optimize_prior_variance <- function(optimize_V, data, betahat, shat2, prior_weights,
                                    alpha = NULL, post_mean2 = NULL,
                                    V_init = NULL, check_null_threshold = 0,
                                    use_servin_stephens = FALSE, beta_1 = NULL, n = NULL) {
  V <- V_init
  if (optimize_V != "simple") {
    if (optimize_V == "optim") {
      lV <- optim(
        par = log(max(c(betahat^2 - shat2, 1), na.rm = TRUE)),
        fn = neg_loglik_logscale, data = data, betahat = betahat, shat2 = shat2,
        prior_weights = prior_weights, method = "Brent", lower = -30,
        upper = 15
      )$par
      ## if the estimated one is worse than current one, don't change it.
      if (neg_loglik_logscale(data, lV, betahat = betahat, shat2 = shat2, prior_weights = prior_weights) >
        neg_loglik_logscale(data, log(V),
          betahat = betahat,
          shat2 = shat2, prior_weights = prior_weights
        )) {
        lV <- log(V)
      }
      V <- exp(lV)
    } else if (optimize_V == "uniroot") {
      V <- est_V_uniroot(data, betahat, shat2, prior_weights)
    } else if (optimize_V == "EM") {
      if (use_servin_stephens) {
        # Servin-Stephens EM update
        V <- sqrt(sum(alpha * (betahat^2 + (beta_1/(n - 2)) + shat2)))
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
  if (loglik(data, 0, betahat, shat2, prior_weights)$lbf_model +
    check_null_threshold >= loglik(data, V, betahat, shat2, prior_weights)$lbf_model) {
    V <- 0
  }

  return(V)
}

# Optimize prior variance for RSS
optimize_prior_variance_rss <- function(optimize_V, data, z, SinvRj, RjSinvRj, shat2,
                                        prior_weights, alpha, post_mean2,
                                        V_init, check_null_threshold) {
  V <- V_init

  if (optimize_V != "simple") {
    if (optimize_V == "optim") {
      # Compute initial value: max((z^T Sigma^{-1} R_j)^2 - 1/RjSinvRj, 1e-6)
      init_vals <- sapply(1:length(z), function(j) sum(SinvRj[, j] * z)^2) - (1 / RjSinvRj)
      init_val <- max(c(init_vals, 1e-6), na.rm = TRUE)

      # Optimize log(V) using Brent method
      lV <- optim(
        par = log(init_val),
        fn = neg_loglik_logscale,
        data = data,
        z = z,
        SinvRj = SinvRj,
        RjSinvRj = RjSinvRj,
        shat2 = shat2,
        prior_weights = prior_weights,
        method = "Brent",
        lower = -30,
        upper = 15
      )$par

      # Check if new estimate improves likelihood
      if (neg_loglik_logscale(data, lV, z, SinvRj, RjSinvRj, shat2, prior_weights) >
        neg_loglik_logscale(data, log(V), z, SinvRj, RjSinvRj, shat2, prior_weights)) {
        lV <- log(V)
      }
      V <- exp(lV)
    } else if (optimize_V == "EM") {
      V <- sum(alpha * post_mean2)
    } else {
      stop("Invalid option for optimize_V method")
    }
  }

  # Check if V=0 gives better likelihood (null check)
  if (loglik(data, 0, z, SinvRj, RjSinvRj, shat2, prior_weights)$lbf_model + check_null_threshold >=
    loglik(data, V, z, SinvRj, RjSinvRj, shat2, prior_weights)$lbf_model) {
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

# Estimate prior variance.
est_V_uniroot <- function(data, betahat, shat2, prior_weights) {
  # Define loglikelihood and gradient as function of lV:=log(V)
  # to improve numerical optimization
  neg_loglik_grad_logscale <- function(lV, data, betahat, shat2, prior_weights) {
    -exp(lV) * loglik(data, exp(lV), betahat, shat2, prior_weights)$gradient
  }

  V.u <- uniroot(neg_loglik_grad_logscale, c(-10, 10),
    extendInt = "upX",
    data = data, betahat = betahat, shat2 = shat2, prior_weights = prior_weights
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
