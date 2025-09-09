#' Single Effect Regression
#'
#' Performs single effect regression (SER) for a single effect in the SuSiE model.
#' Uses S3 method dispatch based on data class for computation.
#'
#' @param data Data object with class determining computation method (individual, ss, or rss_lambda)
#' @param model Current SuSiE model containing alpha, mu, mu2, V, sigma2, and other parameters
#' @param l Integer index of the effect being updated (1 to L)
#' @param residual_variance The residual variance (sigma^2)
#' @param optimize_V Method for optimizing prior variance: "none", "optim", "EM", or "simple"
#' @param check_null_threshold Threshold for setting V to zero for numerical stability
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
           residual_variance = NULL,
           optimize_V = c("none", "optim", "EM", "simple"),
           check_null_threshold = 0) {

    # Match prior variance optimization argument
    optimize_V <- match.arg(optimize_V)

    # Store Prior Variance Value for the lth Effect
    V <- model$V[l]

    # Compute SER statistics (betahat, shat2, initial value for prior variance optimization)
    ser_stats <- compute_ser_statistics(data, model, residual_variance, l)

    # Optimize Prior Variance of lth effect
    if (optimize_V != "EM" && optimize_V != "none") {
      V <- optimize_prior_variance(optimize_V, data, model, ser_stats,
        alpha = NULL, post_mean2 = NULL, V, check_null_threshold)
    }

    # Use loglik to compute logged Bayes factors and posterior inclusion probabilities
    loglik_res <- loglik(data, model, V, ser_stats)
    lbf        <- loglik_res$lbf
    alpha      <- loglik_res$alpha
    lbf_model  <- loglik_res$lbf_model

    # Compute posterior moments
    moments    <- calculate_posterior_moments(data, model, V, residual_variance)

    post_mean  <- moments$post_mean
    post_mean2 <- moments$post_mean2
    beta_1     <- moments$beta_1

    # Expectation-maximization prior variance update using posterior moments
    if (optimize_V == "EM") {
      V <- optimize_prior_variance(optimize_V, data, model, ser_stats,
        alpha, post_mean2, V_init = NULL, check_null_threshold,
        data$use_servin_stephens, moments$beta_1, data$n)
    }

    return(list(
      alpha     = alpha,
      mu        = post_mean,
      mu2       = post_mean2,
      lbf       = lbf,
      lbf_model = lbf_model,
      V         = V
    ))
  }

# Prior Variance Optimization for the lth Effect
optimize_prior_variance <- function(optimize_V, data, model, ser_stats,
                                    alpha = NULL, post_mean2 = NULL,
                                    V_init = NULL, check_null_threshold = 0,
                                    use_servin_stephens = FALSE, beta_1 = NULL, n = NULL) {
  V <- V_init
  if (optimize_V != "simple") {
    if (optimize_V == "optim") {
      V_param_opt <- optim(
        par = ser_stats$optim_init, fn = neg_loglik,
        data = data, model = model,
        ser_stats = ser_stats,
        method = "Brent",
        lower = ser_stats$optim_bounds[1],
        upper = ser_stats$optim_bounds[2]
      )$par

      # Convert optimized parameter to V based on scale of optimization
      V_new <- if (ser_stats$optim_scale == "linear") {
        V_param_opt
      } else {
        exp(V_param_opt)
      }

      # Check if new estimate improves likelihood
      V_param_init <- if (ser_stats$optim_scale == "linear") V else log(V)
      if (neg_loglik(data, model, V_param_opt, ser_stats) >
          neg_loglik(data, model, V_param_init, ser_stats)) {
        V_new <- V
      }
      V <- V_new
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
  if (loglik(data, model, 0, ser_stats)$lbf_model +
    check_null_threshold >= loglik(data, model, V, ser_stats)$lbf_model) {
    V <- 0
  }

  return(V)
}


# Single Effect Update for the lth Effect
single_effect_update <- function(data, model, l, optimize_V, check_null_threshold) {

  # Compute Residuals
  model <- compute_residuals(data, model, l)

  res <- single_effect_regression(data, model, l,
                                  residual_variance = model$residual_variance,
                                  optimize_V = optimize_V,
                                  check_null_threshold = check_null_threshold)

  # Store results from SER
  model$alpha[l, ]        <- res$alpha
  model$mu[l, ]           <- res$mu
  model$mu2[l, ]          <- res$mu2
  model$V[l]              <- res$V
  model$lbf[l]            <- res$lbf_model
  model$lbf_variable[l, ] <- res$lbf

  # Update KL-divergence
  model$KL[l] <- compute_kl(data, model, l)

  # Update fitted values
  model <- update_fitted_values(data, model, l)

  return(model)
}
