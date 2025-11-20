# =============================================================================
# SINGLE EFFECT REGRESSION
#
# Performs single effect regression for the lth effect in the SuSiE model.
# Computes posterior moments, log Bayes factors, and optimizes prior variance.
# =============================================================================
#'
#' @param data Data object (individual, ss, or rss_lambda)
#' @param params Validated params object
#' @param model Current SuSiE model object
#' @param l Effect index being updated
#'
#' @return List with alpha, mu, mu2, lbf, lbf_model, and V for the lth effect
#'
#' @keywords internal
#' @noRd
single_effect_regression <- function(data, params, model, l) {

    # Store Prior Variance Value for the lth Effect
    V <- model$V[l]

    # Compute SER statistics (betahat, shat2, initial value for prior variance optimization)
    ser_stats <- compute_ser_statistics(data, params, model, l)

    # Optimize Prior Variance of lth effect
    if (params$estimate_prior_method != "EM" && params$estimate_prior_method != "none") {
      V <- optimize_prior_variance(data, params, model, ser_stats,
        V_init = V)
    }

    # Compute logged Bayes factors and posterior inclusion probabilities
    loglik_res <- loglik(data, params, model, V, ser_stats)

    # Compute posterior moments
    moments    <- calculate_posterior_moments(data, params, model, V)

    # Expectation-maximization prior variance update using posterior moments
    if (params$estimate_prior_method == "EM") {
      V <- optimize_prior_variance(data, params, model, ser_stats,
        loglik_res$alpha, moments, V_init = NULL)
    }

    return(list(
      alpha     = loglik_res$alpha,
      mu        = moments$post_mean,
      mu2       = moments$post_mean2,
      lbf       = loglik_res$lbf,
      lbf_model = loglik_res$lbf_model,
      V         = V
    ))
  }

# =============================================================================
# PRIOR VARIANCE OPTIMIZATION
#
# Optimizes prior variance for single effects using different methods.
# Handles optim, EM, simple methods and null threshold checking.
# =============================================================================
#'
#' @param data Data object
#' @param params Validated params object
#' @param model Current SuSiE model object
#' @param ser_stats SER statistics and optimization parameters from compute_ser_statistics
#' @param moments Posterior moments from calculate_posterior_moments
#' @param ... Additional method-specific parameters
#'
#' @return Optimized prior variance (scalar)
#'
#' @keywords internal
#' @noRd
optimize_prior_variance <- function(data, params, model, ser_stats,
                                    alpha = NULL, moments = NULL,
                                    V_init = NULL) {
  V <- V_init
  if (params$estimate_prior_method != "simple") {
    if (params$estimate_prior_method == "optim") {
      V_param_opt <- optim(
        par = ser_stats$optim_init,
        fn = function(V_param) neg_loglik(data, params, model, V_param, ser_stats),
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
      if (neg_loglik(data, params, model, V_param_opt, ser_stats) >
          neg_loglik(data, params, model, V_param_init, ser_stats)) {
        V_new <- V
      }
      V <- V_new
    } else if (params$estimate_prior_method == "EM") {
      if (params$use_servin_stephens) {
        # Servin-Stephens EM update
        V <- sqrt(sum(alpha * (ser_stats$betahat^2 + moments$beta_1/(data$n - 2) + ser_stats$shat2)))
      } else {
        # Standard EM update
        V <- sum(alpha * moments$post_mean2)
      }
    } else {
      stop("Invalid option for estimate_prior_method: ", params$estimate_prior_method)
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
  if (loglik(data, params, model, 0, ser_stats)$lbf_model +
    params$check_null_threshold >= loglik(data, params, model, V, ser_stats)$lbf_model) {
    V <- 0
  }

  return(V)
}

# =============================================================================
# SINGLE EFFECT UPDATE
#
# High-level function that updates one effect in the SuSiE model.
# Coordinates residual computation, SER, KL divergence, and fitted value updates.
# =============================================================================
#'
#' @param data Data object (individual, ss, or rss_lambda)
#' @param params Validated params object
#' @param model Current SuSiE model object
#' @param l Effect index being updated
#'
#' @return Updated SuSiE model object with new parameters for effect l
#'
#' @keywords internal
#' @noRd
single_effect_update <- function(data, params, model, l) {

  # Compute Residuals
  model <- compute_residuals(data, params, model, l)

  # Run Single Effect Regression
  res <- single_effect_regression(data, params, model, l)

  # Store results from SER
  model$alpha[l, ]        <- res$alpha
  model$mu[l, ]           <- res$mu
  model$mu2[l, ]          <- res$mu2
  model$V[l]              <- res$V
  model$lbf[l]            <- res$lbf_model
  model$lbf_variable[l, ] <- res$lbf
  model$KL[l]             <- compute_kl(data, params, model, l)

  # Update fitted values
  model <- update_fitted_values(data, params, model, l)

  return(model)
}
