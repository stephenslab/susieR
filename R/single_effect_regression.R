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
#' @return Updated model with alpha, mu, mu2, lbf, lbf_variable, V, and KL stored for the lth effect
#'
#' @keywords internal
#' @noRd
single_effect_regression <- function(data, params, model, l) {

    # Fixed mixture prior path: evaluate BFs on a pre-specified variance grid
    # with given mixture weights, bypassing scalar V optimization entirely.
    # Activated by estimate_prior_method = "fixed_mixture" with non-NULL
    # prior_variance_grid and mixture_weights in params.
    if (params$estimate_prior_method == "fixed_mixture") {
      ser_stats <- compute_ser_statistics(data, params, model, l)
      model <- loglik_mixture(data, params, model, ser_stats, l)
      model <- calculate_posterior_moments_mixture(data, params, model, l)
      model <- compute_kl(data, params, model, l)
      # Store effective V as posterior-weighted grid mean (for diagnostics)
      V_eff <- sum(params$mixture_weights * params$prior_variance_grid)
      model <- set_prior_variance_l(model, l, V_eff)
      return(model)
    }

    # Standard scalar V path (unchanged)

    # Store Prior Variance Value for the lth Effect
    V <- get_prior_variance_l(model, l)

    # Compute SER statistics (betahat, shat2, initial value for prior variance optimization)
    ser_stats <- compute_ser_statistics(data, params, model, l)

    # Optimize Prior Variance of lth effect
    if (params$estimate_prior_method != "EM" && params$estimate_prior_method != "none") {
      V <- optimize_prior_variance(data, params, model, ser_stats,
        V_init = V)
    }

    # Compute logged Bayes factors and posterior inclusion probabilities
    model <- loglik(data, params, model, V, ser_stats, l)

    # Compute posterior moments
    model <- calculate_posterior_moments(data, params, model, V, l)

    # Compute KL divergence
    model <- compute_kl(data, params, model, l)

    # Expectation-maximization prior variance update using posterior moments
    if (params$estimate_prior_method == "EM") {
      V <- optimize_prior_variance(data, params, model, ser_stats,
                                    get_alpha_l(model, l),
                                    get_posterior_moments_l(model, l),
                                    V_init = V)
    }

    # Store prior variance
    model <- set_prior_variance_l(model, l, V)

    return(model)
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
    } else if (params$estimate_prior_method == "uniroot") {
      # Root-finding on the gradient of neg_loglik (on the optimization scale)
      neg_loglik_fn <- function(V_param) neg_loglik(data, params, model, V_param, ser_stats)
      neg_loglik_grad <- function(V_param) {
        h <- max(abs(V_param) * 1e-4, 1e-8)
        (neg_loglik_fn(V_param + h) - neg_loglik_fn(V_param - h)) / (2 * h)
      }

      V_root <- tryCatch(
        uniroot(neg_loglik_grad,
                interval = c(ser_stats$optim_bounds[1], ser_stats$optim_bounds[2]),
                extendInt = "yes",
                tol = .Machine$double.eps^0.25)$root,
        error = function(e) {
          # Fallback: if uniroot fails (no sign change), use initial value
          if (ser_stats$optim_scale == "linear") V else log(V)
        }
      )

      V_new <- if (ser_stats$optim_scale == "linear") V_root else exp(V_root)

      # Check if new estimate improves likelihood
      V_param_init <- if (ser_stats$optim_scale == "linear") V else log(V)
      if (neg_loglik(data, params, model, V_root, ser_stats) >
          neg_loglik(data, params, model, V_param_init, ser_stats)) {
        V_new <- V
      }
      V <- V_new
    } else if (params$estimate_prior_method == "EM") {
      V <- em_update_prior_variance(data, params, model, alpha, moments, V_init)
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
  #
  # For EM, skip this check: the null check would zero V without
  # recomputing the posterior, creating an inconsistent (q, V) pair
  # that can decrease the ELBO. Null effects are handled by
  # trim_null_effects() after convergence instead.
  # see https://github.com/stephenslab/mvsusieR/issues/26
  if (params$estimate_prior_method != "EM" &&
      params$estimate_prior_method != "none") {
    if (loglik(data, params, model, 0, ser_stats) +
      params$check_null_threshold >= loglik(data, params, model, V, ser_stats)) {
      V <- 0
    }
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
  model <- single_effect_regression(data, params, model, l)

  # Update fitted values
  model <- update_fitted_values(data, params, model, l)

  return(model)
}
