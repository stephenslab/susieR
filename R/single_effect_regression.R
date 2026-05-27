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
    if (!params$estimate_prior_method %in%
        c("optim", "uniroot", "simple", "none", "EM", "fixed_mixture")) {
      stop("Invalid option for estimate_prior_method: ",
           params$estimate_prior_method)
    }

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

    # Two S3 hook slots: pre/post loglik. Defaults dispatch on
    # `params$estimate_prior_method`; downstream classes override.

    V <- get_prior_variance_l(model, l)
    ser_stats <- compute_ser_statistics(data, params, model, l)

    out <- pre_loglik_prior_hook(data, params, model, ser_stats,
                                 l = l, V_init = V)
    V     <- out$V
    model <- out$model

    model <- loglik(data, params, model, V, ser_stats, l)
    model <- calculate_posterior_moments(data, params, model, V, l)
    model <- compute_kl(data, params, model, l)

    out <- post_loglik_prior_hook(data, params, model, ser_stats,
                                  l = l, V_init = V)
    V     <- out$V
    model <- out$model

    model <- set_prior_variance_l(model, l, V)
    model
  }

#' Gaussian single-effect log Bayes factors
#'
#' Plain helper, not S3: once a backend has reduced its residual problem to
#' `betahat`, `shat2`, and scalar prior variance `V`, the Gaussian ABF is common.
#' The Wakefield form is algebraically equivalent to the normal-density ratio
#' for positive `shat2`, and avoids Inf - Inf in degenerate null-weight columns.
#'
#' @keywords internal
#' @noRd
gaussian_ser_lbf <- function(betahat, shat2, V) {
  shat2_safe <- pmax(shat2, .Machine$double.eps)
  lbf <- -0.5 * log(1 + V / shat2_safe) +
    0.5 * betahat^2 * V / (shat2_safe * (V + shat2_safe))
  lbf[!is.finite(betahat) | !is.finite(shat2)] <- 0
  lbf
}

#' Apply single-effect log Bayes factors to a model
#'
#' Shared writeback for scalar-prior SER methods. If `l` is NULL, returns only
#' the model-level log Bayes factor for prior-variance optimization.
#'
#' @keywords internal
#' @noRd
apply_ser_lbf <- function(model, lbf, shat2, l = NULL) {
  stable_res  <- lbf_stabilization(lbf, model$pi, shat2)
  weights_res <- compute_posterior_weights(stable_res$lpo)

  if (is.null(l))
    return(list(model = model, lbf_model = weights_res$lbf_model,
                weights = weights_res, stable = stable_res))

  model$alpha[l, ]        <- weights_res$alpha
  model$lbf[l]            <- weights_res$lbf_model
  model$lbf_variable[l, ] <- stable_res$lbf
  list(model = model, lbf_model = weights_res$lbf_model,
       weights = weights_res, stable = stable_res)
}

#' Gaussian posterior moments for one SER
#'
#' Plain helper for the common conjugate normal update.
#'
#' @keywords internal
#' @noRd
gaussian_ser_moments <- function(betahat, shat2, V) {
  post_var  <- V * shat2 / (V + shat2)
  post_mean <- post_var / shat2 * betahat
  no_info <- !is.finite(betahat) | !is.finite(shat2)
  post_var[no_info] <- 0
  post_mean[no_info] <- 0
  list(post_mean = post_mean, post_mean2 = post_var + post_mean^2)
}

#' Gaussian posterior expected log-likelihood for one SER
#' @keywords internal
#' @noRd
gaussian_ser_posterior_e_loglik <- function(alpha, mu, mu2,
                                            betahat, shat2) {
  finite_info <- is.finite(shat2)
  Eb  <- alpha * mu
  Eb2 <- alpha * mu2
  -0.5 * sum((-2 * Eb[finite_info] * betahat[finite_info] +
                Eb2[finite_info]) / shat2[finite_info])
}

#' Store SER posterior moments for effect l
#' @keywords internal
#' @noRd
store_ser_moments <- function(model, l, moments) {
  model$mu[l, ]  <- moments$post_mean
  model$mu2[l, ] <- moments$post_mean2
  model
}

#' Pre-loglik prior-update hook
#'
#' S3 generic, called between SER stats and `loglik`. Default
#' routes to `optimize_prior_variance` for `optim` / `uniroot` /
#' `simple`. Returns `list(V, model)`.
#'
#' @export
#' @keywords internal
pre_loglik_prior_hook <- function(data, params, model, ser_stats,
                                  l, V_init) {
  UseMethod("pre_loglik_prior_hook")
}

#' @export
#' @keywords internal
pre_loglik_prior_hook.default <- function(data, params, model, ser_stats,
                                          l, V_init) {
  if (params$estimate_prior_method %in% c("optim", "uniroot", "simple")) {
    return(optimize_prior_variance(data, params, model, ser_stats,
                                   l = l, V_init = V_init))
  }
  list(V = V_init, model = model)
}

#' Post-loglik prior-update hook
#'
#' S3 generic, called after `loglik` / posterior moments / KL.
#' Default routes to `optimize_prior_variance` for `EM`. Returns
#' `list(V, model)`.
#'
#' @export
#' @keywords internal
post_loglik_prior_hook <- function(data, params, model, ser_stats,
                                   l, V_init) {
  UseMethod("post_loglik_prior_hook")
}

#' @export
#' @keywords internal
post_loglik_prior_hook.default <- function(data, params, model, ser_stats,
                                           l, V_init) {
  if (identical(params$estimate_prior_method, "EM")) {
    return(optimize_prior_variance(
      data, params, model, ser_stats,
      l       = l,
      alpha   = get_alpha_l(model, l),
      moments = get_posterior_moments_l(model, l),
      V_init  = V_init))
  }
  list(V = V_init, model = model)
}

# =============================================================================
# PRIOR VARIANCE OPTIMIZATION
#
# Optimizes prior variance for single effects using different methods.
# Handles optim, EM, simple methods and null threshold checking.
# =============================================================================

#' Per-effect prior variance update (S3 generic)
#'
#' Dispatched on the data class so downstream packages with non-scalar
#' prior structures (e.g., mfsusieR's adaptive mixture-of-normals
#' prior, future cross-modality priors) can run a per-effect prior
#' update step here while reusing the surrounding SER scaffolding.
#'
#' The default path implements the standard susieR scalar-V
#' optimization (`optim` Brent / `uniroot` / `EM` / `simple` /
#' `none`) plus the post-optimization null-threshold check.
#'
#' @param data Data object (e.g., `individual`, `ss`, `rss_lambda`,
#'   or a downstream class such as `mv_individual`, `mf_individual`).
#' @param params Validated params object.
#' @param model Current SuSiE model object.
#' @param ser_stats SER statistics and optimization parameters from
#'   `compute_ser_statistics`.
#' @param l Index of the effect being updated. Used by downstream
#'   methods that need per-effect state (e.g., the EM mixture-weight
#'   path); the default method uses it only for diagnostic purposes.
#' @param alpha Per-SNP posterior weights for effect `l`, supplied by
#'   the EM path (`get_alpha_l(model, l)`); `NULL` on the pre-loglik
#'   call.
#' @param moments Posterior moments for effect `l`, supplied by the
#'   EM path (`get_posterior_moments_l(model, l)`); `NULL` on the
#'   pre-loglik call.
#' @param V_init Initial value for the prior variance scalar.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{`V`}{numeric scalar, the optimized prior variance for
#'     effect `l`.}
#'   \item{`model`}{the (possibly mutated) model object. The default
#'     method leaves `model` unchanged; downstream methods may write
#'     prior-state updates here (e.g., mixture-weight vectors).}
#' }
#'
#' @keywords internal
#' @noRd
optimize_prior_variance <- function(data, params, model, ser_stats,
                                    l       = NULL,
                                    alpha   = NULL,
                                    moments = NULL,
                                    V_init  = NULL) {
  UseMethod("optimize_prior_variance")
}

#' Optimize a scalar SER prior variance
#'
#' Shared implementation for one-dimensional scalar-prior updates. Callers
#' supply likelihood closures on the backend's optimization scale.
#'
#' @keywords internal
#' @noRd
optimize_scalar_prior_variance <- function(V_init, estimate_prior_method,
                                           neg_loglik_fn, loglik_fn,
                                           optim_init, optim_bounds,
                                           optim_scale,
                                           check_null_threshold = 0) {
  V <- V_init
  if (estimate_prior_method == "optim") {
    V_param_opt <- optim(par = optim_init, fn = neg_loglik_fn,
                         method = "Brent",
                         lower = optim_bounds[1],
                         upper = optim_bounds[2])$par

    V_new <- if (optim_scale == "linear") V_param_opt else exp(V_param_opt)
    V_param_init <- if (optim_scale == "linear") V else log(V)
    if (neg_loglik_fn(V_param_opt) > neg_loglik_fn(V_param_init))
      V_new <- V
    V <- V_new
  } else if (estimate_prior_method == "uniroot") {
    neg_loglik_grad <- function(V_param) {
      h <- max(abs(V_param) * 1e-4, 1e-8)
      (neg_loglik_fn(V_param + h) - neg_loglik_fn(V_param - h)) / (2 * h)
    }

    V_root <- tryCatch(
      uniroot(neg_loglik_grad, interval = optim_bounds,
              extendInt = "yes",
              tol = .Machine$double.eps^0.25)$root,
      error = function(e) {
        if (optim_scale == "linear") V else log(V)
      }
    )

    V_new <- if (optim_scale == "linear") V_root else exp(V_root)
    V_param_init <- if (optim_scale == "linear") V else log(V)
    if (neg_loglik_fn(V_root) > neg_loglik_fn(V_param_init))
      V_new <- V
    V <- V_new
  } else if (!estimate_prior_method %in% c("simple", "none")) {
    stop("Invalid option for estimate_prior_method: ", estimate_prior_method)
  }

  if (estimate_prior_method != "none" &&
      loglik_fn(0) + check_null_threshold >= loglik_fn(V))
    V <- 0

  V
}

#' Default scalar-V prior-variance optimization
#'
#' Backbone implementation of `optimize_prior_variance`. Handles the
#' five `params$estimate_prior_method` cases (`optim`, `uniroot`,
#' `EM`, `simple`, `none`) on a scalar prior variance and runs the
#' post-optimization null-threshold check.
#'
#' @inheritParams optimize_prior_variance
#' @return A named list `list(V = ..., model = model)` (see
#'   `optimize_prior_variance` for the full contract). `model` is
#'   returned unchanged by this default method.
#' @keywords internal
#' @noRd
optimize_prior_variance.default <- function(data, params, model, ser_stats,
                                            l       = NULL,
                                            alpha   = NULL,
                                            moments = NULL,
                                            V_init  = NULL) {
  if (params$estimate_prior_method == "EM") {
    V <- em_update_prior_variance(data, params, model, alpha, moments, V_init)
  } else {
    V <- optimize_scalar_prior_variance(
      V_init = V_init,
      estimate_prior_method = params$estimate_prior_method,
      neg_loglik_fn = function(V_param)
        neg_loglik(data, params, model, V_param, ser_stats),
      loglik_fn = function(V_val)
        loglik(data, params, model, V_val, ser_stats),
      optim_init = ser_stats$optim_init,
      optim_bounds = ser_stats$optim_bounds,
      optim_scale = ser_stats$optim_scale,
      check_null_threshold = params$check_null_threshold)
  }

  list(V = V, model = model)
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
