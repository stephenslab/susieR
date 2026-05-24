# =============================================================================
# SINGLE-EFFECT REGRESSION FROM SUMMARY STATISTICS
#
# Lightweight public interface for the one-effect model. This intentionally
# avoids the summary-statistics constructors because they require an LD matrix
# or reference factor; the single-effect model only needs marginal summary
# statistics.
# =============================================================================

#' Single-effect regression from summary statistics
#'
#' @description
#' Fit a single-effect regression (SER) model directly from marginal summary
#' statistics. Unlike \code{\link{susie_rss}}, this interface does not take an
#' LD matrix or reference factor matrix, and it never constructs a diagonal LD
#' matrix internally. The model has exactly one single effect.
#'
#' @param z A vector of z-scores. Provide either \code{z}, or both
#'   \code{bhat} and \code{shat}, but not both.
#' @param bhat Alternative summary statistics: marginal effect estimates.
#' @param shat Standard errors for \code{bhat}; may be a scalar or a vector
#'   with the same length as \code{bhat}.
#' @param n Optional sample size. When supplied, z-scores are adjusted using
#'   the same PVE adjustment as \code{\link{susie_rss}}.
#' @param var_y Optional sample variance of the outcome. When supplied with
#'   \code{n}, this sets the working outcome scale. With \code{bhat} and
#'   \code{shat}, coefficients returned by \code{\link{coef.susie}} are on the
#'   original \code{bhat} scale.
#' @param prior_variance Prior variance on the z-score scale when \code{n} is
#'   not supplied. Default \code{50}, matching the z-only RSS scale.
#' @param scaled_prior_variance Prior variance divided by the working outcome
#'   variance when \code{n} is supplied. Default \code{0.2}.
#' @param prior_weights A vector of prior probabilities, one per variable.
#'   The weights are normalized internally.
#' @param null_weight Prior probability of no effect. If nonzero, a null
#'   variable is appended with zero effect.
#' @param estimate_prior_method Method used for the single prior variance.
#'   \code{"optim"} optimizes the one-effect marginal likelihood;
#'   \code{"simple"} uses the supplied prior variance and applies the
#'   null-threshold check.
#' @param check_null_threshold If the log-likelihood at prior variance zero is
#'   within this amount of the selected prior variance, set the prior variance
#'   to zero.
#' @param prior_tol Prior-variance tolerance used by \code{\link{susie_get_pip}}.
#' @param coverage Coverage for credible-set construction. If \code{NULL},
#'   credible sets are not constructed and no no-LD hint is emitted.
#' @param ethres Effective-number threshold passed to
#'   \code{\link{susie_get_cs_attainable}}.
#'
#' @return A \code{"susie"} object with one single effect.
#'
#' @seealso \code{\link{susie_rss}}, \code{\link{susie_get_cs_attainable}},
#'   \code{\link{susie_get_cs}}
#'
#' @export
susie_ser <- function(z = NULL, bhat = NULL, shat = NULL, n = NULL,
                      var_y = NULL,
                      prior_variance = 50,
                      scaled_prior_variance = 0.2,
                      prior_weights = NULL,
                      null_weight = 0,
                      estimate_prior_method = c("optim", "simple"),
                      check_null_threshold = 0,
                      prior_tol = 1e-9,
                      coverage = 0.95,
                      ethres = NULL) {
  estimate_prior_method <- match.arg(estimate_prior_method)

  input <- prepare_ser_summary(z, bhat, shat, n, var_y,
                               prior_variance, scaled_prior_variance)
  p <- length(input$betahat)

  if (!is.numeric(check_null_threshold) ||
      length(check_null_threshold) != 1L ||
      !is.finite(check_null_threshold) ||
      check_null_threshold < 0)
    stop("check_null_threshold must be a single nonnegative finite numeric.")
  if (!is.numeric(prior_tol) || length(prior_tol) != 1L ||
      !is.finite(prior_tol) || prior_tol < 0)
    stop("prior_tol must be a single nonnegative finite numeric.")
  if (!is.null(coverage) &&
      (!is.numeric(coverage) || length(coverage) != 1L ||
       !is.finite(coverage) || coverage <= 0 || coverage >= 1))
    stop("coverage must be NULL or a single number between 0 and 1.")

  nw <- normalize_null_weight(null_weight, prior_weights, p)
  prior_weights <- normalize_prior_weights(nw$prior_weights,
                                            p + as.integer(nw$add_null))

  variable_names <- input$variable_names
  betahat <- input$betahat
  shat2 <- input$shat2
  scale_factors <- input$scale_factors
  if (nw$add_null) {
    betahat <- c(betahat, 0)
    shat2 <- c(shat2, Inf)
    scale_factors <- c(scale_factors, 1)
    variable_names <- c(variable_names, "null")
  }

  model <- initialize_ser_model(prior_weights, input$V_init,
                                input$sigma2, nw$null_weight)
  V <- model$V[1]

  lbf_model <- function(V_val) {
    lbf <- gaussian_ser_lbf(betahat, shat2, V_val)
    apply_ser_lbf(model, lbf, shat2, l = NULL)$lbf_model
  }

  V <- optimize_scalar_prior_variance(
    V_init = V,
    estimate_prior_method = estimate_prior_method,
    neg_loglik_fn = function(logV) -lbf_model(exp(logV)),
    loglik_fn = lbf_model,
    optim_init = log(max(c(betahat^2 - shat2, 1), na.rm = TRUE)),
    optim_bounds = c(-30, 15),
    optim_scale = "log",
    check_null_threshold = check_null_threshold)

  lbf <- gaussian_ser_lbf(betahat, shat2, V)
  model <- apply_ser_lbf(model, lbf, shat2, l = 1L)$model
  model <- store_ser_moments(model, 1L,
                             gaussian_ser_moments(betahat, shat2, V))
  model$V[1] <- V

  e_loglik <- gaussian_ser_posterior_e_loglik(
    model$alpha[1, ], model$mu[1, ], model$mu2[1, ], betahat, shat2)
  model$KL[1] <- -model$lbf[1] + e_loglik
  model$elbo <- model$lbf[1]

  model$niter <- 1L
  model$converged <- TRUE
  model$X_column_scale_factors <- scale_factors
  model$intercept <- NA_real_
  model$fitted <- NULL
  model$sets <- NULL
  model$z <- input$z
  model$pve_adjustment <- input$pve_adjustment
  model$null_index <- if (nw$add_null) length(betahat) else 0L
  class(model) <- "susie"

  model$pip <- susie_get_pip(model, prior_tol = prior_tol)

  if (!is.null(coverage)) {
    model$sets <- suppressMessages(
      susie_get_cs_attainable(model, coverage = coverage, ethres = ethres)
    )
    warning_message(
      "susie_ser() uses a one-effect credible-set model as in Maller ",
      "et al. 2012. Because no LD reference is supplied, credible sets ",
      "are reported without purity filtering; susie_ser() applies the ",
      "attainable-coverage filter from SparsePro (Zhang et al. 2023). ",
      "Set coverage = NULL to skip CS construction and silence this hint, ",
      "or call susie_get_cs(fit, Xcorr = R, ...) afterward to apply a ",
      "purity filter.",
      style = "hint"
    )
  }

  model <- name_ser_output(model, variable_names)
  return(model)
}

#' @keywords internal
#' @noRd
prepare_ser_summary <- function(z, bhat, shat, n, var_y,
                                prior_variance, scaled_prior_variance) {
  if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 1)
      stop("n must be a single number greater than 1.")
  }
  if (!is.null(var_y)) {
    if (!is.numeric(var_y) || length(var_y) != 1L ||
        !is.finite(var_y) || var_y <= 0)
      stop("var_y must be a single positive finite numeric.")
  }
  if (!is.numeric(prior_variance) || length(prior_variance) != 1L ||
      !is.finite(prior_variance) || prior_variance <= 0)
    stop("prior_variance must be a single positive finite numeric.")
  if (!is.numeric(scaled_prior_variance) ||
      length(scaled_prior_variance) != 1L ||
      !is.finite(scaled_prior_variance) || scaled_prior_variance <= 0)
    stop("scaled_prior_variance must be a single positive finite numeric.")

  summary_input <- normalize_summary_stats_input(z = z, bhat = bhat,
                                                 shat = shat)
  z <- summary_input$z
  bhat <- summary_input$bhat
  shat <- summary_input$shat
  variable_names <- summary_input$variable_names
  if (!is.null(bhat) && (any(!is.finite(bhat)) || any(!is.finite(shat))))
    stop("bhat and shat must be finite.")
  if (any(is.infinite(z)))
    stop("z contains infinite values.")
  z[is.na(z)] <- 0

  z_unadjusted <- z
  pve_adjustment <- NULL
  if (!is.null(n)) {
    pve <- apply_pve_adjustment(z, n)
    z <- pve$z
    pve_adjustment <- pve$adjustment
  }

  working <- summary_stats_working_quantities(
    z = z, n = n, shat = shat, var_y = var_y,
    pve_adjustment = pve_adjustment, prior_variance = prior_variance,
    scaled_prior_variance = scaled_prior_variance)

  if (is.null(variable_names))
    variable_names <- rep("", length(z))

  list(betahat = working$ser_betahat, shat2 = working$ser_shat2,
       sigma2 = working$ser_sigma2, V_init = working$ser_V_init,
       scale_factors = working$ser_scale_factors, z = z_unadjusted,
       pve_adjustment = pve_adjustment,
       variable_names = variable_names)
}

#' @keywords internal
#' @noRd
initialize_ser_model <- function(prior_weights, V, sigma2, null_weight) {
  p <- length(prior_weights)
  list(alpha        = matrix(prior_weights, 1L, p),
       mu           = matrix(0, 1L, p),
       mu2          = matrix(0, 1L, p),
       V            = V,
       KL           = NA_real_,
       lbf          = NA_real_,
       lbf_variable = matrix(NA_real_, 1L, p),
       sigma2       = sigma2,
       pi           = prior_weights,
       null_weight  = null_weight)
}

#' @keywords internal
#' @noRd
name_ser_output <- function(model, variable_names) {
  variable_names[is.na(variable_names)] <- ""
  if (length(variable_names) != ncol(model$alpha))
    return(model)

  if (model$null_index > 0L) {
    pip_names <- variable_names[-model$null_index]
  } else {
    pip_names <- variable_names
  }
  if (all(!nzchar(pip_names)))
    return(model)

  colnames(model$alpha) <- variable_names
  colnames(model$mu) <- variable_names
  colnames(model$mu2) <- variable_names
  colnames(model$lbf_variable) <- variable_names
  names(model$pi) <- variable_names
  names(model$X_column_scale_factors) <- variable_names
  names(model$pip) <- pip_names
  names(model$z) <- pip_names
  return(model)
}
