# RSS R-reference mismatch handling.
#
# Single home for code that targets the discrepancy between the
# supplied R reference and the target population. Active on the SS
# / ss_mixture dispatches; the rss_lambda dispatch (lambda > 0) does
# NOT use any of this (entry-level errors block lambda > 0 with
# R_finite or R_mismatch != "none").
#
#   * 1-D estimator for the variance component lambda_bias
#     (estimate_lambda_bias)
#   * per-variable inflation factor used inside the SER step
#     (compute_shat2_inflation)
#   * model-state storage helper for per-variable inflation
#     (apply_inflation_state)
#   * SER-protected initialization for the recommended EB path
#     (initialize_R_mismatch)
#   * per-iteration region-level fit (fit_R_mismatch)
#   * residual R-mismatch QC diagnostic Q_art (always with R_mismatch)
#   * Bayes-factor attenuation diagnostic for CS reliability
#
# Storage convention on the model:
#   model$lambda_bias    scalar set once per iteration by fit_R_mismatch
#   model$B_corrected    1 / (1/B + lambda_bias)
#   model$shat2_inflation per-variable inflation vector of length p,
#                        consumed by the SER step.

# =============================================================================
# BAYES-FACTOR ATTENUATION DIAGNOSTIC
# =============================================================================

# Store the log-BF attenuation caused by the current R-uncertainty inflation.
# The adjusted log BF is the one used by the SER update; the nominal log BF
# recomputes the same SER with the uninflated standard errors.
#' @keywords internal
record_R_bf_attenuation <- function(model, ser_stats, lbf_adjusted, V, l) {
  if (is.null(l) || is.null(model$shat2_inflation))
    return(model)
  infl <- model$shat2_inflation
  if (length(infl) != length(ser_stats$shat2) ||
      !any(is.finite(infl) & infl > 1))
    return(model)

  shat2_nominal <- ser_stats$shat2 / infl
  lbf_nominal <- gaussian_ser_lbf(ser_stats$betahat, shat2_nominal, V)
  delta <- lbf_nominal - lbf_adjusted
  atten <- pmax(delta, 0)
  atten[!is.finite(atten)] <- NA_real_

  if (is.null(model$R_bf_attenuation)) {
    model$R_bf_attenuation <- matrix(NA_real_, nrow(model$alpha),
                                     length(atten))
  }
  model$R_bf_attenuation[l, ] <- atten
  model
}

# Summarize final BF attenuation after credible sets are available.
#' @keywords internal
summarize_R_bf_attenuation <- function(model, threshold = log(20)) {
  D <- model$R_bf_attenuation
  if (is.null(D) || is.null(model$R_finite_diagnostics))
    return(model)

  label_atten <- function(x) {
    out <- rep(NA_character_, length(x))
    out[is.finite(x) & x < log(3)] <- "stable"
    out[is.finite(x) & x >= log(3) & x < log(20)] <- "mildly_sensitive"
    out[is.finite(x) & x >= log(20) & x < log(150)] <- "sensitive"
    out[is.finite(x) & x >= log(150)] <- "highly_sensitive"
    out
  }

  variant_max <- apply(D, 2, function(x) {
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  })
  variant_component <- apply(D, 2, function(x) {
    if (all(is.na(x))) NA_integer_ else which.max(replace(x, is.na(x), -Inf))
  })

  cs_max <- cs_weighted <- numeric(0)
  cs_top_variable <- integer(0)
  cs_top_attenuation <- numeric(0)
  cs_label <- character(0)
  if (!is.null(model$sets$cs) && length(model$sets$cs) > 0) {
    cs_index <- model$sets$cs_index
    if (is.null(cs_index))
      cs_index <- as.integer(sub("^L", "", names(model$sets$cs)))
    cs_max <- cs_weighted <- cs_top_attenuation <-
      rep(NA_real_, length(model$sets$cs))
    cs_top_variable <- rep(NA_integer_, length(model$sets$cs))
    for (i in seq_along(model$sets$cs)) {
      l <- cs_index[i]
      vars <- model$sets$cs[[i]]
      if (!is.na(l) && l >= 1 && l <= nrow(D) && length(vars) > 0) {
        vals <- D[l, vars]
        cs_max[i] <- if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE)
        cs_weighted[i] <- sum(model$alpha[l, vars] * vals, na.rm = TRUE)
        if (any(is.finite(vals))) {
          j <- which.max(replace(vals, !is.finite(vals), -Inf))
          cs_top_variable[i] <- vars[j]
          cs_top_attenuation[i] <- vals[j]
        }
      }
    }
    cs_label <- label_atten(cs_max)
    names(cs_max) <- names(cs_weighted) <- names(cs_label) <-
      names(cs_top_variable) <- names(cs_top_attenuation) <- names(model$sets$cs)
  }

  flag <- any(is.finite(cs_max) & cs_max >= threshold)

  d <- model$R_finite_diagnostics
  d$bf_attenuation <- list(
    variant_max = variant_max,
    variant_component = variant_component,
    variant_label = label_atten(variant_max),
    cs_max = cs_max,
    cs_weighted = cs_weighted,
    cs_label = cs_label,
    cs_top_variable = cs_top_variable,
    cs_top_attenuation = cs_top_attenuation,
    threshold = threshold
  )
  d$R_sensitivity_flag <- flag
  model$R_finite_diagnostics <- d
  model
}

# =============================================================================
# FINITE-REFERENCE SETUP AND DIAGNOSTICS
# =============================================================================

# Resolve R_finite into an explicit reference sample size B. FALSE is an
# explicit "off" setting; NULL means unspecified. R_finite = TRUE is only
# meaningful when the reference factor X is available; for precomputed R,
# the caller must provide B explicitly.
#' @keywords internal
resolve_R_finite <- function(R_finite, X = NULL, is_multi_panel = FALSE) {
  if (is.null(R_finite))
    return(NULL)
  if (identical(R_finite, FALSE))
    return(NULL)
  if (isTRUE(R_finite)) {
    if (is.null(X))
      stop("R_finite = TRUE requires X input. ",
           "When using a precomputed R matrix, provide a positive number ",
           "specifying the reference sample size B instead.")
    if (is_multi_panel)
      return(vapply(X, nrow, integer(1)))
    return(nrow(X))
  }
  if (!is.numeric(R_finite) || any(!is.finite(R_finite)) ||
      any(R_finite <= 0)) {
    stop("R_finite must be NULL, FALSE, TRUE, or positive numeric value(s).")
  }
  if (is_multi_panel) {
    K <- if (is.null(X)) length(R_finite) else length(X)
    if (length(R_finite) == 1)
      return(rep(as.numeric(R_finite), K))
    if (length(R_finite) == K)
      return(as.numeric(R_finite))
    stop("For multi-panel input, R_finite must be FALSE, TRUE, a single ",
         "positive number, or one positive number per panel.")
  }
  if (length(R_finite) == 1)
    return(as.numeric(R_finite))
  stop("R_finite must be NULL, FALSE, TRUE, or a single positive number.")
}

# Compute finite-reference R diagnostics (debiased Frobenius norm,
# effective rank, r/B ratio, per-variant diagonal deviation from 1).
# Used by the summary-statistics constructors.
#
# @param X Factor matrix (B x p), or NULL.
# @param R Precomputed R matrix (p x p), or NULL.
# @param B Reference panel sample size.
# @param p Number of variants.
# @param x_is_standardized If TRUE, X has been standardized so X'X = R_hat
#   directly (no normalization). If FALSE, R_hat = X'X/B so the Frobenius
#   norm needs a /B^2 correction.
# @return List with B, p, R_frob_sq_debiased, effective_rank, r_over_B,
#   Rhat_diag_deviation.
#' @keywords internal
compute_R_finite_diagnostics <- function(X = NULL, R = NULL, B, p,
                                         x_is_standardized = FALSE) {
  if (!is.null(X)) {
    A <- tcrossprod(X)           # B x B Gram matrix
    R_frob_sq <- sum(A * A)      # ||XX'||_F^2 = ||X'X||_F^2
    if (!x_is_standardized)
      R_frob_sq <- R_frob_sq / nrow(X)^2
    Rhat_diag <- colSums(X^2)
    if (!x_is_standardized)
      Rhat_diag <- Rhat_diag / nrow(X)
  } else if (!is.null(R)) {
    R_frob_sq <- sum(R * R)
    Rhat_diag <- diag(R)
  } else {
    R_frob_sq <- p               # identity fallback
    Rhat_diag <- rep(1, p)
  }

  # Debiased Frobenius norm (Ledoit-Wolf unbiased estimator). In the
  # B = Inf limit there is no finite-reference debiasing term.
  R_frob_sq_db <- if (is.infinite(B)) R_frob_sq else
                    (B * R_frob_sq - p^2) / (B + 1)
  eff_rank <- p^2 / max(R_frob_sq_db, 1)

  list(
    B = B,
    p = p,
    R_frob_sq_debiased = R_frob_sq_db,
    effective_rank = eff_rank,
    r_over_B = eff_rank / B,
    Rhat_diag_deviation = abs(Rhat_diag - 1)
  )
}

# Resolve the current mixture weights for ss_mixture objects. During
# initialization model$omega may not yet exist, so fall back to the constructor's
# initial omega rather than an unrelated uniform mixture.
#' @keywords internal
#' @noRd
get_mixture_omega <- function(data, model) {
  if (!inherits(data, "ss_mixture") || is.null(data$K))
    return(NULL)
  omega <- if (!is.null(model[["omega"]])) model[["omega"]] else data[["omega_init"]]
  if (is.null(omega))
    omega <- rep(1 / data$K, data$K)
  omega <- as.numeric(omega)
  if (length(omega) != data$K || any(!is.finite(omega)))
    omega <- rep(1 / data$K, data$K)
  omega[omega < 0] <- 0
  if (sum(omega) <= 0)
    omega <- rep(1 / data$K, data$K)
  omega / sum(omega)
}

# Current finite reference size for ss_mixture uses the same omega as the
# current LD mixture, B_eff(omega) = 1 / sum_k omega_k^2 / B_k.
#' @keywords internal
#' @noRd
get_current_R_finite_B <- function(data, model) {
  if (inherits(data, "ss_mixture") && !is.null(data$B_list)) {
    omega <- get_mixture_omega(data, model)
    return(1 / sum(omega^2 / data$B_list))
  }
  if (!is.null(model$R_finite_B)) model$R_finite_B else data$R_finite_B
}

# =============================================================================
# 1-D ESTIMATOR FOR lambda_bias
# =============================================================================

# Estimate R-bias variance beyond any supplied finite-reference uncertainty.
# With R_finite_B = Inf this is the B^{-1} = 0 limit, so lambda_bias is the
# total continuous R-mismatch variance component.
# Likelihood on the z-score residual scale,
#   tau_j^2 = sigma2 + (1/R_finite_B + lambda_bias) * s_j,
# either by MLE or by MAP with a half-Cauchy prior on u = sqrt(lambda_bias).
# The Fisher-information boundary SE,
#   SE_0 = sqrt(1 / {0.5 * sum((s_j / tau_j0)^2)}),
#   tau_j0 = sigma2 + s_j / R_finite_B,
# defines a data-driven floor: estimates below 0.1 * SE_0 are zeroed.
# This both suppresses Brent boundary noise and replaces ad-hoc display
# thresholds with one rule; "none" short-circuits before optimization.
#' @keywords internal
estimate_lambda_bias <- function(r, s, sigma2, R_finite_B, method,
                                 R_mismatch_method = "mle") {
  if (is.null(method) || method == "none")
    return(0)
  R_mismatch_method <- match.arg(R_mismatch_method,
                                 c("mle", "map"))
  keep <- is.finite(r) & is.finite(s) & s > .Machine$double.eps
  if (!any(keep) || !is.finite(sigma2) || sigma2 <= .Machine$double.eps)
    return(0)

  cache <- list(r2 = r[keep]^2, s = s[keep])
  cache$base <- sigma2 + cache$s / R_finite_B
  pos <- (cache$r2 - cache$base) / cache$s
  pos <- pos[is.finite(pos) & pos > 0]
  fisher_info0 <- 0.5 * sum((cache$s / cache$base)^2)
  se_boundary <- if (fisher_info0 > 0) sqrt(1 / fisher_info0) else Inf

  if (R_mismatch_method == "mle") {
    upper_lambda <- max(c(1, 100 / R_finite_B, 10 * pos), na.rm = TRUE)
    nll <- function(lambda_bias) {
      tau <- cache$base + lambda_bias * cache$s
      0.5 * sum(log(tau) + cache$r2 / tau)
    }
    while (is.finite(upper_lambda) &&
           upper_lambda < .Machine$double.xmax / 10 &&
           nll(upper_lambda) < nll(upper_lambda / 2)) {
      upper_lambda <- upper_lambda * 10
    }
    opt <- optimize(nll, interval = c(0, upper_lambda), tol = 1e-8)
    if (nll(0) <= opt$objective)
      return(0)
    return(if (opt$minimum < 0.1 * se_boundary) 0 else opt$minimum)
  }

  prior_scale <- sqrt(max(1 / R_finite_B, 1 / 10000))
  upper_lambda <- max(c(1, 100 / R_finite_B, 100 * prior_scale^2,
                        10 * pos), na.rm = TRUE)
  upper_u <- sqrt(upper_lambda)

  nll <- function(u) {
    lambda_bias <- u^2
    tau <- cache$base + lambda_bias * cache$s
    0.5 * sum(log(tau) + cache$r2 / tau) + log1p((u / prior_scale)^2)
  }
  lambda_hat <- optimize(nll, interval = c(0, upper_u), tol = 1e-8)$minimum^2

  if (lambda_hat < 0.1 * se_boundary) 0 else lambda_hat
}

# =============================================================================
# PER-VARIABLE INFLATION
# =============================================================================

# SS-path per-variable inflation factor tau_j^2 / sigma2 with
#   tau_j^2 = sigma2 + (1/R_finite_B + lambda_bias) * (eta_j^2 + v_g),
#   eta_j^2 = XtXr_without_l[j]^2 / (n-1)   (z-score scale)
#   v_g     = sum(b_minus_l * XtXr_without_l).
# Reads the region-level scalar lambda_bias from model (set once per
# iteration by fit_R_mismatch) and applies it to the slot-specific xi_l.
# Returns NULL when no inflation applies, otherwise a list with the
# per-variable inflation vector.
#' @keywords internal
compute_shat2_inflation <- function(data, model, XtXr_without_l, b_minus_l, r) {
  R_finite_B <- get_current_R_finite_B(data, model)
  if (is.null(R_finite_B) ||
      model$sigma2 <= .Machine$double.eps) {
    return(NULL)
  }
  v_g     <- max(sum(b_minus_l * XtXr_without_l), 0)
  eta2    <- XtXr_without_l^2 / (data$n - 1)
  s <- eta2 + v_g
  lambda_bias <- if (is.null(model$lambda_bias)) 0 else model$lambda_bias
  infl <- 1 + (1 / R_finite_B + lambda_bias) * s / model$sigma2
  list(infl = infl)
}

# =============================================================================
# MODEL-STATE STORAGE FOR PER-VARIABLE INFLATION
# =============================================================================

# Unpack the inflation list from compute_shat2_inflation into the model.
# Sets model$shat2_inflation to the per-variant inflation vector.
#' @keywords internal
apply_inflation_state <- function(model, infl_state) {
  if (is.null(infl_state)) {
    model$shat2_inflation <- NULL
    return(model)
  }
  model$shat2_inflation <- infl_state$infl
  model
}

# =============================================================================
# PER-SWEEP REGION-LEVEL fit_R_mismatch
# =============================================================================

#' Fit the region-level lambda_bias from a fitted residual.
#'
#' Math (see archive/ld_mismatch_generativemodel.tex):
#'   beta_bar    = colSums(slot_weight * alpha * mu)        (full posterior mean, betahat scale)
#'   XtXr_full   = X'X * beta_bar = (n-1) * R * beta_bar
#'   r_fit       = (data$Xty - XtXr_full) / sqrt(n-1)       (z-scale fitted residual)
#'   eta_fit_j^2 = XtXr_full[j]^2 / (n-1)                   (z-scale per-variant signal)
#'   v_g_fit     = sum(beta_bar * XtXr_full)                (= beta_bar_z' R beta_bar_z)
#'   xi_fit_j    = eta_fit_j^2 + v_g_fit
#' Estimate lambda_bias on the working likelihood
#'   r_fit_j ~ N(0, sigma2 + (1/B + lambda) * xi_fit_j)
#' using the estimator selected by R_mismatch_method.
#'
#' The same lambda_bias fit is followed by the Q_art diagnostic; see
#' compute_Q_art. The diagnostic does not change lambda_bias or the SER
#' likelihood.
#'
#' @keywords internal
#' @noRd
compute_R_mismatch_state <- function(data, params, model, phase = "sweep") {
  R_mismatch <- if (!is.null(params$R_mismatch)) params$R_mismatch else "none"
  if (R_mismatch == "none") return(model)
  R_finite_B <- get_current_R_finite_B(data, model)
  if (is.null(R_finite_B) || !is.finite(model$sigma2) ||
      model$sigma2 <= .Machine$double.eps)
    return(model)
  if (!inherits(data, c("ss", "ss_mixture"))) return(model)

  sw <- if (!is.null(model$slot_weights)) model$slot_weights else
          rep(1, nrow(model$alpha))
  b_full    <- colSums(sw * model$alpha * model$mu)
  XtXr_full <- if (!is.null(model$XtXr))
                 model$XtXr
               else compute_Rv(data, b_full)
  nm1 <- if (!is.null(data$nm1)) data$nm1 else (data$n - 1)
  if (!is.finite(nm1) || nm1 <= 0) return(model)

  r_fit_z  <- (data$Xty - XtXr_full) / sqrt(nm1)
  v_g_full <- max(sum(b_full * XtXr_full), 0)
  s_full   <- XtXr_full^2 / nm1 + v_g_full

  R_mismatch_method <- if (!is.null(params$R_mismatch_method))
                         params$R_mismatch_method else "mle"
  model$lambda_bias <- estimate_lambda_bias(r_fit_z, s_full, model$sigma2,
                                            R_finite_B, R_mismatch,
                                            R_mismatch_method)
  model$B_corrected <- 1 / (1 / R_finite_B + model$lambda_bias)
  model$R_mismatch_method <- R_mismatch_method

  eigen_R <- get_R_mismatch_eigen(data, model)
  if (is.null(eigen_R))
    stop("R_mismatch requires data$eigen_R; ",
         "summary_stats_constructor should have cached it.")
  eig_delta_rel <- if (!is.null(params$eig_delta_rel))
                     params$eig_delta_rel else 1e-3
  eig_delta_abs <- if (!is.null(params$eig_delta_abs))
                     params$eig_delta_abs else 0
  art <- compute_Q_art(eigen_R, r_fit_z, eig_delta_rel, eig_delta_abs)
  threshold <- if (!is.null(params$artifact_threshold))
                 params$artifact_threshold else 0.1
  flagged <- isTRUE(art$evaluable) && isTRUE(art$Q_art > threshold)

  model$Q_art              <- art$Q_art
  model$artifact_evaluable <- art$evaluable
  model$artifact_flag      <- flagged
  model$low_eigen_count    <- art$low_eigen_count
  model$low_eigen_fraction <- art$low_eigen_count /
                              length(eigen_R$values)
  model$eig_delta          <- art$eig_delta

  model$mode_label <- if (flagged) "warning" else "normal"

  if (isTRUE(params$track_fit)) {
    keep <- is.finite(r_fit_z) & is.finite(s_full) &
            s_full > .Machine$double.eps
    r2 <- r_fit_z[keep]^2
    s_keep <- s_full[keep]
    trace_row <- list(
      sweep = if (is.null(model$R_mismatch_trace)) 1L else
                length(model$R_mismatch_trace) + 1L,
      phase = phase,
      R_mismatch = R_mismatch,
      R_mismatch_method = R_mismatch_method,
      lambda_bias = model$lambda_bias,
      B_corrected = model$B_corrected,
      B = R_finite_B,
      sigma2 = model$sigma2,
      n_effects = nrow(model$alpha),
      n_nonzero_lbf = if (!is.null(model$lbf))
                        sum(is.finite(model$lbf) & model$lbf > 0)
                      else NA_integer_,
      mean_r2 = if (length(r2) > 0) mean(r2) else NA_real_,
      median_r2 = if (length(r2) > 0) stats::median(r2) else NA_real_,
      max_r2 = if (length(r2) > 0) max(r2) else NA_real_,
      mean_s = if (length(s_keep) > 0) mean(s_keep) else NA_real_,
      median_s = if (length(s_keep) > 0) stats::median(s_keep) else NA_real_,
      max_s = if (length(s_keep) > 0) max(s_keep) else NA_real_,
      cor_r2_s = if (length(r2) > 1 && stats::sd(r2) > 0 &&
                     stats::sd(s_keep) > 0)
                   suppressWarnings(stats::cor(r2, s_keep,
                                               method = "spearman"))
                 else NA_real_,
      Q_art = if (!is.null(model$Q_art)) model$Q_art else NA_real_,
      artifact_flag = if (!is.null(model$artifact_flag))
                        model$artifact_flag else NA,
      low_eigen_count = if (!is.null(model$low_eigen_count))
                          model$low_eigen_count else NA_integer_
    )
    model$R_mismatch_trace[[trace_row$sweep]] <- trace_row
  }

  model
}

#' SER-protected initialization for the R_mismatch EB path.
#'
#' The joint EB/sparse objective is path dependent. Starting lambda_bias at
#' zero can let secondary R-mismatch patterns enter as sparse effects before
#' the variance component is estimated. For R_mismatch = "eb", use the
#' initialization procedure described in Sun et al. (2026+).
#' R_mismatch = "eb_ser_init" keeps the previous initialization rule and only
#' initializes in the B = Inf limit. R_mismatch = "eb_force_init" always uses
#' the raw one-SER initializer; R_mismatch = "eb_no_init" always skips it.
#'
#' @keywords internal
#' @noRd
initialize_R_mismatch <- function(data, params, model) {
  R_mismatch <- if (!is.null(params$R_mismatch)) params$R_mismatch else "none"
  R_finite_B <- if (!is.null(model$R_finite_B)) model$R_finite_B else data$R_finite_B
  should_init <- R_mismatch %in% c("eb", "eb_force_init") ||
                 (R_mismatch == "eb_ser_init" && is.infinite(R_finite_B))
  if (!should_init || !inherits(data, c("ss", "ss_mixture")) ||
      nrow(model$alpha) < 1)
    return(model)

  model <- single_effect_update(data, params, model, 1L)
  model$R_mismatch_ser_model <- make_R_mismatch_ser_model(data, params, model, 1L)
  model <- compute_R_mismatch_state(data, params, model, phase = "init_ser")
  init_coherence <- if (R_mismatch == "eb") compute_ser_ld_coherence(data, model, model$alpha[1, ]) else 1
  model$lambda_bias <- init_coherence * model$lambda_bias
  model$R_mismatch_init <- list(
    method = "ser",
    ld_coherence = init_coherence,
    lambda_bias = model$lambda_bias,
    Q_art = if (!is.null(model$Q_art)) model$Q_art else NA_real_,
    artifact_flag = if (!is.null(model$artifact_flag))
                      model$artifact_flag else NA
  )
  model
}

#' @keywords internal
#' @noRd
compute_ser_ld_coherence <- function(data, model, alpha) {
  alpha <- as.numeric(alpha)
  alpha[!is.finite(alpha) | alpha < 0] <- 0
  if (sum(alpha) <= 0)
    return(1)
  prior <- as.numeric(model$pi)
  if (length(prior) != length(alpha) || any(!is.finite(prior)) ||
      sum(prior) <= 0)
    prior <- rep(1 / length(alpha), length(alpha))
  prior <- prior / sum(prior)
  keep <- alpha > 0 & alpha >= 2 * prior
  if (!any(keep))
    return(1)
  alpha <- alpha[keep]
  alpha <- alpha / sum(alpha)
  nm1 <- if (!is.null(data$nm1)) data$nm1 else (data$n - 1)
  R <- if (inherits(data, "ss_mixture") && !is.null(data$panel_R)) {
         omega <- get_mixture_omega(data, model)
         Reduce("+", Map(function(w, Rk) {
           w * Rk[keep, keep, drop = FALSE]
         }, omega, data$panel_R))
       } else if (!is.null(data$XtX)) {
         data$XtX[keep, keep, drop = FALSE] / nm1
       } else if (!is.null(data$X)) {
         crossprod(data$X[, keep, drop = FALSE]) / nm1
       } else {
         NULL
       }
  if (is.null(R))
    return(1)
  # Hadamard squared LD, not the spectral matrix square used by Q_art.
  pmin(1, pmax(0, as.numeric(crossprod(alpha, (R * R) %*% alpha))))
}

#' @keywords internal
#' @noRd
fit_R_mismatch <- function(data, params, model) {
  compute_R_mismatch_state(data, params, model, phase = "sweep")
}

#' One-effect diagnostic model for R_mismatch initialization.
#'
#' @keywords internal
#' @noRd
make_R_mismatch_ser_model <- function(data, params, model, l = 1L) {
  ser <- list(
    alpha = model$alpha[l, , drop = FALSE],
    mu = model$mu[l, , drop = FALSE],
    mu2 = model$mu2[l, , drop = FALSE],
    lbf = model$lbf[l],
    lbf_variable = model$lbf_variable[l, , drop = FALSE],
    V = model$V[l],
    KL = if (!is.null(model$KL)) model$KL[l] else NA_real_,
    sigma2 = model$sigma2,
    pi = model$pi,
    null_index = model$null_index,
    niter = 1L,
    converged = NA
  )
  class(ser) <- c("susie", "list")
  ser$pip <- as.vector(ser$alpha[1, ])
  if (!is.null(data$z))
    ser$z <- data$z
  ser$sets <- get_cs(data, params, ser)
  ser <- get_variable_names(data, ser)
  ser
}

# Eigen accessor for R-mismatch QC. The ordinary SS path stores data$eigen_R.
# The ss_mixture path can change R through omega, so recover the current
# mixture spectrum from panel_R when omega is available; otherwise fall
# back to the initialized X_meta crossproduct.
#' @keywords internal
get_R_mismatch_eigen <- function(data, model) {
  if (!is.null(model$eigen_R))
    return(model$eigen_R)
  if (!is.null(data$eigen_R) && !inherits(data, "ss_mixture"))
    return(data$eigen_R)
  if (inherits(data, "ss_mixture")) {
    omega <- get_mixture_omega(data, model)
    if (!is.null(omega) && !is.null(data$omega_cache)) {
      eig <- eigen_from_reduced(data$omega_cache, omega, data$K, data$p)
      eig$values <- pmax(eig$values, 0)
      return(eig)
    }
    if (!is.null(omega) && !is.null(data$panel_R)) {
      R_mix <- Reduce("+", Map(function(w, R) w * R, omega, data$panel_R))
      R_mix <- 0.5 * (R_mix + t(R_mix))
      eig <- eigen(R_mix, symmetric = TRUE)
      eig$values <- pmax(eig$values, 0)
      return(eig)
    }
    if (!is.null(data$X)) {
      R_init <- crossprod(data$X) / data$nm1
      R_init <- 0.5 * (R_init + t(R_init))
      eig <- eigen(R_init, symmetric = TRUE)
      eig$values <- pmax(eig$values, 0)
      return(eig)
    }
  }
  data$eigen_R
}

# =============================================================================
# Q_art residual R-bias artifact diagnostic
# =============================================================================

# Fraction of the fitted residual projected onto low-eigenvalue directions of R.
#   delta = max(eig_delta_abs, eig_delta_rel * max(d))
#   A_delta = {k : d_k <= delta}
#   Q_art = sum_{k in A_delta} (v_k' r_fit)^2 / sum(r_fit^2)
# This extends the column-space check used for z or Xty in the Zou et al.
# (2022) RSS likelihood: the original check asks whether the input summary
# vector lies in the non-zero eigenspace of R; Q_art asks the same question
# of the fitted residual after R_mismatch correction. A large Q_art means
# the residual still has projection in directions where the supplied R is
# nearly singular, so fine-mapping results should be treated with caution.
#
# Returns a list with Q_art (in [0, 1]), evaluable (FALSE when no
# low-eigenvalues exist or r_fit has negligible norm),
# low_eigen_count, eig_delta. Q_art is a heuristic proportion, not a
# calibrated test statistic; see archive/ld_mismatch_generativemodel.tex
# Sec. "Detecting residual R-bias artifacts".
#' @keywords internal
compute_Q_art <- function(eigen_R, r_fit, eig_delta_rel = 1e-3,
                          eig_delta_abs = 0,
                          residual_norm_floor = 1e-12) {
  d <- eigen_R$values
  delta  <- max(eig_delta_abs, eig_delta_rel * max(d))
  proj <- low_eigen_projection_fraction(eigen_R, r_fit, delta,
                                        residual_norm_floor)
  if (!proj$evaluable) {
    return(list(Q_art = 0, evaluable = FALSE,
                low_eigen_count = proj$low_eigen_count, eig_delta = delta))
  }
  list(Q_art = proj$fraction, evaluable = TRUE,
       low_eigen_count = proj$low_eigen_count, eig_delta = delta)
}
