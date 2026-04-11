# =============================================================================
# SS MIXTURE PANEL METHODS
#
# Class c("ss_mixture", "ss"). Inherits ALL SER/ELBO from ss path.
# Overrides: 5 methods for omega-aware state management.
#
# Model: y ~ N(X(omega)*beta, sigma2*I), X'X = (n-1)*R(omega)
# mu on betahat scale (same as ss). Omega evaluators use z-score scale.
# =============================================================================

# (n-1)*R(omega)*v using current X_meta (or fallback to data$X)
#' @keywords internal
compute_XtXv_mixture <- function(data, model, v) {
  # Use panel_R for accurate R*v (cov2cor-based, not standardize_X)
  if (!is.null(model$omega) && !is.null(data$panel_R)) {
    Rv <- Reduce("+", Map(function(w, R) w * (R %*% v), model$omega, data$panel_R))
    return(data$nm1 * as.vector(Rv))
  }
  # Fallback: data$X = sqrt(n-1)*X_meta_init
  as.vector(compute_Rv(data, v))
}

# 1. Initialize fitted values
#' @keywords internal
initialize_fitted.ss_mixture <- function(data, mat_init) {
  list(XtXr = as.vector(compute_Rv(data, colSums(mat_init$alpha * mat_init$mu))))
}

# 2. Compute residuals using current R(omega)
#' @keywords internal
compute_residuals.ss_mixture <- function(data, params, model, l, ...) {
  sw_l <- get_slot_weight(model, l)
  bl <- model$alpha[l, ] * model$mu[l, ]
  XtXr_without_l <- model$XtXr - sw_l * compute_XtXv_mixture(data, model, bl)

  model$residuals         <- data$Xty - XtXr_without_l
  model$fitted_without_l  <- XtXr_without_l
  model$residual_variance <- model$sigma2
  model$predictor_weights <- rep(data$nm1, data$p)

  if (!is.null(data$sketch_B)) {
    # Convert to z-score scale for inflation computation
    sqnm1 <- sqrt(data$nm1)
    sw <- if (!is.null(model$slot_weights)) model$slot_weights else rep(1, nrow(model$alpha))
    b_minus_l_z <- (colSums(sw * model$alpha * model$mu) - sw_l * bl) * sqnm1
    Rz_without_l_z <- XtXr_without_l / sqnm1
    model$shat2_inflation <- compute_shat2_inflation_rss(
      data, model, Rz_without_l_z, b_minus_l_z)
  }
  return(model)
}

# 3. Update fitted values + precompute z-score quantities for omega
#' @keywords internal
update_fitted_values.ss_mixture <- function(data, params, model, l, ...) {
  sw_l <- get_slot_weight(model, l)
  bl <- model$alpha[l, ] * model$mu[l, ]
  model$XtXr <- model$fitted_without_l + sw_l * compute_XtXv_mixture(data, model, bl)

  # Convert betahat-scale mu to z-score scale for omega evaluators.
  # Weight by slot_weights (c_hat) when active.
  sqnm1 <- sqrt(data$nm1)
  sw <- if (!is.null(model$slot_weights)) model$slot_weights else rep(1, nrow(model$alpha))
  model$Z           <- sw * model$alpha * model$mu * sqnm1
  model$zbar        <- colSums(model$Z)
  model$diag_postb2 <- colSums(sw * model$alpha * model$mu2 * data$nm1)
  return(model)
}

# 4. Update variance: sigma2 (via default ss chain) + omega M-step
#' @keywords internal
update_model_variance.ss_mixture <- function(data, params, model) {
  # Sigma2: reuse default chain (est_residual_variance + bounds)
  if (isTRUE(params$estimate_residual_variance)) {
    model <- update_model_variance.default(data, params, model)
  }

  # Omega M-step
  if (!is.null(data$K) && data$K > 1 && !isTRUE(model$omega_converged)) {
    omega_cur <- if (!is.null(model$omega)) model$omega else rep(1 / data$K, data$K)

    eval_omega <- NULL
    if (!is.null(data$omega_cache)) {
      cache <- data$omega_cache
      iter_cache <- precompute_omega_iteration(cache, model$zbar,
                                                model$diag_postb2, model$Z)
      eval_omega <- function(w) {
        eval_omega_eloglik_reduced(cache, w, iter_cache,
                                    model$sigma2, 0, data$K, data$p)
      }
    } else if (!is.null(data$panel_R)) {
      eval_omega <- function(w) {
        eval_omega_eloglik_R(data$panel_R, w, data$z, model$zbar,
                              model$diag_postb2, model$Z, model$sigma2,
                              0, data$K, data$p)
      }
    }

    if (!is.null(eval_omega)) {
      opt <- optimize_omega(eval_omega, omega_cur, data$K)
      model$omega <- opt$omega
      # Recompute XtXr with updated R(omega)
      b_bar <- colSums(model$alpha * model$mu)
      model$XtXr <- compute_XtXv_mixture(data, model, b_bar)
      if (opt$converged) model$omega_converged <- TRUE
    }
  }
  return(model)
}

# 5. ER2 using current R(omega), not stale data$X
#' @keywords internal
get_ER2.ss_mixture <- function(data, model) {
  B       <- model$alpha * model$mu
  betabar <- colSums(B)
  postb2  <- model$alpha * model$mu2

  XtX_betabar <- compute_XtXv_mixture(data, model, betabar)
  XB2 <- 0
  for (l in seq_len(nrow(B))) {
    bl <- B[l, ]
    XB2 <- XB2 + sum(bl * compute_XtXv_mixture(data, model, bl))
  }

  data$yty - 2 * sum(betabar * data$Xty) + sum(betabar * XtX_betabar) -
    XB2 + data$nm1 * sum(postb2)
}
