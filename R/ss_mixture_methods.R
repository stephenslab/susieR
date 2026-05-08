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
  if (!is.null(data$panel_R)) {
    omega <- get_mixture_omega(data, model)
    Rv <- Reduce("+", Map(function(w, R) w * (R %*% v), omega, data$panel_R))
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

  if (!is.null(data$R_finite_B) && model$sigma2 > .Machine$double.eps) {
    # Region-level scalar lambda_bias is set by fit_R_mismatch once per
    # IBSS sweep; here we just apply it through the slot-specific
    # xi_l = eta_l^2 + v_g,l on z-scale.
    sw <- if (!is.null(model$slot_weights)) model$slot_weights else
            rep(1, nrow(model$alpha))
    b_minus_l <- colSums(sw * model$alpha * model$mu) - sw_l * bl
    nm1  <- data$nm1
    v_g  <- max(sum(b_minus_l * XtXr_without_l), 0)
    xi_l <- XtXr_without_l^2 / nm1 + v_g
    lambda_bias <- if (is.null(model$lambda_bias)) 0 else model$lambda_bias
    R_finite_B <- get_current_R_finite_B(data, model)
    model$shat2_inflation <- 1 + (1 / R_finite_B + lambda_bias) *
                                  xi_l / model$sigma2
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
    omega_cur <- get_mixture_omega(data, model)

    # Omega-objective ridge: small floor used ONLY inside the Eloglik
    # evaluator to stabilize log|sigma2*A(omega)| near rank-deficient
    # vertices. Without it, small eigenvalues of A(omega) produce a huge
    # -0.5 * log|.| penalty at vertex omegas, pulling the optimizer toward
    # the interior (collapse to ~uniform weights). Matches the behavior
    # of the prev rss_lambda path with auto lambda = 1/(n-1). Does NOT
    # affect the ss-SER update, which still uses lambda = 0 (no FDR
    # inflation in the credible-set inference).
    omega_ridge <- 1 / data$nm1
    eval_omega <- NULL
    if (!is.null(data$omega_cache)) {
      cache <- data$omega_cache
      iter_cache <- precompute_omega_iteration(cache, model$zbar,
                                                model$diag_postb2, model$Z)
      eval_omega <- function(w) {
        eval_omega_eloglik_reduced(cache, w, iter_cache,
                                    model$sigma2, omega_ridge, data$K, data$p)
      }
    } else if (!is.null(data$panel_R)) {
      eval_omega <- function(w) {
        eval_omega_eloglik_R(data$panel_R, w, data$z, model$zbar,
                              model$diag_postb2, model$Z, model$sigma2,
                              omega_ridge, data$K, data$p)
      }
    }

    if (!is.null(eval_omega)) {
      opt <- optimize_omega(eval_omega, omega_cur, data$K)
      model$omega <- opt$omega
      if (!is.null(data$R_finite_B) && !is.null(data$B_list))
        model$R_finite_B <- get_current_R_finite_B(data, model)
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
