# =============================================================================
# FUNDAMENTAL COMPUTATIONS
#
# Basic mathematical utilities and core RSS computations. These functions
# handle fundamental operations like sufficient statistics computation and
# eigenvalue inverse calculations.
#
# Functions: compute_suff_stat, compute_Dinv
# =============================================================================

#' @title Compute sufficient statistics for input to \code{susie_ss}
#'
#' @description Computes the sufficient statistics \eqn{X'X, X'y, y'y}
#'   and \eqn{n} after centering (and possibly standardizing) the
#'   columns of \eqn{X} and centering \eqn{y} to have mean zero. We also
#'   store the column means of \eqn{X} and mean of \eqn{y}.
#'
#' @param X An n by p matrix of covariates.
#'
#' @param y An n vector.
#'
#' @param standardize Logical flag indicating whether to standardize
#'   columns of X to unit variance prior to computing summary data
#'
#' @return A list of sufficient statistics (\code{XtX, Xty, yty, n})
#'   and \code{X_colmeans}, \code{y_mean}.
#'
#' @importFrom methods as
#' @importFrom Matrix colMeans
#' @importFrom Matrix crossprod
#'
#' @examples
#' data(N2finemapping)
#' ss <- compute_suff_stat(N2finemapping$X, N2finemapping$Y[, 1])
#'
#' @export
#'
compute_suff_stat <- function(X, y, standardize = FALSE) {
  y_mean <- mean(y)
  y <- y - y_mean
  n <- nrow(X)
  mu <- colMeans(X)
  s <- compute_colSds(X)
  Xty <- drop(y %*% X)
  XtX <- crossprod(X)
  XtX <- as.matrix(XtX)
  XtX <- XtX - n * tcrossprod(mu)
  if (standardize) {
    XtX <- XtX / s
    XtX <- t(XtX)
    XtX <- XtX / s
    Xty <- Xty / s
  }
  n <- length(y)
  yty <- sum(y^2)
  return(list(
    XtX = XtX, Xty = Xty, yty = yty, n = n,
    y_mean = y_mean, X_colmeans = mu
  ))
}

# Compute inverse eigenvalues for RSS-lambda methods
#' @keywords internal
compute_Dinv <- function(model, data) {
  eigen_R <- get_eigen_R(data, model)
  Dinv <- 1 / (model$sigma2 * eigen_R$values + data$lambda)
  Dinv[is.infinite(Dinv)] <- 0
  return(Dinv)
}

# Accessor for eigen_R: check model first (multi-panel), fall through to data
#' @keywords internal
get_eigen_R <- function(data, model) {
  if (!is.null(model$eigen_R)) model$eigen_R else data$eigen_R
}

# Accessor for Vtz: check model first (multi-panel), fall through to data
#' @keywords internal
get_Vtz <- function(data, model) {
  if (!is.null(model$Vtz)) model$Vtz else data$Vtz
}

# =============================================================================
# RSS MODEL METHODS
#
# Core RSS algorithm functions including parameter estimation and model
# preprocessing. These implement the mathematical framework for RSS-based
# fine-mapping and handle iteration-specific computations.
#
# Functions: estimate_s_rss, precompute_rss_lambda_terms
# =============================================================================

#' @title Estimate s in \code{susie_rss} Model Using Regularized LD
#'
#' @description The estimated s gives information about the
#'   consistency between the z scores and LD matrix. A larger \eqn{s}
#'   means there is a strong inconsistency between z scores and LD
#'   matrix. The \dQuote{null-mle} method obtains mle of \eqn{s} under
#'   \eqn{z | R ~ N(0,(1-s)R + s I)}, \eqn{0 < s < 1}. The
#'   \dQuote{null-partialmle} method obtains mle of \eqn{s} under
#'   \eqn{U^T z | R ~ N(0,s I)}, in which \eqn{U} is a matrix containing
#'   the of eigenvectors that span the null space of R; that is, the
#'   eigenvectors corresponding to zero eigenvalues of R. The estimated
#'   \eqn{s} from \dQuote{null-partialmle} could be greater than 1. The
#'   \dQuote{null-pseudomle} method obtains mle of \eqn{s} under
#'   pseudolikelihood \eqn{L(s) = \prod_{j=1}^{p} p(z_j | z_{-j}, s,
#'   R)}, \eqn{0 < s < 1}.
#'
#' @param z A p-vector of z scores.
#'
#' @param R A p by p symmetric, positive semidefinite correlation
#'   matrix.
#'
#' @param n The sample size. (Optional, but highly recommended.)
#'
#' @param r_tol Tolerance level for eigenvalue check of positive
#'   semidefinite matrix of R.
#'
#' @param method a string specifies the method to estimate \eqn{s}.
#'
#' @return A number between 0 and 1.
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' p <- 1000
#' beta <- rep(0, p)
#' beta[1:4] <- 0.01
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' X <- scale(X, center = TRUE, scale = TRUE)
#' y <- drop(X %*% beta + rnorm(n))
#' input_ss <- compute_suff_stat(X, y, standardize = TRUE)
#' ss <- univariate_regression(X, y)
#' R <- cor(X)
#' attr(R, "eigen") <- eigen(R, symmetric = TRUE)
#' zhat <- with(ss, betahat / sebetahat)
#'
#' # Estimate s using the unadjusted z-scores.
#' s0 <- estimate_s_rss(zhat, R)
#'
#' # Estimate s using the adjusted z-scores.
#' s1 <- estimate_s_rss(zhat, R, n)
#'
#' @importFrom stats dnorm
#' @importFrom stats optim
#'
#' @export
#'
estimate_s_rss <- function(z, R, n, r_tol = 1e-08, method = "null-mle") {
  # Check and process input arguments z, R.
  z[is.na(z)] <- 0
  if (is.null(attr(R, "eigen"))) {
    attr(R, "eigen") <- eigen(R, symmetric = TRUE)
  }
  eigenld <- attr(R, "eigen")
  if (any(eigenld$values < -r_tol)) {
    warning_message(
      "The matrix R is not positive semidefinite. Negative ",
      "eigenvalues are set to zero"
    )
  }
  eigenld$values[eigenld$values < r_tol] <- 0

  # Check input n, and adjust the z-scores if n is provided.
  if (missing(n)) {
    warning_message(
      "Providing the sample size (n), or even a rough estimate of n, ",
      "is highly recommended. Without n, the implicit assumption is ",
      "n is large (Inf) and the effect sizes are small (close to zero)."
    )
  } else if (n <= 1) {
    stop("n must be greater than 1")
  }
  if (!missing(n)) {
    sigma2 <- (n - 1) / (z^2 + n - 2)
    z <- sqrt(sigma2) * z
  }

  if (method == "null-mle") {
    negloglikelihood <- function(s, ztv, d) {
      0.5 * sum(log((1 - s) * d + s)) +
        0.5 * tcrossprod(ztv / ((1 - s) * d + s), ztv)
    }
    s <- optim(0.5,
               fn = negloglikelihood, ztv = crossprod(z, eigenld$vectors),
               d = eigenld$values, method = "Brent", lower = 0, upper = 1
    )$par
  } else if (method == "null-partialmle") {
    colspace <- which(eigenld$values > 0)
    if (length(colspace) == length(z)) {
      s <- 0
    } else {
      znull <- crossprod(eigenld$vectors[, -colspace], z) # U2^T z
      s <- sum(znull^2) / length(znull)
    }
  } else if (method == "null-pseudomle") {
    pseudolikelihood <- function(s, z, eigenld) {
      precision <- eigenld$vectors %*% (t(eigenld$vectors) *
                                          (1 / ((1 - s) * eigenld$values + s)))
      postmean <- rep(0, length(z))
      postvar <- rep(0, length(z))
      for (i in 1:length(z)) {
        postmean[i] <- -(1 / precision[i, i]) * precision[i, -i] %*% z[-i]
        postvar[i] <- 1 / precision[i, i]
      }
      return(-sum(dnorm(z, mean = postmean, sd = sqrt(postvar), log = TRUE)))
    }
    s <- optim(0.5,
               fn = pseudolikelihood, z = z, eigenld = eigenld,
               method = "Brent", lower = 0, upper = 1
    )$par
  } else {
    stop("The method is not implemented")
  }
  return(s)
}

# Precompute RSS lambda terms that change per IBSS iteration
#' @keywords internal
precompute_rss_lambda_terms <- function(data, model) {
  # Precompute quantities that change per IBSS iteration
  model$Z           <- model$alpha * model$mu
  model$zbar        <- colSums(model$Z)
  model$diag_postb2 <- colSums(model$alpha * model$mu2)

  return(model)
}

# Dynamic stochastic LD variance inflation for rss_lambda path (z-score scale)
#' @keywords internal
compute_shat2_inflation_rss <- function(data, model, Rz_without_l, b_minus_l) {
  # Use model-level B_eff (updated by omega) if available, else data-level
  B_eff <- if (!is.null(model$stochastic_ld_B)) model$stochastic_ld_B else data$stochastic_ld_B
  if (is.null(B_eff) || model$sigma2 <= .Machine$double.eps) return(NULL)
  v_g  <- max(sum(b_minus_l * Rz_without_l), 0)
  eta2 <- Rz_without_l^2   # z-score scale: no (n-1) division needed
  1 + (eta2 + v_g) / (B_eff * model$sigma2)
}

# =============================================================================
# MULTI-PANEL LD MIXTURE
#
# Functions for combining K reference LD panels with learnable convex weights.
# R(omega) = sum_k omega_k R_hat_k, with X_meta = [sqrt(omega_1) X_1; ...].
#
# Key functions:
#   form_X_meta             -- form composite X from K panels with weights
#   eigen_from_X            -- SVD-based eigendecomposition from X matrix
#   precompute_omega_cache  -- joint SVD for reduced-basis optimization
#   precompute_omega_iteration -- per-IBSS-iter bilinear forms
#   eval_omega_eloglik_reduced -- O(r^3) Eloglik evaluator (Cholesky)
#   eigen_from_reduced      -- recover full p-dim eigen from reduced basis
#   eval_omega_eloglik_R    -- O(p^3) reference implementation (testing)
#   optimize_omega          -- simplex optimizer (Grid+Brent or Frank-Wolfe)
# =============================================================================

# Form composite X from K panels with weights omega
# Pre-allocates output to avoid K intermediate copies
#' @keywords internal
form_X_meta <- function(X_list, omega) {
  K   <- length(X_list)
  p   <- ncol(X_list[[1]])
  nrs <- vapply(X_list, nrow, integer(1))
  X_meta <- matrix(0, sum(nrs), p)
  offset <- 0L
  for (k in seq_len(K)) {
    rows <- offset + seq_len(nrs[k])
    X_meta[rows, ] <- sqrt(omega[k]) * X_list[[k]]
    offset <- offset + nrs[k]
  }
  X_meta
}

# SVD-based eigendecomposition from X matrix (X'X = R)
#' @keywords internal
eigen_from_X <- function(X, p) {
  sv <- svd(X, nu = 0)
  eigen_values <- pmax(sv$d^2, 0)
  eigen_vectors <- sv$v
  if (ncol(eigen_vectors) < p) {
    eigen_vectors <- cbind(eigen_vectors,
                           matrix(0, p, p - ncol(eigen_vectors)))
    eigen_values <- c(eigen_values, rep(0, p - length(eigen_values)))
  }
  idx <- order(eigen_values, decreasing = TRUE)
  list(values = eigen_values[idx], vectors = eigen_vectors[, idx])
}

# Precompute reduced-basis quantities for fast omega optimization.
#
# For K panels with sketch matrices X_k (B_k x p), projects all panel
# correlations into a joint reduced basis V_s (p x r) where r = rank of
# [X_1; ...; X_K]. Each Brent evaluation then works on r x r matrices
# (Cholesky + backsolves) instead of p x p eigendecompositions.
#
# Returns a list to be stored in data$omega_cache.
#' @keywords internal
precompute_omega_cache <- function(X_list, z, r_tol = 1e-8) {
  X_stack <- do.call(rbind, X_list)
  sv <- svd(X_stack, nu = 0)
  keep <- sv$d > r_tol
  V_s <- sv$v[, keep, drop = FALSE]
  r <- ncol(V_s)

  # Project each panel into reduced basis: A_k = V_s' R_k V_s (r x r)
  A_list <- lapply(X_list, function(Xk) {
    Zk <- Xk %*% V_s
    crossprod(Zk)
  })

  list(
    V_s = V_s,
    r = r,
    A_list = A_list,
    Vsz = as.vector(crossprod(V_s, z)),
    z_norm2 = sum(z^2)
  )
}

# Precompute per-IBSS-iteration quantities for the omega Brent evaluator.
# Called once per IBSS iteration (not per Brent eval).
#' @keywords internal
precompute_omega_iteration <- function(cache, zbar, diag_postb2, Z) {
  Vsz_bar <- as.vector(crossprod(cache$V_s, zbar))
  ZVs <- Z %*% cache$V_s   # L x r
  M_postb2 <- crossprod(cache$V_s * diag_postb2, cache$V_s)  # r x r

  list(Vsz_bar = Vsz_bar, ZVs = ZVs, M_postb2 = M_postb2)
}

# Evaluate Eloglik at a candidate omega using reduced basis + Cholesky.
#
# Uses precomputed A_k*vector products from precompute_omega_iteration,
# so per-eval work is dominated by the r x r Cholesky + backsolves.
#' @keywords internal
eval_omega_eloglik_reduced <- function(cache, omega, iter_cache,
                                        sigma2, lambda, K, p) {
  r <- cache$r

  # Form A(omega) = sum_k omega_k A_k once; reused for all terms
  A_omega <- omega[1] * cache$A_list[[1]]
  for (k in seq_len(K)[-1])
    A_omega <- A_omega + omega[k] * cache$A_list[[k]]

  # S_r = sigma2 * A(omega) + lambda * I_r
  S_r <- sigma2 * A_omega + lambda * diag(r)

  # Cholesky factorization: O(r^3/6)
  L <- chol(S_r)

  # log|S| = 2*sum(log(diag(L))) + (p - r)*log(lambda)
  logdet_term <- -0.5 * (2 * sum(log(diag(L))) + (p - r) * log(lambda))

  # S_r^{-1} Vsz via backsolve: O(r^2)
  Sinv_Vsz <- backsolve(L, backsolve(L, cache$Vsz, transpose = TRUE))
  z_null_norm2 <- max(cache$z_norm2 - sum(cache$Vsz^2), 0)
  zSinvz <- sum(cache$Vsz * Sinv_Vsz) + z_null_norm2 / lambda

  # term2: -2 * zbar' R(omega) S^{-1} z
  RSinvz_r <- A_omega %*% Sinv_Vsz
  term2 <- -2 * sum(iter_cache$Vsz_bar * RSinvz_r)

  # term3: zbar' R(omega) S^{-1} R(omega) zbar
  A_Vsz_bar <- A_omega %*% iter_cache$Vsz_bar
  Sinv_A_Vsz_bar <- backsolve(L, backsolve(L, A_Vsz_bar, transpose = TRUE))
  term3 <- sum(A_Vsz_bar * Sinv_A_Vsz_bar)

  # term4: -tr(Z' R(omega) S^{-1} R(omega) Z) via backsolve
  A_ZVs_t <- A_omega %*% t(iter_cache$ZVs)
  Sinv_A_ZVs_t <- backsolve(L, backsolve(L, A_ZVs_t, transpose = TRUE))
  term4 <- -sum(A_ZVs_t * Sinv_A_ZVs_t)

  # term5: tr(diag(postb2) R(omega) S^{-1} R(omega)) via A_omega M A_omega
  AMA_omega <- A_omega %*% iter_cache$M_postb2 %*% A_omega
  Sinv_AMA <- backsolve(L, backsolve(L, AMA_omega, transpose = TRUE))
  term5 <- sum(diag(Sinv_AMA))

  ER2 <- zSinvz + term2 + term3 + term4 + term5
  -p / 2 * log(2 * pi) + logdet_term - 0.5 * ER2
}

# Recover full eigendecomposition from reduced basis after omega is chosen.
# Called once per IBSS iteration (after Brent converges), not per eval.
#' @keywords internal
eigen_from_reduced <- function(cache, omega, K, p) {
  A_omega <- omega[1] * cache$A_list[[1]]
  for (k in seq_len(K)[-1])
    A_omega <- A_omega + omega[k] * cache$A_list[[k]]

  eig <- eigen(0.5 * (A_omega + t(A_omega)), symmetric = TRUE)
  d <- pmax(eig$values, 0)
  V_full <- cache$V_s %*% eig$vectors

  if (cache$r < p) {
    V_full <- cbind(V_full, matrix(0, p, p - cache$r))
    d <- c(d, rep(0, p - cache$r))
  }

  list(values = d, vectors = V_full)
}

# Naive O(p^3) Eloglik evaluator (used for testing and as reference).
# Forms R(omega), eigendecomposes the p x p matrix, evaluates Eloglik.
#' @keywords internal
eval_omega_eloglik_R <- function(panel_R, omega, z, zbar, diag_postb2, Z,
                                  sigma2, lambda, K, p) {
  # Form R(omega) = sum_k omega_k R_k
  R_omega <- omega[1] * panel_R[[1]]
  for (k in seq_len(K)[-1]) R_omega <- R_omega + omega[k] * panel_R[[k]]
  R_omega <- 0.5 * (R_omega + t(R_omega))

  # Eigendecompose
  eig <- eigen(R_omega, symmetric = TRUE)
  D <- pmax(eig$values, 0)
  V <- eig$vectors

  # Eloglik computation
  Vtz <- crossprod(V, z)
  z_null_norm2 <- max(sum(z^2) - sum(Vtz^2), 0)
  S_diag <- sigma2 * D + lambda
  Dinv <- ifelse(S_diag > 0, 1 / S_diag, 0)
  DinvD2 <- Dinv * D^2

  logdet_term <- -0.5 * sum(log(S_diag))
  zSinvz <- sum(Dinv * Vtz^2)
  if (lambda > 0) zSinvz <- zSinvz + z_null_norm2 / lambda

  RSinvz <- V %*% (Dinv * D * Vtz)
  term2 <- -2 * sum(zbar * RSinvz)

  Vtzbar <- crossprod(V, zbar)
  term3 <- sum(Vtzbar^2 * DinvD2)

  ZV <- Z %*% V
  term4 <- -sum(ZV^2 %*% DinvD2)

  diag_RSinvR <- colSums(t(V)^2 * DinvD2)
  term5 <- sum(diag_RSinvR * diag_postb2)

  ER2 <- zSinvz + term2 + term3 + term4 + term5
  -p / 2 * log(2 * pi) + logdet_term - 0.5 * ER2
}

# Optimize omega on the K-simplex by maximizing eval_fn.
# Uses the Frank-Wolfe conditional gradient algorithm: each iteration
# evaluates the objective at all K simplex vertices to find the steepest
# ascent direction, then performs a line search (Brent's method) toward
# that vertex.  For K=2 a coarse grid warm-start is prepended since the
# simplex is 1D and the grid cost is negligible.
# Returns list(omega, converged) where converged indicates max|delta| < tol.
#' @keywords internal
#' @importFrom stats optimize
optimize_omega <- function(eval_fn, omega_cur, K,
                           tol = .omega_tol) {
  omega   <- omega_cur
  cur_val <- eval_fn(omega)

  # K=2 warm-start: coarse grid over the 1D simplex
  if (K == 2) {
    grid <- seq(0, 1, tol$grid_spacing)
    vals <- vapply(grid, function(w1) eval_fn(c(w1, 1 - w1)), numeric(1))
    best_w1 <- grid[which.max(vals)]
    omega   <- c(best_w1, 1 - best_w1)
    cur_val <- max(vals)
  }

  # Frank-Wolfe: conditional gradient on simplex with Brent line search
  for (fw_iter in seq_len(tol$fw_max_iter)) {
    vertex_vals <- vapply(seq_len(K), function(k) {
      e_k <- rep(0, K); e_k[k] <- 1; eval_fn(e_k)
    }, numeric(1))
    k_star <- which.max(vertex_vals)
    s <- rep(0, K); s[k_star] <- 1
    opt <- optimize(function(gamma) eval_fn((1 - gamma) * omega + gamma * s),
                    interval = c(0, 1), maximum = TRUE)
    if (opt$objective - cur_val < tol$fw_stop) break
    omega   <- (1 - opt$maximum) * omega + opt$maximum * s
    cur_val <- opt$objective
  }

  converged <- max(abs(omega - omega_cur)) < tol$convergence
  list(omega = omega, converged = converged)
}

# =============================================================================
# DIAGNOSTIC & QUALITY CONTROL
#
# Functions for RSS model diagnostics, data quality assessment, and
# validation. These help users assess the compatibility between z-scores
# and LD matrices and identify potential data issues.
#
# Functions: kriging_rss
# =============================================================================

#' @title Compute Distribution of z-scores of Variant j Given Other z-scores, and Detect Possible Allele Switch Issue
#'
#' @description Under the null, the rss model with regularized LD
#'   matrix is \eqn{z|R,s ~ N(0, (1-s)R + s I))}. We use a mixture of
#'   normals to model the conditional distribution of z_j given other z
#'   scores, \eqn{z_j | z_{-j}, R, s ~ \sum_{k=1}^{K} \pi_k
#'   N(-\Omega_{j,-j} z_{-j}/\Omega_{jj}, \sigma_{k}^2/\Omega_{jj})},
#'   \eqn{\Omega = ((1-s)R + sI)^{-1}}, \eqn{\sigma_1, ..., \sigma_k}
#'   is a grid of fixed positive numbers. We estimate the mixture
#'   weights \eqn{\pi}  We detect the possible allele switch issue
#'   using likelihood ratio for each variant.
#'
#' @param z A p-vector of z scores.
#'
#' @param R A p by p symmetric, positive semidefinite correlation
#'   matrix.
#'
#' @param n The sample size. (Optional, but highly recommended.)
#'
#' @param r_tol Tolerance level for eigenvalue check of positive
#'   semidefinite matrix of R.
#'
#' @param s an estimated s from \code{estimate_s_rss}
#'
#' @return a list containing a ggplot2 plot object and a table. The plot
#'   compares observed z score vs the expected value. The possible allele
#'   switched variants are labeled as red points (log LR > 2 and abs(z) > 2).
#'   The table summarizes the conditional distribution for each variant
#'   and the likelihood ratio test. The table has the following columns:
#'   the observed z scores, the conditional expectation, the conditional
#'   variance, the standardized differences between the observed z score
#'   and expected value, the log likelihood ratio statistics.
#'
#' @importFrom stats dnorm
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_abline
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 aes_string
#' @importFrom mixsqp mixsqp
#'
#' @examples
#' # See also the vignette, "Diagnostic for fine-mapping with summary
#' # statistics."
#' set.seed(1)
#' n <- 500
#' p <- 1000
#' beta <- rep(0, p)
#' beta[1:4] <- 0.01
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' X <- scale(X, center = TRUE, scale = TRUE)
#' y <- drop(X %*% beta + rnorm(n))
#' ss <- univariate_regression(X, y)
#' R <- cor(X)
#' attr(R, "eigen") <- eigen(R, symmetric = TRUE)
#' zhat <- with(ss, betahat / sebetahat)
#' cond_dist <- kriging_rss(zhat, R, n = n)
#' cond_dist$plot
#'
#' @export
#'
kriging_rss <- function(z, R, n, r_tol = 1e-08,
                        s = estimate_s_rss(z, R, n, r_tol, method = "null-mle")) {
  # Check and process input arguments z, R.
  z[is.na(z)] <- 0
  if (is.null(attr(R, "eigen"))) {
    attr(R, "eigen") <- eigen(R, symmetric = TRUE)
  }
  eigenld <- attr(R, "eigen")
  if (any(eigenld$values < -r_tol)) {
    warning_message(
      "The matrix R is not positive semidefinite. Negative ",
      "eigenvalues are set to zero."
    )
  }
  eigenld$values[eigenld$values < r_tol] <- 0

  # Check and progress input argument s.
  force(s)
  if (s > 1) {
    warning_message("The given s is greater than 1. We replace it with 0.8.")
    s <- 0.8
  } else if (s < 0) {
    stop("The s must be non-negative")
  }

  # Check input n, and adjust the z-scores if n is provided.
  if ((!missing(n)) && (n <= 1)) {
    stop("n must be greater than 1")
  }
  if (missing(n)) {
    warning_message(
      "Providing the sample size (n), or even a rough estimate of n, ",
      "is highly recommended. Without n, the implicit assumption is ",
      "n is large (Inf) and the effect sizes are small (close to zero)."
    )
  } else {
    sigma2 <- (n - 1) / (z^2 + n - 2)
    z <- sqrt(sigma2) * z
  }

  dinv <- 1 / ((1 - s) * eigenld$values + s)
  dinv[is.infinite(dinv)] <- 0
  precision <- eigenld$vectors %*% (t(eigenld$vectors) * dinv)
  condmean <- rep(0, length(z))
  condvar <- rep(0, length(z))
  for (i in 1:length(z)) {
    condmean[i] <- -(1 / precision[i, i]) * precision[i, -i] %*% z[-i]
    condvar[i] <- 1 / precision[i, i]
  }
  z_std_diff <- (z - condmean) / sqrt(condvar)

  # obtain grid
  a_min <- 0.8
  if (max(z_std_diff^2) < 1) {
    a_max <- 2
  } else {
    a_max <- 2 * sqrt(max(z_std_diff^2))
  }
  npoint <- ceiling(log2(a_max / a_min) / log2(1.05))
  a_grid <- 1.05^(seq(-npoint, 0)) * a_max

  # compute likelihood
  sd_mtx <- outer(sqrt(condvar), a_grid)
  matrix_llik <- dnorm(z - condmean, sd = sd_mtx, log = TRUE)
  lfactors <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors

  # estimate weight
  w <- mixsqp(matrix_llik, log = TRUE, control = list(verbose = FALSE))$x

  # Compute denominators in likelihood ratios.
  logl0mix <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors

  # Compute numerators in likelihood ratios.
  matrix_llik <- dnorm(z + condmean, sd = sd_mtx, log = TRUE)
  lfactors <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors
  logl1mix <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors

  # Compute (log) likelihood ratios.
  logLRmix <- logl1mix - logl0mix

  z <- drop(z)
  z_std_diff <- drop(z_std_diff)
  res <- data.frame(
    z = z,
    condmean = condmean,
    condvar = condvar,
    z_std_diff = z_std_diff,
    logLR = logLRmix
  )
  p <- ggplot(res, aes(y = .data$z, x = .data$condmean)) +
    geom_point() +
    labs(y = "Observed z scores", x = "Expected value") +
    geom_abline(intercept = 0, slope = 1) +
    theme_bw()
  idx <- which(logLRmix > 2 & abs(z) > 2)
  if (length(idx) > 0) {
    p <- p + geom_point(
      data = res[idx, ],
      aes(y = .data$z, x = .data$condmean), col = "red"
    )
  }
  return(list(plot = p, conditional_dist = res))
}
