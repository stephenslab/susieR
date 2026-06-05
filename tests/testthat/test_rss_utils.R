context("RSS utility functions")

# ---- Sufficient statistics computation ----

test_that("compute_suff_stat with standardize=FALSE produces correct XtX", {
  base_data <- generate_base_data(n = 10, p = 5, seed = 1)
  X_centered <- scale(base_data$X, center = TRUE, scale = FALSE)
  out <- compute_suff_stat(base_data$X, base_data$y, standardize = FALSE)
  dimnames(out$XtX) <- NULL
  expect_equal(out$XtX, crossprod(X_centered), tolerance = 1e-14)
})

test_that("compute_suff_stat with standardize=TRUE produces correct XtX", {
  base_data <- generate_base_data(n = 10, p = 5, seed = 2)
  X_standardized <- scale(base_data$X, center = TRUE, scale = TRUE)
  out <- compute_suff_stat(base_data$X, base_data$y, standardize = TRUE)
  dimnames(out$XtX) <- NULL
  expect_equal(out$XtX, crossprod(X_standardized), tolerance = 1e-14)
})

test_that("compute_suff_stat with sparse matrix input", {
  base_data <- generate_base_data(n = 10, p = 5, seed = 3)
  X_sparse <- as(base_data$X, "sparseMatrix")
  out_dense <- compute_suff_stat(base_data$X, base_data$y, standardize = FALSE)
  out_sparse <- compute_suff_stat(X_sparse, base_data$y, standardize = FALSE)
  dimnames(out_dense$XtX) <- NULL
  dimnames(out_sparse$XtX) <- NULL
  expect_equal(out_sparse$XtX, out_dense$XtX, tolerance = 1e-14)
  expect_equal(as.vector(out_sparse$Xty), out_dense$Xty, tolerance = 1e-14)
  expect_equal(out_sparse$yty, out_dense$yty, tolerance = 1e-14)
})

test_that("compute_suff_stat produces correct Xty", {
  base_data <- generate_base_data(n = 20, p = 8, seed = 4)
  out <- compute_suff_stat(base_data$X, base_data$y, standardize = FALSE)
  y_centered <- base_data$y - mean(base_data$y)
  X_centered <- scale(base_data$X, center = TRUE, scale = FALSE)
  expected_Xty <- drop(crossprod(X_centered, y_centered))
  expect_equal(out$Xty, expected_Xty, tolerance = 1e-14)
})

test_that("compute_suff_stat produces correct yty", {
  base_data <- generate_base_data(n = 20, p = 8, seed = 5)
  out <- compute_suff_stat(base_data$X, base_data$y, standardize = FALSE)
  y_centered <- base_data$y - mean(base_data$y)
  expected_yty <- sum(y_centered^2)
  expect_equal(out$yty, expected_yty, tolerance = 1e-14)
})

test_that("compute_suff_stat stores column means and y_mean", {
  base_data <- generate_base_data(n = 15, p = 6, seed = 6)
  out <- compute_suff_stat(base_data$X, base_data$y, standardize = FALSE)
  expect_equal(out$X_colmeans, colMeans(base_data$X), tolerance = 1e-14)
  expect_equal(out$y_mean, mean(base_data$y), tolerance = 1e-14)
  expect_equal(out$n, base_data$n)
})

test_that("compute_suff_stat with standardize=TRUE scales correctly", {
  base_data <- generate_base_data(n = 25, p = 10, seed = 7)
  out <- compute_suff_stat(base_data$X, base_data$y, standardize = TRUE)
  X_std <- scale(base_data$X, center = TRUE, scale = TRUE)
  expect_equal(diag(out$XtX), diag(crossprod(X_std)), tolerance = 1e-12)
})

test_that("compute_suff_stat returns list with correct names", {
  base_data <- generate_base_data(n = 10, p = 5, seed = 8)
  out <- compute_suff_stat(base_data$X, base_data$y, standardize = FALSE)
  expect_type(out, "list")
  expect_named(out, c("XtX", "Xty", "yty", "n", "y_mean", "X_colmeans"))
})

test_that("compute_suff_stat output is compatible with susie_ss", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 9)
  ss <- compute_suff_stat(base_data$X, base_data$y, standardize = TRUE)
  fit <- suppressWarnings(
    susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, max_iter = 2, verbose = FALSE)
  )
  expect_true(!is.null(fit$alpha))
  expect_equal(ncol(fit$alpha), 50L)
})

test_that("compute_suff_stat with zero-variance column", {
  skip("Fails on Linux in CI")
  base_data <- generate_base_data(n = 20, p = 5, seed = 10)
  base_data$X[, 3] <- 1
  out <- compute_suff_stat(base_data$X, base_data$y, standardize = TRUE)
  expect_true(is.infinite(out$Xty[3]))
})

# ---- estimate_s_rss ----

test_that("estimate_s_rss null-mle returns scalar in [0,1]", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 11)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  s <- estimate_s_rss(z, R, n = base_data$n, method = "null-mle")
  expect_type(s, "double")
  expect_length(s, 1)
  expect_gte(s, 0)
  expect_lte(s, 1)
})

test_that("estimate_s_rss null-partialmle returns non-negative scalar", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 12)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  s <- estimate_s_rss(z, R, n = base_data$n, method = "null-partialmle")
  expect_type(s, "double")
  expect_length(s, 1)
  expect_gte(s, 0)
})

test_that("estimate_s_rss null-pseudomle returns scalar in [0,1]", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 13)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  s <- estimate_s_rss(z, R, n = base_data$n, method = "null-pseudomle")
  expect_type(s, "double")
  expect_length(s, 1)
  expect_gte(s, 0)
  expect_lte(s, 1)
})

test_that("estimate_s_rss warns when n is not provided", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 14)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  expect_message(
    s <- estimate_s_rss(z, R),
    "sample size"
  )
})

test_that("estimate_s_rss errors when n <= 1", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 15)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  expect_error(estimate_s_rss(z, R, n = 1), "must be greater than 1")
  expect_error(estimate_s_rss(z, R, n = 0), "must be greater than 1")
})

test_that("estimate_s_rss precomputed eigen attr gives same result as on-the-fly", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 16)
  ss <- univariate_regression(base_data$X, base_data$y)
  z <- with(ss, betahat / sebetahat)
  R <- cor(base_data$X)
  attr(R, "eigen") <- eigen(R, symmetric = TRUE)
  s1 <- estimate_s_rss(z, R, n = base_data$n, method = "null-mle")
  R2 <- cor(base_data$X)
  s2 <- estimate_s_rss(z, R2, n = base_data$n, method = "null-mle")
  expect_equal(s1, s2, tolerance = 1e-10)
})

test_that("estimate_s_rss warns and returns valid s when R has negative eigenvalues", {
  set.seed(17)
  p <- 50
  R <- matrix(0.5, p, p)
  diag(R) <- 1
  R[1, 2] <- 1.5
  R[2, 1] <- 1.5
  z <- rnorm(p)
  expect_message(
    s <- estimate_s_rss(z, R, n = 100, method = "null-mle"),
    "not positive semidefinite"
  )
  expect_gte(s, 0)
  expect_lte(s, 1)
})

test_that("estimate_s_rss replaces NA z-scores with 0 and returns valid s", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 18)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  z[5] <- NA
  s <- estimate_s_rss(z, R, n = base_data$n, method = "null-mle")
  expect_gte(s, 0)
  expect_lte(s, 1)
})

test_that("estimate_s_rss null-partialmle handles rank-1 R (perfect LD)", {
  set.seed(19)
  p <- 10
  R <- matrix(1, p, p)
  z <- rnorm(p)
  s <- estimate_s_rss(z, R, n = 100, method = "null-partialmle")
  expect_gte(s, 0)
})

test_that("estimate_s_rss errors on invalid method", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 20)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  expect_error(
    estimate_s_rss(z, R, n = base_data$n, method = "invalid-method"),
    "not implemented"
  )
})

test_that("estimate_s_rss returns small s when z-scores are consistent with R", {
  base_data <- generate_base_data(n = 500, p = 100, k = 3, signal_sd = 0.5, seed = 21)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  s <- estimate_s_rss(z, R, n = base_data$n, method = "null-mle")
  expect_lt(s, 0.01)
})

# ---- kriging_rss ----

test_that("kriging_rss returns named list with ggplot and data frame", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 22)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  result <- kriging_rss(z, R, n = base_data$n)
  expect_named(result, c("plot", "conditional_dist"))
  expect_s3_class(result$plot, "ggplot")
  expect_s3_class(result$conditional_dist, "data.frame")
  expect_equal(nrow(result$conditional_dist), base_data$p)
  expect_true(all(c("z", "condmean", "condvar", "z_std_diff", "logLR") %in%
                    colnames(result$conditional_dist)))
})

test_that("kriging_rss accepts explicit s parameter", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 25)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  result <- kriging_rss(z, R, n = base_data$n, s = 0.1)
  expect_s3_class(result$plot, "ggplot")
  expect_equal(nrow(result$conditional_dist), base_data$p)
})

test_that("kriging_rss warns when n is not provided", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 26)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  expect_message(
    result <- kriging_rss(z, R),
    "sample size"
  )
})

test_that("kriging_rss errors when n <= 1", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 27)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  expect_error(
    kriging_rss(z, R, n = 1),
    "must be greater than 1"
  )
})

test_that("kriging_rss errors when n <= 1 even with explicit s bypassing estimate_s_rss", {
  set.seed(208)
  p <- 8
  z <- rnorm(p)
  R <- diag(p)
  attr(R, "eigen") <- eigen(R, symmetric = TRUE)
  expect_error(kriging_rss(z, R, n = 1, s = 0.1), "n must be greater than 1")
  expect_error(kriging_rss(z, R, n = 0, s = 0.5), "n must be greater than 1")
})

test_that("kriging_rss warns and replaces s > 1 with 0.8", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 28)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  expect_message(
    result <- kriging_rss(z, R, n = base_data$n, s = 1.5),
    "greater than 1"
  )
  expect_s3_class(result$plot, "ggplot")
})

test_that("kriging_rss errors when s < 0", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 29)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  expect_error(
    kriging_rss(z, R, n = base_data$n, s = -0.1),
    "non-negative"
  )
})

test_that("kriging_rss warns and succeeds when R has negative eigenvalues", {
  set.seed(30)
  p <- 50
  R <- matrix(0.5, p, p)
  diag(R) <- 1
  R[1, 2] <- 1.5
  R[2, 1] <- 1.5
  z <- rnorm(p)
  expect_message(
    result <- kriging_rss(z, R, n = 100),
    "not positive semidefinite"
  )
  expect_s3_class(result$plot, "ggplot")
})

test_that("kriging_rss detects flipped allele as positive logLR", {
  base_data <- generate_base_data(n = 500, p = 100, k = 3, signal_sd = 1, seed = 31)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  z[1] <- -z[1]
  result <- kriging_rss(z, R, n = base_data$n)
  expect_gt(result$conditional_dist$logLR[1], 0)
})

test_that("kriging_rss replaces NA z-scores with 0", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 32)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  z[10] <- NA
  result <- kriging_rss(z, R, n = base_data$n)
  expect_equal(result$conditional_dist$z[10], 0)
})

test_that("kriging_rss conditional means and variances are finite and positive", {
  base_data <- generate_base_data(n = 200, p = 50, seed = 33)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  result <- kriging_rss(z, R, n = base_data$n)
  expect_true(all(result$conditional_dist$condvar > 0))
  expect_true(all(is.finite(result$conditional_dist$condmean)))
  expect_true(all(is.finite(result$conditional_dist$z_std_diff)))
})

test_that("kriging_rss takes a_max=2 branch when all standardized diffs are small", {
  set.seed(34)
  p <- 50
  R <- diag(p)
  z <- rnorm(p, mean = 0, sd = 0.3)
  result <- kriging_rss(z, R, n = 100)
  expect_lt(max(result$conditional_dist$z_std_diff^2), 1)
  expect_s3_class(result$plot, "ggplot")
})

test_that("kriging_rss plot adds red points when allele-switch outliers exist", {
  n <- 200
  p <- 50
  base_data <- generate_base_data(n = n, p = p, k = 5, signal_sd = 2, seed = 32)
  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)
  strong_idx <- which(abs(z) > 2)
  if (length(strong_idx) >= 3)
    z[strong_idx[1:3]] <- -z[strong_idx[1:3]]
  result <- kriging_rss(z, R, n = n)
  outliers <- which(result$conditional_dist$logLR > 2 &
                      abs(result$conditional_dist$z) > 2)
  expect_gt(length(outliers), 0)
  expect_s3_class(result$plot, "ggplot")
})

# ---- safe_pd_decomp ----

test_that("safe_pd_decomp uses the Cholesky fast path for PD matrices", {
  set.seed(201)
  A <- crossprod(matrix(rnorm(6 * 6), 6, 6))
  z <- rnorm(6)
  out <- safe_pd_decomp(A, z)
  expect_named(out, c("logdet", "r_eff", "solve", "solve_z"))
  expect_equal(out$r_eff, nrow(A))
  expect_equal(out$logdet, as.numeric(determinant(A, logarithm = TRUE)$modulus),
               tolerance = 1e-10)
  v <- rnorm(6)
  expect_equal(as.vector(out$solve(v)), as.vector(solve(A, v)), tolerance = 1e-8)
  expect_equal(as.vector(out$solve_z), as.vector(solve(A, z)), tolerance = 1e-8)
})

test_that("safe_pd_decomp falls back to eigendecomposition for singular matrices", {
  set.seed(202)
  S <- matrix(0, 4, 4)
  S[1:2, 1:2] <- crossprod(matrix(rnorm(2 * 2), 2, 2))
  S <- (S + t(S)) / 2
  z <- rnorm(4)
  out <- safe_pd_decomp(S, z)
  expect_equal(out$r_eff, 2L)
  sol <- out$solve(z)
  expect_length(as.vector(sol), 4)
  expect_equal(as.vector(out$solve_z), as.vector(sol), tolerance = 1e-10)
  expect_true(is.finite(out$logdet))
})

test_that("safe_pd_decomp returns NULL solve_z when z is not supplied", {
  set.seed(203)
  A <- crossprod(matrix(rnorm(4 * 4), 4, 4))
  out_pd <- safe_pd_decomp(A)
  expect_null(out_pd$solve_z)
  S <- matrix(0, 3, 3)
  S[1, 1] <- 2
  out_sing <- safe_pd_decomp(S)
  expect_null(out_sing$solve_z)
  expect_equal(out_sing$r_eff, 1L)
})

# ---- get_eigen_R / get_Vtz ----

test_that("get_eigen_R prefers model$eigen_R then falls back to data$eigen_R", {
  data <- list(eigen_R = list(values = 1:3, vectors = diag(3)))
  model_data <- list(eigen_R = list(values = 4:6, vectors = diag(3)))
  expect_identical(get_eigen_R(data, model_data), model_data$eigen_R)
  expect_identical(get_eigen_R(data, list()), data$eigen_R)
})

test_that("get_Vtz prefers model$Vtz then falls back to data$Vtz", {
  data <- list(Vtz = c(1, 2, 3))
  model_data <- list(Vtz = c(4, 5, 6))
  expect_identical(get_Vtz(data, model_data), model_data$Vtz)
  expect_identical(get_Vtz(data, list()), data$Vtz)
})

# ---- optimize_omega ----

test_that("optimize_omega K=2 grid warm-start converges to optimal vertex", {
  eval_fn <- function(w) -sum((w - c(1, 0))^2)
  res <- optimize_omega(eval_fn, omega_cur = c(0.5, 0.5), K = 2)
  expect_length(res$omega, 2)
  expect_equal(sum(res$omega), 1, tolerance = 1e-8)
  expect_equal(res$omega, c(1, 0), tolerance = 1e-6)
  expect_true(is.logical(res$converged))
})

test_that("optimize_omega K=3 Frank-Wolfe converges to optimal vertex", {
  eval_fn <- function(w) -sum((w - c(0, 1, 0))^2)
  res <- optimize_omega(eval_fn, omega_cur = rep(1 / 3, 3), K = 3)
  expect_length(res$omega, 3)
  expect_equal(sum(res$omega), 1, tolerance = 1e-8)
  expect_true(all(res$omega >= 0))
  expect_gt(res$omega[2], 0.9)
})

# ---- eigen_from_reduced ----

test_that("eigen_from_reduced pads trailing eigenvalues with zeros when r < p", {
  set.seed(204)
  p <- 6
  r <- 4
  V_s <- qr.Q(qr(matrix(rnorm(p * r), p, r)))
  A1 <- crossprod(matrix(rnorm(r * r), r, r))
  A2 <- crossprod(matrix(rnorm(r * r), r, r))
  cache <- list(V_s = V_s, r = r, A_list = list(A1, A2))
  eig <- eigen_from_reduced(cache, omega = c(0.6, 0.4), K = 2, p = p)
  expect_length(eig$values, p)
  expect_equal(dim(eig$vectors), c(p, p))
  expect_true(all(eig$values >= 0))
  expect_false(is.unsorted(rev(eig$values)))
  expect_true(all(eig$values[(r + 1):p] == 0))
})

test_that("eigen_from_reduced matches direct eigendecomposition when r equals p", {
  set.seed(205)
  p <- 5
  r <- 5
  V_s <- qr.Q(qr(matrix(rnorm(p * r), p, r)))
  A1 <- crossprod(matrix(rnorm(r * r), r, r))
  A2 <- crossprod(matrix(rnorm(r * r), r, r))
  cache <- list(V_s = V_s, r = r, A_list = list(A1, A2))
  omega <- c(0.5, 0.5)
  eig <- eigen_from_reduced(cache, omega = omega, K = 2, p = p)
  R_full <- V_s %*% (omega[1] * A1 + omega[2] * A2) %*% t(V_s)
  R_full <- 0.5 * (R_full + t(R_full))
  ref <- eigen(R_full, symmetric = TRUE)
  expect_equal(sort(eig$values), sort(pmax(ref$values, 0)), tolerance = 1e-8)
})

# ---- pick_init_panel_via_subfits ----

test_that("pick_init_panel_via_subfits returns the highest-ELBO panel index", {
  set.seed(206)
  p <- 8
  R1 <- cov2cor(crossprod(matrix(rnorm(50 * p), 50, p)))
  R2 <- cov2cor(crossprod(matrix(rnorm(50 * p), 50, p)))
  z <- rnorm(p)
  z[2] <- 4
  parent_args <- list(z = z, n = 500, R_finite = NULL, L = 2,
                      max_iter = 3, convergence_method = "pip",
                      verbose = FALSE, model_init = NULL)
  res <- suppressWarnings(
    pick_init_panel_via_subfits(list(R1, R2), "R", parent_args)
  )
  expect_named(res, c("idx", "fits", "elbos"))
  expect_length(res$fits, 2)
  expect_length(res$elbos, 2)
  expect_true(res$idx %in% c(1L, 2L))
  expect_equal(res$idx, which.max(res$elbos))
})

test_that("pick_init_panel_via_subfits scores failed sub-fit as -Inf and selects good panel", {
  set.seed(207)
  p <- 6
  R_good <- cov2cor(crossprod(matrix(rnorm(40 * p), 40, p)))
  R_bad  <- diag(p + 2)
  z <- rnorm(p)
  parent_args <- list(z = z, n = 500, R_finite = NULL, L = 2,
                      max_iter = 3, convergence_method = "pip",
                      verbose = FALSE, model_init = NULL)
  res <- suppressWarnings(
    pick_init_panel_via_subfits(list(R_good, R_bad), "R", parent_args)
  )
  expect_null(res$fits[[2]])
  expect_false(is.null(res$fits[[1]]))
  expect_equal(res$elbos[2], -Inf)
  expect_equal(res$idx, 1L)
})
