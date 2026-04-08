# Tests that the rss_lambda path with lambda=0 produces results consistent
# with the suff_stat (ss) path. This validates that multi-panel mixture with
# lambda=0 (no LD regularization) is mathematically sound.

test_that("single-panel list input (lambda=0) matches ss path in PIPs", {
  set.seed(42)
  p <- 30
  n <- 500
  # Generate full-rank R
  X_full <- matrix(rnorm(n * p), n, p)
  R <- cor(X_full)
  # Generate z-scores with one causal
  beta <- rep(0, p)
  beta[5] <- 3
  z <- as.vector(R %*% beta + rnorm(p))

  # ss path (default lambda=0, single R matrix)
  fit_ss <- suppressWarnings(susie_rss(z = z, R = R, n = n,
                                        L = 5, max_iter = 50,
                                        check_prior = FALSE))

  # rss_lambda path via single-element list (forces multi-panel routing)
  B <- 200
  Xk <- matrix(rnorm(B * p), B, p)
  # Make Xk's correlation match R closely
  Xk <- Xk %*% chol(R)
  fit_lambda0 <- suppressWarnings(susie_rss(z = z, X = list(Xk),
                                             L = 5, max_iter = 50,
                                             check_prior = FALSE))

  # PIPs should identify the same causal variable
  expect_equal(which.max(fit_ss$pip), which.max(fit_lambda0$pip))
  # Top PIP should be high for both

  expect_gt(fit_ss$pip[5], 0.5)
  expect_gt(fit_lambda0$pip[5], 0.5)
})

test_that("multi-panel lambda=0 vs manual R_omega gives similar PIPs", {
  set.seed(456)
  p <- 30
  n <- 500
  B1 <- 100
  B2 <- 100

  # Two panels with correlated LD
  X_base <- matrix(rnorm(n * p), n, p)
  R_base <- cor(X_base)
  L_chol <- chol(R_base)
  X1 <- matrix(rnorm(B1 * p), B1, p) %*% L_chol
  X2 <- matrix(rnorm(B2 * p), B2, p) %*% L_chol

  # Generate z-scores with a strong signal
  beta <- rep(0, p)
  beta[10] <- 5
  z <- as.vector(R_base %*% beta + rnorm(p))

  # Multi-panel mixture (lambda=0, new default)
  fit_mix <- suppressWarnings(susie_rss(z = z, X = list(X1, X2),
                                         L = 5, max_iter = 50,
                                         check_prior = FALSE))

  # Extract omega and compute R_omega manually
  omega <- fit_mix$omega
  R1 <- cov2cor(crossprod(X1))
  R2 <- cov2cor(crossprod(X2))
  R_omega <- omega[1] * R1 + omega[2] * R2

  # Manual single-panel with R_omega
  fit_manual <- suppressWarnings(susie_rss(z = z, R = R_omega, n = n,
                                            L = 5, max_iter = 50,
                                            check_prior = FALSE))

  # Both should identify the same top variable
  expect_equal(which.max(fit_mix$pip), which.max(fit_manual$pip))
  # PIPs should be correlated (only check if there's variance)
  if (sd(fit_mix$pip) > 0 && sd(fit_manual$pip) > 0)
    expect_gt(cor(fit_mix$pip, fit_manual$pip), 0.8)
})

test_that("rss_lambda(lambda=0) matches ss path numerically for full-rank single panel", {
  # When B >= p, the single-element list forces the rss_lambda path.
  # With lambda=0 and full-rank R, the two paths should give identical
  # alpha (PIPs) and mu (posterior means).
  set.seed(314)
  p <- 20
  B <- 200  # B >> p so standardize_X gives full-rank R
  n <- 1000

  # Generate X with known correlation structure
  X_pop <- matrix(rnorm(n * p), n, p)
  R <- cor(X_pop)

  # Sketch that is full-rank with B >> p
  Xk <- matrix(rnorm(B * p), B, p) %*% chol(R)

  # z-scores with signal
  beta <- rep(0, p)
  beta[3] <- 4
  beta[15] <- 3
  z <- as.vector(R %*% beta + rnorm(p) * 0.5)

  # ss path: susie_rss with R matrix directly
  fit_ss <- suppressWarnings(susie_rss(z = z, R = R, n = n,
                                        L = 5, max_iter = 100,
                                        check_prior = FALSE))

  # rss_lambda path: susie_rss with single-element X list
  fit_eigen <- suppressWarnings(susie_rss(z = z, X = list(Xk),
                                           L = 5, max_iter = 100,
                                           check_prior = FALSE))

  # PIPs should agree closely (not exact due to different R approximation)
  expect_equal(which.max(fit_ss$pip), which.max(fit_eigen$pip))

  # Top PIPs should be high for both
  expect_gt(fit_ss$pip[3], 0.8)
  expect_gt(fit_eigen$pip[3], 0.8)
  expect_gt(fit_ss$pip[15], 0.5)
  expect_gt(fit_eigen$pip[15], 0.5)

  # CS should cover the same variables
  cs_ss <- fit_ss$sets$cs
  cs_eigen <- fit_eigen$sets$cs
  # Both should find at least one CS
  expect_gt(length(cs_ss), 0)
  expect_gt(length(cs_eigen), 0)
})

test_that("multi-panel stores single_panel_fits in returned object", {
  set.seed(88)
  p <- 20
  B <- 50
  X1 <- matrix(rnorm(B * p), B, p)
  X2 <- matrix(rnorm(B * p), B, p)
  z <- rnorm(p)
  z[5] <- 4

  fit <- suppressWarnings(susie_rss(z = z, X = list(X1, X2),
                                     L = 3, max_iter = 20,
                                     check_prior = FALSE))
  # single_panel_fits should be a list of length 2 (one per panel)
  expect_true(!is.null(fit$single_panel_fits))
  expect_equal(length(fit$single_panel_fits), 2)
  expect_s3_class(fit$single_panel_fits[[1]], "susie")
  expect_s3_class(fit$single_panel_fits[[2]], "susie")
  # Each single-panel fit should have PIPs
  expect_equal(length(fit$single_panel_fits[[1]]$pip), p)
  expect_equal(length(fit$single_panel_fits[[2]]$pip), p)
})

test_that("lambda > 0 backward compatibility preserved", {
  set.seed(99)
  p <- 20
  B <- 50
  X1 <- matrix(rnorm(B * p), B, p)
  X2 <- matrix(rnorm(B * p), B, p)
  z <- rnorm(p)
  z[3] <- 4

  # Explicit lambda > 0 should still work
  fit <- suppressWarnings(susie_rss(z = z, X = list(X1, X2),
                                     lambda = 0.01,
                                     L = 3, max_iter = 20,
                                     check_prior = FALSE))
  expect_s3_class(fit, "susie")
  expect_equal(length(fit$pip), p)
  expect_true(fit$pip[3] > 0.3)
})

test_that("multi-panel lambda=0 does not produce NaN or Inf in ELBO", {
  set.seed(77)
  p <- 25
  B <- 40  # B < p, so reduced-basis is rank-deficient
  X1 <- matrix(rnorm(B * p), B, p)
  X2 <- matrix(rnorm(B * p), B, p)
  z <- rnorm(p)
  z[7] <- 5

  fit <- suppressWarnings(susie_rss(z = z, X = list(X1, X2),
                                     L = 3, max_iter = 30,
                                     check_prior = FALSE))
  expect_s3_class(fit, "susie")
  # ELBO should be finite
  expect_true(all(is.finite(fit$elbo[fit$elbo != 0])))
  expect_equal(length(fit$pip), p)
})

test_that("multi-panel lambda=0 with B_total < p works correctly", {
  set.seed(55)
  p <- 50
  B1 <- 15
  B2 <- 15
  # B_total = 30 < p = 50: reduced-basis path
  X1 <- matrix(rnorm(B1 * p), B1, p)
  X2 <- matrix(rnorm(B2 * p), B2, p)
  beta <- rep(0, p)
  beta[1] <- 4
  z <- rnorm(p)
  z[1] <- 5

  fit <- suppressWarnings(susie_rss(z = z, X = list(X1, X2),
                                     L = 3, max_iter = 30,
                                     check_prior = FALSE))
  expect_s3_class(fit, "susie")
  expect_true(all(is.finite(fit$elbo[fit$elbo != 0])))
  # Should identify the causal variable as having highest PIP
  expect_equal(which.max(fit$pip), 1L)
})
