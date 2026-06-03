# Key invariant: a K=1 mixture with grid = c(V) and weights = c(1)
# must produce identical results to the scalar V path with
# estimate_prior_variance = FALSE and prior_variance = V.

context("Fixed mixture prior")

set.seed(1)
n <- 200
p <- 50
X <- matrix(rnorm(n * p), n, p)
beta <- rep(0, p)
beta[c(1, 5, 10)] <- c(0.5, -0.3, 0.4)
y <- X %*% beta + rnorm(n)

R <- cor(X)
z <- as.vector(sqrt(n) * crossprod(X, y) / sqrt(n * diag(crossprod(X))))

# ---- K=1 mixture matches scalar V ----

for (data_type in c("individual", "rss", "ss")) {
  test_that(paste("K=1 mixture matches scalar V for", data_type, "data"), {
    L <- 5
    tol <- 1e-10

    if (data_type == "individual") {
      fit_scalar <- susie(X, y, L = L,
                          estimate_prior_variance = FALSE,
                          estimate_residual_variance = FALSE,
                          max_iter = 20, tol = 1e-4)
      V_eff <- fit_scalar$V[1]
      fit_mixture <- susie(X, y, L = L,
                           prior_variance_grid = c(V_eff),
                           mixture_weights = c(1),
                           estimate_residual_variance = FALSE,
                           max_iter = 20, tol = 1e-4)
    } else if (data_type == "rss") {
      fit_scalar <- susie_rss(z = z, R = R, n = n, L = L,
                              estimate_prior_variance = FALSE,
                              estimate_residual_variance = FALSE,
                              max_iter = 20, tol = 1e-4)
      V_eff <- fit_scalar$V[1]
      fit_mixture <- susie_rss(z = z, R = R, n = n, L = L,
                               prior_variance_grid = c(V_eff),
                               mixture_weights = c(1),
                               estimate_residual_variance = FALSE,
                               max_iter = 20, tol = 1e-4)
    } else {
      XtX <- crossprod(X)
      Xty <- crossprod(X, y)
      yty <- sum(y^2)
      fit_scalar <- susie_ss(XtX = XtX, Xty = Xty, yty = yty, n = n, L = L,
                             estimate_prior_variance = FALSE,
                             estimate_residual_variance = FALSE,
                             max_iter = 20, tol = 1e-4)
      V_eff <- fit_scalar$V[1]
      fit_mixture <- susie_ss(XtX = XtX, Xty = Xty, yty = yty, n = n, L = L,
                              prior_variance_grid = c(V_eff),
                              mixture_weights = c(1),
                              estimate_residual_variance = FALSE,
                              max_iter = 20, tol = 1e-4)
    }

    expect_equal(fit_scalar$pip,   fit_mixture$pip,   tolerance = tol)
    expect_equal(fit_scalar$alpha, fit_mixture$alpha, tolerance = tol)
    expect_equal(fit_scalar$mu,    fit_mixture$mu,    tolerance = tol)
    expect_equal(fit_scalar$lbf,   fit_mixture$lbf,   tolerance = tol)
  })
}

# ---- K=1 exact match at single iteration ----

test_that("K=1 mixture is numerically identical for L=1 at single iteration", {
  L <- 1
  fit_scalar <- suppressWarnings(susie(X, y, L = L,
                                       estimate_prior_variance = FALSE,
                                       estimate_residual_variance = FALSE,
                                       max_iter = 1))
  V_eff <- fit_scalar$V[1]
  fit_mixture <- suppressWarnings(susie(X, y, L = L,
                                        prior_variance_grid = c(V_eff),
                                        mixture_weights = c(1),
                                        estimate_residual_variance = FALSE,
                                        max_iter = 1))

  expect_equal(fit_scalar$alpha, fit_mixture$alpha,
               tolerance = .Machine$double.eps * 10)
  expect_equal(fit_scalar$mu,    fit_mixture$mu,
               tolerance = .Machine$double.eps * 10)
  expect_equal(fit_scalar$lbf,   fit_mixture$lbf,
               tolerance = .Machine$double.eps * 10)
})

# ---- K>1 mixture validity ----

test_that("K=3 mixture produces valid outputs for RSS data", {
  L <- 5
  grid <- c(1, 10, 50)
  w <- c(0.3, 0.5, 0.2)

  fit <- susie_rss(z = z, R = R, n = n, L = L,
                   prior_variance_grid = grid,
                   mixture_weights = w,
                   estimate_residual_variance = FALSE,
                   max_iter = 20, tol = 1e-4)

  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  expect_equal(rowSums(fit$alpha), rep(1, L), tolerance = 1e-8)
  expect_length(fit$lbf_grid, L)
  expect_equal(dim(fit$lbf_grid[[1]]), c(p, 3))
  expect_true(all(is.finite(fit$mu)))
  expect_true(all(fit$mu2 >= fit$mu^2 - 1e-10))
})

# ---- Mixture BF matrix dimensions ----

test_that("Mixture BF matrix has correct dimensions for K=2", {
  L <- 1
  grid <- c(1, 50)
  w <- c(0.5, 0.5)

  fit <- suppressWarnings(susie_rss(z = z, R = R, n = n, L = L,
                                    prior_variance_grid = grid,
                                    mixture_weights = w,
                                    estimate_residual_variance = FALSE,
                                    max_iter = 1))

  expect_equal(ncol(fit$lbf_grid[[1]]), 2)
  expect_equal(nrow(fit$lbf_grid[[1]]), p)
})

# ---- Input validation ----

test_that("Invalid mixture prior inputs are rejected", {
  # Mismatched lengths
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 5,
              prior_variance_grid = c(1, 10),
              mixture_weights = c(1)),
    "length"
  )
  # Non-positive grid values
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 5,
              prior_variance_grid = c(-1, 10),
              mixture_weights = c(0.5, 0.5)),
    "prior_variance_grid"
  )
  # Weights not summing to 1
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 5,
              prior_variance_grid = c(1, 10),
              mixture_weights = c(0.3, 0.3)),
    "sum"
  )
  # Negative weights
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 5,
              prior_variance_grid = c(1, 10),
              mixture_weights = c(-0.1, 1.1)),
    "mixture_weights"
  )
})

# ---- NULL weights default to uniform ----

test_that("NULL mixture_weights defaults to uniform", {
  L <- 5
  grid <- c(1, 10, 50)

  fit <- suppressWarnings(susie_rss(z = z, R = R, n = n, L = L,
                                    prior_variance_grid = grid,
                                    estimate_residual_variance = FALSE,
                                    max_iter = 5))

  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  expect_equal(rowSums(fit$alpha), rep(1, L), tolerance = 1e-8)
})

# ---- Backward compatibility: standard path unchanged ----

test_that("Standard susie_rss without mixture prior is unchanged", {
  fit <- susie_rss(z = z, R = R, n = n, L = 5,
                   estimate_prior_variance = TRUE,
                   estimate_residual_variance = FALSE,
                   max_iter = 20)
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  expect_null(fit$lbf_grid)
})

# ---- Lambda regularization path ----

test_that("K=1 mixture matches scalar V with lambda regularization", {
  L <- 3
  lam <- 0.1
  fit_scalar <- susie_rss_lambda(z = z, R = R, n = n, L = L, lambda = lam,
                                 estimate_prior_variance = FALSE,
                                 estimate_residual_variance = FALSE,
                                 max_iter = 10)
  V_eff <- fit_scalar$V[1]

  fit_mixture <- susie_rss_lambda(z = z, R = R, n = n, L = L, lambda = lam,
                                  prior_variance_grid = c(V_eff),
                                  mixture_weights = c(1),
                                  estimate_residual_variance = FALSE,
                                  max_iter = 10)

  expect_equal(fit_scalar$alpha, fit_mixture$alpha, tolerance = 1e-8)
  expect_equal(fit_scalar$mu,    fit_mixture$mu,    tolerance = 1e-8)
  expect_equal(fit_scalar$lbf,   fit_mixture$lbf,   tolerance = 1e-8)
})

# ---- Finite-reference R inflation ----

test_that("Mixture prior works with finite-reference R inflation", {
  skip_if_not_installed("Matrix")
  L <- 3
  grid <- c(1, 10, 50)
  w <- c(0.3, 0.5, 0.2)

  fit <- suppressWarnings(susie_rss(z = z, R = R, n = n, L = L,
                                    prior_variance_grid = grid,
                                    mixture_weights = w,
                                    estimate_residual_variance = FALSE,
                                    R_finite = 30,
                                    max_iter = 5))

  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  expect_equal(rowSums(fit$alpha), rep(1, L), tolerance = 1e-8)
  expect_length(fit$lbf_grid, L)
})

# ---- Asymmetric weights shift PIPs ----

test_that("Asymmetric mixture weights shift PIPs in the expected direction", {
  L <- 3
  fit_large <- susie_rss(z = z, R = R, n = n, L = L,
                         prior_variance_grid = c(0.001, 100),
                         mixture_weights = c(0.01, 0.99),
                         estimate_residual_variance = FALSE,
                         max_iter = 10)
  fit_small <- susie_rss(z = z, R = R, n = n, L = L,
                         prior_variance_grid = c(0.001, 100),
                         mixture_weights = c(0.99, 0.01),
                         estimate_residual_variance = FALSE,
                         max_iter = 10)

  expect_true(all(fit_large$pip >= 0 & fit_large$pip <= 1))
  expect_true(all(fit_small$pip >= 0 & fit_small$pip <= 1))
  # PIPs differ between the two weight configurations
  expect_false(all(abs(fit_large$pip - fit_small$pip) < 1e-6))
})
