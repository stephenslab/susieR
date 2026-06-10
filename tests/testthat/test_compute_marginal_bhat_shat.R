context("compute_marginal_bhat_shat")

set.seed(42)
N <- 50
J <- 5

# ---- basic contracts ----

test_that("vector Y is treated as a single-column matrix", {
  X <- matrix(rnorm(N * J), N, J)
  X <- scale(X, center = TRUE, scale = FALSE)
  y <- rnorm(N)

  out_vec <- compute_marginal_bhat_shat(X, y)
  out_mat <- compute_marginal_bhat_shat(X, matrix(y, ncol = 1))

  expect_equal(dim(out_vec$Bhat), c(J, 1L))
  expect_equal(dim(out_vec$Shat), c(J, 1L))
  expect_equal(out_vec$Bhat, out_mat$Bhat, tolerance = 0)
  expect_equal(out_vec$Shat, out_mat$Shat, tolerance = 0)
})

test_that("matrix Y returns J x T Bhat and Shat", {
  X <- matrix(rnorm(N * J), N, J)
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- matrix(rnorm(N * 3), N, 3)

  out <- compute_marginal_bhat_shat(X, Y)

  expect_equal(dim(out$Bhat), c(J, 3L))
  expect_equal(dim(out$Shat), c(J, 3L))
})

test_that("predictor_weights override matches default colSums(X^2)", {
  X <- matrix(rnorm(N * J), N, J)
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- matrix(rnorm(N * 4), N, 4)

  out_default  <- compute_marginal_bhat_shat(X, Y)
  out_override <- compute_marginal_bhat_shat(X, Y, predictor_weights = colSums(X^2))

  expect_equal(out_default$Bhat, out_override$Bhat, tolerance = 0)
  expect_equal(out_default$Shat, out_override$Shat, tolerance = 0)
})

test_that("sigma2 supplied sets Shat to sqrt(sigma2 / predictor_weights)", {
  X <- matrix(rnorm(N * J), N, J)
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- matrix(rnorm(N * 2), N, 2)

  out      <- compute_marginal_bhat_shat(X, Y, sigma2 = 0.5)
  expected <- matrix(sqrt(0.5 / colSums(X^2)), nrow = J, ncol = 2)

  expect_equal(out$Shat, expected, tolerance = 0)
})

# ---- numerical correctness ----

test_that("Bhat equals X'Y / colSums(X^2) for centred X", {
  X <- matrix(rnorm(N * J), N, J)
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- matrix(rnorm(N * 3), N, 3)

  out      <- compute_marginal_bhat_shat(X, Y)
  expected <- crossprod(X, Y) / colSums(X^2)

  expect_equal(out$Bhat, expected, tolerance = 0)
})

test_that("Shat (no sigma2) matches per-column residual SD / sqrt(n-1)", {
  X <- matrix(rnorm(N * J), N, J)
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- matrix(rnorm(N * 2), N, 2)

  out  <- compute_marginal_bhat_shat(X, Y)
  Bhat <- crossprod(X, Y) / colSums(X^2)

  expected <- matrix(0, nrow = J, ncol = 2)
  for (t in 1:2) {
    for (j in 1:J) {
      r              <- Y[, t] - X[, j] * Bhat[j, t]
      expected[j, t] <- sqrt(var(r))
    }
  }
  expected <- expected / sqrt(N - 1)

  expect_equal(out$Shat, expected, tolerance = 1e-12)
})

test_that("Shat is always a matrix and all positive for T=1 without sigma2", {
  set.seed(99)
  X  <- matrix(rnorm(N * J), N, J)
  X  <- scale(X, center = TRUE, scale = FALSE)
  Y1 <- matrix(rnorm(N), N, 1)

  out <- compute_marginal_bhat_shat(X, Y1)

  expect_true(is.matrix(out$Shat))
  expect_equal(dim(out$Shat), c(J, 1L))
  expect_equal(out$Shat, abs(out$Shat))  # all positive
  expect_true(all(is.finite(out$Shat)))
})
