devtools::load_all(".")


context("sparse multiplication utilities")

# =============================================================================
# compute_Xb
# =============================================================================

test_that("compute_Xb works with dense matrices (centered and scaled)", {
  set.seed(123)
  n <- 50
  p <- 10

  # Create test data
  X_raw <- matrix(rnorm(n * p), n, p)
  b <- rnorm(p)

  # Standardize X and add attributes
  cm <- colMeans(X_raw)
  csd <- apply(X_raw, 2, sd)
  X_std <- scale(X_raw, center = TRUE, scale = TRUE)

  # Add attributes to raw X
  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  # Compute using function
  result <- compute_Xb(X_raw, b)

  # Compute expected result (naive)
  expected <- as.vector(X_std %*% b)

  expect_equal(result, expected, tolerance = 1e-10)
  expect_length(result, n)
})

test_that("compute_Xb works with sparse matrices (centered and scaled)", {
  set.seed(456)
  n <- 100
  p <- 20

  # Create sparse test data (30% non-zero)
  X_raw <- Matrix::Matrix(rbinom(n * p, 1, 0.3) * rnorm(n * p), n, p, sparse = TRUE)
  b <- rnorm(p)

  # Standardize and add attributes
  cm <- Matrix::colMeans(X_raw)
  csd <- apply(as.matrix(X_raw), 2, sd)

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  # Compute using function
  result <- compute_Xb(X_raw, b)

  # Compute expected result (naive with standardization)
  X_std <- scale(as.matrix(X_raw), center = cm, scale = csd)
  expected <- as.vector(X_std %*% b)

  expect_equal(result, expected, tolerance = 1e-10)
  expect_length(result, n)
})

test_that("compute_Xb works with only centering (no scaling)", {
  set.seed(789)
  n <- 30
  p <- 5

  X_raw <- matrix(rnorm(n * p), n, p)
  b <- rnorm(p)

  cm <- colMeans(X_raw)
  csd <- rep(1, p)  # No scaling

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  result <- compute_Xb(X_raw, b)

  # Expected: centered X times b
  X_centered <- scale(X_raw, center = cm, scale = FALSE)
  expected <- as.vector(X_centered %*% b)

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("compute_Xb works with only scaling (no centering)", {
  set.seed(101)
  n <- 30
  p <- 5

  X_raw <- matrix(rnorm(n * p), n, p)
  b <- rnorm(p)

  cm <- rep(0, p)  # No centering
  csd <- apply(X_raw, 2, sd)

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  result <- compute_Xb(X_raw, b)

  # Expected: scaled X times b
  X_scaled <- scale(X_raw, center = FALSE, scale = csd)
  expected <- as.vector(X_scaled %*% b)

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("compute_Xb handles zero vector b", {
  set.seed(202)
  n <- 20
  p <- 8

  X_raw <- matrix(rnorm(n * p), n, p)
  b <- rep(0, p)

  attr(X_raw, "scaled:center") <- colMeans(X_raw)
  attr(X_raw, "scaled:scale") <- apply(X_raw, 2, sd)

  result <- compute_Xb(X_raw, b)

  expect_equal(result, rep(0, n), tolerance = 1e-10)
})

# =============================================================================
# compute_Xty
# =============================================================================

test_that("compute_Xty works with dense matrices (centered and scaled)", {
  set.seed(303)
  n <- 50
  p <- 10

  X_raw <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  # Standardize X
  cm <- colMeans(X_raw)
  csd <- apply(X_raw, 2, sd)
  X_std <- scale(X_raw, center = cm, scale = csd)

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  result <- compute_Xty(X_raw, y)

  # Expected result
  expected <- as.vector(t(X_std) %*% y)

  expect_equal(result, expected, tolerance = 1e-10)
  expect_length(result, p)
})

test_that("compute_Xty works with sparse matrices (centered and scaled)", {
  set.seed(404)
  n <- 100
  p <- 20

  X_raw <- Matrix::Matrix(rbinom(n * p, 1, 0.3) * rnorm(n * p), n, p, sparse = TRUE)
  y <- rnorm(n)

  cm <- Matrix::colMeans(X_raw)
  csd <- apply(as.matrix(X_raw), 2, sd)

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  result <- compute_Xty(X_raw, y)

  # Expected result
  X_std <- scale(as.matrix(X_raw), center = cm, scale = csd)
  expected <- as.vector(t(X_std) %*% y)

  expect_equal(result, expected, tolerance = 1e-10)
  expect_length(result, p)
})

test_that("compute_Xty works with only centering (no scaling)", {
  set.seed(505)
  n <- 30
  p <- 5

  X_raw <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  cm <- colMeans(X_raw)
  csd <- rep(1, p)

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  result <- compute_Xty(X_raw, y)

  X_centered <- scale(X_raw, center = cm, scale = FALSE)
  expected <- as.vector(t(X_centered) %*% y)

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("compute_Xty handles zero vector y", {
  set.seed(606)
  n <- 20
  p <- 8

  X_raw <- matrix(rnorm(n * p), n, p)
  y <- rep(0, n)

  attr(X_raw, "scaled:center") <- colMeans(X_raw)
  attr(X_raw, "scaled:scale") <- apply(X_raw, 2, sd)

  result <- compute_Xty(X_raw, y)

  expect_equal(result, rep(0, p), tolerance = 1e-10)
})

# =============================================================================
# compute_XtX
# =============================================================================

test_that("compute_XtX works with dense matrices (centered and scaled)", {
  set.seed(707)
  n <- 50
  p <- 10

  X_raw <- matrix(rnorm(n * p), n, p)

  cm <- colMeans(X_raw)
  csd <- apply(X_raw, 2, sd)
  X_std <- scale(X_raw, center = cm, scale = csd)

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd
  attr(X_raw, "d") <- colSums(X_std^2)

  result <- compute_XtX(X_raw)

  # Expected result
  expected <- t(X_std) %*% X_std

  expect_equal(result, expected, tolerance = 1e-10)
  expect_equal(dim(result), c(p, p))
})

test_that("compute_XtX produces symmetric matrix", {
  set.seed(808)
  n <- 40
  p <- 8

  X_raw <- matrix(rnorm(n * p), n, p)

  attr(X_raw, "scaled:center") <- colMeans(X_raw)
  attr(X_raw, "scaled:scale") <- apply(X_raw, 2, sd)

  result <- compute_XtX(X_raw)

  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("compute_XtX is positive semi-definite", {
  set.seed(909)
  n <- 60
  p <- 12

  X_raw <- matrix(rnorm(n * p), n, p)

  attr(X_raw, "scaled:center") <- colMeans(X_raw)
  attr(X_raw, "scaled:scale") <- apply(X_raw, 2, sd)

  result <- compute_XtX(X_raw)

  # Check eigenvalues are non-negative
  eigenvalues <- eigen(result, symmetric = TRUE)$values
  expect_true(all(eigenvalues >= -1e-10))  # Allow for numerical error
})

test_that("compute_XtX works with sparse matrices", {
  set.seed(1010)
  n <- 100
  p <- 20

  X_raw <- Matrix::Matrix(rbinom(n * p, 1, 0.3) * rnorm(n * p), n, p, sparse = TRUE)

  cm <- Matrix::colMeans(X_raw)
  csd <- apply(as.matrix(X_raw), 2, sd)

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  result <- compute_XtX(X_raw)

  # Expected result
  X_std <- scale(as.matrix(X_raw), center = cm, scale = csd)
  expected <- t(X_std) %*% X_std

  expect_equal(as.matrix(result), expected, tolerance = 1e-9)
})

test_that("compute_XtX works with only centering (no scaling)", {
  set.seed(1111)
  n <- 30
  p <- 6

  X_raw <- matrix(rnorm(n * p), n, p)

  cm <- colMeans(X_raw)
  csd <- rep(1, p)

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  result <- compute_XtX(X_raw)

  X_centered <- scale(X_raw, center = cm, scale = FALSE)
  expected <- t(X_centered) %*% X_centered

  expect_equal(result, expected, tolerance = 1e-10)
})

# =============================================================================
# compute_MXt
# =============================================================================

test_that("compute_MXt works with dense matrices (centered and scaled)", {
  set.seed(1212)
  n <- 50
  p <- 10
  L <- 3

  X_raw <- matrix(rnorm(n * p), n, p)
  M <- matrix(rnorm(L * p), L, p)

  cm <- colMeans(X_raw)
  csd <- apply(X_raw, 2, sd)
  X_std <- scale(X_raw, center = cm, scale = csd)

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  result <- compute_MXt(M, X_raw)

  # Expected result
  expected <- M %*% t(X_std)

  expect_equal(result, expected, tolerance = 1e-10)
  expect_equal(dim(result), c(L, n))
})

test_that("compute_MXt works with sparse matrices", {
  set.seed(1313)
  n <- 100
  p <- 20
  L <- 5

  X_raw <- Matrix::Matrix(rbinom(n * p, 1, 0.3) * rnorm(n * p), n, p, sparse = TRUE)
  M <- matrix(rnorm(L * p), L, p)

  cm <- Matrix::colMeans(X_raw)
  csd <- apply(as.matrix(X_raw), 2, sd)

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  result <- compute_MXt(M, X_raw)

  # Expected result
  X_std <- scale(as.matrix(X_raw), center = cm, scale = csd)
  expected <- M %*% t(X_std)

  expect_equal(result, expected, tolerance = 1e-9)
})

test_that("compute_MXt works with single row M", {
  set.seed(1414)
  n <- 40
  p <- 8

  X_raw <- matrix(rnorm(n * p), n, p)
  M <- matrix(rnorm(p), 1, p)  # Single row

  cm <- colMeans(X_raw)
  csd <- apply(X_raw, 2, sd)
  X_std <- scale(X_raw, center = cm, scale = csd)

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  result <- compute_MXt(M, X_raw)

  expected <- M %*% t(X_std)

  expect_equal(result, expected, tolerance = 1e-10)
  expect_equal(dim(result), c(1, n))
})

test_that("compute_MXt handles zero matrix M", {
  set.seed(1515)
  n <- 30
  p <- 6
  L <- 2

  X_raw <- matrix(rnorm(n * p), n, p)
  M <- matrix(0, L, p)

  attr(X_raw, "scaled:center") <- colMeans(X_raw)
  attr(X_raw, "scaled:scale") <- apply(X_raw, 2, sd)

  result <- compute_MXt(M, X_raw)

  expect_equal(result, matrix(0, L, n), tolerance = 1e-10)
})

test_that("compute_MXt is equivalent to row-wise compute_Xb", {
  set.seed(1616)
  n <- 40
  p <- 8
  L <- 4

  X_raw <- matrix(rnorm(n * p), n, p)
  M <- matrix(rnorm(L * p), L, p)

  cm <- colMeans(X_raw)
  csd <- apply(X_raw, 2, sd)

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  # Using compute_MXt
  result_MXt <- compute_MXt(M, X_raw)

  # Using row-wise compute_Xb
  result_Xb <- t(apply(M, 1, function(b) compute_Xb(X_raw, b)))

  expect_equal(result_MXt, result_Xb, tolerance = 1e-10)
})

# =============================================================================
# Edge Cases and Consistency
# =============================================================================

test_that("sparse multiplication functions preserve dimensions correctly", {
  set.seed(1717)
  n <- 25
  p <- 7
  L <- 3

  X_raw <- matrix(rnorm(n * p), n, p)
  b <- rnorm(p)
  y <- rnorm(n)
  M <- matrix(rnorm(L * p), L, p)

  attr(X_raw, "scaled:center") <- colMeans(X_raw)
  attr(X_raw, "scaled:scale") <- apply(X_raw, 2, sd)

  # Test dimensions
  expect_length(compute_Xb(X_raw, b), n)
  expect_length(compute_Xty(X_raw, y), p)
  expect_equal(dim(compute_XtX(X_raw)), c(p, p))
  expect_equal(dim(compute_MXt(M, X_raw)), c(L, n))
})

test_that("all sparse functions handle edge case with p=1", {
  set.seed(1818)
  n <- 50
  p <- 1

  X_raw <- matrix(rnorm(n * p), n, p)
  b <- rnorm(p)
  y <- rnorm(n)
  M <- matrix(rnorm(2 * p), 2, p)

  attr(X_raw, "scaled:center") <- colMeans(X_raw)
  attr(X_raw, "scaled:scale") <- apply(X_raw, 2, sd)

  # Should not error
  expect_length(compute_Xb(X_raw, b), n)
  expect_length(compute_Xty(X_raw, y), p)
  expect_equal(dim(compute_XtX(X_raw)), c(p, p))
  expect_equal(dim(compute_MXt(M, X_raw)), c(2, n))
})

test_that("consistency test: compute_Xty(X,y) should equal t(X) %*% y for standardized X", {
  set.seed(1919)
  n <- 45
  p <- 9

  X_raw <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  cm <- colMeans(X_raw)
  csd <- apply(X_raw, 2, sd)

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  # Using compute_Xty
  result1 <- compute_Xty(X_raw, y)

  # Manually standardize and compute
  X_std <- scale(X_raw, center = cm, scale = csd)
  result2 <- as.vector(t(X_std) %*% y)

  expect_equal(result1, result2, tolerance = 1e-10)
})

test_that("consistency test: compute_Xb and compute_XtX relationship", {
  set.seed(2020)
  n <- 50
  p <- 10

  X_raw <- matrix(rnorm(n * p), n, p)
  b1 <- rnorm(p)
  b2 <- rnorm(p)

  cm <- colMeans(X_raw)
  csd <- apply(X_raw, 2, sd)

  attr(X_raw, "scaled:center") <- cm
  attr(X_raw, "scaled:scale") <- csd

  # Compute Xb1 and Xb2
  Xb1 <- compute_Xb(X_raw, b1)
  Xb2 <- compute_Xb(X_raw, b2)

  # Inner product should equal b1' XtX b2
  XtX <- compute_XtX(X_raw)

  result1 <- sum(Xb1 * Xb2)
  result2 <- sum(b1 * (XtX %*% b2))

  expect_equal(result1, result2, tolerance = 1e-9)
})
