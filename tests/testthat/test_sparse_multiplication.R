context("sparse multiplication utilities")

# ---- helpers ----------------------------------------------------------------

make_dense <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  attr(X, "scaled:center") <- colMeans(X)
  attr(X, "scaled:scale") <- apply(X, 2, sd)
  X
}

make_sparse <- function(n, p, seed) {
  set.seed(seed)
  X <- Matrix::Matrix(rbinom(n * p, 1, 0.3) * rnorm(n * p), n, p, sparse = TRUE)
  attr(X, "scaled:center") <- Matrix::colMeans(X)
  attr(X, "scaled:scale") <- apply(as.matrix(X), 2, sd)
  X
}

# Build a fake trend-filtering X with attributes that compute_Xb/Xty/MXt expect.
make_tf <- function(n, order, seed) {
  set.seed(seed)
  X <- Matrix::sparseMatrix(i = NULL, j = NULL, dims = c(n, n))
  attr(X, "matrix.type") <- "tfmatrix"
  attr(X, "order")        <- order
  # Use compute_colstats logic directly (center=TRUE, scale=TRUE)
  cm  <- susieR:::compute_tf_cm(order, n)
  csd <- susieR:::compute_tf_csd(order, n)
  attr(X, "scaled:center") <- cm
  attr(X, "scaled:scale")  <- csd
  X
}

expected_Xb <- function(X_raw, b) {
  cm  <- attr(X_raw, "scaled:center")
  csd <- attr(X_raw, "scaled:scale")
  as.vector(scale(as.matrix(X_raw), center = cm, scale = csd) %*% b)
}

expected_Xty <- function(X_raw, y) {
  cm  <- attr(X_raw, "scaled:center")
  csd <- attr(X_raw, "scaled:scale")
  as.vector(t(scale(as.matrix(X_raw), center = cm, scale = csd)) %*% y)
}

expected_MXt <- function(M, X_raw) {
  cm  <- attr(X_raw, "scaled:center")
  csd <- attr(X_raw, "scaled:scale")
  M %*% t(scale(as.matrix(X_raw), center = cm, scale = csd))
}

# ---- compute_Xb -------------------------------------------------------------

test_that("compute_Xb: dense and sparse paths give values matching scaled multiply", {
  for (type in c("dense", "sparse")) {
    X <- if (type == "dense") make_dense(50, 10, 1) else make_sparse(100, 20, 2)
    set.seed(3)
    b <- rnorm(ncol(X))
    result   <- compute_Xb(X, b)
    expected <- expected_Xb(X, b)
    expect_equal(result, expected, tolerance = 1e-8,
                 label = paste("compute_Xb", type))
    expect_length(result, nrow(X))
  }
})

test_that("compute_Xb: center-only (scale=1) and scale-only (center=0) paths", {
  set.seed(4)
  n <- 30; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  b <- rnorm(p)

  for (variant in c("center_only", "scale_only")) {
    Xv <- X
    if (variant == "center_only") {
      attr(Xv, "scaled:center") <- colMeans(X)
      attr(Xv, "scaled:scale")  <- rep(1, p)
    } else {
      attr(Xv, "scaled:center") <- rep(0, p)
      attr(Xv, "scaled:scale")  <- apply(X, 2, sd)
    }
    expect_equal(compute_Xb(Xv, b), expected_Xb(Xv, b), tolerance = 1e-8,
                 label = paste("compute_Xb", variant))
  }
})

test_that("compute_Xb: zero b vector returns zero vector", {
  X <- make_dense(20, 8, 5)
  expect_equal(compute_Xb(X, rep(0, ncol(X))), rep(0, nrow(X)), tolerance = 1e-8)
})

test_that("compute_Xb: all-zero column is handled without error", {
  set.seed(6)
  n <- 30; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  X[, 2] <- 0
  attr(X, "scaled:center") <- colMeans(X)
  csd <- apply(X, 2, sd); csd[csd == 0] <- 1
  attr(X, "scaled:scale")  <- csd
  b <- rnorm(p)
  expect_equal(compute_Xb(X, b), expected_Xb(X, b), tolerance = 1e-8)
})

test_that("compute_Xb: single-column matrix", {
  set.seed(7)
  n <- 40; p <- 1
  X <- make_dense(n, p, 7)
  b <- rnorm(p)
  expect_equal(compute_Xb(X, b), expected_Xb(X, b), tolerance = 1e-8)
  expect_length(compute_Xb(X, b), n)
})

test_that("compute_Xb: trend-filtering dispatch calls compute_tf_Xb", {
  n <- 50
  X  <- make_tf(n, order = 0, seed = 8)
  cm  <- attr(X, "scaled:center")
  csd <- attr(X, "scaled:scale")
  set.seed(9)
  b <- rnorm(n)
  result   <- compute_Xb(X, b)
  # Reference: use compute_tf_Xb directly
  expected <- susieR:::compute_tf_Xb(0L, b / csd) - sum(cm * b / csd)
  expect_equal(result, expected, tolerance = 1e-8)
  expect_length(result, n)
})

# ---- compute_Xty ------------------------------------------------------------

test_that("compute_Xty: dense and sparse paths give values matching scaled multiply", {
  for (type in c("dense", "sparse")) {
    X <- if (type == "dense") make_dense(50, 10, 10) else make_sparse(100, 20, 11)
    set.seed(12)
    y <- rnorm(nrow(X))
    result   <- compute_Xty(X, y)
    expected <- expected_Xty(X, y)
    expect_equal(result, expected, tolerance = 1e-8,
                 label = paste("compute_Xty", type))
    expect_length(result, ncol(X))
  }
})

test_that("compute_Xty: center-only and scale-only paths", {
  set.seed(13)
  n <- 30; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  for (variant in c("center_only", "scale_only")) {
    Xv <- X
    if (variant == "center_only") {
      attr(Xv, "scaled:center") <- colMeans(X)
      attr(Xv, "scaled:scale")  <- rep(1, p)
    } else {
      attr(Xv, "scaled:center") <- rep(0, p)
      attr(Xv, "scaled:scale")  <- apply(X, 2, sd)
    }
    expect_equal(compute_Xty(Xv, y), expected_Xty(Xv, y), tolerance = 1e-8,
                 label = paste("compute_Xty", variant))
  }
})

test_that("compute_Xty: zero y vector returns zero vector", {
  X <- make_dense(20, 8, 14)
  expect_equal(compute_Xty(X, rep(0, nrow(X))), rep(0, ncol(X)), tolerance = 1e-8)
})

test_that("compute_Xty: trend-filtering dispatch calls compute_tf_Xty", {
  n <- 50
  X  <- make_tf(n, order = 0, seed = 15)
  cm  <- attr(X, "scaled:center")
  csd <- attr(X, "scaled:scale")
  set.seed(16)
  y <- rnorm(n)
  result   <- compute_Xty(X, y)
  expected <- susieR:::compute_tf_Xty(0L, y) / csd - cm / csd * sum(y)
  expect_equal(result, expected, tolerance = 1e-8)
})

# ---- compute_XtX ------------------------------------------------------------

test_that("compute_XtX: dense path gives symmetric PSD matrix equal to scaled X'X", {
  set.seed(17)
  n <- 50; p <- 10
  X <- make_dense(n, p, 17)
  cm  <- attr(X, "scaled:center")
  csd <- attr(X, "scaled:scale")
  X_std    <- scale(X, center = cm, scale = csd)
  expected <- t(X_std) %*% X_std
  result   <- compute_XtX(X)
  expect_equal(result, expected, tolerance = 1e-8)
  expect_equal(dim(result), c(p, p))
  # symmetry
  expect_equal(result, t(result), tolerance = 1e-8)
  # PSD: eigenvalues non-negative
  expect_true(all(eigen(result, symmetric = TRUE)$values >= -1e-8))
})

test_that("compute_XtX: sparse path gives values matching dense path", {
  X <- make_sparse(100, 20, 18)
  cm  <- attr(X, "scaled:center")
  csd <- attr(X, "scaled:scale")
  X_std    <- scale(as.matrix(X), center = cm, scale = csd)
  expected <- t(X_std) %*% X_std
  result   <- compute_XtX(X)
  expect_equal(as.matrix(result), expected, tolerance = 1e-8)
})

test_that("compute_XtX: center-only path", {
  set.seed(19)
  n <- 30; p <- 6
  X <- matrix(rnorm(n * p), n, p)
  attr(X, "scaled:center") <- colMeans(X)
  attr(X, "scaled:scale")  <- rep(1, p)
  X_c     <- scale(X, center = colMeans(X), scale = FALSE)
  expected <- t(X_c) %*% X_c
  expect_equal(compute_XtX(X), expected, tolerance = 1e-8)
})

test_that("compute_XtX: rejects trend-filtering matrices with informative error", {
  set.seed(20)
  X <- matrix(rnorm(50 * 10), 50, 10)
  attr(X, "scaled:center") <- colMeans(X)
  attr(X, "scaled:scale")  <- apply(X, 2, sd)
  attr(X, "matrix.type")   <- "trend_filtering"
  expect_error(compute_XtX(X),
               "compute_XtX not yet implemented for trend filtering matrices")
})

test_that("compute_XtX: bilinear identity b1'(X'X)b2 == (Xb1).(Xb2)", {
  set.seed(21)
  n <- 50; p <- 10
  X <- make_dense(n, p, 21)
  b1 <- rnorm(p); b2 <- rnorm(p)
  XtX  <- compute_XtX(X)
  Xb1  <- compute_Xb(X, b1)
  Xb2  <- compute_Xb(X, b2)
  expect_equal(sum(Xb1 * Xb2), drop(b1 %*% XtX %*% b2), tolerance = 1e-8)
})

# ---- compute_MXt ------------------------------------------------------------

test_that("compute_MXt: dense and sparse paths give values matching M * (scaled X)'", {
  for (type in c("dense", "sparse")) {
    X <- if (type == "dense") make_dense(50, 10, 22) else make_sparse(100, 20, 23)
    set.seed(24)
    L <- 4
    M <- matrix(rnorm(L * ncol(X)), L, ncol(X))
    result   <- compute_MXt(M, X)
    expected <- expected_MXt(M, X)
    expect_equal(result, expected, tolerance = 1e-8,
                 label = paste("compute_MXt", type))
    expect_equal(dim(result), c(L, nrow(X)))
  }
})

test_that("compute_MXt: single-row M", {
  X <- make_dense(40, 8, 25)
  set.seed(26)
  M <- matrix(rnorm(ncol(X)), 1, ncol(X))
  result <- compute_MXt(M, X)
  expect_equal(result, expected_MXt(M, X), tolerance = 1e-8)
  expect_equal(dim(result), c(1, nrow(X)))
})

test_that("compute_MXt: zero M matrix returns zero matrix", {
  X <- make_dense(30, 6, 27)
  L <- 2
  M <- matrix(0, L, ncol(X))
  expect_equal(compute_MXt(M, X), matrix(0, L, nrow(X)), tolerance = 1e-8)
})

test_that("compute_MXt: result equals row-wise compute_Xb", {
  X <- make_dense(40, 8, 28)
  set.seed(29)
  L <- 4
  M <- matrix(rnorm(L * ncol(X)), L, ncol(X))
  result_MXt <- compute_MXt(M, X)
  result_Xb  <- t(apply(M, 1, function(b) compute_Xb(X, b)))
  expect_equal(result_MXt, result_Xb, tolerance = 1e-8)
})

test_that("compute_MXt: trend-filtering dispatch routes through compute_Xb", {
  n <- 50
  X <- make_tf(n, order = 0, seed = 30)
  set.seed(31)
  L <- 3
  M <- matrix(rnorm(L * n), L, n)
  result_MXt <- compute_MXt(M, X)
  result_Xb  <- t(apply(M, 1, function(b) compute_Xb(X, b)))
  expect_equal(result_MXt, result_Xb, tolerance = 1e-8)
  expect_equal(dim(result_MXt), c(L, n))
})

# ---- cross-function dimension consistency ------------------------------------

test_that("all four functions return correct dimensions for p=1", {
  set.seed(32)
  n <- 50; p <- 1
  X <- make_dense(n, p, 32)
  b <- rnorm(p)
  y <- rnorm(n)
  M <- matrix(rnorm(2 * p), 2, p)
  expect_length(compute_Xb(X, b),    n)
  expect_length(compute_Xty(X, y),   p)
  expect_equal(dim(compute_XtX(X)),  c(p, p))
  expect_equal(dim(compute_MXt(M, X)), c(2, n))
})
