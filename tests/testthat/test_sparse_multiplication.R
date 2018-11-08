context("test_sparse_multiplication.R")

test_that("sparse version sparse_multiplication", with(simulate(sparse = TRUE), {
  set.seed(1)
  L = 10
  M = matrix(rnorm(L*p), L, p)
  scaled.X = safe_colScale(X)
  X.sparse = safe_colScale(X.sparse)
  expect_equal(compute_Xb(X.sparse, b), as.numeric(scaled.X%*%b))
  expect_equal(compute_Xb(scaled.X, b), scaled.X%*%b)
  expect_equal(compute_Xty(X.sparse, y), as.numeric(t(scaled.X)%*%y))
  expect_equal(compute_Xty(scaled.X, y), t(scaled.X)%*%y)
  expect_equal(compute_MXt(M, X.sparse),M%*%t(scaled.X),
               check.attributes =  FALSE,tol = 1e-12)
  expect_equal(compute_MXt(M, scaled.X), M%*%t(scaled.X),
               check.attributes = FALSE,tol = 1e-12)
}))
