context("test_sparse_multiplication.R")

test_that("sparse version sparse_multiplication", with(simulate(sparse = TRUE), {
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(1)
  L = 10
  M = matrix(rnorm(L*p), L, p)
  X.dense = set_X_attributes(X)
  X.standardized = t((t(X.dense) - attr(X.dense, 'scaled:center')) / attr(X.dense, 'scaled:scale'))
  X.sparse = set_X_attributes(X.sparse)

  expect_equal(compute_Xb(X.sparse, b), as.numeric(X.standardized%*%b))
  expect_equal(compute_Xb(X.dense, b), as.numeric(X.standardized%*%b))
  expect_equal(compute_Xty(X.sparse, y), as.numeric(t(X.standardized)%*%y))
  expect_equal(compute_Xty(X.dense, y), as.numeric(t(X.standardized)%*%y))
  expect_equal(compute_MXt(M, X.sparse),M%*%t(X.standardized),
               check.attributes =  FALSE,tol = 1e-12)
  expect_equal(compute_MXt(M, X.dense), M%*%t(X.standardized),
               check.attributes = FALSE,tol = 1e-12)
}))
