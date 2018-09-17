test_that("sparse version sparse_multiplication",{
  simulate(sparse=T)
  L = 10
  set.seed(1)
  M = matrix(rnorm(L*p), L, p)
  
  scaled.X = susieR:::safe_colScale(X)
  expect_equal(susieR:::compute_Xb(X.sparse, b), scaled.X%*%b)
  expect_equal(susieR:::compute_Xb(scaled.X, b), scaled.X%*%b)
  expect_equal(susieR:::compute_Xty(X.sparse, y), t(scaled.X)%*%y)
  expect_equal(susieR:::compute_Xty(scaled.X, y), t(scaled.X)%*%y)
  expect_equal(susieR:::compute_MXt(M, X.sparse), M%*%t(scaled.X))
  expect_equal(susieR:::compute_MXt(M, scaled.X), M%*%t(scaled.X))
})
