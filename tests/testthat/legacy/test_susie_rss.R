context("test_susie_rss.R")

test_that("susie_rss (lambda = 0) agrees with previous version",
          with(simulate(200, 500),{
  skip("susie_rss_lambda is not currently being used")
  original.res  = readRDS('susierssfit_lambda0_res.rds')$s
  
  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat

  res = susie_rss_lambda(z, R, lambda = 0, tol=1E-2, check_z = FALSE)

  skip('susie_rss_lambda is updated.')
  expect_equal(res$converged, TRUE)
  expect_equal_susie_rss(res, original.res, tolerance=1e-4)
}))

test_that("susie_rss (lambda = 1) agrees with previous version",
          with(simulate(200, 500),{
  skip("susie_rss_lambda is not currently being used")
  original.res  = readRDS('susierssfit_lambda1_res.rds')$s

  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat

  res = susie_rss_lambda(z, R, lambda = 1, tol=1E-2, check_z = FALSE)

  skip('susie_rss_lambda is updated.')
  expect_equal(res$converged, TRUE)
  expect_equal_susie_rss(res, original.res, tolerance=1e-4)
}))
