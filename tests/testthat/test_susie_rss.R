context("test_susie_rss.R")

test_that("susie_rss (lambda = 0) agrees with previous version", with(simulate(200, 500), {
  original.res  = readRDS('susierssfit_lambda0_res.rds')$s

  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat

  res = susie_rss_lambda(z, R, lambda = 0, tol=1E-2, check_z = FALSE)

  expect_equal(res$converged, TRUE)
  expect_equal_susie_rss(res, original.res, tol=1E-4)
}))

test_that("susie_rss (lambda = 1) agrees with previous version", with(simulate(200, 500), {
  original.res  = readRDS('susierssfit_lambda1_res.rds')$s

  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat

  res = susie_rss_lambda(z, R, lambda = 1, tol=1E-2, check_z = FALSE)

  expect_equal(res$converged, TRUE)
  expect_equal_susie_rss(res, original.res, tol=1E-4)
}))

test_that("susie_rss agrees with susie_rss_lambda when lambda = 0", with(simulate(200, 500), {

  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat

  res1 = susie_rss(z, R, tol=1E-2, check_z = FALSE)
  res2 = susie_rss_lambda(z, R, lambda = 0, tol=1E-2, check_z = FALSE)

  expect_equal(res1$converged, TRUE)
  expect_equal_susie_rss(res1, res2, tol=1e-05)
}))
