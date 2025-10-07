context("test_get_objective_rss.R")

test_that("get_objective (lambda = 0) agrees with previous version", with(simulate(200, 500), {
  original.res = readRDS('objective_rss_lambda0_res.rds')$s

  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat
  R = set_R_attributes(R, 1e-08)
  attr(R, 'lambda') = 0

  res = get_objective_rss(R, z, s)

  expect_equal(res, original.res, tolerance=1e-4)
}))

test_that("get_objective (lambda = 1) agrees with previous version", with(simulate(200, 500), {
  original.res = readRDS('objective_rss_lambda1_res.rds')$s

  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat
  R = set_R_attributes(R, 1e-08)
  attr(R, 'lambda') = 1

  res = get_objective_rss(R, z, s)

  expect_equal(res, original.res, tolerance=1e-4)
}))
