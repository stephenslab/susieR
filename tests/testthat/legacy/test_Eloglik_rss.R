context("test_Eloglik_rss.R")

test_that("Eloglik_rss (lambda = 0) agrees with previous version", with(simulate(200, 500), {
  original.res = readRDS('Eloglik_rss_lambda0_res.rds')$s

  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat
  R = set_R_attributes(R, 1e-08)
  attr(R, 'lambda') = 0

  res  = Eloglik_rss(s$sigma2, R, z, s)

  expect_equal(res, original.res, tolerance=1e-4)
}))

test_that("Eloglik_rss (lambda = 1) agrees with previous version", with(simulate(200, 500), {
  original.res = readRDS('Eloglik_rss_lambda1_res.rds')$s

  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat
  R = set_R_attributes(R, 1e-08)
  attr(R, 'lambda') = 1

  res  = Eloglik_rss(s$sigma2, R, z, s)

  expect_equal(res, original.res, tolerance=1e-4)
}))
