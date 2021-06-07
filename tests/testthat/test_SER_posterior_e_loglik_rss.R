context("test_SER_posterior_e_loglik_rss.R")

test_that("SER_posterior_e_loglik_rss (lambda = 0) agrees with previous version", with(simulate(200, 500), {
  original.res = readRDS('SER_rss_lambda0_res.rds')$s
  Ez = rep(1, p)
  Ez2 = rep(1, p)
  s2 = s$sigma2

  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat
  R = set_R_attributes(R, 1e-08)
  attr(R, 'lambda') = 0
  Sigma = update_Sigma(R, s2, z)

  res = SER_posterior_e_loglik_rss(R, Sigma, z,Ez,Ez2)

  expect_equal(res, original.res, tolerance=1e-4)
}))

test_that("SER_posterior_e_loglik_rss (lambda = 1) agrees with previous version", with(simulate(200, 500), {
  original.res = readRDS('SER_rss_lambda1_res.rds')$s
  Ez = rep(1, p)
  Ez2 = rep(1, p)
  s2 = s$sigma2

  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat
  R = set_R_attributes(R, 1e-08)
  attr(R, 'lambda') = 1
  Sigma = update_Sigma(R, s2, z)

  res = SER_posterior_e_loglik_rss(R, Sigma, z,Ez,Ez2)

  expect_equal(res, original.res, tolerance=1e-4)
}))
