context("test_Eloglik_rss.R")

test_that("Eloglik_rss agrees with previous version", with(simulate(200, 500), {
  original.res = readRDS('Eloglik_rss_lambda0_res.rds')$s

  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat
  R = set_R_attributes(R, 1e-08)

  res  = Eloglik_rss(R, z, s)

  expect_equal(res, original.res, tol=1E-4)
}))

