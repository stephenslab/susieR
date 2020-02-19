context("test_get_ER2_rss.R")

test_that("get_ER2_rss agrees with previous version", with(simulate(200, 500), {
  original.res = readRDS('ER2_rss_lambda0_res.rds')$s

  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat
  R = set_R_attributes(R, 1e-08)

  res = get_ER2_rss(R, z, s)/s$sigma2

  expect_equal(res, original.res, tol=1E-4)
}))
