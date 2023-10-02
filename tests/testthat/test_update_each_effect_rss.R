context("test_update_each_effect_rss.R")

# test_that("update_each_effect_rss (lambda = 0) agrees with previous version ", with(simulate(200, 500), {
#   original.res = readRDS("vbupdate_rss_lambda0_res.rds")$s
#
#   ss = univariate_regression(X, y)
#   R = cor(X)
#   z = ss$betahat/ss$sebetahat
#   s$Rz = rep(5, nrow(R))
#
#   R = set_R_attributes(R, 1e-08)
#   attr(R, "lambda") = 0
#   Sigma = update_Sigma(R, s$sigma2, z)
#
#   res = update_each_effect_rss(R, z, s, Sigma)
#   expect_equal_susie_rss_update(res, original.res, tolerance = 1e-4)
# }))

# test_that("update_each_effect_rss (lambda = 1) agrees with previous version ", with(simulate(200, 500), {
#   original.res = readRDS("vbupdate_rss_lambda1_res.rds")$s
#
#   ss = univariate_regression(X, y)
#   R = cor(X)
#   z = ss$betahat/ss$sebetahat
#   s$Rz = rep(5, nrow(R))
#
#   R = set_R_attributes(R, 1e-08)
#   attr(R, "lambda") = 1
#   Sigma = update_Sigma(R, s$sigma2, z)
#
#   res = update_each_effect_rss(R, z, s, Sigma)
#   expect_equal_susie_rss_update(res, original.res, tolerance = 1e-4)
# }))
