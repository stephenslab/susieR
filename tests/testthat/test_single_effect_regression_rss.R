context("test_single_effect_regression_rss.R")

# test_that("single_effect_regression_rss (lambda = 0) agrees with previous version", with(simulate(200, 500), {
#   original.res = readRDS('singleReg_rss_lambda0_res.rds')$s
#   # the git commit id to previous result is readRDS('singleReg_rss_res.rds')$git_commit_id
#
#   ss = univariate_regression(X, y)
#   R = cor(X)
#   z = ss$betahat/ss$sebetahat
#   R = set_R_attributes(R, 1e-08)
#   attr(R, 'lambda') = 0
#   Sigma = update_Sigma(R, s$sigma2, z)
#
#   s = single_effect_regression_rss(z, Sigma, optimize_V = "none")
#
#   expect_equal(attr(R, 'd'), rep(1, 500), tolerance=1e-4)
#   expect_equal_SER_suff_stat(s, original.res, tolerance=1e-4)
# }))

# test_that("single_effect_regression_rss (lambda = 1) agrees with previous version", with(simulate(200, 500), {
#   original.res = readRDS('singleReg_rss_lambda1_res.rds')$s
#   # the git commit id to previous result is readRDS('singleReg_rss_res.rds')$git_commit_id
#
#   ss = univariate_regression(X, y)
#   R = cor(X)
#   z = ss$betahat/ss$sebetahat
#   R = set_R_attributes(R, 1e-08)
#   attr(R, 'lambda') = 1
#   Sigma = update_Sigma(R, s$sigma2, z)
#
#   s = single_effect_regression_rss(z, Sigma, optimize_V = "none")
#
#   expect_equal(attr(R, 'd'), rep(1, 500), tolerance=1e-4)
#   expect_equal_SER_suff_stat(s, original.res, tolerance=1e-4)
# }))
