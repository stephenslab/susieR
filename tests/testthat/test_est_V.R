context("test_est_V.R")

test_that("est_V has V nonzero for the debug case", {
  debug = readRDS('est_V_debug.rds')
  betahat = debug$b
  shat2 = debug$s2
  prior_weights = debug$pw
  V_debug = est_V_uniroot(betahat, shat2, prior_weights)

  betahat_V0 = rep(0, 512)
  shat2_V0  = rep(1, 512)
  V_0 = est_V_uniroot(betahat_V0, shat2_V0, prior_weights)
  expect_false(V_debug==0)
  expect_equal(V_0, 0)
})
