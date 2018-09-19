test_that("prior weights specification agrees with default",{
  simulate(200,1000)
  res1 = susie(X, y, estimate_prior_variance = TRUE)
  res2 = susie(X, y, estimate_prior_variance = TRUE,
               prior_weights = rep(1/ncol(X), ncol(X)))
  is_equal_susie(res1,res2)
})
