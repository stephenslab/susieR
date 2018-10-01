context("test_null_weight.R")

test_that("null weight specification agrees with default",{
  simulate(200,1000)
  res1 = susie(cbind(X,0), y, estimate_prior_variance = TRUE)
  res2 = susie(X, y, estimate_prior_variance = TRUE, null_weight = 1/(ncol(X)+1))
  expect_equal_susie(res2,res1)
})
