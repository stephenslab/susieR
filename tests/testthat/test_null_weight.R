context("test_null_weight.R")

test_that("null weight specification agrees with default", with(simulate(200,1000), {
  res1 = susie(cbind(X,0), y, estimate_prior_variance = TRUE)
  res2 = susie(X, y, estimate_prior_variance = TRUE, null_weight = 1/(ncol(X)+1))
  expect_equal_susie(res2,res1)
}))

test_that("SS: null weight specification agrees with default", with(simulate(200,1000), {
  ss = compute_ss(X,y)

  res1 = susie_ss(cbind(rbind(ss$XtX,0),0), c(ss$Xty, 0), yty = ss$yty, n = ss$n,
                  estimate_prior_variance = TRUE, estimate_residual_variance = TRUE)

  res2 = susie_ss(ss$XtX, ss$Xty, yty = ss$yty, n = ss$n,
                  estimate_prior_variance = TRUE, estimate_residual_variance = TRUE,
                  null_weight = 1/(ncol(ss$XtX)+1))
  expect_equal_susie_ss(res2,res1)
}))
