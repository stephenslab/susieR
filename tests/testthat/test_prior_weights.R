context("test_prior_weights.R")

test_that("prior weights specification agrees with default", with(simulate(200,1000), {
  res1 = susie(X, y, estimate_prior_variance = TRUE)
  res2 = susie(X, y, estimate_prior_variance = TRUE,
               prior_weights = rep(1/ncol(X), ncol(X)))
  expect_equal_susie(res1,res2)
}))

test_that("SS: prior weights specification agrees with default", with(simulate(200,1000), {
  ss = compute_ss(X, y)
  res1 = susie_ss(ss$XtX, ss$Xty, yty = ss$yty, n = n, estimate_prior_variance = TRUE)
  res2 = susie_ss(ss$XtX, ss$Xty, yty = ss$yty, n = n, estimate_prior_variance = TRUE,
                  prior_weights = rep(1/ncol(ss$XtX), ncol(ss$XtX)))
  expect_equal_susie_ss(res1,res2)
}))
