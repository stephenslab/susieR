context("test_null_weight.R")

test_that("null weight specification agrees with default", with(simulate(200,1000), {
  res1 = susie(cbind(X,0), y, estimate_prior_variance = TRUE)
  res2 = susie(X, y, estimate_prior_variance = TRUE, null_weight = 1/(ncol(X)+1))
  expect_equal_susie(res2,res1,tolerance = 1e-6)
}))

test_that("Sufficient stat (Xty): null weight specification agrees with default", with(simulate(200,1000), {
  ss = compute_ss(X,y)

  res1 = susie_suff_stat(XtX = cbind(rbind(ss$XtX,0),0), Xty = c(ss$Xty, 0), yty = ss$yty, n = ss$n,
                  estimate_prior_variance = TRUE, estimate_residual_variance = TRUE)

  res2 = susie_suff_stat(XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n,
                  estimate_prior_variance = TRUE, estimate_residual_variance = TRUE,
                  null_weight = 1/(ncol(ss$XtX)+1))
  expect_equal_susie_suff_stat(res2,res1,tolerance = 1e-6)
}))

test_that("RSS: null weight specification agrees with default", with(simulate(200,500), {
  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat
  res1 = susie_rss(z = c(z, 0), R = cbind(rbind(R,0),0), n=n,
                   estimate_prior_variance = TRUE)

  res2 = susie_rss(z = z, R = R,n=n,
                   estimate_prior_variance = TRUE,
                   null_weight = 1/(ncol(R)+1))
  expect_equal_susie_suff_stat(res2,res1,tolerance = 1e-04)
}))
