context("test_prior_weights.R")

test_that("prior weights specification agrees with default", with(simulate(200,1000), {
  res1 = susie(X, y, estimate_prior_variance = TRUE)
  res2 = susie(X, y, estimate_prior_variance = TRUE,
               prior_weights = rep(1/ncol(X), ncol(X)))
  expect_equal_susie(res1,res2,tolerance = 0.1)
}))

test_that("Sufficient stat (Xty): prior weights specification agrees with default", with(simulate(200,1000), {
  ss = compute_ss(X, y)
  res1 = susie_suff_stat(XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = n, estimate_prior_variance = TRUE)
  res2 = susie_suff_stat(XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = n, estimate_prior_variance = TRUE,
                  prior_weights = rep(1/ncol(ss$XtX), ncol(ss$XtX)))
  expect_equal_susie_suff_stat(res1,res2,tolerance = 1e-6)
}))

test_that("Sufficient stat (beta): prior weights specification agrees with default", with(simulate(200,1000), {
  ss = univariate_regression(X, y)
  R = cor(X)
  res1 = susie_rss(bhat = ss$betahat,shat = ss$sebetahat,R = R,
                   n = n, var_y = var(y), estimate_prior_variance = TRUE)
  res2 = susie_rss(bhat = ss$betahat, shat = ss$sebetahat, R = R, n = n, var_y = var(y),
                   estimate_prior_variance = TRUE,
                   prior_weights = rep(1/ncol(R), ncol(R)))
  expect_equal_susie_suff_stat(res1,res2,tolerance = 0.001)
}))

test_that("RSS: prior weights specification agrees with default", with(simulate(200,500), {
  ss = univariate_regression(X, y)
  R = cor(X)
  set.seed(1) # otherwise it might fail because the two results will be slighly different
  res1 = susie_rss(z = ss$betahat/ss$sebetahat, R = R, n=n, estimate_prior_variance = TRUE)
  res2 = susie_rss(z = ss$betahat/ss$sebetahat, R = R, n=n, estimate_prior_variance = TRUE,
                   prior_weights = rep(1/ncol(R), ncol(R)))
  expect_equal_susie_suff_stat(res1,res2,tolerance = 1e-05)
}))
