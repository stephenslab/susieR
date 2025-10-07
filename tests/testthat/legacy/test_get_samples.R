context("test_get_samples.R")
test_that("Sampling from posterior distribution agrees with expected value", with(simulate(200, 500), {
  set.seed(1)
  res1 = susie(X, y)
  post_samples1 = susie_get_posterior_samples(res1, 10000)
  expect_equal(dim(post_samples1$b), c(500, 10000))
  expect_equal(dim(post_samples1$gamma), c(500, 10000))
  expect_equal(rowMeans(post_samples1$gamma), res1$pip)
  expect_equal(rowMeans(post_samples1$b),
               susie_get_posterior_mean(res1), tolerance=0.01)

  ss = univariate_regression(X, y)
  R = cor(X)
  z = ss$betahat/ss$sebetahat
  res2 = susie_rss(z,R, n=n)
  post_samples2 = susie_get_posterior_samples(res2, 10000)
  expect_equal(dim(post_samples2$b), c(500, 10000))
  expect_equal(dim(post_samples2$gamma), c(500, 10000))
  expect_equal(rowMeans(post_samples2$gamma), res2$pip, tolerance = 1e-4)
  expect_equal(rowMeans(post_samples2$b), susie_get_posterior_mean(res2), tolerance=0.01)
}))
