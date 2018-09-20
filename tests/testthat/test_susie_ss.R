context("test_susie_ss.R")

test_that("Results from summary stat vs original data",{
  simulate(200,1000)
  X = safe_colScale(X,center=TRUE, scale = TRUE)

  res = susie(X, y, intercept = FALSE, standardize = FALSE, max_iter = 3,
              estimate_residual_variance=FALSE, estimate_prior_variance = TRUE)

  res2 = susie_ss(t(X)%*%X, c(y %*% X), var_y = sum(y^2)/n,
                  residual_variance = var(y), n = n, max_iter = 3,
                  estimate_prior_variance = TRUE)
  expect_equal(coef(res2), coef(res))
  expect_equal(res2$V, res$V)
  expect_equal(res2$elbo, res$elbo)
})
