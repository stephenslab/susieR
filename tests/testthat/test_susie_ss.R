context("test_susie_ss.R")

test_that("Results from summary stat vs original data",{
  simulate(200,1000)
  ss = compute_ss(X, y)

  res = susie(X, y, intercept = TRUE, standardize = TRUE, max_iter = 2,
              estimate_residual_variance=FALSE, estimate_prior_variance = TRUE)

  res2 = susie_ss(ss$XtX, ss$Xty, var_y = ss$vary,
                  residual_variance = var(y), n = ss$n, max_iter = 2,
                  estimate_prior_variance = TRUE)
  expect_equal(res$alpha, res2$alpha)
  expect_equal(res$mu, res2$mu)
  expect_equal(res$mu2, res2$mu2)
  expect_equal(res$V, res2$V)
  expect_equal(res$elbo, res2$elbo)
  cm = colMeans(X, na.rm = TRUE)
  csd = matrixStats::colSds(X, center = cm)
  csd[csd==0] = 1
  X.cs = t( (t(X) - cm) / csd )
  expect_equal(crossprod(X.cs, res$fitted), res2$Xtfitted)
})
