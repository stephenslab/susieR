context("test_susie_ss_interface.R")

test_that("Results from ss bhat interface vs original data: no standardize",{
  simulate(200,1000)
  ss = univariate_regression(X, y)

  R = cor(X)

  orig = susie(X, y, intercept = TRUE, standardize = FALSE, max_iter = 2,
              estimate_residual_variance=FALSE, estimate_prior_variance = TRUE)

  fit = susie_bhat(bhat = ss$betahat, shat = ss$sebetahat, R = R, R_type = 'Cor',
                   var_y = var(y), n = n, standardize = FALSE,
                   max_iter = 2, estimate_prior_variance = TRUE)

  expect_equal(fit$alpha, orig$alpha)
  expect_equal(fit$mu, orig$mu)
  expect_equal(fit$mu2, orig$mu2)
  expect_equal(fit$V, orig$V)
  cm = colMeans(X, na.rm = TRUE)
  X.c = t(t(X) - cm)
  expect_equal(crossprod(X.c, orig$fitted), fit$Xtfitted)
})

test_that("Results from ss bhat interface vs original data: standardize",{
  simulate(200,1000)
  ss = univariate_regression(X, y)
  R = cor(X)
  cm = colMeans(X, na.rm = TRUE)
  csd = matrixStats::colSds(X, center = cm)
  csd[csd==0] = 1
  X.s = t(t(X)/csd)

  orig = susie(X.s, y/sd(y), intercept = TRUE, standardize = FALSE, max_iter = 2,
               estimate_residual_variance=FALSE, estimate_prior_variance = TRUE)

  fit = susie_bhat(bhat = ss$betahat, shat = ss$sebetahat, R = R, R_type = 'Cor',
                   var_y = 1, n = n, standardize = TRUE,
                   max_iter = 2, estimate_prior_variance = TRUE)

  expect_equal(fit$alpha, orig$alpha)
  expect_equal(fit$mu, orig$mu)
  expect_equal(fit$mu2, orig$mu2)
  expect_equal(fit$V, orig$V)

  X.cs = t((t(X) - cm)/csd)
  expect_equal(crossprod(X.cs, orig$fitted), fit$Xtfitted)
})

