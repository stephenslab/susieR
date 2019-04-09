context("test_susie_ss_interface.R")

test_that("Results from ss bhat interface vs original data: no standardize", with(simulate(200,1000), {
  ss = univariate_regression(X, y)

  R = cor(X)

  expect_warning(orig <- susie(X, y, intercept = TRUE, standardize = FALSE, max_iter = 2,
              estimate_residual_variance=FALSE, estimate_prior_variance = FALSE))

  expect_warning(fit <- susie_bhat(bhat = ss$betahat, shat = ss$sebetahat, R = R,
                   var_y = var(y), n = n, standardize = FALSE,
                   max_iter = 2, estimate_prior_variance = FALSE,
                   estimate_residual_variance = FALSE))

  expect_equal(fit$alpha, orig$alpha)
  expect_equal(fit$mu, orig$mu)
  expect_equal(fit$mu2, orig$mu2)
  expect_equal(fit$V, orig$V)
  X.c = set_X_attributes(X, center = TRUE, scale = FALSE)
  X.c = t((t(X.c) - attr(X.c, 'scaled:center'))/attr(X.c, 'scaled:scale'))
  expect_equal(crossprod(X.c, orig$fitted), fit$Xtfitted)
}))

test_that("Results from ss bhat interface vs original data: standardize", with(simulate(200,1000), {
  ss = univariate_regression(X, y)
  R = cor(X)

  X.s = set_X_attributes(X, center = FALSE, scale = TRUE)
  X.s = t((t(X.s) - attr(X.s, 'scaled:center'))/attr(X.s, 'scaled:scale'))
  X.cs = set_X_attributes(X, center = TRUE, scale = TRUE)
  X.cs = t((t(X.cs) - attr(X.cs, 'scaled:center'))/attr(X.cs, 'scaled:scale'))

  expect_warning(orig <- susie(X.s, y, intercept = TRUE, standardize = TRUE, max_iter = 2,
               estimate_residual_variance=FALSE, estimate_prior_variance = FALSE))

  expect_warning(fit <- susie_bhat(bhat = ss$betahat, shat = ss$sebetahat, R = R,
                   n = n, var_y = var(y), standardize = TRUE,
                   max_iter = 2, estimate_prior_variance = FALSE,
                   estimate_residual_variance = FALSE))

  expect_equal(fit$alpha, orig$alpha)
  expect_equal(fit$mu, orig$mu)
  expect_equal(fit$mu2, orig$mu2)
  expect_equal(fit$V, orig$V)
  expect_equal(crossprod(X.cs, orig$fitted), fit$Xtfitted)
}))

test_that("Results from ss bhat interface: t statistics", with(simulate(200,1000), {
  ss = univariate_regression(X, y)
  R = cor(X)

  X.s = set_X_attributes(X, center = FALSE, scale = TRUE)
  X.s = t((t(X.s) - attr(X.s, 'scaled:center'))/attr(X.s, 'scaled:scale'))
  X.cs = set_X_attributes(X, center = TRUE, scale = TRUE)
  X.cs = t((t(X.cs) - attr(X.cs, 'scaled:center'))/attr(X.cs, 'scaled:scale'))

  expect_warning(orig <- susie(X.s, y/sd(y), intercept = TRUE, standardize = TRUE, max_iter = 2,
               estimate_residual_variance=FALSE, estimate_prior_variance = FALSE))

  expect_warning(fit <- susie_bhat(bhat = ss$betahat/ss$sebetahat, shat = 1, R = R,
                   n = n, standardize = TRUE,
                   max_iter = 2, estimate_prior_variance = FALSE,
                   estimate_residual_variance = FALSE))

  expect_equal(fit$alpha, orig$alpha)
  expect_equal(fit$mu, orig$mu)
  expect_equal(fit$mu2, orig$mu2)
  expect_equal(fit$V, orig$V)
  expect_equal(crossprod(X.cs, orig$fitted), fit$Xtfitted)
}))
