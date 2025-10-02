context("test_susie_beta_se.R")

test_that("Results from sufficient stat vs original data: no standardize", with(simulate(200,1000), {
  ss = univariate_regression(X, y)
  suff_stat = compute_suff_stat(X, y)

  R = cor(X)

  orig <- susie(X, y, intercept = TRUE, standardize = FALSE,
                estimate_residual_variance=FALSE, estimate_prior_variance = FALSE)

  fit_suff <- susie_suff_stat(XtX = suff_stat$XtX,
                              Xty = suff_stat$Xty,
                              yty = suff_stat$yty,
                              n = suff_stat$n,
                              standardize = FALSE,
                              estimate_prior_variance = FALSE,
                              estimate_residual_variance = FALSE)

  fit_rss <- susie_rss(bhat = ss$betahat,
                   shat = ss$sebetahat, R = R,
                   var_y = var(y), n = n, standardize = FALSE,
                   estimate_prior_variance = FALSE,
                   estimate_residual_variance = FALSE)

  expect_equal(fit_suff$alpha, orig$alpha)
  expect_equal(fit_suff$mu, orig$mu)
  expect_equal(fit_suff$mu2, orig$mu2)
  expect_equal(fit_suff$V, orig$V)
  X.c = set_X_attributes(X, center = TRUE, scale = FALSE)
  X.c = t((t(X.c) - attr(X.c, 'scaled:center'))/attr(X.c, 'scaled:scale'))
  expect_equal(crossprod(X.c, orig$Xr), fit_suff$XtXr)

  expect_equal(fit_rss$alpha, orig$alpha)
  expect_equal(fit_rss$mu, orig$mu)
  expect_equal(fit_rss$mu2, orig$mu2)
  expect_equal(fit_rss$V, orig$V)
  expect_equal(crossprod(X.c, orig$Xr), fit_rss$XtXr)
}))

test_that("Results from sufficient stat interface vs original data: standardize", with(simulate(200,1000), {
  ss = univariate_regression(X, y)
  suff_stat = compute_suff_stat(X, y)
  R = cor(X)

  X.cs = set_X_attributes(X, center = TRUE, scale = TRUE)
  X.cs = t((t(X.cs) - attr(X.cs, 'scaled:center'))/attr(X.cs, 'scaled:scale'))

  orig <- susie(X, y, intercept = TRUE, standardize = TRUE,
                estimate_residual_variance=FALSE, estimate_prior_variance = FALSE)

  fit_suff <- susie_suff_stat(XtX = suff_stat$XtX,
                              Xty = suff_stat$Xty,
                              yty = suff_stat$yty,
                              n = suff_stat$n,
                              standardize = TRUE,
                              estimate_prior_variance = FALSE,
                              estimate_residual_variance = FALSE)

  fit_rss <- susie_rss(bhat = ss$betahat, shat = ss$sebetahat, R = R,
                       n = n, var_y = var(y), standardize = TRUE,
                       estimate_prior_variance = FALSE,
                       estimate_residual_variance = FALSE)

  expect_equal(fit_suff$alpha, orig$alpha)
  expect_equal(fit_suff$mu, orig$mu)
  expect_equal(fit_suff$mu2, orig$mu2)
  expect_equal(fit_suff$V, orig$V)
  expect_equal(crossprod(X.cs, orig$Xr), fit_suff$XtXr)

  expect_equal(fit_rss$alpha, orig$alpha)
  expect_equal(fit_rss$mu, orig$mu)
  expect_equal(fit_rss$mu2, orig$mu2)
  expect_equal(fit_rss$V, orig$V)
  expect_equal(crossprod(X.cs, orig$Xr), fit_rss$XtXr)
}))

test_that("Results from sufficient stat interface vs original data on standaridzed input", with(simulate(200,1000), {
  ss = univariate_regression(X, y)
  y = y - mean(y)
  y = y/sd(y)
  suff_stat = compute_suff_stat(X, y, standardize = TRUE)
  R = cor(X)

  mu  = colMeans(X)
  s   = compute_colSds(X)
  X.cs = t((t(X) - mu)/s)

  orig <- susie(X.cs, y, intercept = FALSE, standardize = FALSE,
               estimate_residual_variance=FALSE, estimate_prior_variance = FALSE)

  fit_suff <- susie_suff_stat(XtX = suff_stat$XtX,
                              Xty = suff_stat$Xty,
                              yty = suff_stat$yty,
                              n = suff_stat$n,
                              standardize = FALSE,
                              estimate_prior_variance = FALSE,
                              estimate_residual_variance = FALSE)

  fit_rss <- susie_rss(bhat = ss$betahat, shat = ss$sebetahat, R = R,
                       n = n, standardize = TRUE,
                       estimate_prior_variance = FALSE,
                       estimate_residual_variance = FALSE)

  expect_equal(fit_suff$alpha, orig$alpha)
  expect_equal(fit_suff$mu, orig$mu)
  expect_equal(fit_suff$mu2, orig$mu2)
  expect_equal(fit_suff$V, orig$V)
  expect_equal(crossprod(X.cs, orig$Xr), fit_suff$XtXr)

  expect_equal(fit_rss$alpha, orig$alpha)
  expect_equal(fit_rss$mu, orig$mu)
  expect_equal(fit_rss$mu2, orig$mu2)
  expect_equal(fit_rss$V, orig$V)
  expect_equal(crossprod(X.cs, orig$Xr), fit_rss$XtXr)
}))
