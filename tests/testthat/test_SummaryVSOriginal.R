test_that("Results from summary stat vs original data",{
  set.seed(1)
  n = 100
  p = 100
  beta = rep(0,p)
  beta[1:4] = 1
  X = matrix(rnorm(n*p),nrow=n,ncol=p); y = X %*% beta + rnorm(n)
  mean_y = mean(y); y = y-mean_y
  X = scale(X,center=TRUE, scale = FALSE)

  res = susie(X, y, intercept = FALSE, standardize = FALSE,
              estimate_residual_variance=FALSE, estimate_prior_variance = TRUE)

  res2 = susie_ss(t(X)%*%X, t(X)%*%y, var(y), max_iter = 3,
                  estimate_prior_variance = TRUE)
  expect_equal(coef(res2), coef(res))

  expect_equal(res2$sa2, res$sa2)
})
