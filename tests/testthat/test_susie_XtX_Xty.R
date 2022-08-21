context("test_susie_XtX_Xty.R")

test_that(paste("compute_ss XtX calculations are correct with",
                "standardize = FALSE and with standardize = TRUE"),{
  set.seed(1)
  n <- 10
  p <- 5
  X <- matrix(rnorm(n*p),n,p)
  y <- rnorm(n)
  Y1 <- scale(X,center = TRUE,scale = FALSE)
  Y2 <- scale(X,center = TRUE,scale = TRUE)
  out1 <- compute_suff_stat(X,y,standardize = FALSE)
  out2 <- compute_suff_stat(X,y,standardize = TRUE)
  dimnames(out1$XtX) <- NULL
  dimnames(out2$XtX) <- NULL
  expect_equal(out1$XtX,crossprod(Y1),scale = 1,tolerance = 1e-14)
  expect_equal(out2$XtX,crossprod(Y2),scale = 1,tolerance = 1e-14)

  # Run the same checks again, but with X now being a sparse matrix.
  X    <- as(X,"CsparseMatrix")
  out1 <- compute_suff_stat(X,y,standardize = FALSE)
  out2 <- compute_suff_stat(X,y,standardize = TRUE)
  dimnames(out1$XtX) <- NULL
  dimnames(out2$XtX) <- NULL
  expect_equal(out1$XtX,crossprod(Y1),scale = 1,tolerance = 1e-14)
  expect_equal(out2$XtX,crossprod(Y2),scale = 1,tolerance = 1e-14)
})

test_that("Results from sufficient stat vs original data", with(simulate(200,1000), {
  ss = compute_ss(X, y, standardize = FALSE)

  expect_warning(res <- susie(X, y, intercept = TRUE, standardize = TRUE, max_iter = 2,
                              estimate_residual_variance=FALSE, estimate_prior_variance = FALSE))

  expect_warning(res2 <- susie_suff_stat(XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty,
                                  n = ss$n, standardize = TRUE, max_iter = 2,
                                  estimate_prior_variance =FALSE, estimate_residual_variance = FALSE))

  expect_equal(res$alpha, res2$alpha)
  expect_equal(res$mu, res2$mu)
  expect_equal(res$mu2, res2$mu2)
  expect_equal(res$V, res2$V)
  expect_equal(res$elbo, res2$elbo)
  expect_equal(coef(res)[-1], coef(res2)[-1])

  cm  = colMeans(X, na.rm = TRUE)
  csd = compute_colSds(X)
  csd[csd==0] = 1
  X.cs = t( (t(X) - cm) / csd )
  expect_equal(crossprod(X.cs, res$Xr), res2$XtXr)
}))

test_that("Results from sufficient stat vs original data: estimate residual variance", with(simulate(200,1000), {
  ss = compute_ss(X, y, standardize = FALSE)

  expect_warning(res <- susie(X, y, intercept = TRUE, standardize = TRUE, max_iter = 2,
              estimate_residual_variance=TRUE, estimate_prior_variance = FALSE))

  expect_warning(res2 <- susie_suff_stat(XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty,
                  n = ss$n, standardize = TRUE, max_iter = 2,
                  estimate_prior_variance = FALSE, estimate_residual_variance = TRUE))

  expect_equal(res$alpha, res2$alpha)
  expect_equal(res$mu, res2$mu)
  expect_equal(res$mu2, res2$mu2)
  expect_equal(res$V, res2$V)
  expect_equal(res$elbo, res2$elbo)
  expect_equal(res$sigma2, res2$sigma2)
  expect_equal(coef(res)[-1], coef(res2)[-1])

  cm = colMeans(X, na.rm = TRUE)
  csd = compute_colSds(X)
  csd[csd==0] = 1
  X.cs = t( (t(X) - cm) / csd )
  expect_equal(crossprod(X.cs, res$Xr), res2$XtXr)
}))

test_that("MAF filter works", with(simulate(200,1000), {
  X.maf = cbind(X, 0)
  ss = compute_ss(X, y, standardize = FALSE)
  ss.maf = compute_ss(X.maf, y, standardize = FALSE)
  maf = c(rep(1,1000),0.01)
  res1 <- susie_suff_stat(XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty,
                   n = ss$n, standardize = TRUE,
                   estimate_prior_variance =FALSE, estimate_residual_variance = FALSE)

  res2 <- susie_suff_stat(XtX = ss.maf$XtX, Xty = ss.maf$Xty, yty = ss$yty,
                   n = ss$n, maf_thresh = 0.05, maf = maf, standardize = TRUE,
                   estimate_prior_variance =FALSE, estimate_residual_variance = FALSE)

  expect_equal(res1$alpha, res2$alpha)
  expect_equal(res1$mu, res2$mu)
  expect_equal(res1$mu2, res2$mu2)
  expect_equal(res1$V, res2$V)
  expect_equal(res1$elbo, res2$elbo)
  expect_equal(coef(res1), coef(res2))
  expect_equal(res1$XtXr, res2$XtXr)

}))
