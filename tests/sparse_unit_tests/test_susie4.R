create_sparsity_mat = function(sparsity, n, p){
  nonzero = round(n*p*(1-sparsity))
  nonzero.idx = sample(n*p, nonzero)
  mat = numeric(n*p)
  mat[nonzero.idx] = 1
  mat = matrix(mat, nrow=n, ncol=p)
  return(mat)     
}

test_that("sparse version susie",{
  original.res = readRDS('../original_susie_results/susiefit_original_res4.rds')
  set.seed(1)
  n = 1000
  p = 10000
  beta = rep(0,p)
  beta[1]    = 10 
  beta[300]  = 10
  beta[400]  = 10
  beta[1000] = 10
  X = create_sparsity_mat(0.99,n,p)
  y = c(X %*% beta + rnorm(n))
  X.sparse = as(X,'dgCMatrix')
  
  dense.res = susie(X, y, standardize=FALSE, intercept = FALSE)
  sparse.res = susie(X.sparse, y, standardize=FALSE, intercept = FALSE)
  
  sparse.res$alpha = as.matrix(sparse.res$alpha, p, 1)
  sparse.res$mu = as.matrix(sparse.res$mu, p, 1)
  sparse.res$mu2 = as.matrix(sparse.res$mu2, p, 1)
  sparse.res$Xr = as.matrix(sparse.res$Xr, n, 1)
  sparse.res$fitted = as.matrix(sparse.res$fitted, n, 1)
  
  dense.res$alpha = as.matrix(dense.res$alpha, p, 1)
  dense.res$mu = as.matrix(dense.res$mu, p, 1)
  dense.res$mu2 = as.matrix(dense.res$mu2, p, 1)
  dense.res$Xr = as.matrix(dense.res$Xr, n, 1)
  dense.res$fitted = as.matrix(dense.res$fitted, n, 1)
  
  expect_equal(dense.res$alpha, original.res$alpha)
  expect_equal(dense.res$mu, original.res$mu)
  expect_equal(dense.res$mu2, original.res$mu2)
  expect_equal(dense.res$Xr, original.res$Xr)
  expect_equal(dense.res$KL, original.res$KL)
  expect_equal(dense.res$sigma2, original.res$sigma2)
  expect_equal(dense.res$V, original.res$V)
  expect_equal(dense.res$elbo, original.res$elbo)
  expect_equal(dense.res$niter, original.res$niter)
  expect_equal(dense.res$intercept, original.res$intercept)
  expect_equal(dense.res$fitted, original.res$fitted)
  expect_equal(dense.res$X_column_scale_factors, original.res$X_column_scale_factors)
  
  expect_equal(sparse.res$alpha, original.res$alpha)
  expect_equal(sparse.res$mu, original.res$mu)
  expect_equal(sparse.res$mu2, original.res$mu2)
  expect_equal(sparse.res$Xr, original.res$Xr)
  expect_equal(sparse.res$KL, original.res$KL)
  expect_equal(sparse.res$sigma2, original.res$sigma2)
  expect_equal(sparse.res$V, original.res$V)
  expect_equal(sparse.res$elbo, original.res$elbo)
  expect_equal(sparse.res$niter, original.res$niter)
  expect_equal(sparse.res$intercept, original.res$intercept)
  expect_equal(sparse.res$fitted, original.res$fitted)
  expect_equal(sparse.res$X_column_scale_factors, original.res$X_column_scale_factors)
})
