test_that("sparse version get_ER2",{
  original.res = readRDS('../original_susie_results/ER2_original_res.rds')
  set.seed(1)
  n = 1000
  p = 10000
  beta = rep(0,p)
  beta[1]    = 10 
  beta[300]  = 10
  beta[400]  = 10
  beta[1000] = 10
  X.dense = create_sparsity_mat(0.99,n,p)
  X.sparse = as(X.dense,'dgCMatrix')
  y = c(X.dense %*% beta + rnorm(n))
  L = 10
  residual_variance = 0.8
  scaled_prior_variance = 0.2
  s = list(alpha=matrix(1/p,nrow=L,ncol=p),
           mu=matrix(2,nrow=L,ncol=p),
           mu2=matrix(3,nrow=L,ncol=p),
           Xr=rep(5,n), KL=rep(1.2,L),
           sigma2=residual_variance, V=scaled_prior_variance * as.numeric(var(y)))
  
  scaledX.dense = susieR:::safe_colScale(X.dense)
  scaledX.sparse = susieR:::safe_colScale(X.sparse)
  
  dense.res = susieR:::get_ER2(scaledX.dense, y, s)
  sparse.res = susieR:::get_ER2(scaledX.sparse, y, s)
  
  expect_equal(dense.res, original.res)
  expect_equal(sparse.res, original.res)
  
})
