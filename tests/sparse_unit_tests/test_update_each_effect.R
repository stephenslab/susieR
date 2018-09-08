create_sparsity_mat = function(sparsity, n, p){
  nonzero = round(n*p*(1-sparsity))
  nonzero.idx = sample(n*p, nonzero)
  mat = numeric(n*p)
  mat[nonzero.idx] = 1
  mat = matrix(mat, nrow=n, ncol=p)
  return(mat)     
}

test_that("sparse version update_each_effect",{
  original.res = readRDS('../original_susie_results/vbupdate_original_res.rds')
  set.seed(1)
  n = 1000
  p = 10000
  beta = rep(0,p)
  beta[1]    = 10 
  beta[300]  = 10
  beta[400]  = 10
  beta[1000] = 10
  X.dense = create_sparsity_mat(0.99,n,p)
  y = c(X.dense %*% beta + rnorm(n))
  L = 10
  residual_variance = 0.8
  scaled_prior_variance = 0.2
  s = list(alpha=matrix(1/p,nrow=L,ncol=p),
           mu=matrix(2,nrow=L,ncol=p),
           mu2=matrix(3,nrow=L,ncol=p),
           Xr=rep(5,n), KL=rep(1.2,L),
           sigma2=residual_variance, V=scaled_prior_variance * as.numeric(var(y)))
  X.sparse = as(X.dense,'dgCMatrix')

  
  scaledX.dense = susieR:::safe_colScale(X.dense)
  scaledX.sparse = susieR:::safe_colScale(X.sparse)
  
  dense.res = susieR:::update_each_effect(scaledX.dense,y,s)
  sparse.res = susieR:::update_each_effect(scaledX.sparse,y,s)
  
  sparse.res$alpha = as.matrix(sparse.res$alpha, p, 1)
  sparse.res$mu = as.matrix(sparse.res$mu, p, 1)
  sparse.res$mu2 = as.matrix(sparse.res$mu2, p, 1)
  sparse.res$Xr = as.matrix(sparse.res$Xr, n, 1)
  
  dense.res$alpha = as.matrix(dense.res$alpha, p, 1)
  dense.res$mu = as.matrix(dense.res$mu, p, 1)
  dense.res$mu2 = as.matrix(dense.res$mu2, p, 1)
  dense.res$Xr = as.matrix(dense.res$Xr, n, 1)
  
  expect_equal(dense.res$alpha, original.res$alpha)
  expect_equal(dense.res$mu, original.res$mu)
  expect_equal(dense.res$mu2, original.res$mu2)
  expect_equal(dense.res$Xr, original.res$Xr)
  expect_equal(dense.res$KL, original.res$KL)
  expect_equal(dense.res$sigma2, original.res$sigma2)
  expect_equal(dense.res$V, original.res$V)
  
  expect_equal(sparse.res$alpha, original.res$alpha)
  expect_equal(sparse.res$mu, original.res$mu)
  expect_equal(sparse.res$mu2, original.res$mu2)
  expect_equal(sparse.res$Xr, original.res$Xr)
  expect_equal(sparse.res$KL, original.res$KL)
  expect_equal(sparse.res$sigma2, original.res$sigma2)
  expect_equal(sparse.res$V, original.res$V)
  
})
