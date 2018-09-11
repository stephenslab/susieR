create_sparsity_mat = function(sparsity, n, p){
  nonzero = round(n*p*(1-sparsity))
  nonzero.idx = sample(n*p, nonzero)
  mat = numeric(n*p)
  mat[nonzero.idx] = 1
  mat = matrix(mat, nrow=n, ncol=p)
  return(mat)     
}

test_that("sparse version single_effect_regression",{
  original.res = readRDS('../original_susie_results/singleReg_original_res.rds')
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
  X.sparse = as(X.dense,'dgCMatrix')
  V = scaled_prior_variance
  
  scaledX.dense = susieR:::safe_colScale(X.dense)
  scaledX.sparse = susieR:::safe_colScale(X.sparse)
  
  dense.res = susieR:::single_effect_regression(y,scaledX.dense,V)
  sparse.res = susieR:::single_effect_regression(y,scaledX.sparse,V)
  
  sparse.res$alpha = as.matrix(sparse.res$alpha, p, 1)
  sparse.res$mu = as.matrix(sparse.res$mu, p, 1)
  sparse.res$mu2 = as.matrix(sparse.res$mu2, p, 1)
  sparse.res$lbf = as.matrix(sparse.res$lbf, p, 1)
  
  dense.res$alpha = as.matrix(dense.res$alpha, p, 1)
  dense.res$mu = as.matrix(dense.res$mu, p, 1)
  dense.res$mu2 = as.matrix(dense.res$mu2, p, 1)
  dense.res$lbf = as.matrix(dense.res$lbf, p, 1)
  
  expect_equal(dense.res$alpha, original.res$alpha)
  expect_equal(dense.res$mu, original.res$mu)
  expect_equal(dense.res$mu2, original.res$mu2)
  expect_equal(dense.res$lbf, original.res$lbf)
  expect_equal(dense.res$V, original.res$V)
  expect_equal(dense.res$loglik, original.res$loglik)
  
  expect_equal(sparse.res$alpha, original.res$alpha)
  expect_equal(sparse.res$mu, original.res$mu)
  expect_equal(sparse.res$mu2, original.res$mu2)
  expect_equal(sparse.res$lbf, original.res$lbf)
  expect_equal(sparse.res$V, original.res$V)
  expect_equal(sparse.res$loglik, original.res$loglik)
  
})
