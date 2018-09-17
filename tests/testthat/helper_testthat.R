create_sparsity_mat = function(sparsity, n, p){
  nonzero = round(n*p*(1-sparsity))
  nonzero.idx = sample(n*p, nonzero)
  mat = numeric(n*p)
  mat[nonzero.idx] = 1
  mat = matrix(mat, nrow=n, ncol=p)
  return(mat)
}

simulate = function(n=100, p=200, sparse=F) {
  set.seed(1)
  n = n
  p = p
  beta = rep(0,p)
  beta[1]    = 10
  beta[2]  = 10
  beta[3]  = 10
  beta[4] = 10
  if (sparse) {
    X = create_sparsity_mat(0.99,n,p)
    X.sparse = as(X,'dgCMatrix')
  } else {
    X = matrix(rnorm(n*p,3,4),n,p)
    X.sparse = NA
  }
  y = c(X %*% beta + rnorm(n))
  L = 10
  residual_variance = 0.8
  scaled_prior_variance = 0.2
  s = list(alpha=matrix(1/p,nrow=L,ncol=p),
           mu=matrix(2,nrow=L,ncol=p),
           mu2=matrix(3,nrow=L,ncol=p),
           Xr=rep(5,n), KL=rep(1.2,L),
           sigma2=residual_variance, V=scaled_prior_variance * as.numeric(var(y)))
  attach(list(X=X, X.sparse=X.sparse, s=s, y=y, n=n, p=p, b=beta), warn.conflict=F)
}

load_data = function(filename) {
  readRDS(system.file("inst","datafiles",filename,package = "susieR"))
}

is_equal_s = function(sparse.res, original.res){
  sparse.res$alpha = as.matrix(sparse.res$alpha, p, 1)
  sparse.res$mu = as.matrix(sparse.res$mu, p, 1)
  sparse.res$mu2 = as.matrix(sparse.res$mu2, p, 1)
  sparse.res$Xr = as.matrix(sparse.res$Xr, n, 1)

  expect_equal(sparse.res$alpha, original.res$alpha)
  expect_equal(sparse.res$mu, original.res$mu)
  expect_equal(sparse.res$mu2, original.res$mu2)
  expect_equal(sparse.res$Xr, original.res$Xr)
  expect_equal(sparse.res$KL, original.res$KL)
  expect_equal(sparse.res$sigma2, original.res$sigma2)
  expect_equal(sparse.res$V, original.res$V)
}

is_equal_SER = function(sparse.res, original.res){
  sparse.res$alpha = as.matrix(sparse.res$alpha, p, 1)
  sparse.res$mu = as.matrix(sparse.res$mu, p, 1)
  sparse.res$mu2 = as.matrix(sparse.res$mu2, p, 1)
  sparse.res$lbf = as.matrix(sparse.res$lbf, p, 1)
  
  expect_equal(sparse.res$alpha, original.res$alpha)
  expect_equal(sparse.res$mu, original.res$mu)
  expect_equal(sparse.res$mu2, original.res$mu2)
  expect_equal(sparse.res$lbf, original.res$lbf)
  expect_equal(sparse.res$V, original.res$V)
  expect_equal(sparse.res$loglik, original.res$loglik)
}

is_equal_susie = function(sparse.res, dense.res, original.res){
  sparse.res$alpha = as.matrix(sparse.res$alpha, p, 1)
  sparse.res$mu = as.matrix(sparse.res$mu, p, 1)
  sparse.res$mu2 = as.matrix(sparse.res$mu2, p, 1)
  sparse.res$Xr = as.matrix(sparse.res$Xr, n, 1)
  sparse.res$fitted = as.matrix(sparse.res$fitted, n, 1)
  
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
}
