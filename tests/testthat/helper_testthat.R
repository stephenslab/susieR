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
  beta = rep(0,p)
  beta[1:4] = 10
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
  readRDS(system.file("datafiles",filename,package = "susieR"))
}

is_equal_susie_update = function(new.res, original.res){
  new.res$alpha = as.matrix(new.res$alpha, p, 1)
  new.res$mu = as.matrix(new.res$mu, p, 1)
  new.res$mu2 = as.matrix(new.res$mu2, p, 1)
  new.res$Xr = as.matrix(new.res$Xr, n, 1)

  expect_equal(new.res$alpha, original.res$alpha)
  expect_equal(new.res$mu, original.res$mu)
  expect_equal(new.res$mu2, original.res$mu2)
  expect_equal(new.res$Xr, original.res$Xr)
  expect_equal(new.res$KL, original.res$KL)
  expect_equal(new.res$sigma2, original.res$sigma2)
  expect_equal(new.res$V, original.res$V)
}

is_equal_SER = function(new.res, original.res){
  new.res$alpha = as.matrix(new.res$alpha, p, 1)
  new.res$mu = as.matrix(new.res$mu, p, 1)
  new.res$mu2 = as.matrix(new.res$mu2, p, 1)
  new.res$lbf = as.matrix(new.res$lbf, p, 1)
  
  expect_equal(new.res$alpha, original.res$alpha)
  expect_equal(new.res$mu, original.res$mu)
  expect_equal(new.res$mu2, original.res$mu2)
  expect_equal(new.res$lbf, original.res$lbf)
  expect_equal(new.res$V, original.res$V)
  expect_equal(new.res$loglik, original.res$loglik)
}

is_equal_susie = function(new.res, original.res){
  is_equal_susie_update(new.res, original.res)
  new.res$fitted = as.matrix(new.res$fitted, n, 1)
  expect_equal(new.res$elbo, original.res$elbo)
  expect_equal(new.res$niter, original.res$niter)
  expect_equal(new.res$intercept, original.res$intercept)
  expect_equal(new.res$fitted, original.res$fitted)
  expect_equal(new.res$X_column_scale_factors, original.res$X_column_scale_factors)
}
