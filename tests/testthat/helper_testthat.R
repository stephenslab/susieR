create_sparsity_mat = function(sparsity, n, p){
  nonzero = round(n*p*(1-sparsity))
  nonzero.idx = sample(n*p, nonzero)
  mat = numeric(n*p)
  mat[nonzero.idx] = 1
  mat = matrix(mat, nrow=n, ncol=p)
  return(mat)
}

simulate = function(n=100, p=200, sparse=F) {
  suppressWarnings(RNGversion("3.5.0"))
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
           sigma2=residual_variance,
      V=scaled_prior_variance * as.numeric(var(y)))
  return(list(X=X, X.sparse=X.sparse, s=s, y=y, n=n, p=p, b=beta))
}

simulate_tf = function(order){
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(2)
  n = 50
  D = diag(-1, n)
  for (i in 1:(n-1)){
    D[i, i+1] = 1
  }
  if (order==0) {
    beta = c(rep(0,5),rep(1,5),rep(3,5),rep(-2,5),rep(0,30))
    y = beta + rnorm(n)
    X = solve(D)
  } else if (order==1) {
    beta = numeric(n)
    for (i in 1:n){
      if (i <= 5){
        beta[i] = 0.001*i + 2
      } else if (i <= 15){
        beta[i] = 5*0.001*i + 1.6
      } else{
        beta[i] = 6.1 - 10*0.001*i
      }
    }
    y = beta + rnorm(n)
    X = solve(D%*%D)
  } else if (order==2) {
    beta = numeric(n)
    for (i in 1:n){
      if (i <= 5){
        beta[i] = (0.001*i)^2
      } else if (i <= 35){
        beta[i] = -5*(0.001*i)^2 + 0.06
      } else{
        beta[i] = 3*(0.001*i)^2 - 3.86
      }
    }
    y = beta + rnorm(n)
    X = solve(D%*%D%*%D)
  }
  return(list(X=X, y=y))
}

expect_equal_susie_update = function(new.res, original.res){
  expect_equal(new.res$alpha, original.res$alpha)
  expect_equal(new.res$mu, original.res$mu)
  expect_equal(new.res$mu2, original.res$mu2)
  expect_equal(new.res$Xr, original.res$Xr)
  expect_equal(new.res$KL, original.res$KL)
  expect_equal(new.res$sigma2, original.res$sigma2)
  expect_equal(new.res$V, original.res$V)
}

expect_equal_susie_ss_update = function(new.res, original.res){
  expect_equal(new.res$alpha, original.res$alpha)
  expect_equal(new.res$mu, original.res$mu)
  expect_equal(new.res$mu2, original.res$mu2)
  expect_equal(new.res$XtXr, original.res$XtXr)
  expect_equal(new.res$KL, original.res$KL)
  expect_equal(new.res$sigma2, original.res$sigma2)
  expect_equal(new.res$V, original.res$V)
}

expect_equal_SER = function(new.res, original.res){
  expect_equal(new.res$alpha, original.res$alpha)
  expect_equal(new.res$mu, original.res$mu)
  expect_equal(new.res$mu2, original.res$mu2)
  expect_equal(new.res$lbf, original.res$lbf)
  expect_equal(new.res$V, original.res$V)
  expect_equal(new.res$loglik, original.res$loglik)
}

expect_equal_SER_ss = function(new.res, original.res){
  expect_equal(new.res$alpha, original.res$alpha)
  expect_equal(new.res$mu, original.res$mu)
  expect_equal(new.res$mu2, original.res$mu2)
  expect_equal(new.res$lbf, original.res$lbf)
  expect_equal(new.res$V, original.res$V)
  expect_equal(new.res$lbf_model, original.res$lbf_model)
}

expect_equal_susie = function(new.res, original.res){
  expect_equal_susie_update(new.res, original.res)
  expect_equal(new.res$elbo, original.res$elbo)
  expect_equal(new.res$niter, original.res$niter)
  expect_equal(new.res$intercept, original.res$intercept)
  expect_equal(new.res$fitted, original.res$fitted)
  expect_equal(new.res$X_column_scale_factors, original.res$X_column_scale_factors)
}

expect_equal_susie_ss = function(new.res, original.res){
  expect_equal_susie_ss_update(new.res, original.res)
  expect_equal(new.res$elbo, original.res$elbo)
  expect_equal(new.res$niter, original.res$niter)
  expect_equal(new.res$intercept, original.res$intercept)
  expect_equal(new.res$Xtfitted, original.res$Xtfitted)
}
