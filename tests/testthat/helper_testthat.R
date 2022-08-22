# @title sets three attributes for matrix X
# @param X an n by p data matrix that can be either a trend filtering
#   matrix or a regular dense/sparse matrix
# @param center boolean indicating centered by column means or not
# @param scale boolean indicating scaled by column standard deviations or not
# @return X with three attributes e.g. `attr(X, 'scaled:center') is a
# p vector of column means of X if center=TRUE, a p vector of zeros
# otherwise. 'attr(X, 'scaled:scale') is a p vector of columan standard
# deviations of X if scale=TRUE, a p vector of 1s otherwise. 'attr(X,
# 'd') is a p vector of column sums of X.standardized^2,' where
# X.standardized is the matrix X centered by attr(X, 'scaled:center')
# and scaled by attr(X, 'scaled:scale').
#
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
set_X_attributes = function (X, center = TRUE, scale = TRUE) {
    
  # if X is a trend filtering matrix
  if (!is.null(attr(X,"matrix.type"))) {
    order = attr(X,"order")
    n = ncol(X)
    
    # Set three attributes for X.
    attr(X,"scaled:center") = compute_tf_cm(order,n)
    attr(X,"scaled:scale") = compute_tf_csd(order,n)
    attr(X,"d") = compute_tf_d(order,n,attr(X,"scaled:center"),
                               attr(X,"scaled:scale"),scale,center)
    if (!center)
      attr(X,"scaled:center") = rep(0,n)
    if (!scale)
      attr(X,"scaled:scale") = rep(1,n)
  } else {
      
    # If X is either a dense or sparse ordinary matrix.
    # Get column means.
    cm = colMeans(X,na.rm = TRUE)
    
    # Get column standard deviations.
    csd = compute_colSds(X)
    
    # Set sd = 1 when the column has variance 0.
    csd[csd == 0] = 1
    if (!center)
      cm = rep(0,length = length(cm))
    if (!scale) 
      csd = rep(1,length = length(cm))

    # Ah, this code is very inefficient because the matrix becomes
    # dense!
    # 
    #   X.std = as.matrix(X)
    #   X.std = (t(X.std) - cm)/csd
    #   attr(X,"d") = rowSums(X.std * X.std)
    #
    # Set three attributes for X.
    n = nrow(X)
    d = n*colMeans(X)^2 + (n-1)*compute_colSds(X)^2
    d = (d - n*cm^2)/csd^2
    attr(X,"d") = d
    attr(X,"scaled:center") = cm
    attr(X,"scaled:scale") = csd
  }
  return(X)
}

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
    X.sparse = as(X,"CsparseMatrix")
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
      V=scaled_prior_variance * as.numeric(var(y)),
      lbf_variable = matrix(0,L,p))
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

expect_equal_susie_update = function(new.res, original.res, tolerance = .Machine$double.eps^0.5){
  expect_equal(new.res$alpha, original.res$alpha, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu, original.res$mu, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu2, original.res$mu2, scale = 1, tolerance = tolerance)
  expect_equal(new.res$Xr, original.res$Xr, scale = 1, tolerance = tolerance)
  expect_equal(new.res$KL, original.res$KL, scale = 1, tolerance = tolerance)
  expect_equal(new.res$sigma2, original.res$sigma2, scale = 1, tolerance = tolerance)
  expect_equal(new.res$V, original.res$V, scale = 1, tolerance = tolerance)
}

expect_equal_susie_suff_stat_update = function(new.res, original.res, tolerance = .Machine$double.eps^0.5){
  expect_equal(new.res$alpha, original.res$alpha, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu, original.res$mu, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu2, original.res$mu2, scale = 1, tolerance = tolerance)
  expect_equal(new.res$XtXr, original.res$XtXr, scale = 1, tolerance = tolerance)
  expect_equal(new.res$KL, original.res$KL, scale = 1, tolerance = tolerance)
  expect_equal(new.res$sigma2, original.res$sigma2, scale = 1, tolerance = tolerance)
  expect_equal(new.res$V, original.res$V, scale = 1, tolerance = tolerance)
}

expect_equal_susie_rss_update = function(new.res, original.res, tolerance = .Machine$double.eps^0.5){
  expect_equal(new.res$alpha, original.res$alpha, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu, original.res$mu, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu2, original.res$mu2, scale = 1, tolerance = tolerance)
  expect_equal(new.res$Rz, original.res$Rz, scale = 1, tolerance = tolerance)
  expect_equal(new.res$KL, original.res$KL, scale = 1, tolerance = tolerance)
  expect_equal(new.res$sigma2, original.res$sigma2, scale = 1, tolerance = tolerance)
  expect_equal(new.res$V, original.res$V, scale = 1, tolerance = tolerance)
}

expect_equal_SER = function(new.res, original.res){
  expect_equal(new.res$alpha, original.res$alpha)
  expect_equal(new.res$mu, original.res$mu)
  expect_equal(new.res$mu2, original.res$mu2)
  expect_equal(new.res$lbf, original.res$lbf)
  expect_equal(new.res$V, original.res$V)
  expect_equal(new.res$loglik, original.res$loglik)
}

expect_equal_SER_suff_stat = function(new.res, original.res, tolerance = .Machine$double.eps^0.5){
  expect_equal(new.res$alpha, original.res$alpha, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu, original.res$mu, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu2, original.res$mu2, scale = 1, tolerance = tolerance)
  expect_equal(new.res$lbf, original.res$lbf, scale = 1, tolerance = tolerance)
  expect_equal(new.res$V, original.res$V, scale = 1, tolerance = tolerance)
  expect_equal(new.res$lbf_model, original.res$lbf_model, scale = 1, tolerance = tolerance)
}

expect_equal_susie = function(new.res, original.res, tolerance = .Machine$double.eps^0.5){
  expect_equal_susie_update(new.res, original.res, tolerance = tolerance)
  expect_equal(new.res$elbo, original.res$elbo, scale = 1, tolerance = tolerance)
  expect_equal(new.res$niter, original.res$niter, scale = 1, tolerance = tolerance)
  expect_equal(new.res$intercept, original.res$intercept, scale = 1, tolerance = tolerance)
  expect_equal(new.res$fitted, original.res$fitted, scale = 1, tolerance = tolerance)
  expect_equal(new.res$X_column_scale_factors, original.res$X_column_scale_factors, scale = 1, tolerance = tolerance)
}

expect_equal_susie_suff_stat = function(new.res, original.res, tolerance = .Machine$double.eps^0.5){
  expect_equal_susie_suff_stat_update(new.res, original.res, tolerance = tolerance)
  expect_equal(new.res$elbo, original.res$elbo, scale = 1, tolerance = tolerance)
  expect_equal(new.res$niter, original.res$niter, scale = 1, tolerance = tolerance)
  expect_equal(new.res$intercept, original.res$intercept, scale = 1, tolerance = tolerance)
  expect_equal(new.res$Xtfitted, original.res$Xtfitted, scale = 1, tolerance = tolerance)
}

expect_equal_susie_rss = function(new.res, original.res, tolerance = .Machine$double.eps^0.5){
  expect_equal_susie_rss_update(new.res, original.res, scale = 1, tolerance = tolerance)
  expect_equal(new.res$elbo, original.res$elbo, scale = 1, tolerance = tolerance)
  expect_equal(new.res$niter, original.res$niter, scale = 1, tolerance = tolerance)
  expect_equal(new.res$intercept, original.res$intercept, scale = 1, tolerance = tolerance)
  expect_equal(new.res$Rz, original.res$Rz, scale = 1, tolerance = tolerance)
}
