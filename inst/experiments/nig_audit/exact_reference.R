# Exact enumeration of the shared-sigma2 SuSiE model for SMALL p^L only.
# Components may select the same variable; enumerate labeled assignments.
exact_nig = function(X,y,V,alpha0,beta0) {
  X = scale(X,center=TRUE,scale=TRUE)
  y = y-mean(y)
  n = nrow(X)
  p = ncol(X)
  L = length(V)
  stopifnot(p^L <= 10000,alpha0 > 0,beta0 > 0)
  assignments = as.matrix(expand.grid(rep(list(seq_len(p)),L)))
  logs = apply(assignments,1,function(j) {
    Z = sweep(X[,j,drop=FALSE],2,sqrt(V),"*")
    C = diag(n) + tcrossprod(Z)
    U = chol(C)
    whitened = forwardsolve(t(U),y)
    -sum(log(diag(U))) - (alpha0+n-1)/2*log(beta0+sum(whitened^2))
  })
  weight = exp(logs-max(logs))
  weight = weight/sum(weight)
  pip = vapply(seq_len(p),function(j)
    sum(weight[rowSums(assignments == j) > 0]),0.0)
  list(pip=pip,weight=weight,assignments=assignments)
}

run_exact_comparison = function(reps=100) {
  rows = list()
  for (r in seq_len(reps)) {
    set.seed(91000+r)
    n=20
    p=8
    X = matrix(rnorm(n*p),n,p) %*% chol(0.8^abs(outer(1:p,1:p,"-")))
    y = X[,2]-X[,6]+rnorm(n)
    for (L in c(1,2)) {
      ref = exact_nig(X,y,rep(0.5,L),1,1)
      fit = suppressWarnings(susie(X,y,L=L,small=TRUE,alpha0=1,beta0=1,
        scaled_prior_variance=0.5,estimate_prior_variance=FALSE,
        max_iter=1000,tol_small=1e-8,min_abs_corr=0))
      rows[[length(rows)+1L]] = data.frame(replicate=r,L=L,
        max_pip_difference=max(abs(ref$pip-fit$pip)),converged=fit$converged)
    }
  }
  result = do.call(rbind,rows)
  write.csv(result,"inst/experiments/nig_audit/exact_comparison.csv",row.names=FALSE)
  print(aggregate(max_pip_difference ~ L,result,function(x)
    c(mean=mean(x),median=median(x),maximum=max(x))))
  invisible(result)
}
