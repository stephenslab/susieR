#' @title Compute unscaled Xexpand \%*\% b where Xexpand is the n by np matrix formed by, for each of the p columns of Xrank form the n different stumps (changepoint basis functions), and concatentate
#' @param Xord is an n by p matrix whose jth column containing order(X[,j])
#' @param b an np vector whose first n elements correspond to the first column of Xord, next n correspond to second column, etc.
#' @return an n vector
#' @keywords internal
compute_stumps_Xb = function(Xord,b){
  n = nrow(Xord)
  p = ncol(Xord)
  bmat = matrix(b, nrow=n) # convenient way to split b into p lots of n
  res = rep(0,n)
  for(j in 1:p){ # note could parallelize this
    res[Xord[,j]] = res[Xord[,j]] + compute_tf_Xb(order=0, b=bmat[,j])
  }
  return(res)
}

#' @title Compute unscaled t(Xexpand) \%*\% y here Xexpand is the n by np matrix formed by, for each of the p columns of Xrank form the n different stumps (changepoint basis functions), and concatentate
#' @param Xord is an n by p matrix whose jth column containing order(X[,j])

#' @param y an n vector
#' @return an np vector
#' @keywords internal
compute_stumps_Xty = function(Xord,y){
  res = apply(Xord,2,function(o){compute_tf_Xty(order=0, y=y[o])})
  as.vector(res)
}
