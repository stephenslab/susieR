#susie_trendfilter

#' @title perform trend filtering using SuSiE
#' @param y an n vector
#' @param order a scalar for the order of trend filtering
#' @return SuSiE fit for trend filtering
#' @examples
#' set.seed(1)
#' mu = c(rep(0,5),rep(1,5),rep(3,5),rep(-2,5),rep(0,30))
#' y = mu + rnorm(50)
#' s = susie_trendfilter(y, 0)
#' plot(y,pch=".")
#' lines(mu,col=1,lwd=3)
#' lines(predict(s),col=2,lwd=2)
#' s0 = susie_trendfilter(y, 0, estimate_prior_variance = TRUE)
#' s1 = susie_trendfilter(y, 1, L=20)
#' @export
susie_trendfilter = function(y, order, normalize=TRUE,...){
  n = length(y)
  X <- Matrix::sparseMatrix(i=NULL,j=NULL,dims=c(n,n))
  attr(X, "matrix.type") = "tfmatrix"
  attr(X, "order") = order
  if (normalize) {
    attr(X, "d") <- compute_tf_std_d(order, n)
    attr(X, "scaled:center") <- compute_tf_cm(order, n)
    attr(X, "scaled:scale") <- compute_tf_csd(order, n)
    s = susie(X=X, Y=y, ...)
  } else {
    attr(X, "d") <- compute_tf_d(order, n)
    attr(X, "scaled:center") <- rep(0,n)
    attr(X, "scaled:scale") <- rep(1,n)
    s = susie(X=X, Y=y, standardize=FALSE, intercept=FALSE,...)
  }
  
  return(s)
}
