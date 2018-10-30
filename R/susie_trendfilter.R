#susie_trendfilter

#' @title perform trend filtering using SuSiE
#' @param Y an n vector
#' @param order a scalar for the order of trend filtering
#' @return SuSiE fit for trend filtering
#' @examples
#' y = c(rep(0,5),rep(1,5),rep(3,5),rep(-2,5),rep(0,30)) + rnorm(50)
#' susie_trendfilter(y, 0)
#' susie_trendfilter(y, 1, L=20)
#' susie_trenndfilter(y, 0, estimate_prior_variance = TRUE)
#' @export
susie_trendfilter = function(Y, order, standardize=TRUE,...){
  n = length(Y)
  X = matrix(1,n,n)
  class(X) = "tfmatrix"
  attr(X, "order") = order
  if (standardize) {
    attr(X, "d") <- compute_tf_std_d(order, n)
    attr(X, "scaled:center") <- compute_tf_cm(order, n)
    attr(X, "scaled:scale") <- compute_tf_csd(order, n)
  } else {
    attr(X, "d") <- compute_tf_d(order, n)
    attr(X, "scaled:center") <- rep(0,n)
    attr(X, "scaled:scale") <- rep(1,n)
  }
  s = susie(X=X, Y=Y, ...)
  return(s)
}
