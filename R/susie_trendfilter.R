#susie_trendfilter

#' @title perform trend filtering using SuSiE, particularly designed for solving change points problem.
#'        See aslo [Adaptive piecewise polynomial estimation via trend filtering](https://projecteuclid.org/euclid.aos/1395234979)(Ryan J. Tibshirani, 2014)
#'        Implementation details can be found here \code{system.file("inst","mist","susie_trendfilter_imp_detail.pdf",package = "susieR")}
#' @param y an n vector
#' @param order a scalar for the order of trend filtering, defualt order is 0
#'        order > 0 is not recommended due to local convergence problem.
#' @param standardize boolean indicating whether to standardize
#' @param use_mad boolean indicating whether to use MAD method to estimate residual variance in the SuSiE and initialize from this fit
#' @param ... other parameters to pass to susie call
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
susie_trendfilter = function(y, order=0,standardize=FALSE, use_mad=TRUE,...){
  if (order > 0){
    warning("order>0 is not recommended due to the local convergence problem.")
  }
  n = length(y)
  X <- Matrix::sparseMatrix(i=NULL,j=NULL,dims=c(n,n))
  attr(X, "matrix.type") = "tfmatrix"
  attr(X, "order") = order
  if (use_mad){
    mad = estimate_mad_residual_variance(y)
    s_mad_init = susie(X=X, Y=y, standardize = standardize, estimate_residual_variance = FALSE, residual_variance = mad, ...)
    s = susie(X=X, Y=y, standardize=standardize, s_init=s_mad_init, ...)
  } else {
    s = susie(X=X, Y=y, standardize=standardize, ...)
  }
  return(s)
}

#' @title estimate residual variance using MAD estimator
#' @param y an n vector
#' @return a scalar of estimated residual variance
#' @importFrom wavethresh wd
#' @importFrom wavethresh accessD
#' @keywords internal
estimate_mad_residual_variance = function(y){
  n = length(y)
  y_reflect = c(y, rev(y))
  J = floor(log2(2*n))
  y_reflect = y_reflect[1:2^J]
  y_reflect = c(y_reflect, rev(y_reflect))
  ywd = wd(y_reflect, filter.number=1, family="DaubExPhase")
  wc_d = accessD(ywd, level=J-1)
  est_resid = (median(abs(wc_d))/0.6745)^2
  return(est_resid)
}
