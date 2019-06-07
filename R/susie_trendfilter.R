#susie_trendfilter

#' @title Applies susie to perform trend filtering (especially
#'   changepoint problems), a type of non-parametric regression.
#' 
#' @details Fits the non-parametric Gaussian regression model $y=mu
#' +e$, where the mean $mu$ is modelled as $mu=Xb$ where $X$ is a
#' matrix with columns containing an appropriate basis and b is vector
#' with a (sparse) SuSiE prior. In particular, when order=0, the $j$th
#' column of $X$ is a vector with first $j$ elements equal to 0 and
#' remaining elements equal to 1, so b_j corresponds to the change in
#' the mean of y between indices $j$ and $j+1$.  For background on
#' trend filtering see [Adaptive piecewise polynomial estimation via
#' trend
#' filtering](https://projecteuclid.org/euclid.aos/1395234979)(Ryan
#' J. Tibshirani, 2014) The implementation here exploits the special
#' structure of $X$ whiuch means that the matrix-vector product $X'y$
#' is fast to compute in O(n) computation rather than $O(n^2)$ if $X$
#' were formed explicitly. For implementation details, view the
#' "trendfiltering derivations" vignette by running
#' \code{vignette("trendfiltering_derivations")}.
#' 
#' @param y an n vector of observations that are ordered in time or
#'   space (assumed equally-spaced)
#' 
#' @param order an integer specifying the order of trend filtering. Default order=0, which corresponds
#' to "changepoint" problems (i.e. piecewise constant mu). Although order > 0 is implemented, we do not recommend using it since we
#' find there are often problems with convergence of the algorithm to poor local optima, producing unreliable inferences.
#' @param standardize boolean indicating whether to standardize X variables (basis functions); defaults to FALSE as these basis functions already have a natural scale
#' @param use_mad boolean indicating whether to use ``median absolute deviation" (MAD) method to estimate residual variance.
#' If TRUE then susie is run twice,
#' first by fixing the residual variance to the MAD value and then a second time, initializing from the first fit,
#' but with residual variance estimated the usual way (maximizing the ELBO). We have found this strategy
#' typically improves reliability of results by reducing a tendency to converge to poor local optima of the ELBO.
#' @param ... other parameters to pass to susie
#' @return a "susie" object; see ?susie for details
#' @examples
#' set.seed(1)
#' mu = c(rep(0,50),rep(1,50),rep(3,50),rep(-2,50),rep(0,300))
#' y = mu + rnorm(500)
#' s = susie_trendfilter(y)
#' plot(y)
#' lines(mu,col=1,lwd=3)
#' lines(predict(s),col=2,lwd=2)
#' susie_get_cs(s) # returns credible sets (for indices of y that occur just before changepoints)
#' susie_plot_changepoint(s,y) # produces ggplot with credible sets for changepoints on top of plot
#' @export
susie_trendfilter = function(y, order=0,standardize=FALSE, use_mad=TRUE,...){
  if (order > 0){
    warning("order>0 is not recommended (see ?susie_trendfilter for more explanation).")
  }
  n = length(y)
  X <- Matrix::sparseMatrix(i=NULL,j=NULL,dims=c(n,n))
  attr(X, "matrix.type") = "tfmatrix"
  attr(X, "order") = order
  if (use_mad){
    mad = estimate_mad_residual_variance(y)
    s_mad_init = suppressWarnings(susie(X=X, Y=y, standardize = standardize, estimate_residual_variance = FALSE, residual_variance = mad, ...))
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
