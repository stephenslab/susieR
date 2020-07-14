#' @title Apply susie to trend filtering (especially changepoint
#'   problems), a type of non-parametric regression.
#'
#' @description Fits the non-parametric Gaussian regression model
#'   \eqn{y = mu + e}, where the mean \eqn{mu} is modelled as \eqn{mu =
#'   Xb}, X is a matrix with columns containing an appropriate basis,
#'   and b is vector with a (sparse) SuSiE prior. In particular, when
#'   \code{order = 0}, the jth column of X is a vector with the first j
#'   elements equal to zero, and the remaining elements equal to 1, so
#'   that \eqn{b_j} corresponds to the change in the mean of y between
#'   indices j and j+1. For background on trend filtering, see
#'   Tibshirani (2014). See also the "Trend filtering" vignette,
#'   \code{vignette("trend_filtering")}.
#'
#' @details The implementation exploits the special structure of X,
#'   which means that the matrix-vector product \eqn{X^Ty} is fast to
#'   compute in \eqn{O(n)} computation rather than \eqn{O(n^2)} if X
#'   were formed explicitly. For implementation details, view the
#'   "Implementation of SuSiE trend filtering" vignette by running
#'   \code{vignette("trendfiltering_derivations")}.
#'
#' @param y An n-vector of observations ordered in time or space
#'   (assumed to be equally spaced).
#'
#' @param order An integer specifying the order of trend filtering.
#'   The default, \code{order = 0}, corresponds to "changepoint"
#'   problems (\emph{i.e.}, piecewise constant \eqn{mu}). Although
#'   \code{order > 0} is implemented, we do not recommend its use; in
#'   practice, we have found problems with convergence of the algorithm
#'   to poor local optima, producing unreliable inferences.
#' 
#' @param standardize Logical indicating whether to standardize the X
#'   variables ("basis functions"); defaults to \code{standardize = FALSE}
#'   as these basis functions already have a natural scale.
#' 
#' @param use_mad Logical indicating whether to use the "median
#'   absolute deviation" (MAD) method to the estimate residual
#'   variance. If \code{use_mad = TRUE}, susie is run twice, first by
#'   fixing the residual variance to the MAD value, then a second time,
#'   initialized to the first fit, but with residual variance estimated
#'   the usual way (maximizing the ELBO). We have found this strategy
#'   typically improves reliability of the results by reducing a
#'   tendency to converge to poor local optima of the ELBO.
#' 
#' @param ... Other arguments passed to \code{\link{susie}}.
#' 
#' @return A "susie" fit; see \code{\link{susie}} for details.
#'
#' @references R. J. Tibshirani (2014). Adaptive piecewise polynomial
#' estimation via trend filtering. \emph{Annals of Statistics}
#' \bold{42}, 285-323. \url{https://doi.org/10.1214/13-AOS1189}
#' 
#' @examples
#' set.seed(1)
#' mu = c(rep(0,50),rep(1,50),rep(3,50),rep(-2,50),rep(0,300))
#' y = mu + rnorm(500)
#' s = susie_trendfilter(y)
#' plot(y)
#' lines(mu,col=1,lwd=3)
#' lines(predict(s),col=2,lwd=2)
#'
#' # Returns credible sets (indices of y that occur just before
#' # changepoints).
#' susie_get_cs(s)
#'
#' # Produces plot with credible sets for changepoints.
#' susie_plot_changepoint(s,y) 
#'
#' @importFrom Matrix sparseMatrix
#' 
#' @export
#' 
susie_trendfilter = function (y, order = 0, standardize = FALSE,
                              use_mad = TRUE, ...) {
  if (order > 0)
    warning("order > 0 is not recommended")
  n = length(y)
  X = sparseMatrix(i = NULL,j = NULL,dims = c(n,n))
  attr(X,"matrix.type") = "tfmatrix"
  attr(X,"order") = order
  if (use_mad && !("s_init" %in% names(list(...)))) {
    mad = estimate_mad_residual_variance(y)
    s_mad_init = suppressWarnings(susie(X = X,Y = y,standardize = standardize,
      estimate_residual_variance = FALSE,residual_variance = mad,...))
    s = susie(X = X,Y = y,standardize = standardize, s_init = s_mad_init,...)
  } else
    s = susie(X = X,Y = y,standardize = standardize,...)
  return(s)
}

# @title estimate residual variance using MAD estimator
# @param y an n vector
# @return a scalar of estimated residual variance
estimate_mad_residual_variance = function (y)
  0.5*(median(abs(diff(y))/0.6745)^2)
