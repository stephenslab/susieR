#' @title Compute sufficient statistics for input to \code{susie_suff_stat}
#'
#' @description Computes the sufficient statistics \eqn{X'X, X'y, y'y}
#'   and \eqn{n} after centering (and possibly standardizing) the
#'   columns of \eqn{X} and centering \eqn{y} to have mean zero. We also
#'   store the column means of \eqn{X} and mean of \eqn{y}.
#'
#' @param X An n by p matrix of covariates.
#'
#' @param y An n vector.
#'
#' @param standardize Logical flag indicating whether to standardize
#'   columns of X to unit variance prior to computing summary data
#'
#' @return A list of sufficient statistics (\code{XtX, Xty, yty, n})
#'   and \code{X_colmeans}, \code{y_mean}.
#'
#' @importFrom methods as
#' @importFrom Matrix colMeans
#' @importFrom Matrix crossprod
#'
#' @examples
#' data(N2finemapping)
#' ss = compute_suff_stat(N2finemapping$X, N2finemapping$Y[,1])
#'
#' @export
#'
compute_suff_stat = function(X, y, standardize = FALSE) {
  y_mean = mean(y)
  y   = y - y_mean
  n   = nrow(X)
  mu  = colMeans(X)
  s   = compute_colSds(X)
  Xty = drop(y %*% X)
  XtX = crossprod(X)
  XtX = as.matrix(XtX)
  XtX = XtX - n*tcrossprod(mu)
  if (standardize) {
    XtX = XtX/s
    XtX = t(XtX)
    XtX = XtX/s
    Xty = Xty/s
  }
  n   = length(y)
  yty = sum(y^2)
  return(list(XtX = XtX,Xty = Xty,yty = yty,n = n,
              y_mean = y_mean,X_colmeans = mu))
}

#' @title Compute sufficient statistics for input to \code{susie_suff_stat}
#'
#' @description This is a synonym for \code{compute_suff_stat}
#'   included for historical reasons (deprecated).
#'
#' @param X An n by p matrix of covariates.
#'
#' @param y An n vector.
#'
#' @param standardize Logical flag indicating whether to standardize
#'   columns of X to unit variance prior to computing summary data.
#'
#' @return A list of sufficient statistics (X'X, X'y, y'y and n)
#'
#' @importFrom methods as
#'
#' @examples
#' data(N2finemapping)
#' ss = compute_ss(N2finemapping$X, N2finemapping$Y[,1])
#'
#' @export
#'
compute_ss = compute_suff_stat
