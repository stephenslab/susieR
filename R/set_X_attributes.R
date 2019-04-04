#' @title sets the attributes for the X matrix
#' @param X an n by p data matrix that can be a dense, sparse, or trend filtering matrix
#' @param center boolean indicating mean centering or not
#' @param scale boolean indicating scaled by standard deviation or not
#' @return X with three attributes e.g.
#'         attr(X, 'scaled:center') is a p vector of column means of X
#'         attr(X, 'scaled:scale') is a p vector of column standard deviations of X
#'         attr(X, 'd') is a p vector of column sums of X^2
#' if X is a trend filtering matrix, return the initial X with three attributes;
#' if X is a dense matrix, return the scaled X with three attributes;
#' if X is a sparse matrix, return the initial X with three attributes.

set_X_attributes = function(X,
                            center = TRUE,
                            scale = TRUE) {
  # if X is a trend filtering matrix
  if (!is.null(attr(X, "matrix.type"))) {
    order <- attr(X,"order")
    n <- ncol(X)
    attr(X, "scaled:center") <- compute_tf_cm(order, n)
    attr(X, "scaled:scale") <- compute_tf_csd(order, n)
    attr(X, "d") <-
      compute_tf_d(
        order,
        n,
        attr(X, "scaled:center"),
        attr(X, "scaled:scale"),
        scale,
        center
      )
    if (!center) {
      attr(X, "scaled:center") <- rep(0, n)
    }
    if (!scale) {
      attr(X, "scaled:scale") <- rep(1, n)
    }
  } else {
    # if X is a dense or sparse matrix
    X.dense = as.matrix(X)
    # get column means
    cm = colMeans(X.dense, na.rm = TRUE)
    # get column standard deviations
    if (scale) {
      csd = matrixStats::colSds(X.dense)
      # set sd = 1 when the column has variance 0
      csd[csd == 0] = 1
    } else {
      # just divide by 1 if not
      csd = rep(1, length = length(cm))
    }
    if (!center) {
      # just subtract 0
      cm = rep(0, length = length(cm))
    }
    X.dense = t((t(X.dense) - cm) / csd)
    if (is.matrix(X)) {
      # if X is dense, change X into centered and scaled version as required by inputs
      X = X.dense
    }
    # set attributes for either dense or sparse X
    attr(X, "d") <- Matrix::colSums(X.dense * X.dense)
    attr(X, "scaled:center") <- cm
    attr(X, "scaled:scale") <- csd
  }
  return(X)
}
