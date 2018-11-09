# sets the attributes of the X matrix
# @param X can be sparse or dense
# if X is dense return a scaled dense x with associated attributes: scale:center, scale:scale, and scaled.X
# if X is sparse return sparse itself with associated attributes: scale:scale, scale:center, and scaled.X
set_X_attributes = function(X,
                         center = TRUE,
                         scale = TRUE) {
  # if a trendfiltering matrix
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
  } else { # for matrix X

    if (!is.null(attr(X, 'scaled.X'))) {
      X.dense = attr(X, 'scaled.X')
    } else {
      X.dense = as.matrix(X)
    }
    ################
    # Get the column means
    ################
    cm = colMeans(X.dense, na.rm = TRUE)
    ################
    # Get the column sd
    ################
    if (scale) {
      csd = matrixStats::colSds(X.dense, center = cm)
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
      X = X.dense
    }

    attr(X, "d") <- Matrix::colSums(X.dense * X.dense)
    attr(X, "X2t") <- t(X.dense * X.dense)
    attr(X, "scaled:center") <- cm
    attr(X, "scaled:scale") <- csd

  }


  return(X)
}
