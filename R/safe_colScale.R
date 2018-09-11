# Credit: This is heavily based on code from
# https://www.r-bloggers.com/a-faster-scale-function/
# The only change from that code is its treatment of columns with 0 variance.
# This "safe" version scales those columns by 1 instead of 0.
# @param X can be sparse or dense
# if X is dense return a scaled dense x with associated attributes: scale:center, scale:scale, and scaled.X
# if X is sparse return sparse itself with associated attributes: scale:scale, scale:center, and scaled.X
safe_colScale = function(X,
                    center = TRUE,
                    scale = TRUE,
                    add_attr = TRUE,
                    rows = NULL,
                    cols = NULL) {
  
  if (!is.null(rows) && !is.null(cols)) {
    X <- X[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    X <- X[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    X <- X[, cols, drop = FALSE]
  }
  
  if (!is.null(attr(X, 'scaled:center'))){
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
    csd[csd==0] = 1
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }
  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  X.dense = t( (t(X.dense) - cm) / csd )
  
  if (add_attr) {
    attr(X, "scaled.X") <- X.dense
    attr(X, "d") <- Matrix::colSums(compute_X2(X))
    if (center) attr(X, "scaled:center") <- cm
    if (scale) attr(X, "scaled:scale") <- csd
  }
  return(X)
}
