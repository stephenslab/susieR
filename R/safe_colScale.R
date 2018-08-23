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
  x = as.matrix(X)
  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }

  ################
  # Get the column means
  ################
  cm = colMeans(x, na.rm = TRUE)
  ################
  # Get the column sd
  ################
  if (scale) {
    csd = matrixStats::colSds(x, center = cm)
    csd[csd==0] = 1
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }
  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  x = t( (t(x) - cm) / csd )
  if (is.matrix(X)){
    X = x
  }
  if (add_attr) {
    if (center) {
      attr(X, "scaled:center") <- cm
    }
    if (scale) {
      attr(X, "scaled:scale") <- csd
    }
    if (class(X)=='dgCMatrix') {
      attr(X, "scaled.X") <- x
    }
  }
  return(X)
}
