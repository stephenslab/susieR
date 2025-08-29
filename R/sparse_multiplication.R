# @title Computes standardized.X %*% b using sparse multiplication trick
# @param X an n by p unstandardized matrix with three attributes:
# attr(X,"scaled:center"), attr(X,"scaled:scale") and attr(X,"d")
# @param b a p vector
# @return an n vector
#
#' @importFrom Matrix t
#' @importFrom Matrix tcrossprod
compute_Xb <- function(X, b) {
  cm <- attr(X, "scaled:center")
  csd <- attr(X, "scaled:scale")

  # Scale Xb.
  if (!is.null(attr(X, "matrix.type"))) {

    # When X is a trend filtering matrix.
    scaled.Xb <- compute_tf_Xb(attr(X, "order"), b / csd)
  } else {

    # When X is an ordinary sparse/dense matrix.
    scaled.Xb <- tcrossprod(X, t(b / csd))
  }

  # Center Xb.
  Xb <- scaled.Xb - sum(cm * b / csd)
  return(as.numeric(Xb))
}

# @title Computes t(standardized.X) %*% y using sparse multiplication trick
# @param X an n by p unstandardized matrix with three attributes:
# attr(X,"scaled:center"), attr(X,"scaled:scale") and attr(X,"d")
# @param y an n vector
# @return a p vector
#
#' @importFrom Matrix t
#' @importFrom Matrix crossprod
compute_Xty <- function(X, y) {
  cm <- attr(X, "scaled:center")
  csd <- attr(X, "scaled:scale")
  ytX <- crossprod(y, X)

  # Scale Xty.
  if (!is.null(attr(X, "matrix.type"))) {

    # When X is a trend filtering matrix.
    scaled.Xty <- compute_tf_Xty(attr(X, "order"), y) / csd
  } else {

    # When X is an ordinary sparse/dense matrix.
    scaled.Xty <- t(ytX / csd)
  }

  # Center Xty.
  centered.scaled.Xty <- scaled.Xty - cm / csd * sum(y)
  return(as.numeric(centered.scaled.Xty))
}

# @title Computes t(standardized.X) %*% standardized.X using attributes
# @param X an n by p unstandardized matrix with three attributes:
# attr(X,"scaled:center"), attr(X,"scaled:scale") and attr(X,"d")
# @return a p by p matrix representing (scaled X)'(scaled X)
#
#' @importFrom Matrix crossprod
compute_XtX <- function(X) {
  cm <- attr(X, "scaled:center")
  csd <- attr(X, "scaled:scale")
  n <- nrow(X)
  colsums_X <- colSums(X)

  if (!is.null(attr(X, "matrix.type"))) {
    stop("compute_XtX not yet implemented for trend filtering matrices")
  }

  # Compute raw X'X
  XtX_raw <- crossprod(X)

  # Scale columns and rows by 1/csd
  XtX_scaled <- sweep(sweep(XtX_raw, 1, csd, "/"), 2, csd, "/")

  # Adjust for centering
  XtX_centered_scaled <- XtX_scaled -
    n * outer(cm/csd, cm/csd) -
    outer(cm/csd, (colsums_X - n*cm)/csd) -
    outer((colsums_X - n*cm)/csd, cm/csd)

  return(XtX_centered_scaled)
}

# @title Computes M %* %t(standardized.X) using sparse multiplication trick
# @param M a L by p matrix
# @param X an n by p unstandardized matrix with three attributes:
# attr(X,"scaled:center"), attr(X,"scaled:scale") and attr(X,"d")
# @return a L by n matrix
#
#' @importFrom Matrix t
compute_MXt <- function(M, X) {
  cm <- attr(X, "scaled:center")
  csd <- attr(X, "scaled:scale")

  if (!is.null(attr(X, "matrix.type"))) {

    # When X is a trend filtering matrix.
    return(as.matrix(t(apply(M, 1, function(b) compute_Xb(X, b)))))
  } else {

    # When X is an ordinary sparse/dense matrix.
    return(as.matrix(t(X %*% (t(M) / csd)) - drop(M %*% (cm / csd))))
  }

  # This should be the same as
  #
  #   t(apply(M, 1, function(b) compute_Xb(X, b))))
  #
  # as well as
  #
  #   M %*% (t(X)/csd) - drop(tcrossprod(M,t(cm/csd)))
  #
  # but should be more memory-efficient.
}
