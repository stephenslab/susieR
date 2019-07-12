#' @title sets three attributes for matrix X
#' @param X an n by p data matrix that can be either a trend filtering matrix or a regular dense/sparse matrix
#' @param center boolean indicating centered by column means or not
#' @param scale boolean indicating scaled by column standard deviations or not
#' @return X with three attributes e.g.
#'         attr(X, 'scaled:center') is a p vector of column means of X if center=TRUE, a p vector of 0s otherwise.
#'         attr(X, 'scaled:scale') is a p vector of column standard deviations of X if scale=TRUE, a p vector of 1s otherwise.
#'         attr(X, 'd') is a p vector of column sums of X.standardized^2,
#'         where X.standardized is the matrix X centered by attr(X, 'scaled:center') and scaled by attr(X, 'scaled:scale').

set_X_attributes = function(X,
                             center = TRUE,
                             scale = TRUE) {
  # if X is a trend filtering matrix
  if(is.tfmatrix(X)){
    order <- get_order(X)
    n <- get_nrow(X)
    # set three attributes for X
    attr(X, "scaled:center") <- compute_tf_cm(order, n)
    attr(X, "scaled:scale") <- compute_tf_csd(order, n)
    attr(X, "d") <- compute_tf_d(order,n,get_cm(X),get_csd(X),scale,center)
    if (!center) {
      attr(X, "scaled:center") <- rep(0, n)
    }
    if (!scale) {
      attr(X, "scaled:scale") <- rep(1, n)
    }
  } else if(is.tfg_matrix(X)){
    if(center!=TRUE){stop("only center=TRUE implemented for tfg matrix")}
    if(scale!=FALSE){stop("only scale=FALSE implemented for tfg matrix")}
    if(attr(X,"order")!=0){stop("only order=0 implemented for tfg matrix")}
    attr(X, "scaled:center") <- compute_tfg_cm(X)
    attr(X, "scaled:scale") <- rep(1,get_ncol(X))
    attr(X, "d") <- compute_tfg_d(X)
  } else if(is.stumps_matrix(X)){
    n <- get_nrow(X)
    p <- get_ncol(X)
    Xord <- apply(X,2,order)
    X = numeric(0)
    attr(X, "matrix.type") = "stumps_matrix"
    attr(X, "Xord") = Xord #  store ordering  information may want to rethink this later?
    attr(X, "nrow")  = n
    attr(X, "ncol") = p
    # set three attributes for X
    cm = compute_tf_cm(order=0, n)
    csd = compute_tf_csd(order=0, n)
    attr(X, "scaled:center") <- rep(cm,p/n)
    attr(X, "scaled:scale") <- rep(csd,p/n)
    attr(X, "d") <- rep(compute_tf_d(order=0,n,cm,csd,scale,center),p/n)
    if (!center) {
      attr(X, "scaled:center") <- rep(0, p)
    }
    if (!scale) {
      attr(X, "scaled:scale") <- rep(1, p)
    }
  } else {
    # if X is either a dense or sparse ordinary matrix
    # get column means
    cm = Matrix::colMeans(X, na.rm = TRUE)
    # get column standard deviations
    csd = compute_colSds(X)
    # set sd = 1 when the column has variance 0
    csd[csd == 0] = 1
    if (!center) {
      cm = rep(0, length = length(cm))
    }
    if (!scale) {
      csd = rep(1, length = length(cm))
    }
    X.std = (t(X) - cm) / csd
    # set three attributes for X
    attr(X, "d") <- Matrix::rowSums(X.std * X.std)
    attr(X, "scaled:center") <- cm
    attr(X, "scaled:scale") <- csd
  }
  return(X)
}


is.tfmatrix=function(X){
  ifelse(is.null(attr(X, "matrix.type")),FALSE,attr(X,"matrix.type")=="tfmatrix")
}

is.stumps_matrix=function(X){
  ifelse(is.null(attr(X, "matrix.type")),FALSE,attr(X,"matrix.type")=="stumps_matrix")
}

#' @title computes column standard deviations for sparse or dense matrix
#'        replace matrixStats::colSds since this function only takes a dense matrix
#' @param X an n by p matrix of any type, e.g. sparse, dense
#' @return a p vector of column standard deviations
compute_colSds = function(X){
  n = get_nrow(X)
  return(sqrt((colSums(X^2)/n - (colSums(X)/n)^2)*(n/(n-1))))
}
