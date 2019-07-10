# aim of these functions is to implement multiplication
# for "trend filtering general" (tfg)
# where the observed data points may be unequally spaced
# and the locations of changepoints on the x axis can be
# arbitrarily specified. (previous versions will correspond to the
# special case where the )

#' set up a general trend filtering matrix
#' @param t vector of length n specifying locations of data points on x axis
#' @param br vector of length (p-1) specifying break points on x axis (ie where changepoints can occur)
#' Elements of br must increase monotonically.
#' @param order non-negative integer indicating order of trend filtering basis (0 is changepoint basis and is the only case we test and use)
#' @keywords internal
make_tfg_matrix = function(t,br,order=0){
  n = length(t) # number of data points
  p = length(br) + 1 # number of bins specified by breaks
  X <- Matrix::sparseMatrix(i=NULL,j=NULL,dims=c(n,p)) # this is set so that ncol(X) and  nrow(X)  works
  attr(X, "matrix.type") = "tfg_matrix"
  attr(X, "order") = order
  attr(X,"t") <- t
  attr(X,"br") <- br
  attr(X,"order_t") <- order(t)
  attr(X,"t_to_bin") <- .bincode(t,breaks = c(-Inf,br,Inf))
  attr(X,"bin_to_t") <- cumsum(hist(t, breaks = c(-Inf,br,Inf), plot=FALSE)$counts)
  return(X)
}

#' @title Compute unscaled X \%*\% b using the special structure of trend filtering
#' @param X a tfg_matrix created by make_tfg_matrix
#' @param b a p vector of the changes at each change point
#' @return an n vector of the means at each data point
#' @keywords internal
compute_tfg_Xb = function(X,b){
  order = attr(X,"order")
  for(i in 1:(order+1)){
    b = rev(-1*cumsum(rev(b))) # computes mean in each bin
  }
  return(b[attr(X,"t_to_bin")]) #  maps bin means to a mean for each datapoint
}

#' @title Compute t(X) \%*\% y using the special structure of trend filtering
#' @param X a tfg_matrix created by make_tfg_matrix
#' @param y an n vector of data
#' @return a p vector
#' @keywords internal
compute_tfg_Xty = function(X,y){
  order = attr(X,"order")
  y = y[attr(X,"order_t")] # sort y according to increasing t
  for (i in 1:(order+1)){
    y = -1*cumsum(y)
  }
  return(y[attr(X,"bin_to_t")])
}

is.tfg_matrix=function(X){
  ifelse(is.null(attr(X, "matrix.type")),FALSE,attr(X,"matrix.type")=="tfg_matrix")
}

#' @title Compute column mean of the general trend filtering matrix X
#' @param X a general trend filtering matrix
#' @return a p vector of column means
#' @keywords internal
compute_tfg_cm = function(X){
  order = attr(X,"order")
  base = hist(attr(X,"t"), breaks = c(-Inf,attr(X,"br"),Inf), plot=FALSE)$counts
  for (i in 1:(order+1)){
    base = -cumsum(base)
  }
  return(base/nrow(X))
}

#' @title Compute d=diag(Xcent'Xcent) where Xcent is the centered version of the general trend filtering matrix X
#' @param X a general trend filtering matrix
#' @return a p vector d=diag(Xcent'Xcent)
#' @keywords internal
compute_tfg_d = function(X){
  order = attr(X,"order")
  if(order!=0){
    stop("not implemented for order!=0")
  }
  cm = compute_tfg_cm(X) # column means
  return(nrow(X)*(- cm - cm^2))
}

