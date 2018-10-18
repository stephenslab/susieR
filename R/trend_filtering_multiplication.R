#trend filtering multiplication helper functions

#' @title Create D-inverse matrix as input X
#' @param order is the order of trend filtering
#' @param n the length of y
#' @return a D-inverse matrix
#' @importFrom expm %^%
#' @keywords internal
create_Dinv = function(order, n){
  #form a basis D
  D = diag(-1, n)
  for (i in 1:(n-1)){
    D[i, i+1] = 1
  }
  #form D^(k+1)
  # %^% uses the package 'expm'
  D = D %^% (order+1)
  #return inverse of D^(k+1)
  return(solve(D))
}

#' @title Compute unscaled D-inverse \%*\% b using special structure of D-inverse
#' @param order is the order of trend filtering
#' @param b an n=p vector
#' @return an n vector
#' @keywords internal
compute_Dinvb = function(order,b){
  for (i in 1:(order+1)){
    b = rev(-1*cumsum(rev(b)))
  }
  return(b)
}

#' @title Compute unscaled t(D-inverse) \%*\% y using special structure of D-inverse
#' @param order is the order of trend filtering
#' @param y an n vector
#' @return an n vector
#' @keywords internal
compute_Dinvty = function(order,y){
  for (i in 1:(order+1)){
    y = -1*cumsum(y)
  }
  return(y)
}




