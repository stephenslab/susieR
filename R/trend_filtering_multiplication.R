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

#' @title Compute colSums(Dinv * Dinv)
#' @param order is the order of trend filtering
#' @param y an n vector
#' @return an n vector
compute_colSumsDinv2 = function(order, y){
  n = length(y)
  base = rep(-1, n)
  if (order==0) return(cumsum(base^2))
  for (i in 1:order){
    base = cumsum(base)
  }
  return(cumsum(base^2))
}

#' @title Compute sum(t(Dinv * Dinv) * Eb2)
#' @param order is the order of trend filtering
#' @param y an n vector
#' @param Eb2 an n vector
#' @return a scalar
compute_Dinv2tEb2 = function(order, y, Eb2){
  n = length(y)
  base = rep(-1, n)
  if (order==0) return(sum(cumsum(base^2)*Eb2))
  for (i in 1:order){
    base = cumsum(base)
  }
  return(sum(cumsum(base^2)*Eb2))
}

compute_tfcm = function(order, y){
  n = length(y)
  base = rep(1,n)
  for (i in 1:(order+1)){
    base = -cumsum(base)
  }
  return(base/n)
}

compute_tfcsd = functionn(order, y){
  
}




