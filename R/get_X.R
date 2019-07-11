get_csd = function(X){
  return(attr(X,'scaled:scale'))
}

get_cm = function(X){
  return(attr(X,'scaled:center'))
}

get_d = function(X){
  return(attr(X,'d'))
}

get_order = function(X){
  return(attr(X,'order'))
}

get_nrow = function(X){
  return(nrow(X))
}

get_ncol = function(X){
  if(is.stumps_matrix(X)){return(nrow(X)*ncol(X))} # might  want to improve this special case?
  return(ncol(X))
}
