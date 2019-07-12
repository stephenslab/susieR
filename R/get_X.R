get_csd = function(X){
  if(is.list(X)){return(unlist(lapply(X,get_csd)))}
  return(attr(X,'scaled:scale'))
}

get_cm = function(X){
  if(is.list(X)){return(unlist(lapply(X,get_cm)))}
  return(attr(X,'scaled:center'))
}

get_d = function(X){
  if(is.list(X)){return(unlist(lapply(X,get_d)))}
  return(attr(X,'d'))
}

get_order = function(X){
  return(attr(X,'order'))
}

get_nrow = function(X){
  if(is.list(X)){return(get_nrow(X[[1]]))}
  if(is_valid_matrix(X)){return(nrow(X))}
  return(attr(X,"nrow"))
}

get_ncol = function(X){
  if(is.stumps_matrix(X)){return(nrow(X)*ncol(X))} # might  want to improve this special case?
  if(is.list(X)){return(Reduce('+',lapply(X,get_ncol)))}
  if(is_valid_matrix(X)){return(ncol(X))}
  return(attr(X,"ncol"))
}


is_valid_matrix = function(X){
  return((is.double(X) & is.matrix(X)) | inherits(X,"CsparseMatrix"))
}
