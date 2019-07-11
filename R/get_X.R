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
