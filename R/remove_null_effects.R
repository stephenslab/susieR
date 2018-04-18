remove_null_effects = function(s){
  null_indices = (s$sa2==0)
  s$alpha = s$alpha[!null_indices,,drop=FALSE]
  s$mu = s$mu[!null_indices,,drop=FALSE]
  s$mu2 = s$mu2[!null_indices,,drop=FALSE]
  s$sa2 = s$sa2[!null_indices,drop=FALSE]
  return(s)
}

add_null_effect = function(s,sa2=1){
  p = ncol(s$alpha)
  s$alpha = rbind(s$alpha,1/p)
  s$mu = rbind(s$mu,rep(0,p))
  s$mu2 = rbind(s$mu2,rep(0,p))
  s$sa2 = c(s$sa2,sa2)
  return(s)
}
