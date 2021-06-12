remove_null_effects = function (s) {
  null_indices = (s$V == 0)
  s$alpha = s$alpha[!null_indices,,drop = FALSE]
  s$mu    = s$mu[!null_indices,,drop = FALSE]
  s$mu2   = s$mu2[!null_indices,,drop = FALSE]
  s$lbf_variable = s$lbf_variable[!null_indices,,drop = FALSE]
  s$V     = s$V[!null_indices,drop = FALSE]
  return(s)
}

add_null_effect = function (s, V) {
  p = ncol(s$alpha)
  s$alpha = rbind(s$alpha,1/p)
  s$mu    = rbind(s$mu,rep(0,p))
  s$mu2   = rbind(s$mu2,rep(0,p))
  s$lbf_variable = rbind(s$lbf_variable,rep(0,p))
  s$V     = c(s$V,V)
  return(s)
}
