#this is for scaled prior in which effect prior variance is sa * sigma
# s is a susie fit
# a list with elements alpha, mu, sigma2, sa2
elbo = function(X,Y,s){
  L = nrow(s$alpha)
  n = nrow(X)
  p = ncol(X)
  d = colSums(X*X) #note this is currently being computed 2 times - here and in Eloglik; could be made more efficient by avoiding this

  ss <- s$mu2 - s$mu^2 # posterior variance (conditional on inclusion)
  postb2 = s$alpha * s$mu2 # posterior second moment of b
  Ell = Eloglik(X,Y,s)

  sub = s$alpha>0 # this is to avoid taking 0 * log(0) in next line
  KL1 = sum(s$alpha[sub] * log(s$alpha[sub]/(1/p)))

  KL2 = - 0.5* sum(s$alpha * (1 + log(ss)-log(s$sigma2*s$sa2)))
  + 0.5 * sum(postb2)/(s$sigma2*s$sa2)

  return(Ell - KL1 - KL2)
}

Eloglik = function(X,Y,s){
  result =  -(n/2) * log(2*pi* s$sigma2) - (1/(2*s$sigma2)) * get_ER2(X,Y,s)
  return(result)
}


get_ER2 = function(X,Y,s){
  Xr = (s$alpha*s$mu) %*% t(X)
  Xrsum = colSums(Xr)

  d = colSums(X*X)
  postb2 = s$alpha * s$mu2 #posterior second moment

  return(sum((Y-Xrsum)^2) - sum(Xr^2) + sum(d*t(postb2)))
}
