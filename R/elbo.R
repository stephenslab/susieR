#this is for scaled prior in which effect prior variance is sa * sigma
# s is a susie fit
# a list with elements alpha, mu, sigma2, sa2
elbo = function(X,Y,s){
  L = nrow(s$alpha)
  n = nrow(X)
  p = ncol(X)

  Xr = (s$alpha*s$mu) %*% t(X)
  Xrsum = colSums(Xr)

  d = colSums(X*X)
  ss <- s$sa2*s$sigma2/(s$sa2*d + 1)

  postb2 = s$alpha * t(t(s$mu^2) + ss)

  Eloglik =  -(n/2) * log(2*pi* s$sigma2) -
    (1/(2*s$sigma2)) * sum(Y^2) +
    (1/s$sigma2) * sum(Y * Xrsum) -
    (1/(2*s$sigma2)) * sum(Xrsum^2) +
    (1/(2*s$sigma2)) * sum((Xr^2)) -
    (1/(2*s$sigma2)) * sum(d*t(postb2))

  sub = s$alpha>0 # this is to avoid taking 0 * log(0) in next line
  KL1 = sum(s$alpha[sub] * log(s$alpha[sub]/(1/p)))

  KL2 = - 0.5* sum(t(s$alpha) * (1 + log(ss)-log(s$sigma2*s$sa2)))
  + 0.5 * sum(postb2)/(s$sigma2*s$sa2)

  return(Eloglik - KL1 - KL2)
}
