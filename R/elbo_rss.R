#' @title Get objective function from data and susie fit object.
#'
#' @param R p by p corelation matrix
#' @param z length p vector
#' @param s a susie fit object
#' @keywords internal
get_objective_rss = function(R, z, s) {
  return(Eloglik_rss(R,z,s)-sum(s$KL))
}

# @title expected loglikelihood for a susie fit
Eloglik_rss = function(R,z,s){
  d = s$sigma2 * attr(R, 'eigen')$values + attr(R, 'lambda')

  if(attr(R, 'lambda') == 0){
    result =  -(sum(d!=0)/2) * log(2*pi*s$sigma2) -
      0.5 * get_ER2_rss(R,z,s)
  }else{
    result =  -(length(z)/2) * log(2*pi) - 0.5*sum(log(d)) -
      0.5 * get_ER2_rss(R,z,s)
  }
  return(result)
}

# @title expected squared residuals
# @importFrom Matrix diag
get_ER2_rss = function(R,z,s){

  d = s$sigma2 * attr(R, 'eigen')$values + attr(R, 'lambda')

  Dinv = 1/d
  Dinv[is.infinite(Dinv)] = 0
  SinvR = attr(R, 'eigen')$vectors %*% ((Dinv*attr(R, 'eigen')$values) * t(attr(R, 'eigen')$vectors))
  zSinvz = sum(z * (attr(R, 'eigen')$vectors %*% (Dinv * crossprod(attr(R, 'eigen')$vectors, z))))

  Z = s$alpha*s$mu
  if(attr(R, 'lambda') == 0){
    RSinvR = R/s$sigma2
  }else{
    RSinvR = R %*% SinvR
  }
  RZ2 = sum((Z%*%RSinvR) * Z)

  zbar = colSums(Z)

  postb2 = s$alpha * s$mu2 #posterior second moment

  return(zSinvz - 2*sum( (SinvR %*% z) * zbar) + sum(zbar * (RSinvR %*% zbar)) -
           RZ2 + sum(diag(RSinvR)*t(postb2)))
}


#' @title posterior expected loglikelihood for a single effect regression
#' @param R a p by p symmetric and positive semidefinite correlation matrix.
#' @param Sigma residual_var * R + lambda I
#' @param r residuals
#' @param Eb the posterior mean of b (p vector) (alpha * mu)
#' @param Eb2 the posterior second moment of b (p vector) (alpha * mu2)
SER_posterior_e_loglik_rss = function(R, Sigma, r,Ez,Ez2){
  eigenS = attr(Sigma,'eigenS')
  Dinv = 1/(eigenS$values)
  Dinv[is.infinite(Dinv)] = 0
  rR = R %*% r
  SinvEz = eigenS$vectors %*% (Dinv * crossprod(eigenS$vectors, Ez))

  - (0.5) * (- 2*sum(rR*SinvEz) + sum(attr(Sigma, 'RjSinvRj')*as.vector(Ez2)))
}
