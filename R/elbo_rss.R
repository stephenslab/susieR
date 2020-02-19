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
  result =  -(length(attr(R, 'eigen')$values)/2) * log(2*pi*s$sigma2) -
    (1/(2*s$sigma2)) * get_ER2_rss(R,z,s)
  return(result)
}

# @title expected squared residuals
# @importFrom Matrix diag
get_ER2_rss = function(R,z,s){
  Dinv = 1/attr(R, 'eigen')$values
  Utz = crossprod(attr(R, 'eigen')$vectors, z)
  zRinvz = sum( Utz * (Dinv * Utz))

  Z = s$alpha*s$mu
  RZ2 = sum((Z%*%R) * Z)

  zbar = colSums(Z)
  Utzbar = crossprod(attr(R, 'eigen')$vectors, zbar)
  Utz = crossprod(attr(R, 'eigen')$vectors, z)
  postb2 = s$alpha * s$mu2 #posterior second moment

  return(zRinvz - 2*sum( Utzbar * Utz) + sum( (sqrt(attr(R, 'eigen')$values) * Utzbar)^2 ) -
           RZ2 + sum(attr(R, 'd')*t(postb2)))
}


#' @title posterior expected loglikelihood for a single effect regression
#' @param R a p by p symmetric and positive semidefinite correlation matrix.
#' @param r residuals
#' @param s2 residual variance
#' @param Eb the posterior mean of b (p vector) (alpha * mu)
#' @param Eb2 the posterior second moment of b (p vector) (alpha * mu2)
#' @keywords internal
SER_posterior_e_loglik_rss = function(R, r, s2, Ez,Ez2){
  eigenR = attr(R,'eigen')

  rtU = crossprod(attr(R, 'eigen')$vectors, r)
  Utz = crossprod(attr(R, 'eigen')$vectors, Ez)

  - (1/(2*s2)) * (- 2*sum(rtU * Utz) + sum(attr(R, 'd')*as.vector(Ez2)))
}
