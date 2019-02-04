#' @title Get objective function from data and susie fit object.
#'
#' @param XtX a p by p matrix, X'X
#' @param Xty a p vector, X'Y,
#' @param s a susie fit object
#' @keywords internal
get_objective_z = function(R, z, s) {
  return(Eloglik_z(R,z,s)-sum(s$KL))
}

# @title expected loglikelihood for a susie fit
Eloglik_z = function(R,z,s){
  p = length(z)
  result =  -(p/2)*log(2*pi*s$sigma2) - 0.5 * log(attr(R, 'det')) - (1/(2*s$sigma2)) * get_ER2_z(R,z,s)
  return(result)
}

# @title expected squared residuals
# @importFrom Matrix diag
get_ER2_z = function(R,z,s){
  B = s$alpha*s$mu
  XB2 = sum((B%*%R) * B)

  betabar = colSums(B)
  d = attr(R, "d")
  postb2 = s$alpha * s$mu2 #posterior second moment

  return(attr(R, 'ztRinvz') - 2*sum(betabar * z) + sum(betabar * (R %*% betabar)) -
           XB2 + sum(d*t(postb2)))
}

