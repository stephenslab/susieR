#' @title Get objective function from data and susie fit object.
#'
#' @param XtX a p by p matrix, X'X
#' @param Xty a p vector, X'Y,
#' @param s a susie fit object
#' @param yty a scaler, Y'Y, where Y is centered to have mean 0
#' @param n sample size
#' @keywords internal
get_objective_ss = function(XtX, Xty, s, yty, n) {
  return(Eloglik_ss(XtX,Xty,s,yty, n)-sum(s$KL))
}

# @title expected loglikelihood for a susie fit
Eloglik_ss = function(XtX,Xty,s, yty, n){
  result =  -(n/2) * log(2*pi* s$sigma2) - (1/(2*s$sigma2)) * get_ER2_ss(XtX,Xty,s,yty)
  return(result)
}

# @title expected squared residuals
# @importFrom Matrix diag
get_ER2_ss = function(XtX,Xty,s,yty){
  B = s$alpha*s$mu
  XB2 = sum((B%*%XtX) * B)

  betabar = colSums(B)
  d = attr(XtX, "d")
  postb2 = s$alpha * s$mu2 #posterior second moment

  return(yty - 2*sum(betabar * Xty) + sum(betabar * (XtX %*% betabar)) -
           XB2 + sum(d*t(postb2)))
}


# @title posterior expected loglikelihood for a single effect regression
# @param dXtX a p vector of diagonal elements of XtX
# @param Xty a p vector
# @param s2 the residual variance
# @param Eb the posterior mean of b (p vector) (alpha * mu)
# @param Eb2 the posterior second moment of b (p vector) (alpha * mu2)
SER_posterior_e_loglik_ss = function(dXtX,Xty,s2,Eb,Eb2){
  - (0.5/s2) * (- 2*sum(Eb*Xty) + sum(dXtX*as.vector(Eb2)))
}
