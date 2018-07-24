#' @title Get objective function from data and susie fit object.
#'
#' @param XtX a p by p matrix, X'X
#' @param XtY a p vector, X'Y,
#' @param s a susie fit object
#' @param var_y the (sample) variance of the vector Y
#' @param n sample size
#'
#' @export
#'
susie_get_objective_ss = function(XtX, XtY, s, var_y, n) {
  return(Eloglik_ss(XtX,XtY,s,var_y, n)-sum(s$KL))
}

#' @title expected loglikelihood for a susie fit
Eloglik_ss = function(XtX,XtY,s, var_y, n){
  p = ncol(XtX)
  result =  -(n/2) * log(2*pi* s$sigma2) - (1/(2*s$sigma2)) * get_ER2_ss(XtX,XtY,s,var_y,n)
  return(result)
}

# expected squared residuals
get_ER2_ss = function(XtX,XtY,s,var_y,n){
  B = s$alpha*s$mu
  XB2 = sum(t(B) * XtX%*%t(B))

  betabar = colSums(B)
  d = diag(XtX)
  postb2 = s$alpha * s$mu2 #posterior second moment

  return(var_y*n - 2*sum(betabar * XtY) + sum(betabar * (XtX %*% betabar)) -
           XB2 + sum(d*t(postb2)))
}


#' @title posterior expected loglikelihood for a single effect regression
#' @param X an n by p matrix of covariates
#' @param Y an n vector of regression outcome
#' @param s2 the residual variance
#' @param Eb the posterior mean of b (p vector) (alpha * mu)
#' @param Eb2 the posterior second moment of b (p vector) (alpha * mu2)
SER_posterior_e_loglik_ss = function(XtX,XtY,s2,Eb,Eb2){
  d = diag(XtX)
  - (0.5/s2) * (- 2*sum(Eb*XtY) + sum(d*as.vector(Eb2)))
}
