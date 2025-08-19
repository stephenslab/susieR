# @title Get objective function from data and susie fit object.
# @param data A flash data object.
# @param f A flash fit object.
# @keywords internal
get_objective = function (X, Y, s) {
  return(Eloglik(X,Y,s) - sum(s$KL))
}

# Expected loglikelihood for a susie fit.
Eloglik = function (X, Y, s) {
  n = nrow(X)
  return(-(n/2) * log(2*pi*s$sigma2) - (1/(2*s$sigma2)) * get_ER2(X,Y,s))
}

# Expected squared residuals.
# s$Xr is column sum of Xr_L
get_ER2 = function (X, Y, s) {
  Xr_L = compute_MXt(s$alpha * s$mu,X) # L by N matrix
  postb2 = s$alpha * s$mu2 # Posterior second moment.
  return(sum((Y - s$Xr)^2) - sum(Xr_L^2) + sum(attr(X,"d") * t(postb2)))
}

# @title posterior expected loglikelihood for a single effect regression
# @param X an n by p matrix of covariates
# @param Y an n vector of regression outcome
# @param s2 the residual variance
# @param Eb the posterior mean of b (p vector) (alpha * mu)
# @param Eb2 the posterior second moment of b (p vector) (alpha * mu2)
SER_posterior_e_loglik = function (X, Y, s2, Eb, Eb2) {
  n = nrow(X)
  return(-0.5*n*log(2*pi*s2) - 0.5/s2*(sum(Y*Y)
                                       - 2*sum(Y*compute_Xb(X,Eb))
                                       + sum(attr(X,"d") * Eb2)))
}




kl_L1_SS_SER <- function(alpha, m, s2, pi, tau2,
                               sigma2 = NULL,
                               E_inv_sigma2 = NULL, E_log_sigma2 = NULL) {
  eps <- .Machine$double.eps
  v <- pmax(s2 - m^2, eps)  # posterior variance when gamma=1

  # KL for gamma
  KL_gamma <- sum(alpha * (log(pmax(alpha, eps)) - log(pmax(pi, eps))))

  if (!is.null(sigma2)) {
    # Known sigma^2
    KL_norm <- 0.5 * sum(alpha * ( log((sigma2 * tau2) / v) + (v + m^2) / (sigma2 * tau2) - 1 ))
  } else {
    # Integrated sigma^2: need expectations
    if (is.null(E_inv_sigma2) || is.null(E_log_sigma2))
      stop("Provide either sigma2, or both E_inv_sigma2 and E_log_sigma2.")
    KL_norm <- 0.5 * sum(alpha * ( E_log_sigma2 + log(tau2) - log(v) +
                                     E_inv_sigma2 * (v + m^2) / tau2 - 1 ))
  }

  list(KL = as.numeric(KL_gamma + KL_norm),
       components = list(KL_gamma = KL_gamma, KL_beta = KL_norm))
}


objective_L1_SER= function(X,y,s){
 out=  Eloglik(X,y,s)-kl_L1_SS_SER (alpha = s$alpha,
                                     m=s$mu, s2= s$mu^2 -s$mu2,
                                     pi=rep( 1/ncol(s$alpha),ncol(s$alpha)) ,
                                     tau2= s$sigma2,
                                     sigma2 = s$sigma2
  )$KL
 return(out)
}
