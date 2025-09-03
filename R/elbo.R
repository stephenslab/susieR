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

# Compute the ELBO (which is also the marginal log-likelihood since
# the ELBO is exact). The inputs are: ll0, the log-likelihood under
# the "null" model; and lbf, the log-Bayes factors.
compute_elbo <- function (ll0, lbf) {
  x <- max(lbf)
  bf <- exp(lbf - x)
  return(ll0 + x + log(mean(bf)))
}

# Compute the log-normalizing factor for the IG(a,b) distribution.
inv_gamma_factor <- function (a, b)
  a*log(b) - lgamma(a)

# Compute the log-likelihood under the "null" model under the
# normal-inverse-gamma (NIG) prior.
compute_null_loglik_NIG <- function (n, yy, a0, b0) {
  return(-n*log(2*pi)/2 +
           inv_gamma_factor(a0/2,b0/2) -
           inv_gamma_factor((a0 + n)/2,(b0 + yy)/2))
}


objective_L1_SER= function(X,y,s){

  lbf=s$lbf_variable[1,]

  a0= max(s$alpha0,0.01)
  b0=max(s$beta0,0.01)
  yy=sum(y*y)
  n =length(y)


  ll0=  compute_null_loglik_NIG(n, yy, a0, b0)
  out=  compute_elbo (ll0, lbf)
 return(out)
}
