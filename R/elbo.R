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
                         a0, b0, a1, b1) {
  eps <- .Machine$double.eps
  v <- pmax(s2 - m^2, eps)

  # KL for gamma
  KL_gamma <- sum(alpha * (log(pmax(alpha, eps)) - log(pmax(pi, eps))))

  # Expectations under q(sigma^2) ~ IG(a1,b1)
  E_log_sigma2  <- digamma(a1) - log(b1)
  E_inv_sigma2  <- a1 / b1

  # KL for beta given sigma^2
  KL_beta <- 0.5 * sum(alpha * (
    E_log_sigma2 + log(tau2) - log(v) +
      E_inv_sigma2 * (v + m^2) / tau2 - 1
  ))

  # KL between IG posterior and IG prior
  KL_sigma2 <- lgamma(a0) - lgamma(a1) +
    a0 * log(b1/b0) +
    (a1 - a0) * digamma(a1) - a1 + (a1 * b0) / b1

  KL_total <- KL_gamma + KL_beta + KL_sigma2

  list(KL = as.numeric(KL_total),
       components = list(KL_gamma = KL_gamma,
                         KL_beta = KL_beta,
                         KL_sigma2 = KL_sigma2))
}

objective_L1_SER= function(X,y,s){
 out=  Eloglik(X,y,s)-kl_L1_SS_SER (alpha = s$alpha,
                                     m=s$mu, s2= s$mu^2 -s$mu2,
                                     pi=rep( 1/ncol(s$alpha),ncol(s$alpha)) ,
                                     tau2= 1/s$V,
                                     #sigma2 = s$sigma2,
                                    a0=1,b0=1,
                                    a1=nrow(X)/2+1,b1= nrow(X)*s$sigma2+1
  )$KL
 return(out)
}
