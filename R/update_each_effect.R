#' @title update each effect once
#' @param X an n by p matrix of covariantes
#' @param Y an n vector of data
#' @param s a SuSiE fit
#' @param estimate_prior_variance boolean indicating whether to estimate prior variance
#' @param colSum of X^2
update_each_effect <- function (X, Y, s, estimate_prior_variance=FALSE, trendfiltering, order) {

  # Repeat for each effect to update
  L = nrow(s$alpha)

  if(L>0){
    for (l in 1:L){
      # remove lth effect from fitted values
      s$Xr = s$Xr - compute_Xb(X, (s$alpha[l,] * s$mu[l,]), trendfiltering, order)
      #compute residuals
      R = Y - s$Xr

      res <- single_effect_regression(R,X,s$V[l],s$sigma2,s$pi,
                                      estimate_prior_variance, trendfiltering, order)

      # Update the variational estimate of the posterior mean.
      s$mu[l,] <- res$mu
      s$alpha[l,] <- res$alpha
      s$mu2[l,] <- res$mu2
      s$V[l] <- res$V
      s$lbf[l] <- res$lbf_model
      s$KL[l] <- -res$loglik + SER_posterior_e_loglik(X,R,s$sigma2,res$alpha*res$mu,res$alpha*res$mu2,trendfiltering,order)

      s$Xr <- s$Xr + compute_Xb(X, (s$alpha[l,] * s$mu[l,]), trendfiltering, order)
    }
  }

  return(s)
}


