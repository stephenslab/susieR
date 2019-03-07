#' @title update each effect once
#' @param X an n by p matrix of covariantes
#' @param Y an n vector of data
#' @param s a SuSiE fit
#' @param estimate_prior_variance boolean indicating whether to estimate prior variance
#'
update_each_effect <- function (X, Y, s, estimate_prior_variance=FALSE,
                                estimate_prior_method="optim") {
  if(estimate_prior_variance==FALSE) estimate_prior_method="none"

  # Repeat for each effect to update
  L = nrow(s$alpha)

  if(L>0){
    for (l in 1:L){
      # remove lth effect from fitted values
      s$Xr = s$Xr - compute_Xb(X, (s$alpha[l,] * s$mu[l,]))
      #compute residuals
      R = Y - s$Xr

      res <- single_effect_regression(R,X,s$V[l],s$sigma2,s$pi,
                                      estimate_prior_method)

      # Update the variational estimate of the posterior mean.
      s$mu[l,] <- res$mu
      s$alpha[l,] <- res$alpha
      s$mu2[l,] <- res$mu2
      s$V[l] <- res$V
      s$lbf[l] <- res$lbf_model
      s$KL[l] <- -res$loglik + SER_posterior_e_loglik(X,R,s$sigma2,res$alpha*res$mu,res$alpha*res$mu2)

      s$Xr <- s$Xr + compute_Xb(X, (s$alpha[l,] * s$mu[l,]))
    }
  }

  return(s)
}


