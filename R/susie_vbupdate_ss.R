#' @title update each effect once
#' @param XtX a p by p matrix, t(X)X
#' @param XtY a p vector
#' @param s_init a list with elements sigma2, V, alpha, mu, Xr
#' @param estimate_prior_variance boolean indicating whether to estimate prior variance
update_each_effect_ss <- function (XtX, XtY, s_init, estimate_prior_variance=FALSE) {

  # Repeat for each effect to update
  s = s_init
  L = nrow(s$alpha)

  if(L>0){
    for (l in 1:L){
      # remove lth effect from fitted values
      s$XtXr = s$XtXr - XtX %*% (s$alpha[l,] * s$mu[l,])

      #compute residuals
      XtR = XtY - s$XtXr

      res = single_effect_regression_ss(XtR,diag(XtX),s$V[l],s$sigma2,estimate_prior_variance)

      # Update the variational estimate of the posterior mean.
      s$mu[l,] <- res$mu
      s$alpha[l,] <- res$alpha
      s$mu2[l,] <- res$mu2
      s$V[l] <- res$V
      s$KL[l] <- -res$logBF + SER_posterior_e_loglik_ss(XtX,XtR,s$sigma2,res$alpha*res$mu,res$alpha*res$mu2)

      s$XtXr <- s$XtXr + XtX %*% (s$alpha[l,]*s$mu[l,])
    }
  }

  return(s)
}


