#' @title update each effect once
#' @param XtX a p by p matrix, t(X)X
#' @param Xty a p vector
#' @param s_init a list with elements sigma2, V, alpha, mu, Xr
#' @param estimate_prior_variance boolean indicating whether to estimate prior variance
#' @importFrom Matrix diag
update_each_effect_ss <- function (XtX, Xty, s_init, estimate_prior_variance=FALSE) {

  # Repeat for each effect to update
  s = s_init
  L = nrow(s$alpha)
  dXtX = diag(XtX)

  if(L>0){
    for (l in 1:L){
      # remove lth effect from fitted values
      s$XtXr = s$XtXr - XtX %*% (s$alpha[l,] * s$mu[l,])

      #compute residuals
      XtR = Xty - s$XtXr
      res = single_effect_regression_ss(as.matrix(XtR),dXtX,s$V[l],s$sigma2,estimate_prior_variance)
      # Update the variational estimate of the posterior mean.
      s$mu[l,] <- res$mu
      s$alpha[l,] <- res$alpha
      s$mu2[l,] <- res$mu2
      s$V[l] <- res$V
      s$KL[l] <- -res$logBF + SER_posterior_e_loglik_ss(dXtX,XtR,s$sigma2,res$alpha*res$mu,res$alpha*res$mu2)

      s$XtXr <- s$XtXr + XtX %*% (s$alpha[l,]*s$mu[l,])
    }
  }
  s$XtXr = unname(as.matrix(s$XtXr))
  return(s)
}


