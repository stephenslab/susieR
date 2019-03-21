#' @title update each effect once
#' @param R a p by p symmetric and positive semidefinite correlation matrix.
#' @param z a p vector
#' @param s_init a list with elements sigma2, V, alpha, mu, Xr
#' @param Sigma positive definite, R + theta I
#' @param estimate_prior_variance boolean indicating whether to estimate prior variance
#' @importFrom Matrix diag
update_each_effect_rss <- function (R, z, s_init, Sigma, estimate_prior_variance=FALSE,estimate_prior_method="optim") {

  if(estimate_prior_variance==FALSE) estimate_prior_method="none"
  # Repeat for each effect to update
  s = s_init
  L = nrow(s$alpha)

  if(L>0){
    for (l in 1:L){
      # remove lth effect from fitted values
      s$Rz = s$Rz - R %*% (s$alpha[l,] * s$mu[l,])

      #compute residuals
      r = z - s$Rz
      res = single_effect_regression_rss(as.vector(r),Sigma,s$V[l], s$sigma2,s$pi,estimate_prior_method)
      # Update the variational estimate of the posterior mean.
      s$mu[l,] <- res$mu
      s$alpha[l,] <- res$alpha
      s$mu2[l,] <- res$mu2
      s$V[l] <- res$V
      s$lbf[l] <- res$lbf_model
      s$KL[l] <- -res$lbf_model + SER_posterior_e_loglik_rss(R, Sigma,r, res$alpha*res$mu,res$alpha*res$mu2)

      s$Rz <- s$Rz + R %*% (s$alpha[l,]*s$mu[l,])
    }
  }
  s$Rz = unname(as.matrix(s$Rz))
  return(s)
}


