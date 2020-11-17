# @title update each effect once
# @param XtX a p by p matrix, X'X
# @param Xty a p vector
# @param s_init a list with elements sigma2, V, alpha, mu, Xr
# @param estimate_prior_variance boolean indicating whether to
#   estimate prior variance
# @param estimate_prior_method The method used for estimating prior
#   variance, 'optim' or 'EM'.
# @param check_null_threshold float a threshold on the log scale to
#   compare likelihood between current estimate and zero the null
# 
#' @importFrom Matrix diag
update_each_effect_ss = function (XtX, Xty, s_init,
                                  estimate_prior_variance = FALSE,
                                  estimate_prior_method = "optim",
                                  check_null_threshold = 0) {
  if (!estimate_prior_variance)
    estimate_prior_method = "none"
  
  # Repeat for each effect to update.
  s = s_init
  L = nrow(s$alpha)
  if (L > 0) {
    for (l in 1:L) {
        
      # Remove lth effect from fitted values.
      s$XtXr = s$XtXr - XtX %*% (s$alpha[l,] * s$mu[l,])

      # Compute residuals.
      XtR = Xty - s$XtXr
      res = single_effect_regression_ss(as.matrix(XtR),attr(XtX,"d"),s$V[l],
              s$sigma2,s$pi,estimate_prior_method,check_null_threshold)
      
      # Update the variational estimate of the posterior mean.
      s$mu[l,]    = res$mu
      s$alpha[l,] = res$alpha
      s$mu2[l,]   = res$mu2
      s$V[l]      = res$V
      s$lbf[l]    = res$lbf_model
      s$lbf_variable[l,] = res$lbf
      s$KL[l]     = -res$lbf_model +
        SER_posterior_e_loglik_ss(attr(XtX,"d"),XtR,s$sigma2,
                                  res$alpha * res$mu,res$alpha * res$mu2)

      s$XtXr = s$XtXr + XtX %*% (s$alpha[l,] * s$mu[l,])
    }
  }
  s$XtXr = unname(as.matrix(s$XtXr))
  return(s)
}
