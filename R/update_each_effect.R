# @title update each effect once
# @param X an n by p matrix of regressor variables
# @param y an n vector of response variable
# @param s a SuSiE fit
# @param estimate_prior_variance boolean indicating whether to
#   estimate prior variance
# @param check_null_threshold float a threshold on the log scale to
#   compare likelihood between current estimate and zero the null
# @param small Logical. Useful when fitting susie on data with a limited sample size.
#     If set to TRUE, susie is fitted using single-effect regression with the Servin and Stephens prior
#     instead of the default Gaussian prior. This improves the calibration of credible sets.
#     Default is FALSE.
#   @param alpha  numerical parameter for the NIG prior when using Servin
#   and Stephens SER
#
#   @param beta  numerical parameter for the NIG prior when using Servin
#   and Stephens SER
update_each_effect = function (X, y, s, estimate_prior_variance = FALSE,
                               estimate_prior_method = "optim",
                               check_null_threshold,
                               small=FALSE,
                               alpha0=0,
                               beta0=0) {
  if (!estimate_prior_variance)
    estimate_prior_method = "none"

  # Repeat for each effect to update.
  L = nrow(s$alpha)
  if (L > 0)
    for (l in 1:L) {

      # Remove lth effect from fitted values.
      s$Xr = s$Xr - compute_Xb(X,s$alpha[l,] * s$mu[l,])

      # Compute residuals.
      R = y - s$Xr

      res = single_effect_regression(R,X,s$V[l],s$sigma2,s$pi,
                                     estimate_prior_method,
                                     check_null_threshold,
                                     small = small,
                                     alpha0 = alpha0,
                                     beta0  = beta0)
      # Update the variational estimate of the posterior mean.
      s$mu[l,]    = res$mu
      s$alpha[l,] = res$alpha
      s$mu2[l,]   = res$mu2
      s$V[l]      = res$V
      s$lbf[l]    = res$lbf_model
      s$lbf_variable[l,] = res$lbf
      s$KL[l]     = -res$loglik +
        SER_posterior_e_loglik(X,R,s$sigma2,res$alpha * res$mu,
                               res$alpha * res$mu2)

      s$Xr = s$Xr + compute_Xb(X,s$alpha[l,] * s$mu[l,])
    }
  return(s)
}


