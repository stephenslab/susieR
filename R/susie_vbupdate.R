#' @title update each effect once
#' @param X an n by p matrix of covariantes
#' @param Y an n vector of data
#' @param s_init a list with elements sigma2, sa2, alpha, mu, Xr
update_each_effect <- function (X, Y, s_init) {

  # Repeat for each effect to update
  s = s_init
  L = nrow(s$alpha)
  for (l in 1:L){
    # remove lth effect from fitted values
    s$Xr = s$Xr - X %*% (s$alpha[l,] * s$mu[l,])

    #compute residuals
    R = Y - s$Xr

    res = single_effect_regression(R,X,s$sa2[l],s$sigma2)

    # Update the variational estimate of the posterior mean.
    s$mu[l,] <- res$mu
    s$alpha[l,] <- res$alpha
    s$mu2[l,] <- res$mu2

    s$Xr <- s$Xr + X %*% (s$alpha[l,]*s$mu[l,])
  }

  return(s)
}


