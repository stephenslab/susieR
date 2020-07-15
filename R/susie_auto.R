#' @rdname susie
#' 
#' @details susie_auto is an attempt to automate reliable running of
#' susie even on hard problems. Implements a three-stage strategy for
#' each L: first fit susie with very small residual error; next,
#' estimate residual error; finally, estimate the prior variance. If
#' the last step estimates some prior variances to be zero,
#' stop. Otherwise, double L, and repeat.  Initial runs are performed
#' with relaxed tolerance; the final run is performed using the
#' default susie tolerance.
#'
#' @param L_init The initial value of L.
#' 
#' @param L_max The maximum value of L to consider.
#' 
#' @param init_tol The tolerance to passed to susie during early runs
#'   (set large to shorten the initial runs).
#' 
#' @param tol Convergence tolerance.
#'
#' @param \dots Additional arguments passed to \code{\link{susie}}.
#'
#' @seealso \code{\link{susie}}
#' 
#' @importFrom stats sd
#'
#' @export
#' 
susie_auto = function (X, Y, L_init = 1, L_max = 512, verbose = FALSE,
                       init_tol = 1, standardize = TRUE, intercept = TRUE,
                       max_iter = 100,tol = 1e-2,...) {
  L = L_init
  if(verbose)
    message(paste0("Trying L=",L))
  s.0 = susie(X,Y,L = L,residual_variance = 0.01*sd(Y)^2,tol = init_tol,
              scaled_prior_variance = 1,estimate_residual_variance = FALSE,
              estimate_prior_variance = FALSE,standardize = standardize,
              intercept = intercept,max_iter = max_iter,...)
  s.1 = susie(X,Y,s_init = s.0,tol = init_tol,
              estimate_residual_variance = TRUE,
              estimate_prior_variance = FALSE,
              standardize = standardize,intercept = intercept,
              max_iter = max_iter,...)
  s.2 = susie(X,Y,s_init = s.1,tol = init_tol,
              estimate_residual_variance = TRUE,
              estimate_prior_variance = TRUE,
              standardize = standardize,intercept = intercept,
              max_iter = max_iter,...)

  # We call it converged; i.e., L is big enough, if there are any prior
  # variances set to zero.
  converged = !all(s.2$V > 0) 
  while (!converged & (L <= L_max)) {
    for(i in 1:L) {
      s.2 = add_null_effect(s.2,1) # Add in L more effects.
      s.2$sigma2 = 0.01*sd(Y)^2    # Set residual variance to be small
                                   # again for next iteration.
    }
    L = 2*L
    if(verbose)
      message(paste0("Trying L=",L))
    s.0 = susie(X,Y,s_init = s.2,tol = init_tol,
                estimate_residual_variance = FALSE,
                estimate_prior_variance = FALSE,
                standardize = standardize,intercept = intercept,
                max_iter=max_iter,...)
    s.1 = susie(X,Y,s_init = s.0,tol = init_tol,
                estimate_residual_variance = TRUE,
                estimate_prior_variance = FALSE,
                standardize = standardize,intercept = intercept,
                max_iter = max_iter,...)
    s.2 = susie(X,Y,s_init = s.1,tol = init_tol,
                estimate_residual_variance = TRUE,
                estimate_prior_variance = TRUE,
                standardize = standardize,intercept = intercept,
                max_iter = max_iter,...)

    # We call it converged; i.e., L is big enough, if there are any prior
    # variances set to zero.
    converged = !all(s.2$V > 0) 
  }
  
  # Final run at default tolerance to improve fit.
  s.2 = susie(X,Y,s_init = s.2,estimate_residual_variance = TRUE,
              estimate_prior_variance = TRUE,tol = tol,
              standardize = standardize,intercept = intercept,
              max_iter = max_iter,...)
  return(s.2)
}
