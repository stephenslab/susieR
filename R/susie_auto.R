#' @title An attempt to automate reliable running of susie even on hard problems
#' @param X an n by p matrix of covariates
#' @param Y an n vector
#' @param L_init the initial value of L to consider
#' @param Lmax the maximum value of L to consider
#' @param init_tol the tolerance to pass to susie during early runs (set big to run faster)
#' @param standardize logical flag for whether to standardize columns of X to unit variance prior to fitting.
#' Note that `prior_variance` specifies the prior on the coefficients of X after standardization (if performed).
#' If you do not standardize you may need
#' to think carefully about specifying
#' `prior_variance`. Whatever the value of standardize, the coefficients (returned from `coef`) are for X on the original input scale.
#' Any column of X that has zero variance is not standardized, but left as is.
#' @param intercept Should intercept be fitted (default=TRUE) or set to zero (FALSE)
#' @param max_iter maximum number of iterations to perform
#' @param tol convergence tolerance
#' @param verbose if true outputs some progress messages
#' @details Currently runs a 3-stage strategy for each L: first fit susie with very small residual error,
#' then estimate residual error, then estimate prior variance. If the last step estimates
#' some prior variances to be 0 then stop. Otherwise double L and repeat.
#' Initial runs are done with lax tolerance (init_tol); final run done with default tolerance.
#' @export
susie_auto = function(X,Y,L_init=1,L_max=512,verbose=FALSE,init_tol=1,
                      standardize=TRUE,intercept=TRUE,max_iter=100,tol=1e-2){
  L=L_init
  if(verbose){
    message(paste0("Trying L=",L))
  }
  s.0 = susie(X,Y,L=L, residual_variance = 0.01*sd(Y)^2, tol= init_tol, scaled_prior_variance=1,
              estimate_residual_variance = FALSE, estimate_prior_variance = FALSE,
              standardize=standardize,intercept=intercept,max_iter=max_iter)
  s.1 = susie(X,Y,s_init=s.0,tol=init_tol,estimate_residual_variance=TRUE, estimate_prior_variance = FALSE,
              standardize=standardize,intercept=intercept,max_iter=max_iter)
  s.2 = susie(X,Y,s_init=s.1,tol=init_tol,estimate_residual_variance=TRUE, estimate_prior_variance = TRUE,
              standardize=standardize,intercept=intercept,max_iter=max_iter)
  converged = !all(s.2$V>0) # we call it converged, ie L is big enough, if there are any prior variances set to 0

  while(!converged & (L<=L_max)){
    for(i in 1:L){
      s.2 = add_null_effect(s.2) # add in L more effects
      s.2$sigma2 = 0.01*sd(Y)^2 # set residual variance to be small again for next iteration
    }
    L=2*L
    if(verbose){
      message(paste0("Trying L=",L))
    }
    s.0 = susie(X,Y,s_init=s.2,tol= init_tol,estimate_residual_variance = FALSE, estimate_prior_variance = FALSE,
                standardize=standardize,intercept=intercept,max_iter=max_iter)
    s.1 = susie(X,Y,s_init=s.0,tol= init_tol,estimate_residual_variance = TRUE, estimate_prior_variance = FALSE,
                standardize=standardize,intercept=intercept,max_iter=max_iter)
    s.2 = susie(X,Y,s_init=s.1,tol= init_tol,estimate_residual_variance = TRUE, estimate_prior_variance = TRUE,
                standardize=standardize,intercept=intercept,max_iter=max_iter)
    converged = !all(s.2$V>0) # we call it converged, ie L is big enough, if there are any prior variances set to 0
  }
  #final run at default tolerance to improve fit
  s.2 = susie(X,Y,s_init=s.2,estimate_residual_variance = TRUE, estimate_prior_variance = TRUE, tol = tol,
                standardize=standardize,intercept=intercept,max_iter=max_iter)

  return(s.2)
}
