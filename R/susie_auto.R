#' @title An attempt to automate reliable running of susie even on hard problems
#' @param X an n by p matrix of covariates
#' @param Y an n vector
#' @param L_init the initial value of L to consider
#' @param Lmax the maximum value of L to consider
#' @details Currently runs a 3-stage strategy for each L: first fit susie with very small residual error,
#' then estimate residual error, then estimate prior variance. If the last step estimates
#' some prior variances to be 0 then stop. Otherwise double L and repeat.
#' @export
susie_auto = function(X,Y,L_init = 1, L_max= 512, verbose=FALSE){
  L=L_init
  converged = FALSE
  s.0 = susie(X,Y,L=L, residual_variance = 0.01, prior_variance=1, estimate_residual_variance = FALSE, estimate_prior_variance = FALSE)
  while(!converged & (L<=L_max)){
    for(i in 1:L){
      s.0 = add_null_effect(s.0) # add in L more effects
    }
    L=2*L
    if(verbose){
      message(paste0("Trying L=",L))
    }
    s.1 = susie(X,Y,s_init=s.0,estimate_residual_variance = FALSE, estimate_prior_variance = FALSE)
    s.2 = susie(X,Y,s_init=s.1,estimate_residual_variance = TRUE, estimate_prior_variance = FALSE)
    s.0 = susie(X,Y,s_init=s.2,estimate_residual_variance = TRUE, estimate_prior_variance = TRUE)
    converged = !all(s.0$sa2>0) # we call it converged, ie L is big enough, if there are any prior variances set to 0
  }
  return(s.0)
}
