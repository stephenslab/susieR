#' @title Apply susie with z scores and correlation matrix R
#' @param z a p vector of z scores.
#' @param R a p by p symmetric and positive semidefinite matrix
#' @param R_type the type of the matrix R, one of 'Cor', 'XtX' or 'Cov'
#' @param n sample size
#' @param var_y the (sample) variance of the y
#' @param L maximum number of non-zero effects
#' @param scaled_prior_variance the scaled prior variance (vector of length L, or scalar. In latter case gets repeated L times )
#' @param estimate_prior_variance indicates whether to estimate prior (currently not recommended as not working as well)
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param null_weight probability of no effect, for each single effect model
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% confidence interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a confidence set.
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @param ... further arguments to be passed to \code{\link{susie_ss}}
#' @return a susie fit
#'
#' @export
susie_z = function(z, R, R_type = c('Cor', 'Cov', 'XtX'), n, var_y = 1,
                   L=10, scaled_prior_variance=0.2, estimate_prior_variance = TRUE,
                   prior_weights = NULL, null_weight = NULL,
                   coverage=0.95, min_abs_corr=0.5,
                   verbose=FALSE, track_fit = FALSE, ...){
  if(nrow(R) != length(z)){
    stop('The dimension of R does not agree with length of z.')
  }
  R_type = match.arg(R_type)
  if(R_type != 'XtX'){
    XtX = (n-1)*R
  }else{
    XtX = R
  }
  susie_ss(XtX = XtX, Xty = z*sqrt(diag(XtX)*var_y), var_y = var_y, n = n, L = L,
           scaled_prior_variance = scaled_prior_variance,
           estimate_prior_variance = estimate_prior_variance,
           prior_weights = prior_weights, null_weight = null_weight,
           coverage=coverage, min_abs_corr=min_abs_corr,
           verbose=verbose, track_fit = track_fit, ...)
}

#' @title Apply susie with bhat, shat and correlation matrix R
#' @param bhat a p vector of estimated effects.
#' @param shat a p vector of corresponding standard errors.
#' @param R a p by p correlation matrix of X
#' @param R_type the type of the matrix R, one of 'Cor', 'XtX' or 'Cov'
#' @param n sample size
#' @param var_y the (sample) variance of the y
#' @param L maximum number of non-zero effects
#' @param scaled_prior_variance the scaled prior variance (vector of length L, or scalar. In latter case gets repeated L times )
#' @param estimate_prior_variance indicates whether to estimate prior (currently not recommended as not working as well)
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param null_weight probability of no effect, for each single effect model
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% confidence interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a confidence set.
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @param ... further arguments to be passed to \code{\link{susie_ss}}
#' @return a susie fit
#'
#' @export
susie_bhat = function(bhat, shat, R, R_type = c('Cor', 'Cov', 'XtX'), n, var_y = 1, L=10, standardize = TRUE,
                    scaled_prior_variance=0.2,
                    estimate_prior_variance = FALSE,
                    prior_weights = NULL, null_weight = NULL,
                    coverage=0.95, min_abs_corr=0.5,
                    verbose=FALSE, track_fit = FALSE, ...){
  if(length(bhat) != length(shat)){
    stop('The length of bhat does not agree with length of shat.')
  }
  if(nrow(R) != length(bhat)){
    stop('The dimension of R does not agree with length of bhat.')
  }

  R_type = match.arg(R_type)
  if(R_type != 'Cor'){
    R = cov2cor(R)
  }

  that = bhat/shat

  R2 = that^2/(that^2 + n-2)
  sigma2 = var_y*(n-1)*(1-R2)/(n-2)

  if(standardize){
    XtX = (n-1)*R
    Xty = sqrt(sigma2) * sqrt(n-1) * that
  }else{
    XtXdiag = sigma2/(shat^2)
    Xty = bhat * XtXdiag
    XtX = t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
  }

  susie_ss(XtX = XtX, Xty = Xty, var_y = var_y, n = n, L = L,
           scaled_prior_variance = scaled_prior_variance,
           estimate_prior_variance = estimate_prior_variance,
           prior_weights = prior_weights, null_weight = null_weight,
           coverage=coverage, min_abs_corr=min_abs_corr,
           verbose=verbose, track_fit = track_fit, ...)
}
