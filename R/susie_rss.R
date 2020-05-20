#' @title Bayesian sum of single-effect (SuSiE) linear regression using z scores
#' @details Performs sum of single-effect (SuSiE) linear regression with z scores.
#' The summary data required are the p by p correlation matrix R, the p vector z. The summary stats should come from the same individuals.
#' This function fits the regression model z = sum_l Rb_l + e, where e is N(0,residual_variance * R) and the
#' sum_l b_l is a p vector of effects to be estimated.
#' The assumption is that each b_l has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0,var=prior_variance).
#' @param z a p vector of z scores.
#' @param R a p by p symmetric and positive semidefinite correlation matrix. If it is from a reference panel,
#' we recommend a modification on correlation matrix with parameter `z_ld_weight`.
#' @param maf minor allele frequency; to be used along with `maf_thresh` to filter input summary statistics
#' @param maf_thresh variants having MAF smaller than this threshold will be filtered out
#' @param z_ld_weight the weight assigned to the z in the ld matrix, the ld matrix used in model is cov2cor((1-w) R + w zz').
#' We recomment setting `z_ld_weight` as 1/number of samples in reference panel.
#' @param L maximum number of non-zero effects
#' @param prior_variance the prior variance (vector of length L, or scalar. In latter case gets repeated L times )
#' @param residual_variance the residual variance, a scaler between 0 and 1
#' @param r_tol tolerance level for eigen value check of positive semidefinite matrix of R.
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param null_weight probability of no effect, for each single effect model
#' @param estimate_residual_variance indicates whether to estimate residual variance
#' @param estimate_prior_variance indicates whether to estimate prior
#' @param estimate_prior_method The method used for estimating prior
#' variance. "simple" method only compares the loglikelihood between
#' using specified prior variance and using zero, and chose the one that
#' gives larger loglikelihood.
#' @param check_null_threshold when prior variance is estimated, compare the
#' estimate with the null and set prior variance to null (zero) unless the log-likelihood
#' using the estimate is larger than that of null by this threshold. For example,
#' you can set it to 0.1 to nudge the estimate towards zero. Default is 0. Notice that setting it to non-zero
#' may lead to decreasing ELBO in some cases.
#' @param prior_tol when prior variance is estimated, compare the estimated value to this tol at the end of
#' the analysis and exclude a single effect from PIP computation if the estimated prior variance is smaller than it.
#' @param max_iter maximum number of iterations to perform
#' @param s_init a previous susie fit with which to initialize
#' @param intercept_value a value to assign to the intercept (since the intercept cannot be estimated from centered summary data). This
#' value will be used by coef.susie() to assign an intercept value, for consistency with the non-summary-statistic version of this function \code{susie}.
#' Set to NULL if you want coef.susie() not to include an intercept term (and so only return a p vector).
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% credible interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a credible set.
#' @param tol convergence tolerance based on elbo
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @param check_R check whether R is positive semidefinite
#' @param check_z check whether z in space spanned by the non-zero eigenvectors of R
#'
#' @return a susie fit, which is a list with some or all of the following elements\cr
#' \item{alpha}{an L by p matrix of posterior inclusion probabilites}
#' \item{mu}{an L by p matrix of posterior means (conditional on inclusion)}
#' \item{mu2}{an L by p matrix of posterior second moments (conditional on inclusion)}
#' \item{Rr}{an p vector of t(X) times fitted values, the fitted values equal to X times colSums(alpha*mu))}
#' \item{sigma2}{residual variance}
#' \item{V}{prior variance}
#'
#' @export
susie_rss = function(z, R, maf=NULL, maf_thresh=0, z_ld_weight=0,
                     L=10, prior_variance=50, residual_variance=NULL,
                     r_tol=1E-08,
                     prior_weights=NULL, null_weight=NULL,
                     estimate_residual_variance=TRUE,
                     estimate_prior_variance=TRUE,
                     estimate_prior_method=c("optim", "EM", "simple"),
                     check_null_threshold=0, prior_tol = 1E-9,
                     max_iter=100, s_init=list(), intercept_value=0,
                     coverage=0.95, min_abs_corr=0.5,
                     tol=1E-03, verbose=FALSE, track_fit=FALSE, check_R=TRUE, check_z=TRUE){
  # Check input R.
  if(nrow(R) != length(z)) {
    stop(paste0('The dimension of correlation matrix (', nrow(R), ' by ', ncol(R), ') does not agree with expected (', length(z), ' by ', length(z), ')'))
  }
  if(!is_symmetric_matrix(R)){
    stop('R is not a symmetric matrix.')
  }
  if (!(is.double(R) & is.matrix(R)) & !inherits(R,"CsparseMatrix"))
    stop("Input R must be a double-precision matrix, or a sparse matrix.")

  # MAF filter
  if(!is.null(maf)){
    if(length(maf) != length(z)){
      stop(paste0('The length of maf does not agree with expected ', length(z)))
    }
    id = which(maf > maf_thresh)
    R = R[id, id]
    z = z[id]
  }

  if(any(is.infinite(z))){
    stop('z contains infinite value.')
  }
  # Check NA in R
  if(any(is.na(R))){
    stop('R matrix contains missing values.')
  }
  # replace NA in z with 0
  if (any(is.na(z))){
    warning('NA values in z-scores are replaced with 0.')
    z[is.na(z)] = 0
  }

  # modification of R
  if(z_ld_weight > 0){
    R = muffled_cov2cor((1-z_ld_weight)* R + z_ld_weight * tcrossprod(z))
    R = (R + t(R))/2
    check_z = FALSE
  }

  if (is.numeric(null_weight) && null_weight == 0) null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight<0 || null_weight>=1)
      stop('Null weight must be between 0 and 1')
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(R)*(1-null_weight), ncol(R)), null_weight)
    else
      prior_weights = c(prior_weights * (1-null_weight), null_weight)
    R = cbind(rbind(R, 0),0)
    z = c(z, 0)
  }

  if(!is.null(residual_variance) && (residual_variance>1 | residual_variance<0)){
    stop('Residual variance should be a scaler between 0 and 1.')
  }
  if(is.null(residual_variance)){
    residual_variance = 1
  }
  p = ncol(R)
  ## eigen decomposition for R, fileter on eigenvalues
  attr(R, 'eigen') = eigen(R, symmetric = T)
  if (check_R && any(attr(R, 'eigen')$values < -r_tol)) {
    stop(paste0('The correlation matrix (', nrow(R), ' by ', ncol(R), 'is not a positive semidefinite matrix.
                You can bypass this by "check_R = FALSE" which instead sets negative eigenvalues to 0 to allow for continued computations.'))
  }

  ## check whether z in space spanned by the non-zero eigenvectors of R
  if(check_z) {
    proj = check_projection(R, z)
    if(proj$status == FALSE) {
      warning('Input z does not lie in the space of non-zero eigenvectors of R. The result is thus not reliable.
              Please refer to https://github.com/stephenslab/susieR/issues/91 for a possible solution.')
    } else {
      write('Input z is in space spanned by the non-zero eigenvectors of R.
            You can safely set "check_z = FALSE" when you rerun the analysis, to save computation.', stderr())
    }
  }
  R = set_R_attributes(R, r_tol)

  ## transform data
  X = t(attr(R, 'eigen')$vectors[, attr(R, 'eigen')$values!=0]) * attr(R, 'eigen')$values[attr(R, 'eigen')$values!=0]^(0.5)
  Y = (t(attr(R, 'eigen')$vectors[, attr(R, 'eigen')$values!=0]) * attr(R, 'eigen')$values[attr(R, 'eigen')$values!=0]^(-0.5)) %*% z

  s = susie(X, Y, L=L, scaled_prior_variance = prior_variance/var(Y),
            residual_variance=residual_variance,
            prior_weights=prior_weights, null_weight=NULL,
            standardize=FALSE, intercept=FALSE,
            estimate_residual_variance=estimate_residual_variance,
            estimate_prior_variance=estimate_prior_variance,
            estimate_prior_method=estimate_prior_method,
            check_null_threshold=check_null_threshold, prior_tol=prior_tol,
            residual_variance_upperbound=1,
            s_init=s_init, coverage=coverage, min_abs_corr=min_abs_corr,
            compute_univariate_zscore=FALSE,
            na.rm=FALSE, max_iter=max_iter,tol=tol,
            verbose=verbose, track_fit=track_fit)
  s$Rz = crossprod(X, s$Xr)
  s$fitted = s$Rz
  return(s)
}

#' @title Bayesian sum of single-effect (SuSiE) linear regression using z scores with lambda (internal test)
#' @details Performs sum of single-effect (SuSiE) linear regression with z scores.
#' The summary data required are the p by p correlation matrix R, the p vector z.
#' The summary stats should come from the same individuals.
#' This function fits the regression model z = sum_l Rb_l + e, where e is N(0,residual_variance * R + lambda I) and the
#' sum_l b_l is a p vector of effects to be estimated.
#' The assumption is that each b_l has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0,var=prior_variance).
#' @param z a p vector of z scores.
#' @param R a p by p symmetric and positive semidefinite correlation matrix.
#' @param maf minor allele frequency; to be used along with `maf_thresh` to filter input summary statistics
#' @param maf_thresh variants having MAF smaller than this threshold will be filtered out
#' @param L maximum number of non-zero effects
#' @param lambda fudge factor
#' @param prior_variance the prior variance (vector of length L, or scalar. In latter case gets repeated L times )
#' @param residual_variance the residual variance, a scaler between 0 and 1
#' @param r_tol tolerance level for eigen value check of positive semidefinite matrix of R.
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param null_weight probability of no effect, for each single effect model
#' @param estimate_residual_variance indicates whether to estimate residual variance
#' @param estimate_prior_variance indicates whether to estimate prior
#' @param estimate_prior_method The method used for estimating prior
#' variance. "simple" method only compares the loglikelihood between
#' using specified prior variance and using zero, and chose the one that
#' gives larger loglikelihood.
#' @param check_null_threshold when prior variance is estimated, compare the
#' estimate with the null and set prior variance to null (zero) unless the log-likelihood
#' using the estimate is larger than that of null by this threshold. For example,
#' you can set it to 0.1 to nudge the estimate towards zero. Default is 0. Notice that setting it to non-zero
#' may lead to decreasing ELBO in some cases.
#' @param prior_tol when prior variance is estimated, compare the estimated value to this tol at the end of
#' the analysis and exclude a single effect from PIP computation if the estimated prior variance is smaller than it
#' @param max_iter maximum number of iterations to perform
#' @param s_init a previous susie fit with which to initialize
#' @param intercept_value a value to assign to the intercept (since the intercept cannot be estimated from centered summary data). This
#' value will be used by coef.susie() to assign an intercept value, for consistency with the non-summary-statistic version of this function \code{susie}.
#' Set to NULL if you want coef.susie() not to include an intercept term (and so only return a p vector).
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% credible interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a credible set.
#' @param tol convergence tolerance based on elbo
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @param check_R check whether R is positive semidefinite
#' @param check_z check whether z in space spanned by the non-zero eigenvectors of R
#'
#' @return a susie fit, which is a list with some or all of the following elements\cr
#' \item{alpha}{an L by p matrix of posterior inclusion probabilites}
#' \item{mu}{an L by p matrix of posterior means (conditional on inclusion)}
#' \item{mu2}{an L by p matrix of posterior second moments (conditional on inclusion)}
#' \item{Rr}{an p vector of t(X) times fitted values, the fitted values equal to X times colSums(alpha*mu))}
#' \item{sigma2}{residual variance}
#' \item{V}{prior variance}
#'
susie_rss_lambda = function(z, R, maf=NULL, maf_thresh=0,
                            L=10, lambda=0,
                            prior_variance=50, residual_variance=NULL,
                            r_tol = 1e-08, prior_weights=NULL, null_weight=NULL,
                            estimate_residual_variance=TRUE,
                            estimate_prior_variance=TRUE,
                            estimate_prior_method=c("optim","EM","simple"),
                            check_null_threshold=0, prior_tol = 1E-9,
                            max_iter=100, s_init=NULL, intercept_value=0,
                            coverage=0.95, min_abs_corr=0.5,
                            tol=1e-3, verbose=FALSE, track_fit = FALSE,
                            check_R=TRUE, check_z=TRUE) {
  # Check input R.
  if(nrow(R) != length(z)) {
    stop(paste0('The dimension of correlation matrix (', nrow(R), ' by ', ncol(R), ') does not agree with expected (', length(z), ' by ', length(z), ')'))
  }
  if(!is_symmetric_matrix(R)){
    stop('R is not a symmetric matrix.')
  }
  if (!(is.double(R) & is.matrix(R)) & !inherits(R,"CsparseMatrix"))
    stop("Input R must be a double-precision matrix, or a sparse matrix.")

  # MAF filter
  if(!is.null(maf)){
    if(length(maf) != length(z)){
      stop(paste0('The length of maf does not agree with expected ', length(z)))
    }
    id = which(maf > maf_thresh)
    R = R[id, id]
    z = z[id]
  }

  if(any(is.infinite(z))){
    stop('z contains infinite value.')
  }
  # Check NA in R
  if(any(is.na(R))){
    stop('R matrix contains missing values.')
  }
  # replace NA in z with 0
  if (any(is.na(z))){
    warning('NA values in z-scores are replaced with 0.')
    z[is.na(z)] = 0
  }
  if (is.numeric(null_weight) && null_weight == 0) null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight<0 || null_weight>=1)
      stop('Null weight must be between 0 and 1')
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(R)*(1-null_weight), ncol(R)), null_weight)
    else
      prior_weights = c(prior_weights * (1-null_weight), null_weight)
    R = cbind(rbind(R, 0),0)
    z = c(z, 0)
  }

  p = ncol(R)
  ## eigen decomposition for R, fileter on eigenvalues
  attr(R, 'eigen') = eigen(R, symmetric = T)
  if (check_R && any(attr(R, 'eigen')$values < -r_tol)) {
    stop(paste0('The correlation matrix (', nrow(R), ' by ', ncol(R), 'is not a positive semidefinite matrix.
                You can bypass this by "check_R = FALSE" which instead sets negative eigenvalues to 0 to allow for continued computations.'))
  }

  ## check whether z in space spanned by the non-zero eigenvectors of R
  if(check_z) {
    proj = check_projection(R, z)
    if(proj$status == FALSE) {
      warning('Input z does not lie in the space of non-zero eigenvectors of R. The result is thus not reliable.
              Please refer to https://github.com/stephenslab/susieR/issues/91 for a possible solution.')
    } else {
      write('Input z is in space spanned by the non-zero eigenvectors of R.
            You can safely set "check_z = FALSE" when you rerun the analysis, to save computation.', stderr())
    }
  }
  R = set_R_attributes(R, r_tol)
  # initialize susie fit
  s = init_setup_rss(p,L,prior_variance,residual_variance,prior_weights,null_weight)
  if (!missing(s_init)) {
    s = modifyList(s, s_init)
    s = init_finalize_rss(s, R=R)
  } else {
    s = init_finalize_rss(s)
  }

  estimate_prior_method <- match.arg(estimate_prior_method)
  # intialize elbo to NA
  elbo = rep(NA,max_iter+1)
  elbo[1] = -Inf;
  tracking = list()

  attr(R, 'lambda') = lambda
  Sigma = update_Sigma(R, s$sigma2, z) # sigma2*R + lambda I

  for(i in 1:max_iter){
    if (track_fit)
      tracking[[i]] = susie_slim(s)
    s = update_each_effect_rss(R, z, s, Sigma, estimate_prior_variance,estimate_prior_method,check_null_threshold)

    if(verbose){
      print(paste0("before estimate sigma2 objective:",get_objective_rss(R, z, s)))
    }
    # compute objective before updating residual variance
    # because part of the objective s$kl has already been computed
    # under the residual variance before the update
    elbo[i+1] = get_objective_rss(R, z, s)
    if((elbo[i+1]-elbo[i])<tol) {
      s$converged = TRUE
      break;
    }
    if(estimate_residual_variance){

      if(lambda == 0){
        tmp = s
        tmp$sigma2 = 1
        est_sigma2 = (1/sum(attr(R, 'eigen')$values!=0))* get_ER2_rss(R,z,tmp)
        if(est_sigma2 < 0){
          stop('Estimating residual variance failed: the estimated value is negative')
        }
        if(est_sigma2 > 1){
          est_sigma2 = 1
        }
      }else{
        sigma2 = seq(0.1,1,by=0.1)
        tmp = s
        obj = numeric(length(sigma2))
        for(j in 1:length(sigma2)){
          tmp$sigma2 = sigma2[j]
          obj[j] = Eloglik_rss(R, z, tmp)
        }
        est_sigma2 = sigma2[which.max(obj)]
      }

      s$sigma2 = est_sigma2

      if(verbose){
        print(paste0("after estimate sigma2 objective:", get_objective_rss(R, z, s)))
      }
      Sigma = update_Sigma(R, s$sigma2, z)
    }
  }
  elbo = elbo[2:(i+1)] # Remove first (infinite) entry, and trailing NAs.
  s$elbo <- elbo
  s$niter <- i

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in", max_iter, "iterations!"))
    s$converged = FALSE
  }

  s$intercept = intercept_value
  s$fitted = s$Rz

  s$X_column_scale_factors = attr(R,"scaled:scale")

  if (track_fit)
    s$trace = tracking

  ## SuSiE CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    R = muffled_cov2cor(R)
    s$sets = susie_get_cs(s, coverage=coverage, Xcorr=R, min_abs_corr=min_abs_corr)
    s$pip = susie_get_pip(s, prune_by_cs = FALSE, prior_tol = prior_tol)
  }

  return(s)
}

update_Sigma = function(R, sigma2, z){
  Sigma = sigma2*R + attr(R, 'lambda') *diag(length(z))
  eigenS = attr(R, 'eigen')
  eigenS$values = sigma2*eigenS$values + attr(R, 'lambda')

  Dinv = 1/(eigenS$values)
  Dinv[is.infinite(Dinv)] = 0
  attr(Sigma, 'eigenS') = eigenS

  # Sigma^(-1) R_j = U (sigma2 D + lambda)^(-1) D U^T e_j
  attr(Sigma, 'SinvRj') = eigenS$vectors %*% (Dinv*attr(R, 'eigen')$values * t(eigenS$vectors))

  if(attr(R, 'lambda')==0){
    attr(Sigma, 'RjSinvRj') = attr(R, 'd')/sigma2
  }else{
    tmp = t(eigenS$vectors)
    attr(Sigma, 'RjSinvRj') = colSums(tmp * (Dinv*(attr(R, 'eigen')$values^2) * tmp))
  }

  return(Sigma)
}

