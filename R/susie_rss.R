#' @title Sum of Single Effects (SuSiE) Regression using summary statistics
#'
#' @description \code{susie_rss} performs sum of single-effect linear
#'   regression with z scores; all posterior calculations are for
#'   z-scores. This function fits the regression model \eqn{z = \sum_l
#'   R*b_l + e}, where e is \eqn{N(0,R)} and \eqn{\sum_l b_l} is a
#'   p-vector of effects to be estimated. The required summary data are
#'   the p by p correlation matrix, \code{R}, and the p-vector \code{z}.
#'
#' @param z A p-vector of z scores.
#'
#' @param R A p by p symmetric, positive semidefinite correlation
#' matrix.
#'
#' @param maf Minor allele frequency; to be used along with
#'   \code{maf_thresh} to filter input summary statistics.
#'
#' @param maf_thresh Variants having a minor allele frequency smaller
#'   than this threshold are not used.
#'
#' @param z_ld_weight This feature is not recommended. The weights
#'   assigned to the z scores in the LD matrix. If \code{z_ld_weight >
#'   0}, the LD matrix used in the model is \code{cov2cor((1-w)*R +
#'   w*tcrossprod(z))}, where \code{w = z_ld_weight}.
#'
#' @param L Number of components (nonzero coefficients) in the susie
#'   regression model. If L is larger than the number of covariates, p,
#'   L is set to p.
#'
#' @param prior_variance The prior variance. It is either a scalar or
#'   a vector of length L.
#'
#' @param residual_variance Variance of the residual.
#'   If it is not specified, we set it to 1.
#'
#' @param prior_weights A vector of length p, in which each entry
#'   gives the prior probability that SNP j has non-zero effect.
#'
#' @param null_weight Prior probability of no effect (a number between
#'   0 and 1, and cannot be exactly 1).
#'
#' @param estimate_residual_variance The residual variance is
#'   fixed to the value supplied by \code{residual_variance}. We don't
#'   estimate residual variance in susie_rss.
#'
#' @param estimate_prior_variance If \code{estimate_prior_variance =
#'   TRUE}, the prior variance is estimated (this is a separate
#'   parameter for each of the L effects). If provided,
#'   \code{prior_variance} is then used as an initial value for
#'   the optimization. When \code{estimate_prior_variance = FALSE}, the
#'   prior variance for each of the L effects is determined by the
#'   value supplied to \code{prior_variance}.
#'
#' @param estimate_prior_method The method used for estimating prior
#'   variance. When \code{estimate_prior_method = "simple"} is used, the
#'   likelihood at the specified prior variance is compared to the
#'   likelihood at a variance of zero, and the setting with the larger
#'   likelihood is retained.
#'
#' @param check_null_threshold When the prior variance is estimated,
#'   compare the estimate with the null, and set the prior variance to
#'   zero unless the log-likelihood using the estimate is larger by this
#'   threshold amount. For example, if you set
#'   \code{check_null_threshold = 0.1}, this will "nudge" the estimate
#'   towards zero when the difference in log-likelihoods is small. A
#'   note of caution that setting this to a value greater than zero may
#'   lead the IBSS fitting procedure to occasionally decrease the ELBO.
#'
#' @param prior_tol When the prior variance is estimated, compare the
#'   estimated value to \code{prior_tol} at the end of the computation,
#'   and exclude a single effect from PIP computation if the estimated
#'   prior variance is smaller than this tolerance value.
#'
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @param s_init A previous susie fit with which to initialize.
#'
#' @param intercept_value The intercept. (The intercept cannot be
#'   estimated from centered summary data.) This setting will be used by
#'   \code{coef} to assign an intercept value, mainly for consistency
#'   with \code{susie}. Set to \code{NULL} if you want \code{coef} not
#'   to include an intercept term (and so only a p-vector is returned).
#'
#' @param coverage A number between 0 and 1 specifying the
#'   \dQuote{coverage} of the estimated confidence sets.
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
#'
#' @param tol A small, non-negative number specifying the convergence
#'   tolerance for the IBSS fitting procedure. The fitting procedure
#'   will halt when the difference in the variational lower bound, or
#'   \dQuote{ELBO} (the objective function to be maximized), is
#'   less than \code{tol}.
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#'   and a summary of the optimization settings, are printed to the
#'   console.
#'
#' @param track_fit If \code{track_fit = TRUE}, \code{trace}
#'   is also returned containing detailed information about the
#'   estimates at each iteration of the IBSS fitting procedure.
#'
#' @param check_R If \code{check_R = TRUE}, check that \code{R} is
#'   positive semidefinite.
#'
#' @param r_tol Tolerance level for eigenvalue check of positive
#'   semidefinite matrix of R.
#'
#' @param refine If \code{refine = TRUE}, we use a procedure to help
#'   SuSiE get out of local optimum.
#'
#' @return A \code{"susie"} object with some or all of the following
#'   elements:
#'
#' \item{alpha}{An L by p matrix of posterior inclusion probabilites.}
#'
#' \item{mu}{An L by p matrix of posterior means, conditional on
#'   inclusion.}
#'
#' \item{mu2}{An L by p matrix of posterior second moments,
#'   conditional on inclusion.}
#'
#' \item{lbf}{log-Bayes Factor for each single effect.}
#'
#' \item{lbf_variable}{log-Bayes Factor for each variable and single effect.}
#'
#' \item{intercept}{Fixed Intercept.}
#'
#' \item{sigma2}{Fixed Residual variance.}
#'
#' \item{V}{Prior variance of the non-zero elements of b, equal to
#'   \code{scaled_prior_variance * var(Y)}.}
#'
#' \item{elbo}{The value of the variational lower bound, or
#'   \dQuote{ELBO} (objective function to be maximized), achieved at
#'   each iteration of the IBSS fitting procedure.}
#'
#' \item{fitted}{Vector of length n containing the fitted values of
#'   the outcome.}
#'
#' \item{sets}{Credible sets estimated from model fit; see
#'   \code{\link{susie_get_cs}} for details.}
#'
#' \item{pip}{A vector of length p giving the (marginal) posterior
#'   inclusion probabilities for all p covariates.}
#'
#' \item{niter}{Number of IBSS iterations that were performed.}
#'
#' \item{converged}{\code{TRUE} or \code{FALSE} indicating whether
#'   the IBSS converged to a solution within the chosen tolerance
#'   level.}
#'
#' \item{Rr}{An p-vector of \code{t(X)} times fitted values, \code{X
#'   \%*\% colSums(alpha*mu)}.}
#'
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow = n,ncol = p)
#' X = scale(X,center = TRUE,scale = TRUE)
#' y = drop(X %*% beta + rnorm(n))
#'
#' input_ss <- compute_ss(X,y,standardize = TRUE)
#' ss   <- univariate_regression(X,y)
#' R    <- with(input_ss,cov2cor(XtX))
#' zhat <- with(ss,betahat/sebetahat)
#' res  <- susie_rss(zhat,R,L = 10)
#'
#' @export
#'
susie_rss = function (z, R, maf = NULL, maf_thresh = 0, z_ld_weight = 0,
                      L = 10, prior_variance = 50, residual_variance = NULL,
                      prior_weights = NULL, null_weight = NULL,
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      estimate_prior_method = c("optim", "EM", "simple"),
                      check_null_threshold = 0, prior_tol = 1e-9,
                      max_iter = 100, s_init = NULL, intercept_value = 0,
                      coverage = 0.95, min_abs_corr = 0.5,
                      tol = 1e-03, verbose = FALSE, track_fit = FALSE,
                      check_R = FALSE, r_tol = 1e-08, refine = FALSE) {

  # Check input R.
  if (nrow(R) != length(z))
    stop(paste0("The dimension of correlation matrix (", nrow(R)," by ",
                ncol(R),") does not agree with expected (",length(z)," by ",
                length(z),")"))
  if (!is_symmetric_matrix(R))
    stop("R is not a symmetric matrix")
  if (!(is.double(R) & is.matrix(R)) & !inherits(R,"CsparseMatrix"))
    stop("Input R must be a double-precision matrix, or a sparse matrix")

  # MAF filter.
  if (!is.null(maf)) {
    if (length(maf) != length(z))
      stop(paste0("The length of maf does not agree with expected ",length(z)))
    id = which(maf > maf_thresh)
    R = R[id,id]
    z = z[id]
  }

  if (any(is.infinite(z)))
    stop("z contains infinite value")

  # Check for NAs in R.
  if (any(is.na(R)))
    stop("R matrix contains missing values")

  # Replace NAs in z with zeros.
  if (any(is.na(z))) {
    warning("NA values in z-scores are replaced with 0")
    z[is.na(z)] = 0
  }

  if (check_R && any(attr(R,"eigen")$values < -r_tol)){
    semi_pd = check_semi_pd(R,r_tol)
    if (!semi_pd$status){
      stop(paste0("The correlation matrix (",nrow(R)," by ",ncol(R),
                  ") is not a positive semidefinite matrix. The smallest ",
                  "eigenvalue is ",min(semi_pd$eigenvalues),"."))
    }
  }

  # Modify R as needed.
  if (z_ld_weight > 0) {
    warning('From version 0.11.0, the non-zero z_ld_weight is no longer recommended.')
    R = muffled_cov2cor((1-z_ld_weight)*R + z_ld_weight * tcrossprod(z))
    R = (R + t(R))/2
  }

  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(R) * (1 - null_weight),ncol(R)),null_weight)
    else
      prior_weights = c(prior_weights * (1 - null_weight),null_weight)
    R = cbind(rbind(R,0),0)
    z = c(z,0)
  }

  if(estimate_residual_variance){
    warning("SuSiE-RSS no longer estimates residual variance, since we found it didn't help.")
    estimate_residual_variance = FALSE
  }

  if (!is.null(residual_variance) &&
      (residual_variance > 1 | residual_variance < 0))
    stop("Residual variance should be a scalar between 0 and 1")
  if (is.null(residual_variance))
    residual_variance = 1

  s = susie_suff_stat(XtX = R, Xty = z, n = length(z), yty = length(z)-1,
                      L = L, scaled_prior_variance = prior_variance,
                      residual_variance = residual_variance,
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = estimate_prior_variance,
                      estimate_prior_method = estimate_prior_method,
                      check_null_threshold = check_null_threshold, prior_tol = prior_tol,
                      r_tol = r_tol, prior_weights = prior_weights,
                      null_weight = NULL, standardize = FALSE,
                      max_iter = max_iter, s_init = s_init, intercept_value = intercept_value,
                      coverage = coverage, min_abs_corr = min_abs_corr,
                      tol = tol, verbose = verbose, track_fit = track_fit, check_input = FALSE,
                      refine = refine)
  s$fitted = s$Xtfitted
  s$Rr = s$XtXr
  s$Xtfitted = s$XtXr = NULL
  return(s)
}

# Performs sum of single-effect (SuSiE) linear regression with z
# scores (with lambda). The summary data required are the p by p
# correlation matrix R, the p vector z. The summary stats should come
# from the same individuals. This function fits the regression model z
# = sum_l Rb_l + e, where e is N(0,residual_variance * R + lambda I)
# and the sum_l b_l is a p vector of effects to be estimated. The
# assumption is that each b_l has exactly one non-zero element, with
# all elements equally likely to be non-zero. The prior on the
# non-zero element is N(0,var = prior_variance).
#
#' @importFrom stats optimize
susie_rss_lambda = function(z, R, maf = NULL, maf_thresh = 0,
                            L = 10, lambda = 0,
                            prior_variance = 50, residual_variance = NULL,
                            r_tol = 1e-08, prior_weights = NULL,
                            null_weight = NULL,
                            estimate_residual_variance = TRUE,
                            estimate_prior_variance = TRUE,
                            estimate_prior_method = c("optim", "EM", "simple"),
                            check_null_threshold = 0, prior_tol = 1e-9,
                            max_iter = 100, s_init = NULL, intercept_value = 0,
                            coverage = 0.95, min_abs_corr = 0.5,
                            tol = 1e-3, verbose = FALSE, track_fit = FALSE,
                            check_R = TRUE, check_z = FALSE) {

  # Check input R.
  if (nrow(R) != length(z))
    stop(paste0("The dimension of correlation matrix (",nrow(R)," by ",
                ncol(R),") does not agree with expected (",length(z)," by ",
                length(z),")"))
  if (!is_symmetric_matrix(R))
    stop("R is not a symmetric matrix")
  if (!(is.double(R) & is.matrix(R)) & !inherits(R,"CsparseMatrix"))
    stop("Input R must be a double-precision matrix or a sparse matrix")

  # MAF filter.
  if (!is.null(maf)) {
    if (length(maf) != length(z))
      stop(paste0("The length of maf does not agree with expected ",length(z)))
    id = which(maf > maf_thresh)
    R = R[id,id]
    z = z[id]
  }

  if (any(is.infinite(z)))
    stop("z contains infinite values")

  # Check for NAs in R.
  if (any(is.na(R)))
    stop("R matrix contains missing values")

  # Replace NAs in z with zero.
  if (any(is.na(z))) {
    warning("NA values in z-scores are replaced with 0")
    z[is.na(z)] = 0
  }

  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(R)*(1-null_weight),ncol(R)),null_weight)
    else
      prior_weights = c(prior_weights * (1-null_weight),null_weight)
    R = cbind(rbind(R,0),0)
    z = c(z,0)
  }

  # Eigen decomposition for R, fileter on eigenvalues.
  p = ncol(R)
  attr(R,"eigen") = eigen(R,symmetric = TRUE)
  if (check_R && any(attr(R,"eigen")$values < -r_tol))
    stop(paste0("The correlation matrix (",nrow(R)," by ",ncol(R),
                "is not a positive semidefinite matrix. ",
                "The smallest eigenvalue is ",min(attr(R,"eigen")$values),
                ". You can bypass this by \"check_R = FALSE\" which instead ",
                "sets negative eigenvalues to 0 to allow for continued ",
                "computations."))

  # Check whether z in space spanned by the non-zero eigenvectors of R.
  if (check_z) {
    proj = check_projection(R,z)
    if (!proj$status)
      warning("Input z does not lie in the space of non-zero eigenvectors ",
              "of R.")
    else
      message("Input z is in space spanned by the non-zero eigenvectors of ",
              "R.")
  }
  R = set_R_attributes(R,r_tol)

  if (lambda == 'estimate'){
    colspace = which(attr(R,"eigen")$values > 0)
    if(length(colspace) == length(z)){
      lambda = 0
    }else{
      znull = crossprod(attr(R,"eigen")$vectors[,-colspace], z) # U2^T z
      lambda = sum(znull^2)/length(znull)
    }
  }

  # Initialize susie fit.
  s = init_setup_rss(p,L,prior_variance,residual_variance,prior_weights,
                     null_weight)
  if (!missing(s_init) && !is.null(s_init)) {
    if (!inherits(s_init,"susie"))
      stop("s_init should be a susie object")
    if (max(s_init$alpha) > 1 || min(s_init$alpha) < 0)
      stop("s_init$alpha has invalid values outside range [0,1]; please ",
           "check your input")
    # First, remove effects with s_init$V = 0
    s_init = susie_prune_single_effects(s_init, verbose=FALSE)
    # Then prune or expand
    s_init = susie_prune_single_effects(s_init, L, s$V, verbose)
    s = modifyList(s,s_init)
    s = init_finalize_rss(s,R = R)
  } else
    s = init_finalize_rss(s)

  s$sigma2 = s$sigma2 - lambda
  estimate_prior_method = match.arg(estimate_prior_method)

  # Intialize elbo to NA.
  elbo = rep(NA,max_iter+1)
  elbo[1] = -Inf;
  tracking = list()

  attr(R,"lambda") = lambda
  Sigma = update_Sigma(R,s$sigma2,z)  # sigma2*R + lambda*I

  for (i in 1:max_iter) {
    if (track_fit)
      tracking[[i]] = susie_slim(s)
    s = update_each_effect_rss(R,z,s,Sigma,estimate_prior_variance,
                               estimate_prior_method,check_null_threshold)
    if (verbose)
      print(paste0("before estimate sigma2 objective:",
                   get_objective_rss(R,z,s)))

    # Compute objective before updating residual variance because part
    # of the objective s$kl has already been computed under the
    # residual variance, before the update.
    elbo[i+1] = get_objective_rss(R,z,s)
    if ((elbo[i+1] - elbo[i]) < tol) {
      s$converged = TRUE
      break
    }
    if (estimate_residual_variance) {
      if (lambda == 0) {
        est_sigma2 = (1/sum(attr(R,"eigen")$values != 0))*get_ER2_rss(1,R,z,s)
        if (est_sigma2 < 0)
          stop("Estimating residual variance failed: the estimated value ",
               "is negative")
        if (est_sigma2 > 1)
          est_sigma2 = 1
      } else {
        est_sigma2 = optimize(Eloglik_rss, interval = c(1e-4, 1-lambda),
                              R = R, z = z, s = s, maximum = TRUE)$maximum
        if(Eloglik_rss(est_sigma2, R, z, s) < Eloglik_rss(1-lambda, R, z, s)){
          est_sigma2 = 1-lambda
        }
      }
      s$sigma2 = est_sigma2

      if (verbose)
        print(paste0("after estimate sigma2 objective:",
                     get_objective_rss(R,z,s)))
      Sigma = update_Sigma(R,s$sigma2,z)
    }
  }

  # Remove first (infinite) entry, and trailing NAs.
  elbo = elbo[2:(i+1)]
  s$elbo = elbo
  s$niter = i
  s$lambda = lambda

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!"))
    s$converged = FALSE
  }

  s$intercept = intercept_value
  s$fitted = s$Rz

  s$X_column_scale_factors = attr(R,"scaled:scale")

  if (track_fit)
    s$trace = tracking

  # SuSiE CS and PIP.
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    R = muffled_cov2cor(R)
    s$sets = susie_get_cs(s,coverage = coverage,Xcorr = R,
                          min_abs_corr = min_abs_corr)
    s$pip = susie_get_pip(s,prune_by_cs = FALSE,prior_tol = prior_tol)
  }
  if (!is.null(names(z))) {
    variable_names = names(z)
    if (!is.null(null_weight))
      variable_names = c("null",variable_names)
    colnames(s$alpha) = variable_names
    colnames(s$mu) = variable_names
    colnames(s$mu2) = variable_names
    colnames(s$lbf_variable) = variable_names
    names(s$pip) = variable_names
  }

  return(s)
}

update_Sigma = function (R, sigma2, z) {
  Sigma = sigma2*R + attr(R,"lambda") * diag(length(z))
  eigenS = attr(R,"eigen")
  eigenS$values = sigma2 * eigenS$values + attr(R,"lambda")

  Dinv = 1/(eigenS$values)
  Dinv[is.infinite(Dinv)] = 0
  attr(Sigma,"eigenS") = eigenS

  # Sigma^(-1) R_j = U (sigma2 D + lambda)^(-1) D U^T e_j
  attr(Sigma,"SinvRj") = eigenS$vectors %*% (Dinv * attr(R,"eigen")$values *
                           t(eigenS$vectors))

  if (attr(R,"lambda") == 0)
    attr(Sigma,"RjSinvRj") = attr(R,"d")/sigma2
  else {
    tmp = t(eigenS$vectors)
    attr(Sigma,"RjSinvRj") =
      colSums(tmp * (Dinv*(attr(R,"eigen")$values^2) * tmp))
  }

  return(Sigma)
}

