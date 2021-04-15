#' @rdname susie
#'
#' @param bhat A p-vector of estimated effects.
#'
#' @param shat A p-vector of standard errors.
#'
#' @param R A p by p symmetric, positive semidefinite matrix. It
#'   can be \eqn{X'X}, the covariance matrix \eqn{X'X/(n-1)}, or a
#'   correlation matrix. It should be estimated from the same samples
#'   used to compute \code{bhat} and \code{shat}. Using an out-of-sample
#'   matrix may produce unreliable results.
#'
#' @param n The sample size.
#'
#' @param var_y The sample variance of y, defined as \eqn{y'y/(n-1)}.
#'   When the sample variance cannot be provided, the coefficients
#'   (returned from \code{coef}) are computed on the "standardized" X, y
#'   scale.
#'
#' @param XtX A p by p matrix \eqn{X'X} in which the columns of X
#'   are centered to have mean zero.
#'
#' @param Xty A p-vector \eqn{X'y} in which y and the columns of X are
#'   centered to have mean zero.
#'
#' @param yty A scalar \eqn{y'y} in which y is centered to have mean
#'   zero.
#'
#' @param maf Minor allele frequency; to be used along with
#'   \code{maf_thresh} to filter input summary statistics.
#'
#' @param maf_thresh Variants having a minor allele frequency smaller
#'   than this threshold are not used.
#'
#' @param r_tol Tolerance level for eigenvalue check of positive
#'   semidefinite matrix of R.
#'
#' @param intercept_value The intercept. (The intercept cannot be
#'   estimated from centered summary data.) This setting will be used by
#'   \code{coef} to assign an intercept value, mainly for consistency
#'   with \code{susie}. Set to \code{NULL} if you want \code{coef} not
#'   to include an intercept term (and so only a p-vector is returned).
#'
#' @param check_input If \code{check_input = TRUE},
#'   \code{susie_suff_stat} performs additional checks on \code{XtX} and
#'   \code{Xty}. The checks are: (1) check that \code{XtX} is positive
#'   semidefinite; (2) check that \code{Xty} is in the space spanned by
#'   the non-zero eigenvectors of \code{XtX}.
#'
#' @export
#'
susie_suff_stat = function (bhat, shat, R, n, var_y, XtX, Xty, yty,
                            maf = NULL, maf_thresh = 0, L = 10,
                            scaled_prior_variance = 0.2,
                            residual_variance = NULL,
                            estimate_residual_variance = TRUE,
                            estimate_prior_variance = TRUE,
                            estimate_prior_method = c("optim","EM","simple"),
                            check_null_threshold = 0, prior_tol = 1e-9,
                            r_tol = 1e-08, prior_weights = NULL,
                            null_weight = NULL, standardize = TRUE,
                            max_iter = 100, s_init = NULL,
                            intercept_value = 0, coverage = 0.95,
                            min_abs_corr = 0.5, tol = 1e-3, verbose = FALSE,
                            track_fit = FALSE, check_input = FALSE, refine = FALSE) {

  # Process input estimate_prior_method.
  estimate_prior_method = match.arg(estimate_prior_method)

  if (missing(n))
    stop("n must be provided")

  # Check sufficient statistics.
  missing_bhat = c(missing(bhat), missing(shat), missing(R))
  missing_XtX = c(missing(XtX), missing(Xty), missing(yty))

  if (all(missing_bhat) & all(missing_XtX))
    stop("Please provide either all of bhat, shat, R, n, var_y or all of ",
         "XtX, Xty, yty, n")
  if (any(missing_bhat) & any(missing_XtX))
    stop("Please provide either all of bhat, shat, R, n, var_y or all of ",
         "XtX, Xty, yty, n")
  if (all(missing_bhat) & any(missing_XtX))
    stop("Please provide all of XtX, Xty, yty, n")
  if (all(missing_XtX) & any(missing_bhat))
    stop("Please provide all of bhat, shat, R, n, var_y")
  if ((!any(missing_XtX)) & (!all(missing_bhat)))
    warning("Only using information from XtX, Xty, yty, n")
  if (!any(missing_bhat)) {
    if (!all(missing_XtX))
      warning("Only using information from bhat, shat, R, n, var_y")

    # Compute XtX, Xty, yty from bhat, shat, R, n, var_y.
    if (length(shat) == 1)
      shat = rep(shat,length(bhat))
    if (length(bhat) != length(shat))
      stop("The length of bhat does not agree with length of shat")
    if (anyNA(bhat) || anyNA(shat))
      stop("The input summary statistics have missing values")
    if (any(shat == 0))
      stop("shat contains one or more zeros")

    that = bhat/shat
    that[is.na(that)] = 0
    R2 = that^2/(that^2 + n-2)
    sigma2 = (n-1)*(1-R2)/(n-2)

    # Convert any input R to correlation matrix.
    # If R has 0 colums and rows, cov2cor produces NaN and warning.
    X0 = diag(R) == 0
    R = muffled_cov2cor(R)

    # Change the columns and rows with NaN to 0.
    if (sum(X0) > 0)
      R[X0,] = R[,X0] = 0

    if (missing(var_y)) {
      XtX = (n-1)*R
      Xty = sqrt(sigma2) * sqrt(n-1) * that
      var_y = 1
    } else {
      XtXdiag = var_y * sigma2/(shat^2)
      Xty = that * var_y * sigma2/shat
      XtX = t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
    }
    yty = var_y * (n-1)
  }

  # Check input XtX.
  if (ncol(XtX) != length(Xty))
    stop(paste0("The dimension of XtX (",nrow(XtX)," by ",ncol(XtX),
                ") does not agree with expected (",length(Xty)," by ",
                length(Xty),")"))
  if (!is_symmetric_matrix(XtX))
    stop("XtX is not a symmetric matrix")

  # MAF filter.
  if (!is.null(maf)) {
    if (length(maf) != length(Xty))
      stop(paste("The length of maf does not agree with expected",length(Xty)))
    id = which(maf > maf_thresh)
    XtX = XtX[id,id]
    Xty = Xty[id]
  }

  if (any(is.infinite(Xty)))
    stop("Xty contains infinite values")
  if (!(is.double(XtX) & is.matrix(XtX)) & !inherits(XtX,"CsparseMatrix"))
    stop("Input X must be a double-precision matrix, or a sparse matrix")
  if (any(is.na(XtX)))
    stop("XtX matrix contains NAs")

  if (check_input) {

    # Check whether XtX is positive semidefinite.
    semi_pd = check_semi_pd(XtX,r_tol)
    if (!semi_pd$status)
      stop("XtX is not a positive semidefinite matrix")

    # Check whether Xty in space spanned by the non-zero eigenvectors of XtX
    proj = check_projection(semi_pd$matrix,Xty)
    if (!proj$status)
      warning("Xty does not lie in the space of the non-zero eigenvectors ",
              "of XtX")
  }

  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (is.null(prior_weights))
      prior_weights = c(rep(1/ncol(XtX)*(1-null_weight),ncol(XtX)),null_weight)
    else
      prior_weights = c(prior_weights*(1 - null_weight),null_weight)
    XtX = cbind(rbind(XtX,0),0)
    Xty = c(Xty,0)
  }

  p = ncol(XtX)

  if (standardize) {
    dXtX = diag(XtX)
    csd = sqrt(dXtX/(n-1))
    csd[csd == 0] = 1
    XtX = (1/csd) * t((1/csd) * XtX)
    Xty = (1/csd) * Xty
  } else
    csd = rep(1, length = p)
  attr(XtX,"d") = diag(XtX)
  attr(XtX,"scaled:scale") = csd

  # Initialize susie fit.
  s = init_setup(0,p,L,scaled_prior_variance,residual_variance,prior_weights,
                 null_weight,yty/(n-1),standardize)
  s$Xr = NULL
  s$XtXr = rep(0,p)

  if (!missing(s_init)&& !is.null(s_init)) {
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
    s = init_finalize(s,X = XtX)
    s$XtXr = s$Xr
    s$Xr = NULL
  } else
    s = init_finalize(s)

  # Initialize elbo to NA.
  elbo = rep(as.numeric(NA),max_iter + 1)
  elbo[1] = -Inf;
  tracking = list()

  for (i in 1:max_iter) {
    if (track_fit)
      tracking[[i]] = susie_slim(s)
    s = update_each_effect_ss(XtX,Xty,s,estimate_prior_variance,
                              estimate_prior_method,check_null_threshold)

    if (verbose)
      print(paste0("objective:",get_objective_ss(XtX,Xty,s,yty,n)))

    # Compute objective before updating residual variance because part
    # of the objective s$kl has already been computed under the
    # residual variance before the update.
    elbo[i+1] = get_objective_ss(XtX,Xty,s,yty,n)
    if ((elbo[i+1] - elbo[i]) < tol) {
      s$converged = TRUE
      break
    }
    if (estimate_residual_variance) {
      est_sigma2 = estimate_residual_variance_ss(XtX,Xty,s,yty,n)
      if (est_sigma2 < 0)
        stop("Estimating residual variance failed: the estimated value ",
             "is negative")
      s$sigma2 = est_sigma2
      if (verbose)
        print(paste0("objective:",get_objective_ss(XtX,Xty,s,yty,n)))
    }
  }
  elbo = elbo[2:(i+1)] # Remove first (infinite) entry, and trailing NAs.
  s$elbo = elbo
  s$niter = i

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!"))
    s$converged = FALSE
  }

  s$intercept = intercept_value
  s$Xtfitted = s$XtXr

  s$X_column_scale_factors = attr(XtX,"scaled:scale")

  if (track_fit)
    s$trace = tracking

  # SuSiE CS and PIP.
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    R = muffled_cov2cor(XtX)
    s$sets = susie_get_cs(s,coverage = coverage,Xcorr = R,
                          min_abs_corr = min_abs_corr)
    s$pip = susie_get_pip(s,prune_by_cs = FALSE,prior_tol = prior_tol)
  }

  if(refine){
    if(!is.null(null_weight) && null_weight!=0){
      ## if null_weight is specified
      ## we remove the extra 0 column
      XtX = XtX[1:(ncol(XtX)-1), 1:(ncol(XtX)-1)]
      Xty = Xty[1:(ncol(XtX)-1)]
    }
    conti = TRUE
    while(conti){
      m = list()
      for(cs in 1:length(s$sets$cs)){
        if(!missing(s_init) && !is.null(s_init)){
          warning('The given s_init is not used in refinement.')
        }
        pw = rep(1, ncol(XtX))
        pw[s$sets$cs[[cs]]] = 0
        s2 = susie_suff_stat(XtX = XtX, Xty = Xty, yty = yty, n = n, L = L,
                             prior_weights = pw, s_init = NULL,
                             scaled_prior_variance = scaled_prior_variance,
                             residual_variance = residual_variance,
                             estimate_residual_variance = estimate_residual_variance,
                             estimate_prior_variance = estimate_prior_variance,
                             estimate_prior_method = estimate_prior_method,
                             check_null_threshold = check_null_threshold, prior_tol = prior_tol,
                             r_tol = r_tol, max_iter = max_iter,
                             null_weight = NULL, standardize = standardize,
                             intercept_value = intercept_value, coverage = coverage,
                             min_abs_corr = min_abs_corr, tol = tol, verbose = FALSE,
                             track_fit = FALSE, check_input = FALSE, refine = FALSE)
        sinit2 = s2[c('alpha', 'mu', 'mu2')]
        class(sinit2) = 'susie'
        s3 = susie_suff_stat(XtX = XtX, Xty = Xty, yty = yty, n = n, L = L,
                             prior_weights = NULL, s_init = sinit2,
                             scaled_prior_variance = scaled_prior_variance,
                             residual_variance = residual_variance,
                             estimate_residual_variance = estimate_residual_variance,
                             estimate_prior_variance = estimate_prior_variance,
                             estimate_prior_method = estimate_prior_method,
                             check_null_threshold = check_null_threshold, prior_tol = prior_tol,
                             r_tol = r_tol, max_iter = max_iter,
                             null_weight = NULL, standardize = standardize,
                             intercept_value = intercept_value, coverage = coverage,
                             min_abs_corr = min_abs_corr, tol = tol, verbose = FALSE,
                             track_fit = FALSE, check_input = FALSE, refine = FALSE)
        m = c(m, list(s3))
      }
      elbo = sapply(m, function(x) susie_get_objective(x))
      if((max(elbo) - susie_get_objective(s)) <= 0){
        conti=FALSE
      }else{
        s = m[[which.max(elbo)]]
      }
    }
  }

  return(s)
}

# @title Check whether A is positive semidefinite
# @param A a symmetric matrix
# @return a list of result:
# \item{matrix}{The matrix with eigen decomposition}
# \item{status}{whether A is positive semidefinite}
# \item{eigenvalues}{eigenvalues of A truncated by r_tol}
check_semi_pd = function (A, tol) {
  attr(A,"eigen") = eigen(A,symmetric = TRUE)
  eigenvalues = attr(A,"eigen")$values
  eigenvalues[abs(eigenvalues) < tol] = 0
  return(list(matrix = A,
              status = !any(eigenvalues < 0),
              eigenvalues = eigenvalues))
}

# @title Check whether b is in space spanned by the non-zero eigenvectors
#   of A
# @param A a p by p matrix
# @param b a length p vector
# @return a list of result:
# \item{status}{whether b in space spanned by the non-zero
#  eigenvectors of A}
# \item{msg}{msg gives the difference between the projected b and b if
#   status is FALSE}
check_projection = function (A, b) {
  if (is.null(attr(A,"eigen")))
    attr(A,"eigen") = eigen(A,symmetric = TRUE)
  B = attr(A,"eigen")$vectors[,attr(A,"eigen")$values > .Machine$double.eps]
  msg = all.equal(as.vector(B %*% crossprod(B,b)),as.vector(b),check.names = FALSE)
  if (!is.character(msg))
    return(list(status = TRUE,msg = NA))
  else
    return(list(status = FALSE,msg = msg))
}
