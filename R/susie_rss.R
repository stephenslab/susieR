#' @title Sum of Single Effects (SuSiE) Regression using summary statistics
#'
#' @description \code{susie_rss} performs variable selection under a
#'   sparse Bayesian multiple linear regression of \eqn{Y} on \eqn{X}
#'   using only the z-scores from standard univariate regression
#'   of \eqn{Y} on each column of \eqn{X}, and an estimate \eqn{R} of
#'   the correlation matrix between columns of \eqn{X}. It does this by
#'   combining the "RSS likelihood" from Zhu and Stephens (2017) with
#'   the Sum of Single Effects" model from Wang et al (2020).
#'
#'
#' @details In some applications, particularly genetic applications,
#' it is desired to fit a regression model (\eqn{Y = X\tilde{b} + E}
#' say, which we refer to as "the original regression model" or ORM)
#' without access to the actual values of \eqn{Y} and \eqn{X}, but
#' given only some summary statistics. \code{susie_rss} assumes the
#' availability of \eqn{z} scores from standard univariate regression
#' of \eqn{Y} on each column of \eqn{X}, and an estimate \eqn{R} of
#' the correlation matrix between columns of \eqn{X} (\eqn{R} is
#' sometimes called the LD matrix in genetic applications). See Zhu
#' and Stephens (2017), and references therein, for further
#' background.
#'
#' The \code{susie_rss} function is based on the model (2.10) from
#' Zhu and Stephens, \eqn{z | R, b ~ N(Rb,R)} where \eqn{b} is a
#' vector of length p representing the effects to be estimated. The
#' effect \eqn{b_j} is simply a multiple of the coefficient
#' \eqn{\tilde{b}_j} in the ORM, and so \eqn{b_j} is non-zero if and
#' only if \eqn{\tilde{b}_j} is non-zero. In this sense the variable
#' selection problem in this model is the same as the variable
#' selection problem in the ORM, and so the credible sets and PIPs
#' computed by \code{susie_rss} can be interpreted as credible sets
#' and PIPs for the ORM. However, converting posterior estimates of
#' \eqn{b_j} to estimates of \eqn{\tilde{b}_j} would require
#' computation of the scaling factor (not done here).
#'
#' More precisely, \code{susie_rss} assumes the log-likelihood for
#' \eqn{b} is \eqn{l(b; z,R) = -0.5(b'Rb - 2z'b)}, which is equivalent
#' to model (2.10) from Zhu and Stephens if \eqn{R} is invertible, but
#' does not require \eqn{R} to be invertible. It combines this
#' likelihood with the \dQuote{susie prior} which assumes that \eqn{b
#' = \sum_{l=1}^L b_l} where each \eqn{b_l} is a vector of length p
#' with exactly one non-zero element; see \code{\link{susie}} and Wang
#' et al (2020) for details.
#'
#' In practice, this is accomplished by calling \code{susie_suff_stat}
#' with \code{XtX = R} and \code{Xty = z}, and fixing
#' \code{residual_variance = 1}. (Values for \code{n} and \code{yty}
#' are also required by \code{susie_suff_stat}. They do not affect
#' inference when the residual variance is fixed, but they do affect
#' the interpretation of \code{scaled_prior_variance}; we set
#' \code{n=2, yty=1} so that \eqn{var(y) = yty/(n-1) = 1}.) Additional
#' arguments to be passed to \code{\link{susie_suff_stat}} can be
#' provided via \code{...}.
#'
#' @param z A p-vector of z scores.
#'
#' @param R A p by p symmetric, positive semidefinite correlation
#' matrix.
#'
#' @param z_ld_weight This parameter is included for backwards
#'   compatibility with previous versions of the function, but it is no
#'   longer recommended to use a non-zero value. If \code{z_ld_weight
#'   > 0}, the matrix R used in the model is adjusted to be
#'   \code{cov2cor((1-w)*R + w*tcrossprod(z))}, where \code{w =
#'   z_ld_weight}.
#'
#' @param L Maximum number of components (nonzero coefficients) in the
#'   susie regression model. If L is larger than the number of
#'   covariates, p, L is set to p.
#'
#' @param prior_variance The prior variance(s) for the non-zero
#'   element of \eqn{b_l}. It is either a scalar, or a vector of length
#'   L. When \code{estimate_prior_variance = TRUE} (highly recommended)
#'   this simply provides an initial value for the prior variance, and
#'   the default value of 50 is simply intended to be a large initial
#'   value.
#'
#' @param estimate_prior_variance If \code{estimate_prior_variance =
#'   TRUE}, which is highly recommended, the prior variance is estimated
#'   (this is a separate parameter for each of the L effects). If
#'   provided, \code{prior_variance} is then used as an initial value
#'   for the optimization. When \code{estimate_prior_variance = FALSE}
#'   (not recommended) the prior variance for each of the L effects is
#'   determined by the value supplied to \code{prior_variance}.
#'
#' @param \dots Other parameters to be passed to
#' \code{\link{susie_suff_stat}}.
#'
#' @return A \code{"susie"} object with the following
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
#' \item{V}{Prior variance of the non-zero elements of b.}
#'
#' \item{elbo}{The value of the variational lower bound, or
#'   \dQuote{ELBO} (objective function to be maximized), achieved at
#'   each iteration of the IBSS fitting procedure.}
#'
#' \item{Rr}{A vector of length p, equal to \code{R \%*\% colSums(alpha*mu)}.}
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
#' @references
#'
#' G. Wang, A. Sarkar, P. Carbonetto and M. Stephens (2020). A simple
#'   new approach to variable selection in regression, with application
#'   to genetic fine-mapping. \emph{Journal of the Royal Statistical
#'   Society, Series B} \bold{82} \doi{10.1101/501114}.
#'
#' X. Zhu and M. Stephens (2017). Bayesian large-scale multiple regression
#'   with summary statistics from genome-wide association studies.
#'   \emph{Annals of Applied Statistics} \bold{11}, 1561-1592.
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
#' input_ss = compute_suff_stat(X,y,standardize = TRUE)
#' ss   = univariate_regression(X,y)
#' R    = with(input_ss,cov2cor(XtX))
#' zhat = with(ss,betahat/sebetahat)
#' res  = susie_rss(zhat,R,L = 10)
#'
#' @export
#'
susie_rss = function (z, R, z_ld_weight = 0, L=10, prior_variance = 50,
                      estimate_prior_variance=TRUE, ...) {

  # Check input R.
  if (nrow(R) != length(z))
    stop(paste0("The dimension of correlation matrix (", nrow(R)," by ",
                ncol(R),") does not agree with expected (",length(z)," by ",
                length(z),")"))
  if (!is_symmetric_matrix(R))
    stop("R is not a symmetric matrix")
  if (!(is.double(R) & is.matrix(R)) & !inherits(R,"CsparseMatrix"))
    stop("Input R must be a double-precision matrix, or a sparse matrix")

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

  # Modify R by z_ld_weight; this modification was designed to ensure
  # the column space of R contained z, but susie_suff_stat does not
  # require this, and it is no longer recommended
  if (z_ld_weight > 0) {
    warning("As of version 0.11.0, use of non-zero z_ld_weight is no longer ",
            "recommended")
    R = muffled_cov2cor((1-z_ld_weight)*R + z_ld_weight * tcrossprod(z))
    R = (R + t(R))/2
  }

  # The choice of n=2, yty=1 is arbitrary except in that it ensures
  # var(y) = yty/(n-1) = 1 and because of this scaled_prior_variance =
  # prior_variance.
  s = susie_suff_stat(XtX = R,Xty = z,n = 2,yty = 1,L = L,
                      scaled_prior_variance = prior_variance,
                      estimate_prior_variance = estimate_prior_variance,
                      residual_variance = 1,
                      estimate_residual_variance = FALSE,standardize = FALSE,
                      ...)

  s$Rr = s$XtXr
  s$XtXr = NULL
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

  # Eigen decomposition for R, filter on eigenvalues.
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

  if (lambda == "estimate"){
    colspace = which(attr(R,"eigen")$values > 0)
    if (length(colspace) == length(z))
      lambda = 0
    else {
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

