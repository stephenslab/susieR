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
#' @param R A p by p correlation matrix.
#'
#' @param z_ld_weight This parameter is included for backwards
#'   compatibility with previous versions of the function, but it is no
#'   longer recommended to use a non-zero value. If \code{z_ld_weight
#'   > 0}, the matrix R used in the model is adjusted to be
#'   \code{cov2cor((1-w)*R + w*tcrossprod(z))}, where \code{w =
#'   z_ld_weight}.
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
#' @param check_prior When \code{check_prior = TRUE}, it checks if the
#'   estimated prior variance becomes unreasonably large (comparing with
#'   10 * max(abs(z))^2).
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
#' G. Wang, A. Sarkar, P. Carbonetto and M. Stephens (2020). A simple
#'   new approach to variable selection in regression, with application
#'   to genetic fine-mapping. \emph{Journal of the Royal Statistical
#'   Society, Series B} \bold{82}, 1273-1300 \doi{10.1101/501114}.
#'
#'   Y. Zou, P. Carbonetto, G. Wang and M. Stephens (2021).
#'   Fine-mapping from summary data with the \dQuote{Sum of Single Effects}
#'   model. \emph{bioRxiv} \doi{10.1101/2021.11.03.467167}.
#' 
#'   X. Zhu and M. Stephens (2017). Bayesian large-scale multiple
#'   regression with summary statistics from genome-wide association
#'   studies. \emph{Annals of Applied Statistics} \bold{11}, 1561-1592.
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
#' res  = susie_rss(zhat,R)
#'
#' # Toy example illustrating behaviour susie_rss when the z-scores
#' # are mostly consistent with a non-invertible correlation matrix.
#' # Here the CS should contain both variables, and two PIPs should
#' # be nearly the same.
#' z = c(6,6.01)
#' R = matrix(1,2,2)
#' fit = susie_rss(z,R)
#' print(fit$sets$cs)
#' print(fit$pip)
#'
#' # In this second toy example, the only difference is that one
#' # z-score is much larger than the other. Here we expect that the
#' # second PIP will be much larger than the first.
#' z = c(6,7)
#' R = matrix(1,2,2)
#' fit = susie_rss(z,R)
#' print(fit$sets$cs)
#' print(fit$pip)
#'
#' @export
#'
susie_rss = function (z, R, z_ld_weight = 0, prior_variance = 50,
                      estimate_prior_variance=TRUE, check_prior=TRUE, ...) {

  # Check input R.
  if (nrow(R) != length(z))
    stop(paste0("The dimension of correlation matrix (", nrow(R)," by ",
                ncol(R),") does not agree with expected (",length(z)," by ",
                length(z),")"))

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
  s = susie_suff_stat(XtX = R,Xty = z,n = 2,yty = 1,
                      scaled_prior_variance = prior_variance,
                      estimate_prior_variance = estimate_prior_variance,
                      residual_variance = 1,estimate_residual_variance = FALSE,
                      standardize = FALSE,check_prior = check_prior,...)

  s$Rr = s$XtXr
  s$XtXr = NULL
  return(s)
}

