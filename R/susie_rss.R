#' @title Sum of Single Effects (SuSiE) Regression using Summary Statistics
#'
#' @description \code{susie_rss} performs variable selection under a
#'   sparse Bayesian multiple linear regression of \eqn{Y} on \eqn{X}
#'   using the z-scores from standard univariate regression
#'   of \eqn{Y} on each column of \eqn{X}, an estimate, \eqn{R}, of
#'   the correlation matrix for the columns of \eqn{X}, and optionally,
#'   \emph{but strongly recommended}, the sample size n. See
#'   \dQuote{Details} for other ways to call \code{susie_rss}
#'
#' @details In some applications, particularly genetic applications,
#' it is desired to fit a regression model (\eqn{Y = Xb + E} say,
#' which we refer to as "the original regression model" or ORM)
#' without access to the actual values of \eqn{Y} and \eqn{X}, but
#' given only some summary statistics. \code{susie_rss} assumes
#' availability of z-scores from standard univariate regression of
#' \eqn{Y} on each column of \eqn{X}, and an estimate, \eqn{R}, of the
#' correlation matrix for the columns of \eqn{X} (in genetic
#' applications \eqn{R} is sometimes called the \dQuote{LD matrix}).
#'
#' With the inputs \code{z}, \code{R} and sample size \code{n},
#' \code{susie_rss} computes PVE-adjusted z-scores \code{z_tilde}, and
#' calls \code{susie_suff_stat} with \code{XtX = (n-1)R}, \code{Xty =
#' } \eqn{\sqrt{n-1} z_tilde}, \code{yty = n-1}, \code{n = n}. The
#' output effect estimates are on the scale of \eqn{b} in the ORM with
#' \emph{standardized} \eqn{X} and \eqn{y}. When the LD matrix
#' \code{R} and the z-scores \code{z} are computed using the same
#' matrix \eqn{X}, the results from \code{susie_rss} are same as, or
#' very similar to, \code{susie} with \emph{standardized} \eqn{X} and
#' \eqn{y}.
#'
#' Alternatively, if the user provides \code{n}, \code{bhat} (the
#' univariate OLS estimates from regressing \eqn{y} on each column of
#' \eqn{X}), \code{shat} (the standard errors from these OLS
#' regressions), the in-sample correlation matrix \eqn{R =
#' cov2cor(crossprod(X))}, and the variance of \eqn{y}, the results
#' from \code{susie_rss} are same as \code{susie} with \eqn{X} and
#' \eqn{y}. The effect estimates are on the same scale as the
#' coefficients \eqn{b} in the ORM with \eqn{X} and \eqn{y}.
#'
#' In rare cases in which the sample size, \eqn{n}, is unknown,
#' \code{susie_rss} calls \code{susie_suff_stat} with \code{XtX = R}
#' and \code{Xty = z}, and with \code{residual_variance = 1}. The
#' underlying assumption of performing the analysis in this way is
#' that the sample size is large (\emph{i.e.}, infinity), and/or the
#' effects are small. More formally, this combines the log-likelihood
#' for the noncentrality parameters, \eqn{\tilde{b} = \sqrt{n} b},
#' \deqn{L(\tilde{b}; z, R) = -(\tilde{b}'R\tilde{b} -
#' 2z'\tilde{b})/2,} with the \dQuote{susie prior} on
#' \eqn{\tilde{b}}; see \code{\link{susie}} and Wang \emph{et al}
#' (2020) for details. In this case, the effect estimates returned by
#' \code{susie_rss} are on the noncentrality parameter scale.
#'
#' The \code{estimate_residual_variance} setting is \code{FALSE} by
#' default, which is recommended when the LD matrix is estimated from
#' a reference panel. When the LD matrix \code{R} and the summary
#' statistics \code{z} (or \code{bhat}, \code{shat}) are computed
#' using the same matrix \eqn{X}, we recommend setting
#' \code{estimate_residual_variance = TRUE}.
#'
#' @param z p-vector of z-scores.
#'
#' @param R p x p correlation matrix.
#'
#' @param n The sample size.
#'
#' @param bhat Alternative summary data giving the estimated effects
#'   (a vector of length p). This, together with \code{shat}, may be
#'   provided instead of \code{z}.
#'
#' @param shat Alternative summary data giving the standard errors of
#'   the estimated effects (a vector of length p). This, together with
#'   \code{bhat}, may be provided instead of \code{z}.
#'
#' @param var_y The sample variance of y, defined as \eqn{y'y/(n-1)}.
#'   When the sample variance is not provided, the coefficients
#'   (returned from \code{coef}) are computed on the
#'   \dQuote{standardized} X, y scale.
#'
#' @param z_ld_weight This parameter is included for backwards
#'   compatibility with previous versions of the function, but it is no
#'   longer recommended to set this to a non-zero value. When
#'   \code{z_ld_weight > 0}, the matrix \code{R} is adjusted to be
#'   \code{cov2cor((1-w)*R + w*tcrossprod(z))}, where \code{w =
#'   z_ld_weight}.
#'
#' @param prior_variance The prior variance(s) for the non-zero
#'   noncentrality parameterss \eqn{\tilde{b}_l}. It is either a scalar,
#'   or a vector of length L. When the \code{susie_suff_stat} option
#'   \code{estimate_prior_variance} is set to \code{TRUE} (which is
#'   highly recommended) this simply provides an initial value for the
#'   prior variance. The default value of 50 is simply intended to be a
#'   large initial value. Note this setting is only relevant when
#'   \code{n} is unknown. If \code{n} is known, the relevant option is
#'   \code{scaled_prior_variance} in \code{\link{susie_suff_stat}}.
#'
#' @param estimate_residual_variance The default is FALSE, the
#'   residual variance is fixed to 1 or variance of y. If the in-sample
#'   LD matrix is provided, we recommend setting
#'   \code{estimate_residual_variance = TRUE}.
#'
#' @param check_prior When \code{check_prior = TRUE}, it checks if the
#'   estimated prior variance becomes unreasonably large (comparing with
#'   100 * max(abs(z))^2).
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
#' new approach to variable selection in regression, with application
#' to genetic fine-mapping. \emph{Journal of the Royal Statistical
#' Society, Series B} \bold{82}, 1273-1300 \doi{10.1101/501114}.
#'
#' Y. Zou, P. Carbonetto, G. Wang, G and M. Stephens
#' (2022). Fine-mapping from summary data with the \dQuote{Sum of
#' Single Effects} model. \emph{PLoS Genetics} \bold{18},
#' e1010299. \doi{10.1371/journal.pgen.1010299}.
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
#' res  = susie_rss(zhat,R, n=n)
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
susie_rss = function (z, R, n, bhat, shat, var_y,
                      z_ld_weight = 0,
                      estimate_residual_variance = FALSE,
                      prior_variance = 50,
                      check_prior = TRUE, ...) {

  if (estimate_residual_variance)
    warning_message("For estimate_residual_variance = TRUE, please check ",
                    "that R is the \"in-sample\" LD matrix; that is, the ",
                    "correlation matrix obtained using the exact same data ",
                    "matrix X that was used for the other summary ",
                    "statistics. Also note, when covariates are included in ",
                    "the univariate regressions that produced the summary ",
                    "statistics, also consider removing these effects from ",
                    "X before computing R.",style = "hint")

  # Check input R.
  if (missing(z))
    p <- length(bhat)
  else
    p <- length(z)
  if (nrow(R) != p)
      stop(paste0("The dimension of R (",nrow(R)," x ",ncol(R),") does not ",
                  "agree with expected (",p," x ",p,")"))

  # Check input n.
  if (!missing(n))
    if (n <= 1)
      stop("n must be greater than 1")

  # Check inputs z, bhat and shat. Note that bhat is no longer used
  # after this step.
  if (sum(c(missing(z),missing(bhat) || missing(shat))) != 1)
    stop("Please provide either z or (bhat, shat), but not both")
  if (missing(z)) {
    if (length(shat) == 1)
      shat = rep(shat,length(bhat))
    if (length(bhat) != length(shat))
      stop("The lengths of bhat and shat do not agree")
    if (anyNA(bhat) || anyNA(shat))
      stop("bhat, shat cannot have missing values")
    if (any(shat <= 0))
      stop("shat cannot have zero or negative elements")
    z = bhat/shat
  }
  if (length(z) < 1)
    stop("Input vector z should have at least one element")
  z[is.na(z)] = 0

  # When n is provided, compute the PVE-adjusted z-scores.
  if (!missing(n)) {
    adj = (n-1)/(z^2 + n - 2)
    z   = sqrt(adj) * z
  }

  # Modify R by z_ld_weight; this modification was designed to ensure
  # the column space of R contained z, but susie_suff_stat does not
  # require this, and is no longer recommended.
  if (z_ld_weight > 0) {
    warning_message("As of version 0.11.0, use of non-zero z_ld_weight is no longer ",
            "recommended")
    R = muffled_cov2cor((1-z_ld_weight)*R + z_ld_weight*tcrossprod(z))
    R = (R + t(R))/2
  }

  # Call susie_suff_stat. We call susie_suff_stat in two different
  # ways depending on whether n is provided.
  if (missing(n)) {

    # The sample size (n) is not provided, so use unadjusted z-scores.
    # The choice of n=2, yty=1 is mostly arbitrary except in that it
    # ensures var(y) = yty/(n-1) = 1, and because of this
    # scaled_prior_variance = prior_variance.
    warning_message("Providing the sample size (n), or even a rough estimate of n, ",
            "is highly recommended. Without n, the implicit assumption is ",
            "n is large (Inf) and the effect sizes are small (close to zero).")
    s = susie_suff_stat(XtX = R,Xty = z,n = 2,yty = 1,
                        scaled_prior_variance = prior_variance,
                        estimate_residual_variance = estimate_residual_variance,
                        standardize = FALSE,check_prior = check_prior,...)
  } else {

    # The sample size (n) is provided, so use PVE-adjusted z-scores.
    if (!missing(shat) & !missing(var_y)) {

      # var_y, shat (and bhat) are provided, so the effects are on the
      # *original scale*.
      XtXdiag = var_y * adj/(shat^2)
      XtX = t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
      XtX = (XtX + t(XtX))/2
      Xty = z * sqrt(adj) * var_y / shat
    } else {

      # The effects are on the *standardized* X, y scale.
      XtX = (n-1)*R
      Xty = sqrt(n-1)*z
      var_y = 1
    }
    s = susie_suff_stat(XtX = XtX,Xty = Xty,n = n,yty = (n-1)*var_y,
                        estimate_residual_variance = estimate_residual_variance,
                        check_prior = check_prior,...)
  }
  return(s)
}

