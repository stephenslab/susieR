#' @title Attempt at Automating SuSiE for Hard Problems
#'
#' @description \code{susie_auto} is an attempt to automate reliable
#'   running of susie even on hard problems. It implements a three-stage
#'   strategy for each L: first, fit susie with very small residual
#'   error; next, estimate residual error; finally, estimate the prior
#'   variance. If the last step estimates some prior variances to be
#'   zero, stop. Otherwise, double L, and repeat. Initial runs are
#'   performed with relaxed tolerance; the final run is performed using
#'   the default susie tolerance.
#'
#' @param X An n by p matrix of covariates.
#'
#' @param y The observed responses, a vector of length n.
#'
#' @param L_init The initial value of L.
#'
#' @param L_max The largest value of L to consider.
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#'   and a summary of the optimization settings, are printed to the
#'   console.
#'
#' @param init_tol The tolerance to passed to \code{susie} during
#'   early runs (set large to shorten the initial runs).
#'
#' @param standardize If \code{standardize = TRUE}, standardize the
#'   columns of X to unit variance prior to fitting. Note that
#'   \code{scaled_prior_variance} specifies the prior on the
#'   coefficients of X \emph{after} standardization (if it is
#'   performed). If you do not standardize, you may need to think more
#'   carefully about specifying \code{scaled_prior_variance}. Whatever
#'   your choice, the coefficients returned by \code{coef} are given for
#'   \code{X} on the original input scale. Any column of \code{X} that
#'   has zero variance is not standardized.
#'
#' @param intercept If \code{intercept = TRUE}, the intercept is
#'   fitted; it \code{intercept = FALSE}, the intercept is set to
#'   zero. Setting \code{intercept = FALSE} is generally not
#'   recommended.
#'
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @param tol A small, non-negative number specifying the convergence
#'   tolerance for the IBSS fitting procedure. The fitting procedure
#'   will halt when the difference in the variational lower bound, or
#'   \dQuote{ELBO} (the objective function to be maximized), is
#'   less than \code{tol}.
#'
#' @param \dots Additional arguments passed to \code{\link{susie}}.
#'
#' @return See \code{\link{susie}} for a description of return values.
#'
#' @seealso \code{\link{susie}}
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
#' res = susie_auto(X,y)
#' plot(beta,coef(res)[-1])
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#' plot(y,predict(res))
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#'
#' @importFrom stats sd
#'
#' @export
#'
susie_auto = function (X, y, L_init = 1, L_max = 512, verbose = FALSE,
                       init_tol = 1, standardize = TRUE, intercept = TRUE,
                       max_iter = 100,tol = 1e-2, ...) {
  L = L_init
  if (verbose)
    message(paste0("Trying L=",L))
  s.0 = susie(X,y,L = L,residual_variance = 0.01*sd(y)^2,tol = init_tol,
              scaled_prior_variance = 1,estimate_residual_variance = FALSE,
              estimate_prior_variance = FALSE,standardize = standardize,
              intercept = intercept,max_iter = max_iter,...)
  s.1 = susie(X,y,s_init = s.0,tol = init_tol,
              estimate_residual_variance = TRUE,
              estimate_prior_variance = FALSE,
              standardize = standardize,intercept = intercept,
              max_iter = max_iter,...)
  s.2 = susie(X,y,s_init = s.1,tol = init_tol,
              estimate_residual_variance = TRUE,
              estimate_prior_variance = TRUE,
              standardize = standardize,intercept = intercept,
              max_iter = max_iter,...)

  # We call it converged---i.e., L is "big enough"---if there are any
  # prior variances set to zero.
  converged = !all(s.2$V > 0)
  while (!converged & (L <= L_max)) {
    for (i in 1:L) {
      s.2 = add_null_effect(s.2,1) # Add in L more effects.
      s.2$sigma2 = 0.01*sd(y)^2    # Set residual variance to be small
                                   # again for next iteration.
    }
    L = 2*L
    if (verbose)
      message(paste0("Trying L=",L))
    s.0 = susie(X,y,s_init = s.2,tol = init_tol,
                estimate_residual_variance = FALSE,
                estimate_prior_variance = FALSE,
                standardize = standardize,intercept = intercept,
                max_iter = max_iter,...)
    s.1 = susie(X,y,s_init = s.0,tol = init_tol,
                estimate_residual_variance = TRUE,
                estimate_prior_variance = FALSE,
                standardize = standardize,intercept = intercept,
                max_iter = max_iter,...)
    s.2 = susie(X,y,s_init = s.1,tol = init_tol,
                estimate_residual_variance = TRUE,
                estimate_prior_variance = TRUE,
                standardize = standardize,intercept = intercept,
                max_iter = max_iter,...)

    # We call it converged---i.e., L is "big enough"---if there are
    # any prior variances set to zero.
    converged = !all(s.2$V > 0)
  }

  # Final run at default tolerance to improve fit.
  s.2 = susie(X,y,s_init = s.2,estimate_residual_variance = TRUE,
              estimate_prior_variance = TRUE,tol = tol,
              standardize = standardize,intercept = intercept,
              max_iter = max_iter,...)
  return(s.2)
}
