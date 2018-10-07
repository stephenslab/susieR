#' @title Bayesian sum of single-effect (susie) linear regression using summary stat
#' @details Performs sum of single-effect (susie) linear regression of Y on X when
#' only summary statistics are available. The summary data required are
#' the p by p matrix X'X, the p vector X'Y, and the sample variance of Y, or (1/n)Y'Y. The summary stats should come from the same individuals.
#' Both the columns of X and the vector Y
#' should be centered to have mean 0 before
#' computing these summary statistics; you may also want to scale each column of X and Y to have variance 1 (see examples).
#' This function fits the regression model Y = sum_l Xb_l + e, where elements of e are iid N(0,var=residual_variance) and the
#' sum_l b_l is a p vector of effects to be estimated.
#' The assumption is that each b_l has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0,var=scaled_prior_variance*var(Y)).
#' @param XtX a p by p matrix, X'X, where columns of X are centered to have mean 0
#' @param Xty a p vector, X'Y, where columns of X are centered and Y is centered to have mean 0
#' @param var_y the (sample) variance of the vector Y
#' @param n sample size
#' @param L maximum number of non-zero effects
#' @param scaled_prior_variance the scaled prior variance (vector of length L, or scalar. In latter case gets repeated L times )
#' @param residual_variance the residual variance (defaults to variance of Y)
#' @param estimate_prior_variance indicates whether to estimate prior (currently not recommended as not working as well)
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param max_iter maximum number of iterations to perform
#' @param s_init a previous susie fit with which to initialize
#' @param intercept_value a value to assign to the intercept (since the intercept cannot be estimated from centered summary data). This
#' value will be used by coef.susie() to assign an intercept value, for consistency with the non-summary-statistic version of this function \code{susie}.
#' Set to NULL if you want coef.susie() not to include an intercept term (and so only return a p vector).
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% confidence interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a confidence set.
#' @param tol convergence tolerance based on alpha
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @return a susie fit, which is a list with some or all of the following elements\cr
#' \item{alpha}{an L by p matrix of posterior inclusion probabilites}
#' \item{mu}{an L by p matrix of posterior means (conditional on inclusion)}
#' \item{mu2}{an L by p matrix of posterior second moments (conditional on inclusion)}
#' \item{XtXr}{an p vector of t(X) times fitted values, the fitted values equal to X times colSums(alpha*mu))}
#' \item{sigma2}{residual variance}
#' \item{V}{prior variance}
#'
#' @examples
#' set.seed(1)
#' n    <- 1000
#' p    <- 1000
#' beta <- rep(0,p)
#' beta[1:4] <- 1
#' X        <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' y        <- c(X %*% beta + rnorm(n))
#' input_ss <- compute_ss(X,y,standardize = TRUE)
#' res      <- with(input_ss,susie_ss(XtX,Xty,vary,n))
#' coef(res)
#'
#' @export
susie_ss = function(XtX,Xty,var_y = 1, n, L=10,scaled_prior_variance=0.2,residual_variance=NULL,estimate_prior_variance = FALSE,
                    prior_weights = NULL, max_iter=100,s_init = NULL, intercept_value=0,
                    coverage=0.95, min_abs_corr=0.5,
                    tol=1e-4, verbose=FALSE, track_fit = FALSE){
  # Check input XtX.
  if (!(is.double(XtX) & is.matrix(XtX)) & !inherits(XtX,"CsparseMatrix"))
    stop("Input X must be a double-precision matrix, or a sparse matrix.")
  p = ncol(XtX)

  # initialize susie fit
  if(!is.null(s_init)){
    if(!missing(L) || !missing(scaled_prior_variance) || !missing(residual_variance))
      stop("if provide s_init then L, scaled_prior_variance and residual_variance must not be provided")
    s = s_init
  } else {

    if(is.null(residual_variance)){
      residual_variance = var_y
    }
    residual_variance= as.numeric(residual_variance) #avoid problems with dimension if entered as matrix


    if(length(scaled_prior_variance)==1){
      scaled_prior_variance = rep(scaled_prior_variance,L)
    }

    if(is.null(prior_weights)){
      prior_weights = rep(1/p, p)
    }

    # Check inputs sigma and sa.
    if (length(residual_variance) != 1)
      stop("Inputs residual_variance must be scalar")
    # Check inputs sigma and sa.
    if (length(scaled_prior_variance) != L)
      stop("Inputs scaled_prior_variance must be of length 1 or L")

    #initialize susie fit
    s = list(alpha=matrix(1/p,nrow=L,ncol=p), mu=matrix(0,nrow=L,ncol=p),
             mu2 = matrix(0,nrow=L,ncol=p), XtXr=rep(0,p), sigma2= residual_variance,
             V= scaled_prior_variance*var_y, KL = rep(NA,L),
             pi = prior_weights)
    class(s) <- "susie"
  }

  #intialize elbo to NA
  elbo = rep(NA,max_iter+1)
  elbo[1] = -Inf;
  tracking = list()

  alpha_new = s$alpha
  for(i in 1:max_iter){
    if (track_fit)
      tracking[[i]] = s
    alpha_old = alpha_new
    s = update_each_effect_ss(XtX, Xty, s, estimate_prior_variance)
    alpha_new = s$alpha

    elbo[i+1] = susie_get_objective_ss(XtX, Xty, s, var_y, n)

    if(verbose){
      print(paste0("objective:",elbo[i+1]))
    }
    if(max(abs(alpha_new - alpha_old)) < tol) break;

    # if((elbo[i+1]-elbo[i])<tol) break;
  }
  elbo = elbo[1:(i+1)] #remove trailing NAs
  s$elbo <- elbo
  s$niter <- i

  s$intercept = intercept_value
  s$Xtfitted = s$XtXr

  s$X_column_scale_factors = 1

  if (track_fit)
    s$trace = tracking

  ## SuSiE CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = susie_get_CS(s, coverage=coverage, Xcorr=cov2cor(XtX), min_abs_corr=min_abs_corr)
    s$pip = susie_get_PIP(s,s$sets$cs_index)
  }

  return(s)
}
