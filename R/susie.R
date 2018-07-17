#' @title Bayesian sum of single-effect (susie) linear regression of Y on X
#' @details Performs sum of single-effect (susie) linear regression of Y on X.
#' That is, this function
#' fits the regression model Y= sum_l Xb_l + e, where elements of e are iid N(0,s2) and the
#' sum_l b_l is a p vector of effects to be estimated.
#' The assumption is that each b_l has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0,var=sa2*s2).
#' @param X an n by p matrix of covariates
#' @param Y an n vector
#' @param L maximum number of non-zero effects
#' @param prior_variance the scaled prior variance (vector of length L, or scalar. In latter case gets repeated L times )
#' @param residual_variance the residual variance (defaults to variance of Y)
#' @param standardize logical flag for whether to standardize columns of X to unit variance prior to fitting.
#' Note that `prior_variance` specifies the prior on the coefficients of X after standardization (if performed).
#' If you do not standardize you may need
#' to think carefully about specifying
#' `prior_variance`. Whatever the value of standardize, the coefficients (returned from `coef`) are for X on the original input scale.
#' Any column of X that has zero variance is not standardized, but left as is.
#' @param intercept Should intercept be fitted (default=TRUE) or set to zero (FALSE)
#' @param max_iter maximum number of iterations to perform
#' @param tol convergence tolerance
#' @param estimate_residual_variance indicates whether to estimate residual variance
#' @param estimate_prior_variance indicates whether to estimate prior (currently not recommended as not working as well)
#' @param s_init a previous susie fit with which to initialize
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @return a susie fit, which is a list with some or all of the following elements\cr
#' \item{alpha}{an L by p matrix of posterior inclusion probabilites}
#' \item{mu}{an L by p matrix of posterior means (conditional on inclusion)}
#' \item{mu2}{an L by p matrix of posterior second moments (conditional on inclusion)}
#' \item{Xr}{an n vector of fitted values, equal to X times colSums(alpha*mu))}
#' \item{sigma2}{residual variance}
#' \item{sa2}{scaled prior variance; ie prior variance is sigma2*sa2}
#' \item{elbo}{vector of values of elbo achieved (objective function)}
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1] = 1
#' beta[2] = 1
#' beta[3] = 1
#' beta[4] = 1
#' X = matrix(rnorm(n*p),nrow=n,ncol=p)
#' y = X %*% beta + rnorm(n)
#' res =susie(X,y,L=10)
#' coef(res)
#' plot(y,predict(res))
#' @export
susie = function(X,Y,L=10,prior_variance=0.2,residual_variance=NULL,standardize=TRUE,intercept=TRUE,max_iter=100,tol=1e-2,estimate_residual_variance=TRUE,estimate_prior_variance = FALSE, s_init = NULL, verbose=FALSE, track_fit=FALSE){
  # Check input X.
  if (!is.double(X) || !is.matrix(X))
    stop("Input X must be a double-precision matrix")
  p = ncol(X)
  n = nrow(X)
  mean_y = mean(Y)

  if(intercept){ # center Y and X
    Y = Y-mean_y
    X = safe_colScale(X,center=TRUE, scale = FALSE)
  } else {
    attr(X,"scaled:center")=rep(0,p)
  }

  if(standardize){
    X = safe_colScale(X,center=FALSE, scale=TRUE)
  } else {
    attr(X,"scaled:scale")=rep(1,p)
  }


  # initialize susie fit
  if(!is.null(s_init)){
    if(!missing(L) || !missing(prior_variance) || !missing(residual_variance))
      stop("if provide s_init then L, sa2 and sigma2 must not be provided")
    keys = c('alpha', 'mu', 'mu2', 'sa2')
    if(!all(keys %in% names(s_init)))
      stop(paste("s_init requires all of the following attributes:", paste(keys, collapse = ', ')))
    if (!all(dim(s_init$mu) == dim(s_init$mu2)))
      stop("dimension of mu and mu2 in s_init do not match")
    if (!all(dim(s_init$mu) == dim(s_init$alpha)))
      stop("dimension of mu and alpha in s_init do not match")
    if (dim(s_init$alpha)[1] != length(s_init$sa2))
      stop("sa2 must have length of nrow of alpha in s_init")
    if (is.null(s_init$Xr)) s_init$Xr = X%*%colSums(s_init$mu*s_init$alpha)
    if (is.null(s_init$sigma2)) s_init$sigma2 = var(Y)
    # reset KL
    s_init$KL = rep(NA, nrow(s_init$alpha))
    s = s_init
  } else {

    if(is.null(residual_variance)){
      residual_variance=var(Y)
    }
    residual_variance= as.numeric(residual_variance) #avoid problems with dimension if entered as matrix


    if(length(prior_variance)==1){
      prior_variance = rep(prior_variance,L)
    }

    # Check inputs sigma and sa.
    if (length(residual_variance) != 1)
      stop("Inputs residual_variance must be scalar")
    # Check inputs sigma and sa.
    if (length(prior_variance) != L)
      stop("Inputs prior_variance must be of length 1 or L")

    # initialize susie fit
    s = list(alpha=matrix(1/p,nrow=L,ncol=p),
             mu=matrix(0,nrow=L,ncol=p),
             mu2=matrix(0,nrow=L,ncol=p),
             Xr=rep(0,n), KL=rep(NA,L),
             sigma2=residual_variance, sa2=prior_variance)
  }
  class(s) = "susie"

  #intialize elbo to NA
  elbo = rep(NA,max_iter+1)
  elbo[1] = -Inf;
  tracking = list()

  for(i in 1:max_iter){
    #s = add_null_effect(s,0)
    if (track_fit)
      tracking[[i]] = s
    s = update_each_effect(X, Y, s, estimate_prior_variance)
    if(verbose){
        print(paste0("objective:",susie_get_objective(X,Y,s)))
    }
    if(estimate_residual_variance){
      new_sigma2 = estimate_residual_variance(X,Y,s)
      #s$sa2 = (s$sa2*s$sigma2)/new_sigma2 # this is so prior variance does not change with update
      s$sigma2 = new_sigma2
      if(verbose){
        print(paste0("objective:",susie_get_objective(X,Y,s)))
      }
    }
    #s = remove_null_effects(s)

    elbo[i+1] = susie_get_objective(X,Y,s)
    if((elbo[i+1]-elbo[i])<tol) break;
  }
  elbo = elbo[1:(i+1)] #remove trailing NAs
  s$elbo <- elbo
  s$niter <- i

  if(intercept){
    s$intercept = mean_y - sum(attr(X,"scaled:center")* (colSums(s$alpha*s$mu)/attr(X,"scaled:scale")))# estimate intercept (unshrunk)
    s$fitted = s$Xr + mean_y
  } else {
    s$intercept = 0
    s$fitted = s$Xr
  }

  s$X_column_scale_factors = attr(X,"scaled:scale")
  if (track_fit)
    s$trace = tracking

  return(s)
}
