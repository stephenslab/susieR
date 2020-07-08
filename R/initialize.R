#' @title Initialize a susie object using regression coefficients
#'
#' @param coef_index An L-vector for indices of nonzero effects.
#' 
#' @param coef_value An L-vector for initial estimated beta values.
#' 
#' @param p A scalar the number of variables in the data.
#' 
#' @return A list (alpha, mu, mu2) to be used by \code{susie}.
#'
#' @examples
#' # Add example(s) here.
#' 
#' @export
#' 
susie_init_coef = function (coef_index, coef_value, p) {
  L = length(coef_index)
  if (L <= 0)
    stop("Need at least one non-zero effect")
  if (!all(coef_value != 0))
    stop("Input coef_value must be non-zero for all its elements")
  if (L != length(coef_value))
    stop("Inputs coef_index and coef_value must of the same length")
  if (max(coef_index) > p)
    stop("Input coef_index exceeds the boundary of p")
  alpha = matrix(0,nrow = L,ncol = p)
  mu = matrix(0,nrow = L,ncol = p)
  for(i in 1:L){
    alpha[i,coef_index[i]] = 1
    mu[i,coef_index[i]] = coef_value[i]
  }
  return(list(alpha = alpha,mu = mu,mu2 = mu*mu))
}

# Set default susie initialization.
init_setup = function (n, p, L, scaled_prior_variance, residual_variance,
                       prior_weights, null_weight, varY, standardize) {
  if (!is.numeric(scaled_prior_variance) || scaled_prior_variance < 0)
    stop("Scaled prior variance should be positive number")
  if (scaled_prior_variance > 1 && standardize)
    stop("Scaled prior variance should be no greater than 1 when ",
         "standardize = TRUE")
  if(is.null(residual_variance))
    residual_variance = varY
  if(is.null(prior_weights))
    prior_weights = rep(1/p,p)
  else
    prior_weights = prior_weights / sum(prior_weights)
  if(length(prior_weights) != p)
    stop("Prior weights must have length p")
  if (p < L)
    L = p
  s = list(alpha  = matrix(1/p,nrow = L,ncol = p),
           mu     = matrix(0,nrow = L,ncol = p),
           mu2    = matrix(0,nrow = L,ncol = p),
           Xr     = rep(0,n),
           KL     = rep(as.numeric(NA),L),
           lbf    = rep(as.numeric(NA),L),
           sigma2 = residual_variance,
           V      = scaled_prior_variance*varY,
           pi     = prior_weights)
  if (is.null(null_weight))
    s$null_index = 0
  else
    s$null_index = p
  class(s) = "susie"
  return(s)
}

# Update a susie fit object in order to initialize susie model.
init_finalize = function (s, X = NULL, Xr = NULL) {
  if(length(s$V) == 1)
    s$V = rep(s$V,nrow(s$alpha))
  
  # Check sigma2.
  if (!is.numeric(s$sigma2))
    stop("Input residual variance sigma2 must be numeric")
  
  # Avoid problems with dimension if input is a 1 x 1 matrix.
  s$sigma2 = as.numeric(s$sigma2)
  if (length(s$sigma2) != 1)
    stop("Input residual variance sigma2 must be a scalar")
  if (s$sigma2 <= 0)
    stop("Residual variance sigma2 must be positive (is your var(Y) zero?)")
  
  # check prior variance
  if (!is.numeric(s$V))
    stop("Input prior variance must be numeric")
  if (!all(s$V >= 0))
    stop("prior variance must be non-negative")
  if (!all(dim(s$mu) == dim(s$mu2)))
    stop("dimension of mu and mu2 in input object do not match")
  if (!all(dim(s$mu) == dim(s$alpha)))
    stop("dimension of mu and alpha in input object do not match")
  if (nrow(s$alpha) != length(s$V))
    stop("Input prior variance V must have length of nrow of alpha in ",
         "input object")
  
  # Update Xr.
  if (!missing(Xr))
    s$Xr = Xr
  if (!missing(X))
    s$Xr = compute_Xb(X,colSums(s$mu * s$alpha))
  
  # Reset KL and lbf.
  s$KL = rep(as.numeric(NA),nrow(s$alpha))
  s$lbf = rep(as.numeric(NA),nrow(s$alpha))
  class(s) = "susie"
  return(s)
}
