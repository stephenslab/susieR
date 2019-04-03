# @title Set default susie initialization
init_setup_rss = function(p, L, prior_variance, residual_variance, prior_weights, null_weight, varY) {
  if (!is.numeric(prior_variance) || prior_variance < 0)
    stop("Scaled prior variance should be positive number.")
  if(is.null(residual_variance))
    residual_variance = varY
  if(is.null(prior_weights))
    prior_weights = rep(1/p, p)
  if(length(prior_weights) != p)
    stop("Prior weights must have length p.")
  if (p < L) L = p
  s = list(alpha=matrix(1/p,nrow=L,ncol=p),
           mu=matrix(0,nrow=L,ncol=p),
           mu2=matrix(0,nrow=L,ncol=p),
           Rz=rep(0,p), KL=rep(NA,L),
           lbf=rep(NA,L),
           sigma2=residual_variance,
           V=prior_variance * varY,
           pi=prior_weights)
  if (is.null(null_weight)) s$null_index = 0
  else s$null_index = p
  class(s) = 'susie'
  return(s)
}

# @title Update a susie fit object in order to initialize susie model.
init_finalize_rss = function(s, R=NULL, Rz=NULL) {
  if(length(s$V) == 1)
    s$V = rep(s$V, nrow(s$alpha))
  ## check sigma2
  if (!is.numeric(s$sigma2))
    stop("Input residual variance `sigma2` must be numeric")
  ## avoid problems with dimension if input is a 1 by 1 matrix
  s$sigma2 = as.numeric(s$sigma2)
  if (length(s$sigma2) != 1)
    stop("Input residual variance `sigma2` must be a scalar")
  if (s$sigma2 <= 0)
    stop("residual variance `sigma2` must be positive (is your var(Y) zero?)")
  ## check prior variance
  if (!is.numeric(s$V))
    stop("Input prior variance must be numeric")
  if (!all(s$V >= 0))
    stop("prior variance must be non-negative")
  if (!all(dim(s$mu) == dim(s$mu2)))
    stop("dimension of mu and mu2 in input object do not match")
  if (!all(dim(s$mu) == dim(s$alpha)))
    stop("dimension of mu and alpha in input object do not match")
  if (nrow(s$alpha) != length(s$V))
    stop("Input prior variance V must have length of nrow of alpha in input object")
  ## update Xr
  if (!missing(Rz))
    s$Rz = Rz
  if (!missing(R))
    s$Rz = compute_Xb(R, colSums(s$mu*s$alpha))
  ## reset KL and lbf
  s$KL = rep(NA, nrow(s$alpha))
  s$lbf = rep(NA, nrow(s$alpha))
  class(s) = 'susie'
  return(s)
}
