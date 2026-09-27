# Fixed coding probabilities are inputs to this fitter. An external coding-
# prior learning loop can refit with a new delta_prior and use alpha_delta.
.slide_validate_prior <- function(grid, prior) {
  if (is.null(prior)) return(NULL)
  if (!is.numeric(grid) || !is.null(dim(grid)) || !length(grid) ||
      any(!is.finite(grid)) || any(abs(grid)>1) ||
      any(diff(grid)<=0) || !any(grid==0))
    stop("delta_grid must be a strictly increasing numeric vector in [-1,1] containing zero.")
  if (!is.numeric(prior) || !is.null(dim(prior)) || length(prior)!=length(grid) ||
      any(!is.finite(prior)) || any(prior<0) || !any(prior>0))
    stop("delta_prior must have one finite nonnegative weight per grid value and positive total mass.")
  # Scale first so that large but finite supplied weights cannot overflow.
  prior <- prior/max(prior)
  unname(prior/sum(prior))
}

.slide_grid_weights <- function(data) {
  weights <- matrix(rep(data$delta_prior,each=data$p),data$p)
  forced <- data$forced | seq_len(data$p)>data$input_p
  weights[forced,] <- 0
  weights[forced,which(data$delta_grid==0)] <- 1
  weights
}

.slide_prior_ser <- function(data,model,V) {
  ans <- slide_prior_ser_native(as.double(data$xx),as.double(data$xh),
    as.double(data$hh),as.double(model$residuals),as.double(model$hresiduals),
    as.double(V),as.double(model$sigma2),as.double(data$delta_grid),
    as.double(data$delta_prior),data$forced | seq_len(data$p)>data$input_p)
  colnames(ans$summary) <- c("delta","lbf","mu","mu2","mu_delta",
                            "mu2_delta","mu2_delta2")
  ans
}

.slide_prior_initialize <- function(data,model) {
  L <- nrow(model$alpha); K <- length(data$delta_grid)
  model$delta_grid <- data$delta_grid
  model$delta_prior <- data$delta_prior
  for (field in c("mu_delta","mu2_delta","mu2_delta2"))
    model[[field]] <- matrix(0,L,data$p)
  for (field in c("alpha_delta","mu_grid","mu2_grid"))
    model[[field]] <- array(0,c(L,data$p,K))
  weights <- .slide_grid_weights(data)
  for (l in seq_len(L)) {
    model$alpha_delta[l,,] <- weights*model$alpha[l,]
    model$delta[l,] <- drop(weights %*% data$delta_grid)
  }
  model
}

.slide_prior_warm_start <- function(data,params,model,warm) {
  matrices <- c("alpha","mu","mu2","delta","mu_delta","mu2_delta","mu2_delta2")
  arrays <- c("alpha_delta","mu_grid","mu2_grid")
  valid <- function(x,shape) is.numeric(x) && identical(dim(x),shape) && all(is.finite(x))
  if (!inherits(warm,"susie_slide") ||
      !identical(warm$delta_grid,data$delta_grid) ||
      !all(vapply(matrices,function(f) valid(warm[[f]],dim(model$alpha)),logical(1))) ||
      !all(vapply(arrays,function(f) valid(warm[[f]],dim(model$alpha_delta)),logical(1))) ||
      !isTRUE(all.equal(unname(warm$X_column_scale_factors),unname(data$scale))) ||
      !identical(unname(warm$delta_forced),unname(data$forced)) ||
      !identical(colnames(warm$alpha),colnames(data$X)) ||
      length(warm$V)!=nrow(model$alpha) || any(!is.finite(warm$V)) || any(warm$V<0) ||
      any(warm$alpha_delta<0) || any(warm$alpha<0) ||
      any(abs(rowSums(warm$alpha)-1)>1e-8) ||
      any(abs(apply(warm$alpha_delta,c(1,2),sum)-warm$alpha)>1e-8))
    stop("A prior-slider warm start requires matching grid, SNPs, dimensions, scales and count filter, with valid posterior arrays.")
  for (field in c(matrices,arrays,"V")) model[[field]] <- warm[[field]]
  if (isTRUE(params$estimate_residual_variance)) {
    if (length(warm$sigma2)!=1 || !is.finite(warm$sigma2) || warm$sigma2<=0)
      stop("Invalid residual variance in prior-slider warm start.")
    model$sigma2 <- warm$sigma2
  }
  # Keep the newly supplied SNP and slider priors, not the warm fit's priors.
  for (l in seq_len(nrow(model$alpha)))
    model$component_fitted[,l] <- .slide_product(data,model$alpha[l,]*model$mu[l,],
                                               model$alpha[l,]*model$mu_delta[l,])
  model$Xr <- rowSums(model$component_fitted)
  model
}

.slide_mu_delta <- function(model) {
  if (!is.null(model$delta_prior)) model$mu_delta else model$mu*model$delta
}

.slide_second_moment <- function(data,model,l=NULL) {
  if(!is.null(l)) {
    if(is.null(data$delta_prior)) return(model$mu2[l,]*model$slide_s[l,])
    return(pmax(model$mu2[l,]*data$xx+2*model$mu2_delta[l,]*data$xh+
                  model$mu2_delta2[l,]*data$hh,0))
  }
  if (is.null(data$delta_prior)) return(model$mu2*model$slide_s)
  pmax(sweep(model$mu2,2,data$xx,"*")+
           2*sweep(model$mu2_delta,2,data$xh,"*")+
           sweep(model$mu2_delta2,2,data$hh,"*"),0)
}

# The additive coefficient posterior is a mixture, so Gaussian approximation
# from only mu/mu2 is inappropriate for sampling and sign probabilities.
.slide_prior_lfsr <- function(fit) {
  prob <- fit$alpha_delta
  positive <- prob
  variance <- pmax(0,fit$mu2_grid-fit$mu_grid^2)
  positive[] <- stats::pnorm(fit$mu_grid,sd=sqrt(variance))
  positive[variance==0] <- as.numeric(fit$mu_grid[variance==0]>=0)
  mass <- apply(prob*positive,c(1,2),sum)
  1-rowSums(pmax(mass,fit$alpha-mass))
}

.slide_prior_samples <- function(fit,num_samples) {
  if (!is.numeric(num_samples) || length(num_samples)!=1 ||
      !is.finite(num_samples) || num_samples<1 || num_samples!=floor(num_samples))
    stop("num_samples must be a positive integer.")
  p <- ncol(fit$alpha); K <- length(fit$delta_grid)
  b <- h <- matrix(0,p,num_samples)
  for (l in which(fit$V>1e-9)) {
    prob <- matrix(fit$alpha_delta[l,,],p,K)
    chosen <- sample.int(p*K,num_samples,replace=TRUE,prob=as.vector(prob))
    j <- (chosen-1L)%%p+1L
    k <- (chosen-1L)%/%p+1L
    means <- matrix(fit$mu_grid[l,,],p,K)[chosen]
    seconds <- matrix(fit$mu2_grid[l,,],p,K)[chosen]
    beta <- stats::rnorm(num_samples,means,sqrt(pmax(0,seconds-means^2)))/
      fit$X_column_scale_factors[j]
    if (fit$null_index>0) beta[j==fit$null_index] <- 0
    idx <- cbind(j,seq_len(num_samples))
    b[idx] <- b[idx]+beta
    h[idx] <- h[idx]+beta*fit$delta_grid[k]
  }
  list(b=b,gamma=1*((b!=0)|(h!=0)),b_heterozygote=h)
}
