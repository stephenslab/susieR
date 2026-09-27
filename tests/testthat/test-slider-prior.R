test_that("finite-slider SER agrees with direct Gaussian mixture calculations", {
  set.seed(881)
  X <- genotypes(42,4)
  y <- .6*X[,2]-.4*(X[,2]==1)+rnorm(nrow(X),sd=.8)
  grid <- c(-1,-.5,0,.5,1); w <- c(.1,.2,.4,0,.3)
  pi <- c(.1,.2,0,.7); V <- .6; sigma2 <- .64
  for(intercept in c(FALSE,TRUE)) for(standardize in c(FALSE,TRUE)) {
    fit <- quiet(susieRSlidePrior::susie,X,y,L=1,delta_grid=grid,delta_prior=w,
      prior_weights=pi,scaled_prior_variance=V/var(y),residual_variance=sigma2,
      estimate_prior_variance=FALSE,estimate_residual_variance=FALSE,
      min_obs=0,intercept=intercept,standardize=standardize,coverage=NULL,tol=1e-10)
    r <- y-if(intercept) mean(y) else 0
    scale <- if(standardize) apply(X,2,sd) else rep(1,ncol(X))
    bf <- m <- v <- matrix(0,ncol(X),length(grid))
    for(j in seq_len(ncol(X))) for(k in seq_along(grid)) {
      z <- X[,j]+grid[k]*(X[,j]==1)
      z <- (z-if(intercept) mean(z) else 0)/scale[j]
      covariance <- diag(sigma2,length(y))+V*tcrossprod(z)
      inverse_r <- solve(covariance,r)
      bf[j,k] <- -.5*as.numeric(determinant(covariance,logarithm=TRUE)$modulus)+
        length(y)/2*log(sigma2)-.5*sum(r*inverse_r)+sum(r*r)/(2*sigma2)
      m[j,k] <- V*sum(z*inverse_r)
      v[j,k] <- V-V^2*sum(z*solve(covariance,z))
    }
    evidence <- rowSums(exp(bf)*rep(w,each=ncol(X)))
    rho <- exp(bf)*outer(pi,w); rho <- rho/sum(rho)
    conditional <- exp(bf)*rep(w,each=ncol(X))/evidence
    expect_equal(unname(fit$alpha_delta[1,,]),rho,tolerance=1e-10)
    expect_equal(unname(fit$lbf_variable[1,]),log(evidence),tolerance=1e-10)
    expect_equal(fit$lbf,log(sum(pi*evidence)),tolerance=1e-10)
    expect_equal(unname(fit$mu[1,]),rowSums(conditional*m),tolerance=1e-10)
    expect_equal(unname(fit$mu2[1,]),rowSums(conditional*(v+m*m)),tolerance=1e-10)
    expect_equal(unname(fit$mu_delta[1,]),rowSums(conditional*m*rep(grid,each=ncol(X))),
                 tolerance=1e-10)
    expect_equal(fit$delta_prior,w,tolerance=1e-15)
    expect_true(all(fit$alpha_delta[,3,]==0))
    expect_true(all(fit$alpha_delta[,,4]==0))
    expect_equal(predict(fit,newx=X),fit$fitted,tolerance=1e-10)
  }
})

test_that("17-point IBSS equals an explicitly expanded Gaussian design", {
  set.seed(1234); X <- genotypes(220,5)
  y <- .8*(X[,1]-.5*(X[,1]==1))-.6*(X[,4]+.5*(X[,4]==1))+rnorm(nrow(X),sd=.8)
  grid <- seq(-1,1,length.out=17); p <- ncol(X)
  args <- list(L=2,scaled_prior_variance=.5/var(y),residual_variance=.64,
    estimate_prior_variance=FALSE,estimate_residual_variance=FALSE,max_iter=300,
    tol=1e-10,coverage=NULL)
  fit <- quiet(do.call,susieRSlidePrior::susie,c(list(X=X,y=y),args))
  Z <- do.call(cbind,lapply(grid,function(d) {
    raw <- X+d*(X==1)
    sweep(sweep(raw,2,colMeans(raw)),2,apply(X,2,sd),"/")
  }))
  ref <- quiet(do.call,susie_additive,c(list(X=Z,y=y-mean(y),intercept=FALSE,
                                         standardize=FALSE),args))
  expect_equal(unname(fit$fitted),unname(ref$fitted+mean(y)),tolerance=1e-8)
  for(l in 1:2) {
    expect_equal(as.vector(fit$alpha_delta[l,,]),unname(ref$alpha[l,]),tolerance=1e-8)
    expect_equal(as.vector(fit$mu_grid[l,,]),unname(ref$mu[l,]),tolerance=1e-8)
    expect_equal(as.vector(fit$mu2_grid[l,,]),unname(ref$mu2[l,]),tolerance=1e-8)
  }
  expect_equal(unname(apply(fit$alpha_delta,c(1,2),sum)),unname(fit$alpha),tolerance=1e-12)
  expect_equal(unname(fit$delta_weights),unname(apply(fit$alpha_delta,c(1,3),sum)))
  expect_equal(fit$delta_prior,rep(1/17,17))
  expect_true(all(diff(fit$elbo)>-1e-7))

  # Independently sum second moments using all expanded columns.
  means <- sapply(1:2,function(l) drop(Z %*% as.vector(fit$alpha_delta[l,,]*fit$mu_grid[l,,])))
  second <- sum(vapply(1:2,function(l)
    sum(as.vector(fit$alpha_delta[l,,]*fit$mu2_grid[l,,])*colSums(Z^2)),numeric(1)))
  erss <- sum((y-mean(y)-rowSums(means))^2)+second-sum(means^2)
  expect_equal(fit$expected_squared_residuals,erss,tolerance=1e-9)
  KL <- vapply(1:2,function(l) {
    rho <- as.vector(fit$alpha_delta[l,,]); m <- as.vector(fit$mu_grid[l,,])
    m2 <- as.vector(fit$mu2_grid[l,,]); variance <- m2-m*m
    active <- rho>0
    sum(rho[active]*log(rho[active]/(1/(p*17))))+
      .5*sum(rho[active]*(m2[active]/fit$V[l]-1+log(fit$V[l]/variance[active])))
  },numeric(1))
  expect_equal(fit$KL,KL,tolerance=1e-8)
  expect_equal(tail(fit$elbo,1),-length(y)/2*log(2*pi*fit$sigma2)-erss/(2*fit$sigma2)-sum(KL),
               tolerance=1e-8)
  expect_gt(max(abs(fit$mu_delta-fit$mu*fit$delta)),1e-4)
})

test_that("a point slider prior reduces to the corresponding fixed coding", {
  set.seed(400); X <- genotypes(160,6); y <- .7*X[,2]+rnorm(nrow(X))
  for(d in c(-1,0,.5,1)) {
    w <- rep(0,17); w[which(seq(-1,1,length.out=17)==d)] <- 1
    args <- list(X=X,y=y,L=2,estimate_prior_variance=FALSE,
      estimate_residual_variance=FALSE,coverage=NULL,tol=1e-9)
    prior <- quiet(do.call,susieRSlidePrior::susie,c(args,list(delta_prior=w)))
    fixed <- quiet(do.call,susieRSlidePrior::susie,c(args,list(delta=d)))
    for(field in c("alpha","mu","mu2","delta","mu_delta","fitted","sigma2","V"))
      expect_equal(prior[[field]],fixed[[field]],tolerance=2e-7,info=field)
  }
})

test_that("fixed grid probabilities survive variance learning and warm starts", {
  set.seed(19); X <- genotypes(260,8); y <- X[,3]+.5*(X[,3]==1)+rnorm(nrow(X),sd=.6)
  w <- rep(1,17); w[9] <- 8; w <- w/sum(w)
  for(method in c("optim","EM","simple")) {
    fit <- quiet(susieRSlidePrior::susie,X,y,L=2,delta_prior=w,
      estimate_prior_method=method,max_iter=1000,tol=1e-8,coverage=NULL)
    expect_equal(fit$delta_prior,w,tolerance=1e-15)
    expect_true(fit$converged)
    expect_true(all(diff(fit$elbo)>-1e-6))
    expect_equal(fit$sigma2,fit$expected_squared_residuals/length(y),tolerance=1e-4)
    warm <- quiet(susieRSlidePrior::susie,X,y,L=2,delta_prior=w,
      estimate_prior_method=method,model_init=fit,max_iter=1000,tol=1e-8,coverage=NULL)
    expect_equal(warm$fitted,fit$fitted,tolerance=1e-4)
  }
  original <- quiet(susieRSlidePrior::susie,X,y,L=1,estimate_prior_variance=FALSE,
    estimate_residual_variance=FALSE,delta_prior=w,coverage=NULL)
  w2 <- rep(0,17); w2[c(1,9,17)] <- c(.7,.2,.1)
  pi2 <- seq_len(ncol(X))/sum(seq_len(ncol(X)))
  warm <- quiet(susieRSlidePrior::susie,X,y,L=1,estimate_prior_variance=FALSE,
    estimate_residual_variance=FALSE,delta_prior=w2,prior_weights=pi2,
    model_init=original,coverage=NULL)
  cold <- quiet(susieRSlidePrior::susie,X,y,L=1,estimate_prior_variance=FALSE,
    estimate_residual_variance=FALSE,delta_prior=w2,prior_weights=pi2,coverage=NULL)
  expect_equal(warm$delta_prior,w2,tolerance=1e-15)
  expect_equal(warm$pi,pi2,tolerance=1e-15)
  expect_equal(warm$alpha_delta,cold$alpha_delta,tolerance=1e-12)
})

test_that("forced additive SNPs and null effects have separate prior bookkeeping", {
  set.seed(211); X <- genotypes(120,4)
  X[,2] <- rep(c(0,1),60); X[,3] <- 1
  y <- .5*X[,1]+rnorm(nrow(X))
  w <- rep(0,17); w[c(1,17)] <- c(.8,.2)
  fit <- quiet(susieRSlidePrior::susie,X,y,L=2,delta_prior=w,null_weight=.2,
    estimate_prior_variance=FALSE,estimate_residual_variance=FALSE,coverage=NULL)
  forced <- which(fit$delta_forced)
  expect_true(all(fit$alpha_delta[,forced,-9]==0))
  expect_equal(unname(fit$alpha_delta[,forced,9]),unname(fit$alpha[,forced]))
  free <- which(!fit$delta_forced & seq_len(ncol(fit$alpha))<=ncol(X))
  expect_equal(fit$delta_prior_counts,apply(fit$alpha_delta[,free,,drop=FALSE],c(1,3),sum))
  expect_equal(predict(fit,newx=X),fit$fitted,tolerance=1e-10)
  no_free <- quiet(susieRSlidePrior::susie,X[,2,drop=FALSE],y,L=1,delta_prior=w,
    estimate_prior_variance=FALSE,estimate_residual_variance=FALSE,coverage=NULL)
  expect_true(all(no_free$delta_prior_counts==0))
  zero <- quiet(susieRSlidePrior::susie,X,y,L=2,delta_prior=w,scaled_prior_variance=0,
    estimate_prior_variance=FALSE,estimate_residual_variance=FALSE,coverage=NULL)
  expect_true(all(zero$delta_prior_counts==0))
  expect_true(all(zero$mu_delta==0))
  expect_equal(unname(apply(zero$alpha_delta,c(1,2),sum)),unname(zero$alpha))
  expect_equal(predict(zero,newx=X),rep(mean(y),nrow(X)))
})

test_that("grid validation, strong signals and degenerate states are stable", {
  X <- matrix(rep(0:2,each=8),ncol=1); y <- X[,1]+rep(c(-.1,.1),12)
  for(w in list(rep(0,17),rep(-1,17),rep(NA_real_,17),rep(Inf,17),1,matrix(1,1,17)))
    expect_error(susieRSlidePrior::susie(X,y,delta_prior=w),"delta_prior")
  for(grid in list(c(-1,1),c(0,0),c(0,2),c(1,0),c(0,NA_real_)))
    expect_error(susieRSlidePrior::susie(X,y,delta_grid=grid,delta_prior=rep(1,length(grid))),"delta_grid")
  fit <- quiet(susieRSlidePrior::susie,X,y,L=1,delta_prior=rep(1e308,17),
    estimate_prior_variance=FALSE,estimate_residual_variance=FALSE,coverage=NULL)
  expect_equal(fit$delta_prior,rep(1/17,17))
  fit <- quiet(susieRSlidePrior::susie,X,1000*y,L=1,scaled_prior_variance=.1,
    residual_variance=1,estimate_prior_variance=FALSE,estimate_residual_variance=FALSE,
    coverage=NULL)
  expect_true(all(is.finite(fit$alpha_delta)))
  expect_equal(sum(fit$alpha_delta),1)
  # No heterozygotes: all codings coincide, so their posterior equals the prior.
  X <- matrix(rep(c(0,2),each=15),ncol=1); y <- X[,1]+rep(c(-.1,.1),15)
  w <- seq_len(17); w <- w/sum(w)
  fit <- quiet(susieRSlidePrior::susie,X,y,L=1,min_obs=0,delta_prior=w,
    estimate_prior_variance=FALSE,estimate_residual_variance=FALSE,coverage=NULL)
  expect_equal(as.vector(fit$alpha_delta),w,tolerance=1e-12)
})

test_that("finite-grid sampling and lfsr retain the mixture distribution", {
  set.seed(511); X <- genotypes(60,3); y <- .4*X[,1]+rnorm(nrow(X))
  fit <- quiet(susieRSlidePrior::susie,X,y,L=1,estimate_prior_variance=FALSE,
    estimate_residual_variance=FALSE,coverage=NULL)
  rho <- fit$alpha_delta[1,,]; m <- fit$mu_grid[1,,]; m2 <- fit$mu2_grid[1,,]
  positive <- rowSums(rho*pnorm(m,sd=sqrt(pmax(0,m2-m^2))))
  expected <- 1-sum(pmax(positive,fit$alpha[1,]-positive))
  expect_equal(unname(susie_get_lfsr(fit)),unname(expected),tolerance=1e-12)
  draws <- susie_get_posterior_samples(fit,20000)
  expect_named(draws,c("b","gamma","b_heterozygote"))
  expect_equal(rowMeans(draws$b),unname(coef(fit)[-1,1]),tolerance=.012)
  expect_equal(rowMeans(draws$b_heterozygote),unname(coef(fit)[-1,2]),tolerance=.012)
  expect_error(susie_get_posterior_samples(fit,0),"num_samples")
})
