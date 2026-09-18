test_that("compiled evidence and moments agree with a direct Gaussian calculation", {
  set.seed(11)
  x <- rep(0:2,c(9,11,13)); y <- .4+.6*x-.25*(x==1)+rnorm(length(x),sd=.7)
  for(center in c(FALSE,TRUE)) for(scale in c(FALSE,TRUE)) for(d in c(-1,-.6,0,.6,1)) {
    result <- ser_native(x,y,delta=d,intercept=center,standardize=scale)
    z <- x+d*(x==1)
    if(center) {z <- z-mean(z); r <- y-mean(y)} else r <- y
    if(scale) z <- z/sd(x)
    covariance <- diag(.49,length(x))+.8*tcrossprod(z)
    logbf <- -.5*as.numeric(determinant(covariance,logarithm=TRUE)$modulus)+
      length(x)/2*log(.49)-.5*drop(crossprod(r,solve(covariance,r)))+sum(r*r)/(2*.49)
    variance <- 1/(1/.8+sum(z*z)/.49)
    expect_equal(result[1,2],logbf,tolerance=1e-10)
    expect_equal(result[1,3],variance*sum(z*r)/.49,tolerance=1e-10)
    expect_equal(result[1,7],variance,tolerance=1e-10)
  }
})

test_that("bounded cubic solver finds the global maximum including degenerate designs", {
  set.seed(29)
  for(k in 1:70) {
    counts <- sample(0:60,3,replace=TRUE); if(sum(counts)<3) counts[1] <- 3
    x <- rep(0:2,counts); y <- rnorm(3,sd=3)[x+1]+rnorm(length(x))
    V <- exp(runif(1,-7,3)); sigma2 <- exp(runif(1,-5,3))
    ans <- ser_native(x,y,V,sigma2)
    z <- x-mean(x); h <- (x==1)-mean(x==1); r <- y-mean(y)
    objective <- function(d) {
      s <- pmax(0,sum(z*z)+2*d*sum(z*h)+d*d*sum(h*h))
      t <- sum(z*r)+d*sum(h*r)
      .5*(V*t*t/(sigma2*(sigma2+V*s))-log1p(V*s/sigma2))
    }
    grid <- seq(-1,1,length.out=4001); j <- which.max(objective(grid))
    local <- optimize(objective,c(grid[max(1,j-1)],grid[min(length(grid),j+1)]),
                       maximum=TRUE,tol=1e-11)
    expect_gte(ans[1,2]+1e-7,max(objective(grid),local$objective))
    expect_lte(abs(ans[1,1]),1)
  }
  expect_equal(ser_native(rep(1,10),1:10)[1,c(1,2,3,4)],rep(0,4))
  expect_equal(ser_native(rep(0:2,5),1:15,V=0)[1,1:4],rep(0,4))
})

test_that("delta zero matches the actual upstream additive engine", {
  set.seed(22); X <- genotypes(); y <- .9*X[,2]-.7*X[,9]+rnorm(nrow(X),sd=.65)
  for(standardize in c(FALSE,TRUE)) for(intercept in c(FALSE,TRUE))
    for(estimate in c(FALSE,TRUE)) {
      args <- list(X=X,y=y,L=2,standardize=standardize,intercept=intercept,
                   estimate_prior_variance=estimate,estimate_residual_variance=estimate,
                   residual_variance=.65^2,max_iter=300,tol=1e-8)
      base <- quiet(do.call,susie_additive,args)
      slide <- quiet(do.call,susieSlide::susie,c(args,list(delta=0)))
      for(field in c("alpha","mu","mu2","lbf_variable","V","sigma2","fitted","pip","elbo"))
        expect_equal(slide[[field]],base[[field]],tolerance=2e-6,info=field)
      expect_equal(slide$sets,base$sets,tolerance=1e-7)
      expect_true(all(slide$delta==0))
      expect_equal(predict(slide,newx=X),slide$fitted,tolerance=1e-10)
    }
})

test_that("slider predictions, expected residuals, and conditional ELBO agree independently", {
  set.seed(33); X <- genotypes(600,22)
  y <- 1.1*(X[,3]-.7*(X[,3]==1))-.9*(X[,15]+.55*(X[,15]==1))+rnorm(600,sd=.5)
  fit <- quiet(susieSlide::susie,X,y,L=2,residual_variance=.25,
                estimate_prior_variance=FALSE,estimate_residual_variance=FALSE,tol=1e-9)
  F <- matrix(0,nrow(X),2); expected <- 0; KL <- numeric(2)
  for(l in 1:2) {
    Z <- X+sweep((X==1)*1,2,fit$delta[l,],"*")
    Z <- sweep(sweep(Z,2,colMeans(Z)),2,fit$X_column_scale_factors,"/")
    s <- colSums(Z^2); a <- fit$alpha[l,]; m <- fit$mu[l,]
    F[,l] <- drop(Z %*% (a*m))
    expected <- expected+sum(a*fit$mu2[l,]*s)
    v <- fit$V[l]/(1+fit$V[l]*s/fit$sigma2)
    KL[l] <- sum(a*log(a/(fit$pi+sqrt(.Machine$double.eps))))+
      .5*sum(a*((v+m*m)/fit$V[l]-1+log(fit$V[l]/v)))
  }
  erss <- sum((y-mean(y)-rowSums(F))^2)+expected-sum(F^2)
  expect_equal(fit$expected_squared_residuals,erss,tolerance=1e-9)
  expect_equal(fit$KL,KL,tolerance=1e-8)
  expect_equal(tail(fit$elbo,1),-length(y)/2*log(2*pi*fit$sigma2)-erss/(2*fit$sigma2)-sum(KL),tolerance=1e-8)
  expect_true(all(diff(fit$elbo)>-1e-6))
  expect_equal(predict(fit,newx=X),fit$fitted,tolerance=1e-10)
  expect_equal(rowSums(F)+mean(y),unname(fit$fitted),tolerance=1e-10)
  expect_equal(fit$delta_cs$summary,slider_cs_table(fit))
  expect_s3_class(summary(fit),"summary.susie")
  expect_lt(abs(fit$delta[which.max(fit$alpha[,3]),3]+.7),.16)
  expect_lt(abs(fit$delta[which.max(fit$alpha[,15]),15]-.55),.16)
})

test_that("count threshold includes missing classes and takes precedence over fixed delta", {
  for(category in 0:2) for(nrare in c(0,4,5)) {
    counts <- c(15,15,15); counts[category+1] <- nrare
    X <- matrix(rep(0:2,counts),ncol=1)
    y <- seq_len(nrow(X))/nrow(X)+X[,1]
    fit <- quiet(susieSlide::susie,X,y,L=1,delta=1,estimate_prior_variance=FALSE,
                  estimate_residual_variance=FALSE)
    expect_equal(unname(fit$genotype_counts[1,]),counts)
    expect_equal(unname(fit$delta_forced),nrare<5)
    expect_equal(unname(fit$delta[1,1]),if(nrare<5) 0 else 1)
    expect_equal(predict(fit,newx=X),fit$fitted,tolerance=1e-10)
  }
  X <- cbind(rep(0:2,c(20,20,2)),rep(1,42))
  fit <- quiet(susieSlide::susie,X,seq_len(42),L=1,estimate_prior_variance=FALSE)
  expect_true(all(is.finite(fit$mu)))
  expect_true(all(fit$delta==0))
})

test_that("allele reversal, sparse inputs, caching, and null columns are handled", {
  set.seed(52); X <- genotypes(450,12); y <- 1.1*(X[,2]-.6*(X[,2]==1))+rnorm(450,sd=.6)
  args <- list(y=y,L=1,residual_variance=.36,estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE,null_weight=.1,tol=1e-8)
  fit <- quiet(do.call,susieSlide::susie,c(list(X=X),args))
  flip <- quiet(do.call,susieSlide::susie,c(list(X=2-X),args))
  expect_equal(fit$pip,flip$pip,tolerance=1e-8)
  expect_equal(fit$delta,-flip$delta,tolerance=1e-7)
  expect_equal(fit$fitted,flip$fitted,tolerance=1e-9)
  sparse <- quiet(do.call,susieSlide::susie,c(list(X=Matrix::Matrix(X,sparse=TRUE)),args))
  cached <- quiet(do.call,susieSlide::susie,c(list(X=X,cache_heterozygotes=TRUE,chunk_size=3),args))
  expect_equal(sparse$fitted,fit$fitted,tolerance=1e-8)
  expect_equal(cached$fitted,fit$fitted,tolerance=1e-9)
  expect_equal(ncol(fit$delta),ncol(X)+1L)
  expect_equal(ncol(fit$delta_cs$delta),ncol(X))
  expect_identical(colnames(fit$delta_cs$delta),colnames(X))
  expect_length(fit$pip,ncol(X))
  expect_equal(nrow(coef(fit)),ncol(X)+1L)
  expect_equal(predict(fit,newx=X),fit$fitted,tolerance=1e-9)
  expect_equal(nrow(summary(fit)$vars),ncol(X))
})

test_that("learned variances, EM updates, and warm starts converge consistently", {
  set.seed(16); X <- genotypes(600,15); y <- .7*(X[,2]+.6*(X[,2]==1))+rnorm(600,sd=.7)
  for(method in c("optim","EM")) {
    fit <- quiet(susieSlide::susie,X,y,L=1,estimate_prior_method=method,max_iter=1000,tol=1e-9)
    expect_true(fit$converged)
    expect_true(all(diff(fit$elbo)>-1e-6))
    expect_equal(fit$sigma2,fit$expected_squared_residuals/length(y),tolerance=1e-5)
    warm <- quiet(susieSlide::susie,X,y,L=1,estimate_prior_method=method,
                   model_init=fit,max_iter=1000,tol=1e-9)
    expect_equal(warm$fitted,fit$fitted,tolerance=1e-5)
    expect_equal(warm$delta,fit$delta,tolerance=1e-5)
  }
})

test_that("missing outcome counts and invalid/unsupported inputs are explicit", {
  X <- matrix(rep(0:2,each=5),ncol=1); y <- seq_len(15); y[15] <- NA
  fit <- quiet(susieSlide::susie,X,y,L=1,na.rm=TRUE,delta=1,estimate_prior_variance=FALSE)
  expect_equal(unname(fit$genotype_counts[1,]),c(5,5,4))
  expect_true(fit$delta_forced)
  expect_equal(unname(fit$delta[1,1]),0)
  expect_error(susieSlide::susie(X,y),"finite")
  expect_error(susieSlide::susie(X+.1,1:15),"hard-call")
  expect_error(susieSlide::susie(X,1:15,min_obs=-1),"min_obs")
  expect_error(susieSlide::susie(X,1:15,delta=2),"delta")
  expect_error(susieSlide::susie(X,1:15,refine=TRUE),"Unsupported")
  expect_error(susieSlide::susie(X,1:15,estimate_residual_method="NIG"),"Gaussian")
})

test_that("non-singleton credible-set purity uses the fitted slider coordinates", {
  set.seed(882)
  X <- genotypes(400,4)
  X[,2] <- X[,1]
  y <- rnorm(400)
  fit <- quiet(susieSlide::susie,X,y,L=1,delta=c(-1,1,-.4,.4),
    estimate_prior_variance=FALSE,estimate_residual_variance=FALSE,
    scaled_prior_variance=.0001,coverage=.99,min_abs_corr=0,n_purity=-1)
  expect_length(fit$sets$cs,1)
  j <- fit$sets$cs[[1]]
  expect_gt(length(j),1)
  z <- X[,j,drop=FALSE]+sweep((X[,j,drop=FALSE]==1)*1,2,fit$delta[1,j],"*")
  correlations <- abs(cor(z)); correlations <- correlations[upper.tri(correlations)]
  expected <- c(min(correlations),mean(correlations),median(correlations))
  expect_equal(as.numeric(fit$sets$purity[1,]),expected,tolerance=1e-10)
  expect_identical(fit$delta_cs$summary$snp_index,j)
})
