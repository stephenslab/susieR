# Run after compare_models.R, from the package source root.
# Reproduces the core fitting settings inspected in the user's workhorse.R:
# stacked X, I(X>=1), I(X==2); standardize=FALSE; prior-variance EM;
# residual variance estimated; uniform predictor priors; L=10.
# Synthetic Gaussian outcomes are used; GTEx preprocessing is not reproduced.
library(susieSlide)
outdir <- Sys.getenv("SLIDE_COMPARISON_OUT","validation")
quiet <- function(fun,...) suppressMessages(fun(...))
examples <- readRDS(file.path(outdir,"comparison_examples.rds"))
methods <- c("Additive_EM","Slider_EM","Susie_mix_EM")
rows <- list(); fits <- list(); i <- 0L
group_pip <- function(fit,p) {
  a <- fit$alpha[,seq_len(p),drop=FALSE]+fit$alpha[,p+seq_len(p),drop=FALSE]+
    fit$alpha[,2*p+seq_len(p),drop=FALSE]
  -expm1(colSums(log1p(-pmin(a[fit$V>1e-9,,drop=FALSE],1))))
}
stack <- function(X) {
  Z <- cbind(X,(X>=1)*1,(X==2)*1)
  colnames(Z) <- paste(rep(c("A","D","R"),each=ncol(X)),rep(colnames(X),3),sep=":")
  Z
}
for(mode in c("Additive","Recessive","Partially_recessive","Dominant","Partially_dominant","Mixed")) {
  example <- examples[[paste(mode,2,"Slider",sep="_")]]
  X <- example$X; y <- example$y; p <- ncol(X)
  for(method in methods) {
    Z <- if(method=="Susie_mix_EM") stack(X) else X
    # No class would be dropped by the workhorse's >=5 column-sum rule here.
    stopifnot(all(colSums(Z)>=5))
    fun <- if(method=="Slider_EM") susieSlide::susie else susieR::susie
    warns <- character()
    tm <- system.time(fit <- withCallingHandlers(
      quiet(fun,Z,y,L=10,standardize=FALSE,estimate_prior_method="EM",
             min_abs_corr=0,max_iter=1000,tol=1e-6),
      warning=function(w) {warns <<- c(warns,conditionMessage(w)); invokeRestart("muffleWarning")}))
    newx <- if(method=="Susie_mix_EM") stack(example$Xtest) else example$Xtest
    predicted <- predict(fit,newx=newx)
    pip <- if(method=="Susie_mix_EM") group_pip(fit,p) else fit$pip
    i <- i+1L
    rows[[i]] <- data.frame(mode=mode,method=method,seconds=unname(tm["elapsed"]),
      mean_function_mse=mean((predicted-example$mean_test)^2),
      mean_causal_pip=mean(pip[example$causal]),sigma2=fit$sigma2,
      active_components=sum(fit$V>1e-9),converged=fit$converged,niter=fit$niter,
      largest_elbo_drop=min(c(0,diff(fit$elbo))),warnings=paste(warns,collapse="; "))
    fits[[paste(mode,method,sep="_")]] <- fit
  }
  cat("Workhorse comparison:",mode,"\n")
}
workhorse <- do.call(rbind,rows)
write.csv(workhorse,file.path(outdir,"workhorse_comparison.csv"),row.names=FALSE)
saveRDS(fits,file.path(outdir,"workhorse_fits.rds"))

# Small null diagnostic with variance optimization. This is an explicit check
# for overfitting from optimized deltas, not a claim of genome-wide calibration.
null_rows <- list(); i <- 0L
for(rep in 1:30) {
  set.seed(8000+rep)
  X <- matrix(rbinom(300*60,2,.3),300,60); storage.mode(X) <- "double"
  y <- rnorm(300)
  for(method in c("Additive","Slider","Stack_equal")) {
    Z <- if(method=="Stack_equal") stack(X) else X
    fun <- if(method=="Slider") susieSlide::susie else susieR::susie
    warns <- character()
    fit <- withCallingHandlers(quiet(fun,Z,y,L=3,max_iter=1000,tol=1e-6,min_abs_corr=.5),
      warning=function(w) {warns <<- c(warns,conditionMessage(w)); invokeRestart("muffleWarning")})
    pip <- if(method=="Stack_equal") group_pip(fit,ncol(X)) else fit$pip
    i <- i+1L
    null_rows[[i]] <- data.frame(repetition=rep,method=method,max_pip=max(pip),
      sum_pip=sum(pip),snps_above_half=sum(pip>.5),cs_count=length(fit$sets$cs),
      active_components=sum(fit$V>1e-9),converged=fit$converged,
      warnings=paste(warns,collapse="; "))
  }
}
null <- do.call(rbind,null_rows)
write.csv(null,file.path(outdir,"null_diagnostic.csv"),row.names=FALSE)
write.csv(aggregate(null[,c("max_pip","sum_pip","snps_above_half","cs_count","active_components")],
                      by=null["method"],FUN=mean),
          file.path(outdir,"null_summary.csv"),row.names=FALSE)

# Actual one-million-SNP inference pass, with precomputed sufficient statistics.
# No full genotype matrix / IBSS fit / genotype I/O is included in this timing.
set.seed(320)
m <- 1000000L; n <- 1000; sigma <- .6
maf <- runif(m,.03,.5)
n1 <- round(n*2*maf*(1-maf)); n2 <- pmax(1,round(n*maf^2)); n0 <- n-n1-n2
beta <- (runif(m)<.05)*sample(c(-.8,-.3,.3,.8),m,TRUE)
delta <- runif(m,-1,1); sx <- n1+2*n2
sy0 <- sigma*sqrt(n0)*rnorm(m)
sy1 <- n1*beta*(1+delta)+sigma*sqrt(n1)*rnorm(m)
sy2 <- 2*n2*beta+sigma*sqrt(n2)*rnorm(m); sy <- sy0+sy1+sy2
s2 <- (n1+4*n2-sx*sx/n)/(n-1); s <- sqrt(s2)
xx <- rep(n-1,m); xh <- (n1-sx*n1/n)/s2; hh <- (n1-n1*n1/n)/s2
xy <- (sy1+2*sy2-sx*sy/n)/s; hy <- (sy1-n1*sy/n)/s
fixed <- rep(NA_real_,m); fixed[pmin(n0,n1,n2)<5] <- 0
times <- numeric(3); check <- numeric(3)
for(rep in 1:3) {
  timing <- system.time(result <- susieSlide:::slide_ser_native(xx,xh,hh,xy,hy,.5,.36,fixed))
  times[rep] <- timing["elapsed"]; check[rep] <- sum(result[,2])
}
stopifnot(length(unique(check))==1,all(is.finite(result)),all(abs(result[,1])<=1))
write.csv(data.frame(repetition=1:3,snps=m,seconds=times,forced_snps=sum(!is.na(fixed)),checksum=check),
          file.path(outdir,"compiled_million_snp.csv"),row.names=FALSE)
print(workhorse,row.names=FALSE)
cat("Compiled million-SNP elapsed seconds:",times,"\n")
