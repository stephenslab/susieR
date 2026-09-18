# End-to-end in-memory fitting timings; no file I/O. Run at package source root.
# One simulated data set per size, one warm-up, three timed repetitions.
library(susieSlide)
outdir <- Sys.getenv("SLIDE_COMPARISON_OUT","validation")
dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
rows <- list(); i <- 0L
for(p in c(1000L,10000L)) {
  set.seed(909+p)
  X <- matrix(rbinom(500*p,2,.3),500,p)
  storage.mode(X) <- "double"
  y <- .8*(X[,15]-.5*(X[,15]==1))-.9*(X[,p-20]+.5*(X[,p-20]==1))+rnorm(500,sd=.6)
  scale <- apply(X,2,sd)
  # Stacked designs are divided by the original additive genotype SD.
  Z <- sweep(cbind(X,2*(X>=1),2*(X==2)),2,rep(scale,3),"/")
  for(method in c("Additive","Slider","Stack_equal")) {
    call_fit <- function() {
      fun <- if(method=="Slider") susieSlide::susie else susieR::susie
      suppressMessages(fun(if(method=="Stack_equal") Z else X,y,L=2,
        standardize=method!="Stack_equal",scaled_prior_variance=.5/var(y),
        residual_variance=.36,estimate_prior_variance=FALSE,
        estimate_residual_variance=FALSE,max_iter=200,tol=1e-6,coverage=NULL))
    }
    invisible(call_fit())
    for(rep in 1:3) {
      tm <- system.time(fit <- call_fit())
      stopifnot(fit$converged)
      i <- i+1L
      rows[[i]] <- data.frame(n=nrow(X),p=p,method=method,repetition=rep,
        seconds=unname(tm["elapsed"]),niter=fit$niter,converged=fit$converged)
    }
  }
  cat("Completed end-to-end benchmark p =",p,"\n")
}
results <- do.call(rbind,rows)
write.csv(results,file.path(outdir,"fit_benchmark.csv"),row.names=FALSE)
print(aggregate(seconds~n+p+method,results,median),row.names=FALSE)
writeLines(c(capture.output(sessionInfo()),
  "Elapsed fit times include package input preparation and IBSS, but exclude genotype simulation and stacking.",
  "Fixed residual/coefficient variances; L=2; one warm-up and three repetitions; default slider H cache is FALSE."),
  file.path(outdir,"benchmark_session.txt"))
