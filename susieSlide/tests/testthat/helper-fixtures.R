quiet <- function(fun,...) suppressMessages(fun(...))
genotypes <- function(n=350,p=18) {
  X <- matrix(rbinom(n*p,2,.3),n,p)
  storage.mode(X) <- "double"
  colnames(X) <- paste0("rs",seq_len(p))
  X
}
ser_native <- function(x,y,V=.8,sigma2=.49,delta=NA_real_,standardize=FALSE,intercept=TRUE) {
  h <- as.numeric(x==1)
  s <- if(standardize) sd(x) else 1
  if(s==0) s <- 1
  z <- (x-if(intercept) mean(x) else 0)/s
  h <- (h-if(intercept) mean(h) else 0)/s
  r <- y-if(intercept) mean(y) else 0
  .Call("slide_ser",as.double(sum(z*z)),as.double(sum(z*h)),as.double(sum(h*h)),
        as.double(sum(z*r)),as.double(sum(h*r)),as.double(V),as.double(sigma2),
        as.double(delta),PACKAGE="susieSlide")
}
