# Reproduce from the package source root:
#   Rscript inst/examples/compare_models.R
# Output path and repetitions can be changed through environment variables.
# All methods use the same Gaussian coefficient prior in additive-SD units.
library(susieSlide)
stopifnot(requireNamespace("susieR",quietly=TRUE))
outdir <- Sys.getenv("SLIDE_COMPARISON_OUT","validation")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
repetitions <- as.integer(Sys.getenv("SLIDE_COMPARISON_REPS","10"))
stopifnot(repetitions>=1)
quiet <- function(fun,...) suppressMessages(fun(...))

# Two independent haplotypes per individual; adjacent SNPs share alleles
# with probability rho. Blocks are independent; marginal MAF is stationary.
simulate_ld <- function(n,p=80,rho=.8,block=20) {
  H <- matrix(0,2*n,p)
  for(j in seq_len(p)) {
    maf <- c(.25,.35,.3,.4)[ceiling(j/block)]
    fresh <- rbinom(2*n,1,maf)
    H[,j] <- if((j-1)%%block==0) fresh else ifelse(runif(2*n)<rho,H[,j-1],fresh)
  }
  X <- H[seq_len(n),,drop=FALSE]+H[n+seq_len(n),,drop=FALSE]
  storage.mode(X) <- "double"
  colnames(X) <- paste0("rs",seq_len(p))
  X
}
stack_design <- function(X,scales) {
  Z <- cbind(X,2*(X>=1),2*(X==2))
  Z <- sweep(Z,2,rep(scales,3),"/")
  colnames(Z) <- paste(rep(c("A","D","R"),each=ncol(X)),
                       rep(colnames(X),3),sep=":")
  Z
}
fit_stack <- function(X,y,L=2,V=.5,sigma2=.36,learn_weights=FALSE,
                       learn_V=FALSE,max_em=40,em_tol=1e-5) {
  p <- ncol(X); scales <- apply(X,2,sd); scales[scales==0] <- 1
  Z <- stack_design(X,scales)
  weights <- rep(1/3,3); fit <- NULL; converged_em <- !learn_weights
  history <- list()
  for(iter in seq_len(if(learn_weights) max_em else 1)) {
    fit <- quiet(susieR::susie,Z,y,L=L,standardize=FALSE,
      residual_variance=sigma2,estimate_residual_variance=FALSE,
      scaled_prior_variance=V/var(y),estimate_prior_variance=learn_V,
      prior_weights=rep(weights,each=p)/p,model_init=fit,
      max_iter=400,tol=1e-7,min_abs_corr=0)
    history[[iter]] <- c(weights,elbo=tail(fit$elbo,1))
    if(!learn_weights) break
    mass <- colSums(fit$alpha)
    # Variational EM M step: expected coding counts summed over effects/SNPs.
    # Components whose estimated V is zero provide no coding information.
    active <- which(fit$V>1e-9)
    if(!length(active)) {converged_em <- TRUE; break}
    mass <- colSums(fit$alpha[active,,drop=FALSE])
    new_weights <- vapply(0:2,function(k) sum(mass[k*p+seq_len(p)]),numeric(1))
    new_weights <- new_weights/sum(new_weights)
    if(max(abs(new_weights-weights))<em_tol) {converged_em <- TRUE; break}
    # On the last iteration leave weights consistent with the returned fit.
    if(iter<max_em) weights <- new_weights
  }
  group_alpha <- fit$alpha[,seq_len(p),drop=FALSE]+
    fit$alpha[,p+seq_len(p),drop=FALSE]+fit$alpha[,2*p+seq_len(p),drop=FALSE]
  group_pip <- -expm1(colSums(log1p(-pmin(group_alpha[fit$V>1e-9,,drop=FALSE],1))))
  list(fit=fit,alpha=group_alpha,pip=group_pip,scales=scales,weights=weights,
       em_converged=converged_em,em_iter=iter,em_history=do.call(rbind,history))
}
component_sets <- function(alpha,V,coverage=.95) {
  lapply(which(V>1e-9),function(l) {
    j <- order(alpha[l,],decreasing=TRUE)
    k <- which(cumsum(alpha[l,j])>=coverage)[1]
    if(is.na(k)) k <- length(j)
    j[seq_len(k)]
  })
}
score_fit <- function(method,fit,Xtest,ytest,mean_test,causal,delta_true,
                       elapsed,mode,rep,L,learn_V) {
  stacked <- grepl("Stack",method)
  f <- if(stacked) fit$fit else fit
  alpha <- if(stacked) fit$alpha else f$alpha
  pip <- if(stacked) fit$pip else f$pip
  prediction <- if(stacked) predict(f,newx=stack_design(Xtest,fit$scales)) else predict(f,newx=Xtest)
  sets <- component_sets(alpha,f$V)
  lead <- apply(alpha,1,which.max)
  active <- which(f$V>1e-9)
  hit <- causal %in% lead[active]
  covered <- vapply(causal,function(j) any(vapply(sets,function(s) j %in% s,logical(1))),logical(1))
  row <- data.frame(mode=mode,repetition=rep,L=L,learn_prior_variance=learn_V,
    method=method,seconds=elapsed,test_mse=mean((ytest-prediction)^2),
    mean_function_mse=mean((mean_test-prediction)^2),
    mean_causal_pip=mean(pip[causal]),lead_recall=mean(hit),
    causal_in_any_95cs=mean(covered),
    mean_cs_size=if(length(sets)) mean(lengths(sets)) else NA_real_,
    active_components=length(active),converged=f$converged,
    em_converged=if(stacked) fit$em_converged else TRUE,
    em_iterations=if(stacked) fit$em_iter else 0,
    null_snps_pip_above_half=sum(pip[-causal]>.5))
  estimates <- lapply(seq_along(causal),function(k) {
    j <- causal[k]; l <- which.max(alpha[,j])
    data.frame(mode=mode,repetition=rep,L=L,learn_prior_variance=learn_V,method=method,
      snp=j,component=l,delta_true=delta_true[k],
      delta_estimated=if(method=="Slider") f$delta[l,j] else NA_real_,pip=pip[j],
      rank=rank(-pip,ties.method="min")[j])
  })
  list(metrics=row,estimates=do.call(rbind,estimates))
}

modes <- list(Additive=c(0,0),Recessive=c(-1,-1),
               Partially_recessive=c(-.5,-.5),Dominant=c(1,1),
               Partially_dominant=c(.5,.5),Mixed=c(-.6,.6))
metrics <- estimates <- examples <- list(); counter <- 0L
for(rep in seq_len(repetitions)) for(mode in names(modes)) {
  set.seed(10000+100*rep+match(mode,names(modes)))
  X <- simulate_ld(1000); Xtest <- simulate_ld(2000)
  causal <- c(12L,52L); beta <- c(.9,-.8); delta_true <- modes[[mode]]
  signal <- function(X) as.numeric((X[,causal,drop=FALSE]+
    sweep((X[,causal,drop=FALSE]==1)*1,2,delta_true,"*")) %*% beta)
  y <- .3+signal(X)+rnorm(nrow(X),sd=.6)
  mean_test <- .3+signal(Xtest); ytest <- mean_test+rnorm(nrow(Xtest),sd=.6)
  # Main comparison uses exactly two components and fixed matched priors.
  # For replicate 1, also allow extra components and learn their variances.
  for(setting in seq_len(if(rep==1) 2 else 1)) {
    L <- if(setting==1) 2 else 4; learn_V <- setting==2
    for(method in c("Additive","Slider","Stack_equal","Stack_EM")) {
      timing <- system.time({
        if(method=="Additive")
          fit <- quiet(susieR::susie,X,y,L=L,standardize=TRUE,
            residual_variance=.36,estimate_residual_variance=FALSE,
            scaled_prior_variance=.5/var(y),estimate_prior_variance=learn_V,
            max_iter=400,tol=1e-7,min_abs_corr=0)
        else if(method=="Slider")
          fit <- quiet(susieSlide::susie,X,y,L=L,standardize=TRUE,
            residual_variance=.36,estimate_residual_variance=FALSE,
            scaled_prior_variance=.5/var(y),estimate_prior_variance=learn_V,
            max_iter=400,tol=1e-7,min_abs_corr=0)
        else fit <- fit_stack(X,y,L=L,learn_V=learn_V,learn_weights=method=="Stack_EM")
      })
      counter <- counter+1L
      scored <- score_fit(method,fit,Xtest,ytest,mean_test,causal,delta_true,
                           unname(timing["elapsed"]),mode,rep,L,learn_V)
      metrics[[counter]] <- scored$metrics; estimates[[counter]] <- scored$estimates
      if(rep==1) examples[[paste(mode,L,method,sep="_")]] <-
        list(X=X,y=y,Xtest=Xtest,mean_test=mean_test,causal=causal,delta_true=delta_true,fit=fit)
    }
  }
  cat("Completed replicate",rep,"/",repetitions,"mode",mode,"\n")
}
metrics <- do.call(rbind,metrics); estimates <- do.call(rbind,estimates)
write.csv(metrics,file.path(outdir,"comparison_runs.csv"),row.names=FALSE)
write.csv(estimates,file.path(outdir,"causal_estimates.csv"),row.names=FALSE)
saveRDS(examples,file.path(outdir,"comparison_examples.rds"))
main <- metrics[metrics$L==2,]
summary <- aggregate(main[,c("seconds","test_mse","mean_function_mse","mean_causal_pip",
                             "lead_recall","causal_in_any_95cs","mean_cs_size",
                             "null_snps_pip_above_half")],
                      by=main[,c("mode","method")],FUN=mean)
write.csv(summary,file.path(outdir,"comparison_summary.csv"),row.names=FALSE)
writeLines(c(capture.output(sessionInfo()),paste("Replicates:",repetitions),
  "Known residual variance 0.36. Main: L=2, fixed coefficient prior variance 0.5.",
  "All coefficient priors matched on the original additive genotype SD scale.",
  "Stack_EM is a documented variational-EM comparator, not an inspected external susie-mix package.",
  "Coverage metric is membership in any component's 95% SNP-level candidate set, with no purity filter.",
  "These selected non-null examples do not establish PIP/FDR or credible-set calibration."),
  file.path(outdir,"comparison_session.txt"))

colors <- c(Additive="#C15B42",Slider="#087F8C",Stack_equal="#8B7D3C",Stack_EM="#65519B")
png(file.path(outdir,"comparison_plot.png"),width=1500,height=800,res=150)
par(mfrow=c(1,2),mar=c(7,4.5,3,1),las=1)
for(metric in c("mean_function_mse","mean_causal_pip")) {
  heights <- sapply(names(modes),function(mode) {
    z <- summary[summary$mode==mode,]; z[[metric]][match(names(colors),z$method)]
  })
  positions <- barplot(heights,beside=TRUE,col=colors,border=NA,axisnames=FALSE,
    ylab=if(metric=="mean_function_mse") "Test error against true mean" else "Mean causal SNP PIP",
    main=if(metric=="mean_function_mse") "Prediction of the genetic signal" else "Support for the causal variants",
    ylim=if(metric=="mean_causal_pip") c(0,1.13) else NULL)
  text(colMeans(positions),par("usr")[3]-.04*diff(par("usr")[3:4]),
       labels=gsub("_","\n",names(modes)),xpd=NA,adj=c(.5,1),cex=.68)
  if(metric=="mean_function_mse") legend("topleft",legend=gsub("_"," ",names(colors)),fill=colors,bty="n",cex=.8)
}
dev.off()
print(summary,row.names=FALSE)
cat("Nonconverged inner fits:",sum(!metrics$converged),
    "; nonconverged coding EM fits:",sum(!metrics$em_converged),"\n")
