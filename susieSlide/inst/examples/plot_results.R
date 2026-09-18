# Run after both comparison scripts, from the package source directory.
library(susieSlide)
outdir <- Sys.getenv("SLIDE_COMPARISON_OUT","validation")
summary <- read.csv(file.path(outdir,"comparison_summary.csv"))
modes <- c("Additive","Recessive","Partially_recessive","Dominant","Partially_dominant","Mixed")
colors <- c(Additive="#C15B42",Slider="#087F8C",Stack_equal="#8B7D3C",Stack_EM="#65519B")
png(file.path(outdir,"comparison_plot.png"),width=1700,height=900,res=150)
par(mfrow=c(1,2),mar=c(7,4.5,3,1),las=1)
for(metric in c("mean_function_mse","mean_causal_pip")) {
  heights <- sapply(modes,function(mode) {
    z <- summary[summary$mode==mode,]; z[[metric]][match(names(colors),z$method)]
  })
  positions <- barplot(heights,beside=TRUE,col=colors,border=NA,axisnames=FALSE,
    ylab=if(metric=="mean_function_mse") "Test MSE against true mean" else "Mean causal SNP PIP",
    main=if(metric=="mean_function_mse") "Prediction of the genetic signal" else "Support for the causal variants",
    ylim=if(metric=="mean_causal_pip") c(0,1.13) else NULL)
  text(colMeans(positions),par("usr")[3]-.04*diff(par("usr")[3:4]),
       labels=gsub("_","\n",modes),xpd=NA,adj=c(.5,1),cex=.75)
  if(metric=="mean_function_mse") legend("topleft",
    legend=c("Additive","Slider","Stacked, equal weights","Stacked, coding-weight EM"),
    fill=colors,bty="n",cex=.8)
}
dev.off()

examples <- readRDS(file.path(outdir,"comparison_examples.rds"))
png(file.path(outdir,"effect_shapes.png"),width=1800,height=1150,res=160)
par(mfrow=c(2,3),mar=c(4,4,3,1),oma=c(1,1,3,0),las=1)
for(mode in modes) {
  e <- examples[[paste(mode,2,"Slider",sep="_")]]
  j <- e$causal[1]; p <- ncol(e$X); x <- 0:2
  truth <- .9*(x+e$delta_true[1]*(x==1))
  fitted <- sapply(c("Additive","Slider","Stack_equal"),function(method) {
    obj <- examples[[paste(mode,2,method,sep="_")]]$fit
    if(method=="Slider") {
      cf <- coef(obj)
      x*cf[j+1,1]+(x==1)*cf[j+1,2]
    } else if(method=="Additive") x*coef(obj)[j+1] else {
      cf <- coef(obj$fit)[-1]
      (x*cf[j]+2*(x>=1)*cf[p+j]+2*(x==2)*cf[2*p+j])/obj$scales[j]
    }
  })
  plot(x,truth,type="l",lwd=5,col="#C9CED3",xaxt="n",
    xlab="Original genotype",ylab="Effect relative to genotype 0",
    main=gsub("_"," ",mode),ylim=range(c(truth,fitted)))
  axis(1,at=x)
  for(k in seq_len(ncol(fitted))) lines(x,fitted[,k],type="b",pch=c(17,16,15)[k],
    col=colors[k],lwd=2,lty=k,cex=.9)
  if(mode=="Additive") legend("topleft",legend=c("True","Additive","Slider","Stacked (L=2)"),
    col=c("#C9CED3",colors[1:3]),lwd=c(5,2,2,2),lty=c(1,1,2,3),bty="n",cex=.8)
}
mtext("Recovered effect shape: first causal SNP, replicate 1",side=3,outer=TRUE,line=1,font=2,cex=1.2)
mtext("Matched priors and L=2. Stacking can use additional components when L is larger.",
      side=1,outer=TRUE,cex=.8)
dev.off()

main <- read.csv(file.path(outdir,"comparison_runs.csv"))
est <- read.csv(file.path(outdir,"causal_estimates.csv"))
null <- read.csv(file.path(outdir,"null_diagnostic.csv"))
cat("All main/sensitivity fits converged:",all(main$converged & main$em_converged),"\n")
cat("All null fits converged:",all(null$converged),"\n")
print(aggregate(cbind(any_half=as.numeric(snps_above_half>0),active_components)~method,null,mean))
cat("Mean absolute delta error, L=2:\n")
print(aggregate(abs(delta_estimated-delta_true)~mode,est[est$method=="Slider" & est$L==2,],mean))
cat("L=4 sensitivity, single replicate:\n")
print(main[main$L==4,c("mode","method","mean_function_mse","converged")],row.names=FALSE)
