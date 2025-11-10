# Script to summarize the results of running the small_sim.R with
# different versions of susieR.
load("small_sim_out.RData")
cat("n =",n,"\n")
methods <- c("susie","susie_small")

# Summarize coverage, power, purity, and CS sizes.
N <- length(causal_snps)
power           <- c(0,0)
coverage        <- c(0,0)
names(power)    <- methods
names(coverage) <- methods
for (i in 1:N) {
  x1 <- length(intersect(causal_snps[[i]],unlist(res_susie[[i]]$cs)))
  x2 <- length(intersect(causal_snps[[i]],unlist(res_susie_small[[i]]$cs)))
  power["susie"]          <- power["susie"]          + x1
  power["susie_small"]    <- power["susie_small"]    + x2
  coverage["susie"]       <- coverage["susie"]       + x1
  coverage["susie_small"] <- coverage["susie_small"] + x2
}

power <- power / sum(sapply(causal_snps,length))
coverage["susie"] <- 
  coverage["susie"] / sum(sapply(res_susie,function (x) length(x$cs)))
coverage["susie_small"] <- 
  coverage["susie_small"] / sum(sapply(res_susie_small,
                                       function (x) length(x$cs)))
