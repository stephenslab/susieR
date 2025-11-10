# Script to summarize the results of running the small_sim.R with
# different versions of susieR.
#
# TO DO: Compare the running times of the different methods.
# 
library(ggplot2)
library(cowplot)
load("../datafiles/small_sim_out_v0.14.25.RData")
res_susie_small_old <- res_susie_small
causal_snps_old     <- causal_snps
runtimes_old        <- runtimes
load("../datafiles/small_sim_out_v0.14.34.RData")
res_susie_small_new <- res_susie_small
runtimes <- data.frame(susie           = runtimes$susie,
                       susie_small_old = runtimes_old$susie_small,
                       susie_small_new = runtimes$susie_small)
methods <- c("susie","susie_small_old","susie_small_new")

# Summarize coverage, power, purity, and CS sizes.
N <- length(causal_snps)
power           <- c(0,0,0)
coverage        <- c(0,0,0)
names(power)    <- methods
names(coverage) <- methods
for (i in 1:N) {
  x1 <- length(intersect(causal_snps[[i]],unlist(res_susie[[i]]$cs)))
  x2 <- length(intersect(causal_snps_old[[i]],
                         unlist(res_susie_small_old[[i]]$cs)))
  x3 <- length(intersect(causal_snps[[i]],
                         unlist(res_susie_small_new[[i]]$cs)))
  power["susie"]              <- power["susie"]              + x1
  power["susie_small_old"]    <- power["susie_small_old"]    + x2
  power["susie_small_new"]    <- power["susie_small_new"]    + x3
  coverage["susie"]           <- coverage["susie"]           + x1
  coverage["susie_small_old"] <- coverage["susie_small_old"] + x2
  coverage["susie_small_new"] <- coverage["susie_small_new"] + x3
}

power["susie"] <- 
  power["susie"] / sum(sapply(causal_snps,length))
power["susie_small_old"] <- 
  power["susie_small_old"] / sum(sapply(causal_snps_old,length))
power["susie_small_new"] <- 
  power["susie_small_new"] / sum(sapply(causal_snps,length))
coverage["susie"] <- 
  coverage["susie"] / sum(sapply(res_susie,function (x) length(x$cs)))
coverage["susie_small_old"] <- 
  coverage["susie_small_old"] / sum(sapply(res_susie_small_old,
                                           function (x) length(x$cs)))
coverage["susie_small_new"] <- 
  coverage["susie_small_new"] / sum(sapply(res_susie_small_new,
                                           function (x) length(x$cs)))
cat("power:\n")
print(power)
cat("coverage:\n")
print(coverage)
cat("running times:\n")
print(summary(runtimes))

# Plot the distribution of CS sizes.
sizes_susie <- unlist(lapply(res_susie,function (x) sapply(x$cs,length)))
sizes_susie_small_old <- unlist(lapply(res_susie_small_old,
                                       function (x) sapply(x$cs,length)))
sizes_susie_small_new <- unlist(lapply(res_susie_small_new,
                                       function (x) sapply(x$cs,length)))
pdat <- rbind(data.frame(method = "susie",
                         size   = sizes_susie),
              data.frame(method = "susie_small_old",
                         size   = sizes_susie_small_old),
              data.frame(method = "susie_small_new",
                         size   = sizes_susie_small_new))
p <- ggplot(pdat,aes(x = size,fill = method)) +
  geom_histogram(color = "white",position = "dodge") +
  scale_fill_manual(values = c("darkblue","dodgerblue","darkorange")) +
  labs(x = "CS size",y = "number of CSs",fill = "") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "bottom",
        legend.direction = "vertical")
print(p)
