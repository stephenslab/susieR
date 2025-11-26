# Script to summarize the results of running small_sim.R.
library(ggplot2)
library(cowplot)
load("../datafiles/small_sim_out_v0.14.48.RData")
runtimes <- data.frame(susie    = runtimes$susie,
                       susie_ss = runtimes$susie_small)
methods <- c("susie","susie_ss")

# Summarize coverage, power and running times.
N <- length(causal_snps)
power           <- c(0,0)
coverage        <- c(0,0)
names(power)    <- methods
names(coverage) <- methods
V1_true <- NULL
V2_true <- NULL
V1_false <- NULL
V2_false <- NULL
for (i in 1:N) {
  get_tp <- function (cs) {
    if (length(cs) == 0)
      return(NULL)
    else
      return(names(which(sapply(cs, 
        function (x) length(intersect(causal_snps[[i]],x))>0))))
  }
  V1 <- res_susie[[i]]$V
  V2 <- res_susie_small[[i]]$V
  all_cs <- paste0("L",1:10)
  names(V1) <- all_cs
  names(V2) <- all_cs
  x1 <- get_tp(res_susie[[i]]$sets$cs)
  x2 <- get_tp(res_susie_small[[i]]$sets$cs)
  V1_true <- c(V1_true,V1[x1])
  V2_true <- c(V2_true,V2[x2])
  V1_false <- c(V1_false,V1[setdiff(all_cs,x1)])
  V2_false <- c(V2_false,V2[setdiff(all_cs,x2)])
  power["susie"]       <- power["susie"]       + length(x1)
  power["susie_ss"]    <- power["susie_ss"]    + length(x2)
  coverage["susie"]    <- coverage["susie"]    + length(x1)
  coverage["susie_ss"] <- coverage["susie_ss"] + length(x2)
}

num_causal   <- sum(sapply(causal_snps,length))
num_susie    <- sum(sapply(res_susie,function (x) length(x$sets$cs)))
num_susie_ss <- sum(sapply(res_susie_small,function (x) length(x$sets$cs)))
power <- power / num_causal
coverage["susie"]    <- coverage["susie"] / num_susie
coverage["susie_ss"] <- coverage["susie_ss"] / num_susie_ss
cat("power:\n")
print(power)
cat("coverage:\n")
print(coverage)
cat("running times:\n")
print(summary(runtimes))

# Summarize the CS sizes.
get_cs_sizes <- function (res)
  unlist(lapply(res,function (x) sapply(x$sets$cs,length)))
sizes_susie    <- get_cs_sizes(res_susie)
sizes_susie_ss <- get_cs_sizes(res_susie_small)
cat("median CS size:\n")
cat("susie =",median(sizes_susie),"\n")
cat("susie_ss =",median(sizes_susie_ss),"\n")
pdat <- rbind(data.frame(method = "susie",   size = sizes_susie),
              data.frame(method = "susie_ss",size = sizes_susie_ss))
pdat <- subset(pdat,size <= 50)
p1 <- ggplot(pdat,aes(x = size,fill = method)) +
  geom_histogram(color = "white",position = "dodge",bins = 16) +
  scale_fill_manual(values = c("darkblue","dodgerblue","darkorange")) +
  labs(x = "CS size",y = "number of CSs",fill = "") +
  theme_cowplot(font_size = 10) +
  theme(legend.position = "bottom",
        legend.direction = "vertical")
ggsave("small_sim_sizes.pdf",p1,height = 3,width = 3)

# Summarize the prior variances.
pdat <- rbind(data.frame(method = "susie",   causal = TRUE, V = V1_true),
              data.frame(method = "susie",   causal = FALSE,V = V1_false),
              data.frame(method = "susie_ss",causal = TRUE, V = V2_true),
              data.frame(method = "susie_ss",causal = FALSE,V = V2_false))
pdat <- transform(pdat,sigma = sqrt(V))
pdat <- subset(pdat,sigma < 2)
p2 <- ggplot(pdat,aes(x = sigma,color = causal,fill = causal)) +
  facet_grid(rows = vars(method),scales = "free_y") +
  geom_histogram(bins = 24,position = "dodge",linewidth = 0.05) +
  scale_color_manual(values = c("darkblue","orangered")) +
  scale_fill_manual(values = c("darkblue","orangered")) +
  labs(x = "prior s.d.") +
  theme_cowplot(font_size = 10)
ggsave("small_sim_V.pdf",p2,height = 3.5,width = 3)
