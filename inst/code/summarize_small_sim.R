# Script to summarize the results of running the small_sim.R with
# different versions of susieR.
library(ggplot2)
library(cowplot)
load("../datafiles/small_sim_out_v0.14.25.RData")
res_susie_small_old <- res_susie_small
runtimes_old        <- runtimes
load("../datafiles/small_sim_out_v0.14.34.RData")
res_susie_small_new <- res_susie_small
runtimes <- data.frame(susie           = runtimes$susie,
                       susie_small_old = runtimes_old$susie_small,
                       susie_small_new = runtimes$susie_small)
methods <- c("susie","susie_small_old","susie_small_new")

# Summarize coverage, power and running times.
N <- length(causal_snps)
power           <- c(0,0,0)
coverage        <- c(0,0,0)
names(power)    <- methods
names(coverage) <- methods
V1_true <- NULL
V2_true <- NULL
V3_true <- NULL
V1_false <- NULL
V2_false <- NULL
V3_false <- NULL
for (i in 1:N) {
  get_tp <- function (cs) {
    if (length(cs) == 0)
      return(NULL)
    else
      return(names(which(sapply(cs, 
        function (x) length(intersect(causal_snps[[i]],x))>0))))
  }
  V1 <- res_susie[[i]]$V
  V2 <- res_susie_small_old[[i]]$V
  V3 <- res_susie_small_new[[i]]$V
  all_cs <- paste0("L",1:10)
  names(V1) <- all_cs
  names(V2) <- all_cs
  names(V3) <- all_cs
  x1 <- get_tp(res_susie[[i]]$sets$cs)
  x2 <- get_tp(res_susie_small_old[[i]]$sets$cs)
  x3 <- get_tp(res_susie_small_new[[i]]$sets$cs)
  V1_true <- c(V1_true,V1[x1])
  V2_true <- c(V2_true,V2[x2])
  V3_true <- c(V3_true,V3[x3])
  V1_false <- c(V1_false,V1[setdiff(all_cs,x1)])
  V2_false <- c(V2_false,V2[setdiff(all_cs,x2)])
  V3_false <- c(V3_false,V3[setdiff(all_cs,x3)])
  power["susie"]              <- power["susie"]              + length(x1)
  power["susie_small_old"]    <- power["susie_small_old"]    + length(x2)
  power["susie_small_new"]    <- power["susie_small_new"]    + length(x3)
  coverage["susie"]           <- coverage["susie"]           + length(x1)
  coverage["susie_small_old"] <- coverage["susie_small_old"] + length(x2)
  coverage["susie_small_new"] <- coverage["susie_small_new"] + length(x3)
}

num_causal          <- sum(sapply(causal_snps,length))
num_susie           <- sum(sapply(res_susie,function (x) length(x$sets$cs)))
num_susie_small_old <- sum(sapply(res_susie_small_old,
                                  function (x) length(x$sets$cs)))
num_susie_small_new <- sum(sapply(res_susie_small_new,
                                  function (x) length(x$sets$cs)))
power <- power / num_causal
coverage["susie"] <- coverage["susie"] / num_susie
coverage["susie_small_old"] <- 
  coverage["susie_small_old"] / num_susie_small_old
coverage["susie_small_new"] <- 
    coverage["susie_small_new"] / num_susie_small_new
cat("power:\n")
print(power)
cat("coverage:\n")
print(coverage)
cat("running times:\n")
print(summary(runtimes))

# Summarize the CS sizes.
get_cs_sizes <- function (res)
  unlist(lapply(res,function (x) sapply(x$sets$cs,length)))
sizes_susie           <- get_cs_sizes(res_susie)
sizes_susie_small_old <- get_cs_sizes(res_susie_small_old)
sizes_susie_small_new <- get_cs_sizes(res_susie_small_new)
cat("median CS size:\n")
cat("susie =",median(sizes_susie),"\n")
cat("susie_small_old =",median(sizes_susie_small_old),"\n")
cat("susie_small_new =",median(sizes_susie_small_new),"\n")
pdat <- rbind(data.frame(method = "susie",
                         size   = sizes_susie),
              data.frame(method = "susie_small_old",
                         size   = sizes_susie_small_old),
              data.frame(method = "susie_small_new",
                         size   = sizes_susie_small_new))
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
pdat <- rbind(data.frame(method = "susie",causal = TRUE,V = V1_true),
              data.frame(method = "susie",causal = FALSE,V = V1_false),
              data.frame(method = "susie_small_old",causal=TRUE,V=V2_true),
              data.frame(method = "susie_small_old",causal=FALSE,V=V2_false),
              data.frame(method = "susie_small_new",causal=TRUE,V=V3_true),
              data.frame(method = "susie_small_new",causal=FALSE,V=V3_false))
pdat <- transform(pdat,sigma = sqrt(V))
p2 <- ggplot(pdat,aes(x = sigma,fill = causal)) +
  facet_grid(rows = vars(method),scales = "free_y") +
  geom_histogram(color = "white",bins = 36) +
  scale_fill_manual(values = c("darkblue","orangered")) +
  labs(x = "prior s.d.") +
  theme_cowplot(font_size = 10)
ggsave("small_sim_V.pdf",p2,height = 3.5,width = 3)
