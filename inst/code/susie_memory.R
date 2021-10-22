# export MEM_CHECK_INTERVAL=0.01
# python3 monitor_memory.py Rscript susie_memory.R
# max rss_memory (without susie_rss): 0.64 GB
# max rss_memory (with susie_rss): 1.48 GB
library(susieR)
set.seed(1)
p <- 4000
n <- 1000
X <- matrix(rnorm(n*p),n,p)
y <- rnorm(n)
ss <- susieR:::univariate_regression(X,y)
z <- ss$betahat/ss$sebetahat
R <- cor(X)
print(object.size(R),unit = "GB")
# out <- susie_rss(z,R,estimate_prior_variance = FALSE,
#                  min_abs_corr = 0,verbose = TRUE)
# print(sapply(out$sets$cs,length))
