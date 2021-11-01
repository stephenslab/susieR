# export MEM_CHECK_INTERVAL=0.01
# python3 monitor_memory.py Rscript susie_rss_memory.R
#
# NOTES:
#
# - Without any improvements:
#   Size of X: 0.5 GB
#   max rss_memory: 4.15 GB
#
# - With some improvements, the initial steps right before the main
#   loop use 0.86 GB.
#
# - susie_rss right before the CS and PIP calculations uses 1.6 GB.
#
# library(susieR)
devtools::load_all()
# set.seed(1)
# p  <- 8000
# n  <- 1000
# X  <- matrix(rnorm(n*p),n,p)
# y  <- rnorm(n)
# ss <- susieR:::univariate_regression(X,y)
# z  <- ss$betahat/ss$sebetahat
# R  <- cor(X)
# save(list = c("X","y","z","R"),file = "susie_rss_data.RData")
load("susie_rss_data.RData")
cat("Size of R:\n")
print(object.size(R),unit = "GB")
cat("Running susie_rss.\n")
set.seed(1)
out <- susie_rss(z,R,estimate_prior_variance = FALSE,min_abs_corr = 0,
                 check_input = FALSE,refine = FALSE,verbose = TRUE)
print(sapply(out$sets$cs,length))
