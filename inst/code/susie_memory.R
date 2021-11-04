# export MEM_CHECK_INTERVAL=0.01
# python3 monitor_memory.py Rscript susie_memory.R
#
# NOTES:
#
# - Without any improvements:
#   Size of X: 1 GB
#   max rss_memory: 4.70 GB
#
# - With the improvements:
#   Size of X: 1 GB
#   max rss_memory: 3.00 GB
#
# set.seed(1)
# p <- 16000
# n <- 8000
# X <- matrix(rnorm(n*p),n,p)
# X <- scale(X,center = TRUE,scale = TRUE)
# y <- rnorm(n)
# save(list = c("X","y"),file = "susie_data.RData")
# library(susieR)
devtools::load_all()
load("susie_data.RData")
cat("Size of X:\n")
print(object.size(X),unit = "GB")
cat("Running susie.\n")
set.seed(1)
out <- susie(X,y,estimate_prior_variance = FALSE,min_abs_corr = 0,
             verbose = TRUE)
print(sapply(out$sets$cs,length))
