# export MEM_CHECK_INTERVAL=0.01
# python3 monitor_memory.py Rscript susie_memory.R
#
# NOTES:
#
# - Without any improvements:
#   Size of X: 1 GB
#
# - The initial checks use about 2 GB. Most of that memory usage is
#   due to set_X_attributes.
#
# library(susieR)
devtools::load_all()
# set.seed(1)
# p <- 16000
# n <- 8000
# X <- matrix(rnorm(n*p),n,p)
# y <- rnorm(n)
# save(list = c("X","y"),file = "susie_data.RData")
load("susie_data.RData")
cat("Size of X:\n")
print(object.size(X),unit = "GB")
cat("Running susie.\n")
set.seed(1)
out <- susie(X,y,estimate_prior_variance = FALSE,min_abs_corr = 0,
             verbose = TRUE)
print(sapply(out$sets$cs,length))
