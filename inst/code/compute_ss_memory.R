# export MEM_CHECK_INTERVAL=0.01
# python3 monitor_memory.py Rscript compute_ss_memory.R
#
# NOTES:
#
# - Without any improvements:
#   Size of X: 0.3 GB
#   max rss_memory: 1.65 GB
#
# - The original centering and scaling steps require about 1 GB.
#
# - With the improvements:
#   Size of X: 0.3 GB
#   max rss_memory: 0.66 GB
#
# set.seed(1)
# p <- 2000
# n <- 20000
# X <- matrix(rnorm(n*p),n,p)
# y <- rnorm(n)
# save(list = c("X","y"),file = "compute_ss_data.RData")
# library(susieR)
devtools::load_all()
set.seed(1)
load("compute_ss_data.RData")
cat("Size of X:\n")
print(object.size(X),unit = "GB")
cat("Running compute_ss.\n")
out <- compute_ss(X,y,standardize = TRUE)
