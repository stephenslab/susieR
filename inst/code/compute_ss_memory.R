# export MEM_CHECK_INTERVAL=0.01
# python3 monitor_memory.py Rscript compute_ss_memory.R
#
# NOTES:
#
# - *Add notes here.*
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
set.seed(1)
cat("Size of X:\n")
print(object.size(X),unit = "GB")
cat("Running susie.\n")
out <- compute_ss(X,y,standardize = TRUE)
