# export MEM_CHECK_INTERVAL=0.01
# python3 monitor_memory.py Rscript susie_rss_memory.R
#
# NOTES:
#
# - Add notes here.
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
                 verbose = TRUE)
print(sapply(out$sets$cs,length))
