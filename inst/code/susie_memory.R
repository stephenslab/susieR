# export MEM_CHECK_INTERVAL=0.01
# python3 monitor_memory.py Rscript susie_memory.R
# library(susieR)
devtools::load_all()
set.seed(1)
p <- 12000 # 8000
n <- 4000 # 1000
X <- matrix(rnorm(n*p),n,p)
y <- rnorm(n)
# ss <- susieR:::univariate_regression(X,y)
# z <- ss$betahat/ss$sebetahat
# R <- cor(X)
cat("Size of X:\n")
print(object.size(X),unit = "GB")
# cat("Size of R:\n")
# print(object.size(R),unit = "GB")
cat("Running susie_rss.\n")
out <- susie(X,y,estimate_prior_variance = FALSE,min_abs_corr = 0,
             verbose = TRUE)
# out <- susie_rss(z,R,estimate_prior_variance = FALSE,min_abs_corr = 0,
#                  verbose = TRUE)
print(sapply(out$sets$cs,length))

