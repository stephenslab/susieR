library(susieR)
set.seed(1)
p <- 4000
n <- 1000
X <- matrix(rnorm(n*p),n,p)
y <- rnorm(n)
ss <- susieR:::univariate_regression(X,y)
z <- ss$betahat/ss$sebetahat
R <- cor(X)
out <- susie_rss(z,R,estimate_prior_variance = FALSE)
