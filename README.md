# susieR

an R package for "sum of single effects" regression

This is very much work in progress, and the interface will likely change under you! 
If you want to use it, I recommend
you contact me at mstephens@uchicago.edu.


# Quick Start

This fits a sparse linear regression model with up to $L$ non-zero effects.
Generally there is no harm in over-stating $L$ (that is, the method
is pretty robust to overfitting) except that computation will grow as $L$ grows.

Here is a minimal example:
```
devtools::install_github("stephenslab/susieR")
set.seed(1)
n = 1000
p = 1000
beta = rep(0,p)
beta[1] = 1
beta[2] = 1
beta[3] = 1
beta[4] = 1
X = matrix(rnorm(n*p),nrow=n,ncol=p)
y = X %*% beta + rnorm(n)
res =susie(X,y,niter=20,L=10,calc_elbo = TRUE)
coef(res)
plot(y,predict(res))
```
