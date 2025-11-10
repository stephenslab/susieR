set.seed(1)
load("../datafiles/small_sample_test2.RData")
fit <- susie(X,y,L = 1,small = TRUE,alpha0 = 0.01,beta0 = 0.01,
             verbose = TRUE)
print(diff(fit$elbo))
print(fit$V)

