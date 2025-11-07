set.seed(1)

## ----load-data---------------------------------------
library(susieR)
data(data_small)
y <- data_small$y
X <- data_small$X
dim(X)

## ----run-susie, fig.height=5, fig.width=5, message=FALSE----
t0 <- proc.time()
res_susie <- susie(X,y,L = 10,min_abs_corr = 0,verbose = FALSE)
t1 <- proc.time()
print(t1 - t0)
print(res_susie$sets$cs)
print(round(res_susie$V,digits = 4))
print(susie_plot(res_susie,y = "PIP"))

# susie with small = TRUE, L = 1.
out <- susie(X,y,L = 1,small = TRUE,alpha0 = 0.1,beta0 = 0.1,
             verbose = TRUE)
print(out$sets$cs)
print(out$V)

## ----run-susie-small, message=FALSE, warning=FALSE----
t0 <- proc.time()
res_susie_small <- susie(X,y,L = 10,small = TRUE,min_abs_corr = 0,
                         alpha0 = 0.1,beta0 = 0.1,verbose = TRUE)
t1 <- proc.time()
print(t1 - t0)
# print(res_susie_small$sets$cs)
print(round(res_susie_small$V,digits = 4))

## ----------------------------------------------------
print(susie_plot(res_susie_small,y = "PIP"))

