#create an s fit for trend filtering from orginal SuSiE version
devtools::install_github("stephenslab/susieR")
library(susieR)
n = 1000 
D = diag(-1, n)
for (i in 1:(n-1)){
  D[i, i+1] = 1
}
#order 0 fit
set.seed(1)
beta = c(rep(0,100),rep(1,100),rep(3,100),rep(-2,100),rep(0,600))
y = beta + rnorm(n)
D1inv = solve(D)
s0 = susie(D1inv,y,L=10)
saveRDS(s0, "tf_original_s0.rds")

#order 1 fit
set.seed(1)
beta = numeric(n)
for (i in 1:n){
  if (i <= 100){
    beta[i] = 0.001*i + 2
  } else if (i <= 300){
    beta[i] = 5*0.001*i + 1.6
  } else{
    beta[i] = 6.1 - 10*0.001*i
  }
}
y = beta + rnorm(n)
D2inv = solve(D%*%D)
s1 = susie(D2inv,y,L=10)
saveRDS(s1, "tf_original_s1.rds")

#order 2 fit
set.seed(1)
beta = numeric(n)
for (i in 1:n){
  if (i <= 100){
    beta[i] = (0.001*i)^2
  } else if (i <= 700){
    beta[i] = -5*(0.001*i)^2 + 0.06
  } else{
    beta[i] = 3*(0.001*i)^2 - 3.86
  }
}
y = beta + rnorm(n)
D3inv = solve(D%*%D%*%D)
s2 = susie(D3inv,y,L=10)
saveRDS(s2, "tf_original_s2.rds")

