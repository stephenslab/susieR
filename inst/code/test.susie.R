library(Matrix)
library(susieR)

create_sparsity_mat = function(sparsity, n, p){
  nonzero = round(n*p*(1-sparsity))
  nonzero.idx = sample(n*p, nonzero)
  mat = numeric(n*p)
  mat[nonzero.idx] = 1
  mat = matrix(mat, nrow=n, ncol=p)
  return(mat)     
}

set.seed(1)
n = 1000
p = 10000
beta = rep(0,p)
beta[1]    = 10 
beta[300]  = 10
beta[400]  = 10
beta[1000] = 10

X.dense = create_sparsity_mat(0.99,n,p)
X.sparse = as(X.dense,'dgCMatrix')

y = c(X.dense %*% beta + rnorm(n))

set.seed(1)
cat("Running with dense X:\n")
timing <- system.time(susie.dense.fit <- susie(X.dense,y,verbose = TRUE,
                                               tol = 1e-4))
print(timing)
cat("Running with sparse X:\n")
set.seed(1)
timing <- system.time(susie.sparse.fit <- susie(X.sparse,y,verbose = TRUE,
                                                tol = 1e-4))
print(timing)
