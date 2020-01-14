## results from original susie
devtools::install_github("stephenslab/susieR")
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
n = 100
p = 200
beta = rep(0,p)
beta[1]    = 10
beta[2]  = 10
beta[3]  = 10
beta[4] = 10
X.dense = create_sparsity_mat(0.99,n,p)
y = c(X.dense %*% beta + rnorm(n))
L = 10
residual_variance = 0.8
scaled_prior_variance = 0.2
s = list(alpha=matrix(1/p,nrow=L,ncol=p),
         mu=matrix(2,nrow=L,ncol=p),
         mu2=matrix(3,nrow=L,ncol=p),
         Xr=rep(5,n), KL=rep(1.2,L),
         sigma2=residual_variance, V=scaled_prior_variance * as.numeric(var(y)))
X = susieR:::set_X_attributes(X.dense)
Eb = rep(1, p)
Eb2 = rep(1, p)
s2 = residual_variance
V = scaled_prior_variance


objective.original.res = susieR::susie_get_objective(s)
saveRDS(objective.original.res, 'objective_original_res.rds')

Eloglik.original.res = susieR:::Eloglik(X,y,s)
saveRDS(Eloglik.original.res, 'Eloglik_original_res.rds')

ER2.original.res = susieR:::get_ER2(X,y,s)
saveRDS(ER2.original.res, 'ER2_original_res.rds')

SER.original.res = susieR:::SER_posterior_e_loglik(X,y,s2,Eb,Eb2)
saveRDS(SER.original.res, 'SER_original_res.rds')

singleReg.original.res = susieR:::single_effect_regression(y,X,V)
saveRDS(singleReg.original.res, 'singleReg_original_res.rds')

vbupdate.original.res = susieR:::update_each_effect(X, y, s)
saveRDS(vbupdate.original.res, 'vbupdate_original_res.rds')

susiefit.original.res = susie(X.dense,y)
saveRDS(susiefit.original.res, 'susiefit_original_res.rds')

susiefit.original.res2 = susie(X.dense, y, standardize = TRUE, intercept = FALSE)
susiefit.original.res3 = susie(X.dense, y, standardize = FALSE, intercept = TRUE)
susiefit.original.res4 = susie(X.dense, y, standardize = FALSE, intercept = FALSE)
saveRDS(susiefit.original.res2, 'susiefit_original_res2.rds')
saveRDS(susiefit.original.res3, 'susiefit_original_res3.rds')
saveRDS(susiefit.original.res4, 'susiefit_original_res4.rds')
