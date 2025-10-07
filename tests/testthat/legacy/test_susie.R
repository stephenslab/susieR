context("test_susie.R")

# test_that("susie agrees with version 0.3", with(simulate(sparse=T), {
#   original.res  = readRDS('susiefit_original_res.rds')
#   original.res2 = readRDS('susiefit_original_res2.rds')
#   original.res3 = readRDS('susiefit_original_res3.rds')
#   original.res4 = readRDS('susiefit_original_res4.rds')
#
#   original.res$Xr = as.vector(original.res$Xr)
#   original.res2$Xr = as.vector(original.res2$Xr)
#   original.res3$Xr = as.vector(original.res3$Xr)
#   original.res4$Xr = as.vector(original.res4$Xr)
#   dense.res = susie(X, y, tol=1E-2, estimate_prior_variance = FALSE)
#   sparse.res = susie(X.sparse, y, tol=1E-2, estimate_prior_variance = FALSE)
#
#   dense.res2 = susie(X, y, standardize=TRUE, intercept = FALSE, tol=1E-2, estimate_prior_variance = FALSE)
#   sparse.res2 = susie(X.sparse, y, standardize=TRUE, intercept = FALSE,
#                       tol=1E-2, estimate_prior_variance = FALSE)
#
#   dense.res3 = susie(X, y, standardize=FALSE, intercept = TRUE, tol=1E-2, estimate_prior_variance = FALSE)
#   sparse.res3 = susie(X.sparse, y, standardize=FALSE, intercept = TRUE,
#                       tol=1E-2, estimate_prior_variance = FALSE)
#
#   dense.res4 = susie(X, y, standardize=FALSE, intercept = FALSE, tol=1E-2, estimate_prior_variance = FALSE)
#   sparse.res4 = susie(X.sparse, y, standardize=FALSE, intercept = FALSE,
#                       tol=1E-2, estimate_prior_variance = FALSE)
#   expect_equal_susie(sparse.res, original.res)
#   expect_equal_susie(dense.res, original.res)
#   expect_equal_susie(sparse.res2, original.res2)
#   expect_equal_susie(dense.res2, original.res2)
#   expect_equal_susie(sparse.res3, original.res3)
#   expect_equal_susie(dense.res3, original.res3)
#   expect_equal_susie(sparse.res4, original.res4)
#   expect_equal_susie(dense.res4, original.res4)
# }))
