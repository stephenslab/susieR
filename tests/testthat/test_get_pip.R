 context("test_get_pip.R")
 test_that("PIP computation is correct in light of null weight and estimated prior",{
   res = list(alpha = matrix(c(rep(1,10),rep(2,10),rep(3,10))/10, 3, 10, byrow=T), sets=list(cs_index = c(2)), null_index = 10, V = rep(1,3))
   class(res) = 'susie'
   expect_equal(susie_get_pip(res), rep(0.496, 9))
   #
   expect_equal(susie_get_pip(res, prune_by_cs=T), rep(0.2, 9))
   res$V = rep(0,3)
   expect_equal(susie_get_pip(res), rep(0,9))
   res$V = rep(1,3)
   res$sets$cs_index=NULL
   expect_equal(susie_get_pip(res,prune_by_cs=T), rep(0,9))
 })
