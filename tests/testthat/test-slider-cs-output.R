test_that("CS delta output follows reported components and includes nonmember SNPs", {
  # The middle component was filtered out; the remaining CSs were reordered.
  fit <- list(
    delta=rbind(c(-.8,-.6,-.4,0),c(-.2,0,.2,0),c(.4,.6,.8,0)),
    alpha=matrix(.25,3,4,dimnames=list(NULL,c("rsA","rsB","rsC","null"))),
    pip=c(.5,.6,.7),delta_forced=rep(FALSE,4),null_index=4L,
    sets=list(cs=list(L3=2L,L1=c(1L,3L)),cs_index=c(3L,1L)))
  output <- .slider_cs_output(fit)
  expected <- rbind(L3=c(.4,.6,.8),L1=c(-.8,-.6,-.4))
  colnames(expected) <- c("rsA","rsB","rsC")
  expect_named(output,c("summary","delta"))
  expect_identical(output$delta,expected)
  expect_equal(output$summary$cs,c("L3","L1","L1"))
  expect_equal(output$summary$component,c(3L,1L,1L))
  expect_equal(output$summary$snp,c("rsB","rsA","rsC"))
  expect_equal(output$summary$delta,c(.6,-.8,-.4))
  expect_equal(output$summary,slider_cs_table(fit))
})

test_that("CS delta output retains matrix shape for one SNP and no reported CSs", {
  X <- matrix(rep(0:2,each=20),ncol=1)
  y <- X[,1]+rep(c(-.1,.1),30)
  for(null_weight in c(0,.1)) {
    fit <- quiet(susieSlide::susie,X,y,L=1,delta=.25,null_weight=null_weight,
      estimate_prior_variance=FALSE,estimate_residual_variance=FALSE)
    expect_equal(fit$delta_cs$delta,
                 matrix(.25,1,1,dimnames=list("L1","SNP1")))
    expect_equal(fit$delta_cs$summary$delta,.25)

    fit <- quiet(susieSlide::susie,X,y,L=1,coverage=NULL,null_weight=null_weight,
      estimate_prior_variance=FALSE,estimate_residual_variance=FALSE)
    expect_s3_class(fit$delta_cs$summary,"data.frame")
    expect_equal(nrow(fit$delta_cs$summary),0L)
    expect_identical(dim(fit$delta_cs$delta),c(0L,1L))
    expect_identical(colnames(fit$delta_cs$delta),"SNP1")
  }
})
