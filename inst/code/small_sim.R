# Small script to evaluate the NIG prior version of SuSiE in
# simulated data sets.
library(matrixStats)
library(susieR)
susie_version <- packageVersion("susieR")
print(susie_version)
N <- 100
n <- 80
set.seed(1)
geno <- readRDS("../datafiles/Thyroid.FMO2.1Mb.RDS")$X
storage.mode(geno) <- "double"

causal_snps     <- vector("list",N)
res_susie       <- vector("list",N)
res_susie_small <- vector("list",N)

for (iter in 1:N) {
  cat(iter,"")

  # Subsample the genotypes.
  i <- sample(nrow(geno),n)
  X <- geno[i,]

  # Remove SNPs that show no variation in the subset.
  j <- which(colSds(X) > 0)
  X <- X[,j]

  # Simulate b and y.
  p <- ncol(X)
  b <- rep(0,p)
  names(b) <- colnames(X)
  p1   <- sample(3,1)
  j    <- sample(p,p1)
  b[j] <- sample(c(-1,1),p1,replace = TRUE)
  e    <- rnorm(n,sd = 0.1)
  y    <- drop(X %*% b + e)
  causal_snps[[iter]] <- j

  # Run susie with normal prior.
  fit1 <- susie(X,y,L = 10,standardize = FALSE,min_abs_corr = 0.5,
                estimate_prior_method = "EM",small = FALSE,
                verbose = FALSE)
  res_susie[[iter]] <- fit1$sets

  # Run susie with NIG prior. 
  fit2 <- suppressMessages(
            susie(X,y,L = 10,standardize = FALSE,min_abs_corr = 0.5,
                  estimate_prior_method = "EM",small = TRUE,
                  alpha0 = 2,beta0 = 2,verbose = FALSE))
  res_susie_small[[iter]] <- fit2$sets
}
cat("\n")

# Save the results.
save(list = c("n","susie_version","causal_snps","res_susie",
              "res_susie_small"),
     file = "small_sim_out.RData")
