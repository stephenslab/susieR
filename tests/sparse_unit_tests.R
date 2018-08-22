devtools::install_local("~/Desktop/M/M-github-repos/susieR-sparse-v3")
library(testthat)
library(susieR)

test_results <- test_dir("~/Desktop/M/M-github-repos/susieR-sparse-v3/tests/sparse_unit_tests", reporter="summary")
