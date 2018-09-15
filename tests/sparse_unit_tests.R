devtools::install_local("../susieR-sparse-v3") #relative path to susieR repo
library(testthat)
library(susieR)

test_results <- test_dir("../susieR-sparse-v3/tests/sparse_unit_tests", reporter="summary")
