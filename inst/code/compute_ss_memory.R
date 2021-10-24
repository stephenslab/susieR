# export MEM_CHECK_INTERVAL=0.01
# python3 monitor_memory.py Rscript compute_ss_memory.R
#
# NOTES:
#
# - *Add notes here.*
#
# library(susieR)
set.seed(1)
devtools::load_all()
load("susie_rss_data.RData")
cat("Size of X:\n")
print(object.size(X),unit = "GB")
cat("Running compute_ss.\n")
out <- compute_ss(X,y,standardize = TRUE)
