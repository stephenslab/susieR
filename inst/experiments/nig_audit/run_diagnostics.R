# Run from the package root. Optional arguments: package_path output_prefix.
args = commandArgs(TRUE)
package_path = if (length(args)) args[1] else "."
prefix = if (length(args) > 1) args[2] else "inst/experiments/nig_audit/corrected"
pkgload::load_all(package_path,quiet = TRUE)
quiet_fit = function(...) suppressMessages(suppressWarnings(susie(...)))
set.seed(20261110)
n = 30
p = 20
X = matrix(rnorm(n*p),n,p)
X[,2] = 0.8*X[,1] + 0.6*X[,2]
y = 0.7*X[,1] + rnorm(n)
fit = function(X,y,...) quiet_fit(X,y,small = TRUE,L = 1,
                                 min_abs_corr = 0,max_iter = 500,...)
f = fit(X,y,estimate_prior_variance = FALSE)
f_y = fit(X,10*y,estimate_prior_variance = FALSE)
f_x = fit(sweep(X,2,seq(0.2,3,length.out=p),"*"),y,
          estimate_prior_variance = FALSE)
e = fit(X,y)
e_y = fit(X,10*y)
e_x = fit(sweep(X,2,seq(0.2,3,length.out=p),"*"),y)
g = fit(X,y,prior_weights = c(0,rep(1,p-1)))
null = fit(cbind(X,0),y,prior_weights = c(rep(0.01,p),0.8),
           scaled_prior_variance = 0,estimate_prior_variance = FALSE)
out = data.frame(check = c("fixed_y_scale_PIP","fixed_X_units_PIP",
                           "estimated_y_scale_PIP","estimated_X_units_PIP",
                           "fixed_y_scale_mean","fixed_y_scale_second_moment",
                           "zero_prior_PIP","zero_effect_null_probability",
                           "estimated_converged","estimated_y_converged"),
                 value = c(max(abs(f$alpha-f_y$alpha)),
                           max(abs(f$alpha-f_x$alpha)),
                           max(abs(e$alpha-e_y$alpha)),
                           max(abs(e$alpha-e_x$alpha)),
                           max(abs(f$mu-f_y$mu/10)),
                           max(abs(f$mu2-f_y$mu2/100)),g$alpha[1,1],
                           null$alpha[1,p+1],e$converged,e_y$converged))
write.csv(out,paste0(prefix,"_diagnostics.csv"),row.names=FALSE)
print(out)
saveRDS(list(fixed=f,estimated=e),paste0(prefix,"_fits.rds"))
