# Reproducible paired simulation. Run from package root:
# Rscript inst/experiments/nig_audit/run_simulations.R package_path label reps
args = commandArgs(TRUE)
package_path = if (length(args)) args[1] else "."
label = if (length(args) > 1) args[2] else "corrected"
reps = if (length(args) > 2) as.integer(args[3]) else 200L
pkgload::load_all(package_path,quiet=TRUE)
settings = expand.grid(n=c(20,50),scenario=c("single","one_overfit","three","null"),
                       stringsAsFactors=FALSE)
p = 40
Sigma = 0.8^abs(outer(1:p,1:p,"-"))
factor = chol(Sigma)
rows = list()
idx = 0L
for (k in seq_len(nrow(settings))) {
  n = settings$n[k]
  scenario = settings$scenario[k]
  causal_count = switch(scenario,single=1,one_overfit=1,three=3,null=0)
  L = if (scenario == "single") 1 else 5
  for (r in seq_len(reps)) {
    set.seed(400000 + 1000*k + r)
    X = matrix(rnorm(n*p),n,p) %*% factor
    causal = if (causal_count) sample.int(p,causal_count) else integer()
    b = numeric(p)
    if (causal_count) b[causal] = sample(c(-1,1),causal_count,TRUE)
    # Population h2=0.5. Gaussian errors, with no in-sample rescaling.
    noise_var = if (causal_count) drop(crossprod(b,Sigma %*% b)) else 1
    y = drop(X %*% b + rnorm(n,sd=sqrt(noise_var)))
    methods = if (label == "baseline") c("gaussian","baseline_nig") else "corrected_nig"
    for (method in methods) {
      elapsed = system.time({
        fit = tryCatch(suppressMessages(suppressWarnings(susie(X,y,L=L,
          small=method != "gaussian",max_iter=200,min_abs_corr=0.5))),
          error=function(e) e)
      })[["elapsed"]]
      bad = inherits(fit,"error")
      cs = if (bad) list() else fit$sets$cs
      hits = vapply(cs,function(s) any(s %in% causal),TRUE)
      idx = idx+1L
      rows[[idx]] = data.frame(n=n,scenario=scenario,replicate=r,method=method,
        ncs=length(cs),covered=sum(hits),false_cs=sum(!hits),
        any_false=any(!hits),mean_size=if(length(cs)) mean(lengths(cs)) else NA_real_,
        causal_recovered=if(causal_count) sum(causal %in% unlist(cs)) else 0,
        causal_count=causal_count,error=bad,
        converged=if(bad) FALSE else fit$converged,
        iterations=if(bad) NA_integer_ else fit$niter,seconds=elapsed)
    }
  }
  cat("Completed",label,"n",n,scenario,"replicates",reps,"\n")
}
results = do.call(rbind,rows)
prefix = file.path("inst/experiments/nig_audit",label)
write.csv(results,paste0(prefix,"_simulations.csv"),row.names=FALSE)
summaries = lapply(split(results,interaction(results$n,results$scenario,results$method)),
  function(z) {
    coverage = if(sum(z$ncs)) sum(z$covered)/sum(z$ncs) else NA_real_
    data.frame(n=z$n[1],scenario=z$scenario[1],method=z$method[1],
      reps=nrow(z),ncs=sum(z$ncs),coverage=coverage,
      coverage_mcse=if(sum(z$ncs)) sd(z$covered-coverage*z$ncs)/
        mean(z$ncs)/sqrt(nrow(z)) else NA_real_,
      mean_cs=mean(z$ncs),mean_size=mean(z$mean_size,na.rm=TRUE),
      any_false=mean(z$any_false),
      power=if(sum(z$causal_count)) sum(z$causal_recovered)/sum(z$causal_count) else NA_real_,
      errors=sum(z$error),nonconverged=sum(!z$converged),seconds=sum(z$seconds))
  })
summary = do.call(rbind,summaries)
write.csv(summary,paste0(prefix,"_summary.csv"),row.names=FALSE)
print(summary,row.names=FALSE)
capture.output(sessionInfo(),file=paste0(prefix,"_session.txt"))
