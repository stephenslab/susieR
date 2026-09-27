# Compare separately installed susieR, original susieSlide, and finite-prior
# susieRSlidePrior. Workers are sequential: package load/startup is outside timing.
# Run with Rscript; see benchmark_settings.txt for the exact measured scope.
# Required: SLIDER_BENCH_PRIOR_LIB (library containing susieRSlidePrior >= 0.3.0).
# Optional: SLIDER_BENCH_LEGACY_LIB (original library, otherwise default),
# SLIDER_BENCH_MB_LIB, SLIDER_BENCH_OUT, SLIDER_BENCH_REPS (default 5),
# SLIDER_BENCH_PILOT=true (only the smallest data set),
# SLIDER_BENCH_LEGACY_SOURCE (original R/slider.R, to verify its interface).
# Install the finite-prior package with R CMD INSTALL --preclean first;
# pkgload development objects may otherwise retain -O0 compiler flags.
# Pass --summarize to aggregate existing worker CSVs without repeating fits.

Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
argv <- commandArgs(trailingOnly=TRUE)
outdir <- Sys.getenv("SLIDER_BENCH_OUT","validation/slider_prior_benchmark")
prior_lib <- Sys.getenv("SLIDER_BENCH_PRIOR_LIB")
legacy_lib <- Sys.getenv("SLIDER_BENCH_LEGACY_LIB",.libPaths()[1])
mb_lib <- Sys.getenv("SLIDER_BENCH_MB_LIB",prior_lib)
repetitions <- as.integer(Sys.getenv("SLIDER_BENCH_REPS","5"))
pilot <- identical(Sys.getenv("SLIDER_BENCH_PILOT"),"true")
stopifnot(nzchar(prior_lib),is.finite(repetitions),repetitions>=1L)
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
cases <- data.frame(n=c(500L,500L,1000L),p=c(1000L,1000L,3000L),
                    variance=c("fixed","estimated","fixed"))
if(pilot) cases <- cases[cases$n==500L,,drop=FALSE]
cases$L <- 10L
cases$case <- paste0("n",cases$n,"_p",cases$p,"_L",cases$L,"_",cases$variance)
methods <- c("susie","susie_slide","susie_slide_prior")

if(length(argv)==2L && argv[1]=="--worker") {
  method <- argv[2]
  stopifnot(method %in% methods)
  if(method=="susie") {
    ns <- loadNamespace("susieR")
  } else {
    ns <- loadNamespace(if(method=="susie_slide") "susieSlide" else "susieRSlidePrior",
                        lib.loc=if(method=="susie_slide") legacy_lib else prior_lib)
  }
  fit_fun <- get("susie",envir=ns)
  stopifnot(identical("delta_prior" %in% names(formals(fit_fun)),method=="susie_slide_prior"))
  mb_ns <- loadNamespace("microbenchmark",lib.loc=mb_lib)
  mb_fun <- get("microbenchmark",envir=mb_ns)
  package_path <- getNamespaceInfo(ns,"path")
  package_name <- getNamespaceName(ns)
  package_version <- as.character(getNamespaceVersion(ns))
  baseline_source <- Sys.getenv("SLIDER_BENCH_LEGACY_SOURCE")
  baseline_verified <- NA
  if(method=="susie_slide" && nzchar(baseline_source)) {
    baseline <- new.env(parent=baseenv())
    sys.source(baseline_source,envir=baseline)
    baseline_verified <- identical(formals(fit_fun),formals(baseline$susie)) &&
      identical(deparse(body(fit_fun)),deparse(body(baseline$susie)))
    stopifnot(baseline_verified)
  }
  results <- diagnostics <- list()
  cat("Worker:",method,";",package_name,package_version,"\n"); flush.console()
  for(i in seq_len(nrow(cases))) {
    cfg <- cases[i,]
    set.seed(271000L+cfg$n+cfg$p)
    X <- matrix(rbinom(cfg$n*cfg$p,2,.3),cfg$n,cfg$p)
    storage.mode(X) <- "double"
    colnames(X) <- paste0("rs",seq_len(cfg$p))
    causal <- c(15L,as.integer(cfg$p/2),cfg$p-20L)
    delta <- c(-.5,0,.5); beta <- c(.9,-.7,.8)
    signal <- drop((X[,causal,drop=FALSE]+
      sweep((X[,causal,drop=FALSE]==1)*1,2,delta,"*")) %*% beta)
    y <- .3+signal+rnorm(cfg$n,sd=.8)
    estimated <- cfg$variance=="estimated"
    args <- list(X=X,y=y,L=cfg$L,standardize=TRUE,intercept=TRUE,
      scaled_prior_variance=.5/var(y),residual_variance=.64,
      estimate_prior_variance=estimated,estimate_residual_variance=estimated,
      estimate_prior_method="optim",max_iter=300,tol=1e-6,coverage=NULL,verbose=FALSE)
    if(method!="susie") args$min_obs <- 5L
    if(method=="susie_slide_prior") args$delta_prior <- rep(1/17,17)
    run_fit <- function() suppressMessages(do.call(fit_fun,args))
    warm_seconds <- system.time(warm <- run_fit())["elapsed"]
    stopifnot(isTRUE(warm$converged),all(is.finite(warm$fitted)),
              all(is.finite(warm$elbo)))
    if(method=="susie_slide_prior")
      stopifnot(isTRUE(all.equal(warm$delta_prior,rep(1/17,17))))
    diag <- data.frame(case=cfg$case,n=cfg$n,p=cfg$p,L=cfg$L,variance=cfg$variance,
      method=method,package=package_name,version=package_version,
      warm_seconds=unname(warm_seconds),niter=warm$niter,converged=warm$converged,
      sigma2=warm$sigma2,active_effects=sum(warm$V>1e-9),
      forced_snps=if(method=="susie") NA_integer_ else sum(warm$delta_forced),
      returned_fit_mb=as.numeric(object.size(warm))/1024^2,
      final_elbo=tail(warm$elbo,1),baseline_interface_verified=baseline_verified)
    diagnostics[[i]] <- diag
    cat(cfg$case,": warm-up",round(warm_seconds,3),"s;",warm$niter,
        "iterations; timing",repetitions,"repeats\n"); flush.console()
    rm(warm); invisible(gc())
    # The only timed expression is the complete fit. Object setup, simulation,
    # namespace loading, output writing and the warm-up are outside timing.
    mb <- mb_fun(run_fit(),times=repetitions,unit="ms",
                 control=list(order="inorder"),setup=invisible(gc()))
    results[[i]] <- data.frame(case=cfg$case,n=cfg$n,p=cfg$p,L=cfg$L,
      variance=cfg$variance,method=method,repetition=seq_len(nrow(mb)),
      milliseconds=mb$time/1e6,niter=diag$niter)
    saveRDS(mb,file.path(outdir,paste0(method,"_",cfg$case,"_microbenchmark.rds")))
    write.csv(do.call(rbind,results),file.path(outdir,paste0(method,"_raw.csv")),row.names=FALSE)
    write.csv(do.call(rbind,diagnostics),file.path(outdir,paste0(method,"_diagnostics.csv")),row.names=FALSE)
    cat("Median",round(median(mb$time)/1e6,2),"ms\n"); flush.console()
    rm(X,y,args,mb); invisible(gc())
  }
  writeLines(c(paste("Method:",method),paste("Package path:",package_path),
    paste("Baseline interface verified:",baseline_verified),
    paste("OMP_NUM_THREADS:",Sys.getenv("OMP_NUM_THREADS")),
    paste("Processor:",Sys.getenv("PROCESSOR_IDENTIFIER")),
    paste("Logical processors:",Sys.getenv("NUMBER_OF_PROCESSORS")),
    capture.output(sessionInfo())),file.path(outdir,paste0(method,"_session.txt")))
  quit(status=0)
}

stopifnot(length(argv)==0L || identical(argv,"--summarize"))
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
stopifnot(length(script_arg)==1L)
script_path <- normalizePath(sub("^--file=","",script_arg),winslash="/",mustWork=TRUE)
rscript <- file.path(R.home("bin"),"Rscript.exe")
if(!file.exists(rscript)) rscript <- file.path(R.home("bin"),"Rscript")
set.seed(272026)
worker_order <- sample(methods)
if(length(argv)==0L) writeLines(c(
  "Complete individual-level fits; actual separately installed implementations.",
  paste("Sequential worker order:",paste(worker_order,collapse=", ")),
  paste("Repetitions per method and scenario:",repetitions),
  "Genotypes are iid Binomial(2, 0.3), complete 0/1/2 hard calls.",
  "Three causal SNPs: effects (0.9,-0.7,0.8), sliders (-0.5,0,0.5); noise SD 0.8.",
  "Matched L=10; intercept=TRUE; original additive-SD scaling; tol=1e-6; max_iter=300.",
  "Initial coefficient-prior variance 0.5 and noise variance 0.64 for every method.",
  "Fixed: neither Gaussian variance estimated. Estimated: both estimated, prior method optim.",
  "The new slider probabilities remain fixed at 1/17 in both variance settings.",
  "Slider defaults: min_obs=5, chunk_size=1000, cache_heterozygotes=FALSE.",
  "coverage=NULL omits credible-set construction/purity in all methods.",
  "One full fit warm-up per method/case; convergence verified before timing.",
  "microbenchmark uses its default timer warmup; setup=invisible(gc()) is outside each timing.",
  "Includes input preparation, IBSS, fitting output construction; excludes simulation, startup and I/O.",
  "All workers run sequentially with OMP/MKL/OpenBLAS thread environment variables set to 1.",
  "One simulated data set per size; method/variance settings use exactly the same X and y.",
  "Different iteration counts reflect different models. Per-iteration times include amortized setup/output.",
  "Returned-fit memory is object.size, not peak process memory."
),file.path(outdir,"benchmark_settings.txt"))
if(length(argv)==0L) for(method in worker_order) {
  status <- system2(rscript,c("--vanilla",shQuote(script_path),"--worker",method))
  if(status!=0L) stop("Benchmark worker failed: ",method," (",status,").")
}
raw <- do.call(rbind,lapply(methods,function(method)
  read.csv(file.path(outdir,paste0(method,"_raw.csv")))))
diagnostics <- do.call(rbind,lapply(methods,function(method)
  read.csv(file.path(outdir,paste0(method,"_diagnostics.csv")))))
summary <- do.call(rbind,lapply(split(raw,list(raw$case,raw$method),drop=TRUE),function(x) {
  info <- x[1,c("case","n","p","L","variance","method","niter")]
  cbind(info,repetitions=nrow(x),median_ms=median(x$milliseconds),
    q25_ms=unname(quantile(x$milliseconds,.25)),q75_ms=unname(quantile(x$milliseconds,.75)),
    min_ms=min(x$milliseconds),max_ms=max(x$milliseconds),
    median_ms_per_iteration=median(x$milliseconds)/x$niter[1])
}))
summary$ratio_to_susie <- NA_real_
summary$ratio_to_slider <- NA_real_
for(case in unique(summary$case)) {
  idx <- which(summary$case==case)
  additive <- summary$median_ms[idx[summary$method[idx]=="susie"]]
  slider <- summary$median_ms[idx[summary$method[idx]=="susie_slide"]]
  summary$ratio_to_susie[idx] <- summary$median_ms[idx]/additive
  summary$ratio_to_slider[idx] <- summary$median_ms[idx]/slider
}
summary <- summary[order(summary$n,summary$p,summary$variance,summary$method),]
write.csv(raw,file.path(outdir,"benchmark_raw.csv"),row.names=FALSE)
write.csv(diagnostics,file.path(outdir,"benchmark_diagnostics.csv"),row.names=FALSE)
write.csv(summary,file.path(outdir,"benchmark_summary.csv"),row.names=FALSE)
saveRDS(list(raw=raw,diagnostics=diagnostics,summary=summary,cases=cases),
        file.path(outdir,"benchmark_results.rds"))
print(summary[,c("n","p","variance","method","niter","median_ms",
                  "ratio_to_susie","ratio_to_slider")],row.names=FALSE)
