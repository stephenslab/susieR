# Run in a fresh R session from the corrected repository root after installation.
stopifnot(read.dcf('DESCRIPTION')[1,'Package']=='susieSlide',
          !file.exists('susieSlide/DESCRIPTION'))
library(susieSlide)
stopifnot(packageVersion('susieSlide')>='0.2.0',
          !'susieR' %in% loadedNamespaces(),
          !'susieR' %in% names(getNamespaceImports('susieSlide')))
cat('Installed package:',find.package('susieSlide'),'\n')
cat('Version:',as.character(packageVersion('susieSlide')),'\n')
cat('Root package and namespace ownership: PASS\n')
examples <- readRDS('validation/comparison_examples.rds')
modes <- c('Additive','Recessive','Partially_recessive','Dominant','Partially_dominant','Mixed')
rows <- lapply(modes,function(mode) {
  example <- examples[[paste(mode,2,'Slider',sep='_')]]
  fit <- suppressMessages(susieSlide::susie(example$X,example$y,L=2,
    standardize=TRUE,residual_variance=.36,estimate_residual_variance=FALSE,
    scaled_prior_variance=.5/var(example$y),estimate_prior_variance=FALSE,
    max_iter=400,tol=1e-7,min_abs_corr=0))
  fields <- c('alpha','mu','mu2','delta','lbf_variable','V','sigma2','fitted','pip','elbo')
  differences <- vapply(fields,function(field) max(abs(fit[[field]]-example$fit[[field]])),numeric(1))
  stopifnot(fit$converged,max(differences)<1e-8,
            max(abs(predict(fit,newx=example$X)-fit$fitted))<1e-10,
            identical(dim(fit$delta),dim(fit$alpha)),
            isTRUE(all.equal(summary(fit)$delta,fit$delta_cs)))
  data.frame(mode=mode,converged=fit$converged,max_field_difference=max(differences),
             mean_function_mse=mean((predict(fit,newx=example$Xtest)-example$mean_test)^2))
})
result <- do.call(rbind,rows)
stopifnot(!'susieR' %in% loadedNamespaces())
print(result,row.names=FALSE)
cat('All six scenario regressions and slider predictions: PASS\n')
cat('susieR namespace loaded:', 'susieR' %in% loadedNamespaces(),'\n')
write.csv(result,'validation/root-scenario-regression.csv',row.names=FALSE)
