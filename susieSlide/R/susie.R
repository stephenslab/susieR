#' Fit SuSiE with an empirical-Bayes heterozygote slider
#' @export
susie <- function(X, y, L = min(10, ncol(X)), scaled_prior_variance = 0.2,
                  residual_variance = NULL, prior_weights = NULL, null_weight = 0,
                  standardize = TRUE, intercept = TRUE,
                  estimate_residual_variance = TRUE,
                  estimate_prior_variance = TRUE,
                  estimate_prior_method = "optim",
                  min_obs = 5, delta = NULL, chunk_size = 1000L,
                  cache_heterozygotes = FALSE, coverage = 0.95,
                  min_abs_corr = 0.5, max_iter = 100, tol = 1e-3,
                  verbose = FALSE, ...) {
  call <- match.call()
  if ((!is.matrix(X) && !inherits(X,"sparseMatrix")) ||
      nrow(X)<2 || ncol(X)<1) stop("X must be a sample-by-SNP genotype matrix.")
  values <- if (inherits(X,"sparseMatrix")) X@x else X
  if (!is.numeric(values) || any(!is.finite(values)) || any(!values %in% 0:2))
    stop("X must contain finite hard-call genotypes 0, 1, or 2.")
  if (!is.numeric(y) || length(y)!=nrow(X)) stop("y must have one numeric value per row of X.")
  if (!is.numeric(min_obs) || length(min_obs)!=1 || !is.finite(min_obs) ||
      min_obs<0 || min_obs!=floor(min_obs)) stop("min_obs must be a nonnegative integer.")
  if (!is.numeric(chunk_size) || length(chunk_size)!=1 || !is.finite(chunk_size) ||
      chunk_size<1 || chunk_size!=floor(chunk_size)) stop("chunk_size must be a positive integer.")
  for (flag in list(standardize,intercept,cache_heterozygotes,
                    estimate_residual_variance,estimate_prior_variance,verbose))
    if (!is.logical(flag) || length(flag)!=1 || is.na(flag)) stop("Logical options must be TRUE or FALSE.")
  if (!is.numeric(L) || length(L)!=1 || !is.finite(L) || L<1 || L!=floor(L))
    stop("L must be a positive integer.")
  if (!is.numeric(max_iter) || length(max_iter)!=1 || !is.finite(max_iter) ||
      max_iter<1 || max_iter!=floor(max_iter)) stop("max_iter must be a positive integer.")
  if (!is.numeric(tol) || length(tol)!=1 || !is.finite(tol) || tol<=0)
    stop("tol must be positive and finite.")
  estimate_prior_method <- match.arg(estimate_prior_method,c("optim","EM","simple"))
  dots <- list(...)
  supported <- c("check_null_threshold","prior_tol","residual_variance_upperbound",
                 "residual_variance_lowerbound","na.rm","n_purity","median_abs_corr",
                 "estimate_residual_method","track_fit","model_init")
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) ||
                       any(!names(dots) %in% supported)))
    stop("Unsupported slider option(s): ",paste(setdiff(names(dots),supported),collapse=", "),
         ". See ?susieSlide::susie for the supported individual-level interface.")
  if (!is.null(dots$estimate_residual_method) &&
      !dots$estimate_residual_method %in% c("MoM","MLE"))
    stop("The slider currently supports Gaussian residual variance updates (MoM or MLE).")
  if (isTRUE(dots$na.rm)) {
    keep <- !is.na(y); X <- X[keep,,drop=FALSE]; y <- y[keep]
  }
  if (length(y)<2 || any(!is.finite(y)) || !is.finite(var(y)) || var(y)<=0)
    stop("y must be finite with positive variance after removing missing observations.")
  if (!is.null(coverage) && (length(coverage)!=1 || !is.finite(coverage) || coverage<=0 || coverage>1))
    stop("coverage must be NULL or in (0,1].")
  for (value in list(min_abs_corr,dots$median_abs_corr))
    if (!is.null(value) && (length(value)!=1 || !is.finite(value) || value<0 || value>1))
      stop("Purity thresholds must be NULL or in [0,1].")
  p <- ncol(X); L <- min(L,p)
  if (!is.null(delta)) {
    if (!is.numeric(delta) || any(!is.finite(delta)) || any(abs(delta)>1) ||
        (!is.null(dim(delta)) && !identical(dim(delta),c(as.integer(L),as.integer(p)))) ||
        (is.null(dim(delta)) && !length(delta) %in% c(1L,p)))
      stop("delta must be NULL, a value in [-1,1], a length-p vector, or an L-by-p matrix.")
    if (is.null(dim(delta))) delta <- matrix(rep(delta,length.out=p),L,p,byrow=TRUE)
  }
  if (is.matrix(X)) storage.mode(X) <- "double"
  warm <- dots$model_init; dots$model_init <- NULL
  args <- c(list(X=X,y=y,L=L,scaled_prior_variance=scaled_prior_variance,
    residual_variance=residual_variance,prior_weights=prior_weights,null_weight=null_weight,
    standardize=standardize,intercept=intercept,
    estimate_residual_variance=estimate_residual_variance,
    estimate_prior_variance=estimate_prior_variance,estimate_prior_method=estimate_prior_method,
    coverage=coverage,min_abs_corr=min_abs_corr,max_iter=max_iter,tol=tol,
    verbose=verbose,init_only=TRUE),dots)
  setup <- do.call(susieR::susie,args)
  data <- setup$data; params <- setup$params
  # Counts are based on the retained individuals, before any centering/scaling.
  data$chunk_size <- as.integer(min(chunk_size,max(1,floor(64e6/(8*data$n)))))
  data$H <- if(cache_heterozygotes) (data$X==1)*1 else NULL
  n1 <- numeric(data$p)
  for (j in .chunks(data$p,data$chunk_size)) n1[j] <- Matrix::colSums(.H(data,j))
  sx <- Matrix::colSums(data$X)
  n2 <- (sx-n1)/2
  data$counts <- cbind(`0`=data$n-n1-n2,`1`=n1,`2`=n2)
  data$forced <- apply(data$counts,1,min)<min_obs
  data$min_obs <- min_obs
  data$hmean <- if(intercept) n1/data$n else numeric(data$p)
  data$xmean <- attr(data$X,"scaled:center")
  data$scale <- attr(data$X,"scaled:scale")
  data$xx <- pmax(0,attr(data$X,"d"))
  data$xh <- (n1-data$n*data$xmean*data$hmean)/data$scale^2
  data$hh <- pmax(0,(n1-data$n*data$hmean^2)/data$scale^2)
  if (!is.null(delta) && data$p>p) delta <- cbind(delta,0)
  data$fixed_delta <- delta
  data$input_p <- p
  data$warm <- warm
  class(data) <- c("slide_individual",class(data))
  fit <- susieR::susie_workhorse(data,params)
  fit$call <- call
  fit$delta_method <- if(is.null(delta)) "empirical Bayes (plug-in)" else "fixed"
  fit$log_bf_type <- "conditional on fitted delta; beta integrated analytically"
  fit$objective_type <- "ELBO conditional on component-by-SNP deltas"
  fit$delta_cs <- slider_cs_table(fit)
  fit
}

.chunks <- function(p,size) lapply(seq.int(1L,p,by=size),function(i) seq.int(i,min(p,i+size-1L)))
.H <- function(data,j) {
  if (!is.null(data$H)) data$H[,j,drop=FALSE] else (data$X[,j,drop=FALSE]==1)*1
}
.h_score <- function(data,r) {
  ans <- numeric(data$p)
  for(j in .chunks(data$p,data$chunk_size)) ans[j] <- as.numeric(Matrix::crossprod(.H(data,j),r))
  (ans-data$hmean*sum(r))/data$scale
}
.slide_product <- function(data,b,bh) {
  result <- .engine("compute_Xb")(data$X,b)
  active <- which(bh!=0)
  if(length(active)) {
    for(k in .chunks(length(active),data$chunk_size)) {
      j <- active[k]
      result <- result+as.numeric(.H(data,j) %*% (bh[j]/data$scale[j]))
    }
    result <- result-sum(data$hmean*bh/data$scale)
  }
  result
}
