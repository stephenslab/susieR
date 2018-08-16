#' @title Local false sign rate (lfsr) for susie confidence sets
#' @details This computes the average lfsr across SNPs for each l, weighted by the
#' posterior inclusion probability alpha
#' @param res a susie fit, the output of `susieR::susie()`
#' @return an l vector of lfsr for confidence sets
#' @export
susie_get_lfsr = function(res){
  pos_prob = pnorm(0,mean=t(res$mu),sd=sqrt(res$mu2-res$mu^2))
  neg_prob = 1-pos_prob
  1-rowSums(res$alpha*t(pmax(pos_prob,neg_prob)))
}

#find how many variables in the CS
# x is a probability vector
n_in_CS_x = function(x, coverage = 0.9){
  sum(cumsum(sort(x,decreasing = TRUE))<coverage)+1
}

# return binary vector indicating if each point is in CS
# x is a probability vector
in_CS_x = function(x, coverage = 0.9){
  n = n_in_CS_x(x, coverage)
  o = order(x,decreasing=TRUE)
  result = rep(0,length(x))
  result[o[1:n]] = 1
  return(result)
}

# returns an l by p binary matrix
# indicating which variables are in susie confidence sets
in_CS = function(res, coverage = 0.9){
  if (class(res) == "susie")
    res = res$alpha
  t(apply(res,1,function(x) in_CS_x(x, coverage)))
}

n_in_CS = function(res, coverage = 0.9){
  if (class(res) == "susie")
    res = res$alpha
  apply(res,1,function(x) n_in_CS_x(x, coverage))
}

get_purity = function(pos, X, Xcorr, n = 100) {
  if (length(pos) == 1) {
    c(1,1,1)
  } else {
    if (length(pos) > n) pos = sample(pos, n)
    if (is.null(Xcorr)) {
      X_sub = X[,pos]
      if (length(pos) > n) {
        # remove identical columns
        pos_rm = sapply(1:ncol(X_sub), function(i) all( abs(X_sub[,i] - mean(X_sub[,i])) < .Machine$double.eps ^ 0.5 ))
        if (length(pos_rm)) {
          X_sub = X_sub[,-pos_rm]
        }
      }
      value = abs(cor(X_sub))
    } else {
      value = abs(Xcorr[pos, pos])
    }
    c(min(value), mean(value), median(value))
  }
}

#' @title Extract confidence sets from SuSiE model
#' @details It reports indices of variables in each confidence set identified,
#' as well as summaries of correlation between variables within each set.
#' @param res a susie fit, the output of `susieR::susie()`, or simply the posterior
#' inclusion probability matrix alpha.
#' @param X N by P matrix of variables.
#' When provided, correlation between variables will be computed and used to remove
#' confidence sets whose minimum correlation between variables is smaller than `min_abs_corr` (see below).
#' @param Xcorr P by P matrix of correlations between variables.
#' when provided, it will be used to remove confidence sets
#' whose minimum correlation between variables is smaller than `min_abs_corr` (see below).
#' @param coverage coverage of confident sets. Default to 0.9 for 90\% confidence interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a confidence set.
#' Default set to 0.5 to correspond to squared correlation of 0.25,
#' a commonly used threshold for genotype data in genetics studies.
#' @return a list of `cs`, and additionally `purity` and selected `cs_index` if `X` or `Xcorr` was provided.
#' @export
susie_get_CS = function(res,
                        X = NULL, Xcorr = NULL,
                        coverage = 0.9,
                        min_abs_corr = 0.5,
                        dedup = TRUE) {
  if (class(res) == "susie")
    res = res$alpha
  if (!is.null(X) && !is.null(Xcorr)) {
    stop("Only one of X or Xcorr should be specified")
  }
  if (!is.null(X) && ncol(X) != length(res[1,])) {
    stop("Column of X matrix does not match length of res posterior")
  }
  if (!is.null(Xcorr) && ncol(Xcorr) != length(res[1,])) {
    stop("Dimension of Xcorr matrix does not match length of res posterior")
  }
  if (!is.null(Xcorr) && !isSymmetric(Xcorr)) {
    stop("Xcorr matrix must be symmetric")
  }
  # L by P binary matrix
  status = in_CS(res, coverage)
  # an L list of CS positions
  cs = lapply(1:nrow(status), function(i) which(status[i,]!=0))
  cs = cs[lapply(cs, length) > 0]
  # FIXME: see issue 21
  # https://github.com/stephenslab/susieR/issues/21
  if (dedup) cs = cs[!duplicated(cs)]
  # compute and filter by "purity"
  if (is.null(Xcorr) && is.null(X)) {
    return(list(cs=cs))
  } else {
    purity = data.frame(do.call(rbind, lapply(1:length(cs), function(i) get_purity(cs[[i]], X, Xcorr))))
    colnames(purity) = c('min.abs.corr', 'mean.abs.corr', 'median.abs.corr')
    is_pure = which(purity$min.abs.corr > min_abs_corr)
    if (length(is_pure) > 0) {
      cs = cs[is_pure]
      purity = purity[is_pure,]
      row_names = paste0("L", is_pure)
      names(cs) = row_names
      rownames(purity) = row_names
      ## re-order CS list and purity rows based on purity
      ordering = order(purity$min.abs.corr, decreasing=T)
      return(list(cs = cs[ordering], purity = purity[ordering,], cs_index = is_pure[ordering]))
    } else {
      return(list(cs = NULL))
    }
  }
}

#' @title Compute posterior inclusion probability (PIP) for all variables
#' @param res a susie fit, the output of `susieR::susie()`, or simply the posterior
#' inclusion probability matrix alpha.
#' @param include_index index of single effect models to consider when calculating PIP. Default to NULL.
#' @return a vector of posterior inclusion probability.
#' @export
susie_get_PIP = function(res, include_index = NULL) {
  if (class(res) == "susie")
    res = res$alpha
  if (!is.null(include_index)) {
    res = res[include_index,,drop=FALSE]
  }
  return(as.vector(1 - apply(1 - res, 2, prod)))
}

# compute standard error for regression coef
calc_stderr = function(X, residuals) {
    # S = (X'X)^-1 \Sigma
    sqrt(diag(sum(residuals^2) / (nrow(X) - 2) * chol2inv(chol(t(X) %*% X))))
  }

# univariate regression between each column of X and y
# Remove covariates if Z is not NULL
univariate_regression = function(X, y, Z=NULL, centered=FALSE, return_residuals=FALSE) {
  if (!centered) {
    y = y - mean(y)
    X = safe_colScale(X, center=TRUE, scale = FALSE)
  }
  if (!is.null(Z)) {
    if (!centered) Z = safe_colScale(Z, center=TRUE, scale=FALSE)
    y = .lm.fit(Z, y)$residuals
  }
  output = try(do.call(rbind,
                       lapply(1:ncol(X), function(i) {
                         g = .lm.fit(cbind(1, X[,i]), y)
                         return(c(coef(g)[2], calc_stderr(cbind(1, X[,i]), g$residuals)[2]))
                       })),
               silent = TRUE)
  if (class(output) == 'try-error') {
    output = matrix(0,ncol(X),2)
    for (i in 1:ncol(X)) {
      fit = summary(lm(y ~ X[,i]))$coef
      if (nrow(fit) == 2) {
        output[i,] = as.vector(summary(lm(y ~ X[,i]))$coef[2,1:2])
      } else {
        output[i,] = c(0,0)
      }
    }
  }
  if (return_residuals) {
    return(list(betahat = output[,1], sebetahat = output[,2],
                residuals = y))
  } else {
    return(list(betahat = output[,1], sebetahat = output[,2]))
  }
}

# computes z score (t-statistic) for association between each
# column of X and y
calc_z = function(X,y){
  out = univariate_regression(X,y)
  return(out$betahat/out$sebetahat)
}

#' @title Plot per variable summary in SuSiE CSs
#' @param data can be raw data \code{X} and \code{y} as \code{list(X=X,y=y)}, or be a vector of z-scores, or p-values.
#' @param res a susie fit, the output of `susieR::susie()`, or simply the posterior
#' inclusion probability matrix alpha.
#' @param dtype a string indicating the input data type (choices are raw_data (for raw data input), p (for p-value input), z (for z-scores input) and PIP (for PIP input))
#' @param CS the output of `susieR::susie_get_CS()`
#' @param coverage coverage of confident sets. Default to 0.9 for 90\% confidence interval.
#' @param add_bar add horizontal bar to signals in confidence interval.
#' @param b for simulated data, specify b = true effects (highlights in red).
#' @param max_cs the biggest CS to display, based on purity (set max_cs in between 0 and 1) or size (>1).
#' @export
susie_pplot = function(data,res=NULL,dtype='raw_data',CS=NULL,coverage=0.9,add_bar=FALSE,pos=NULL,b=NULL,max_cs=400,...){
  if (dtype=='raw_data') {
    z = calc_z(data$X,data$y)
    zneg = -abs(z)
    p = -log10(pnorm(zneg))
  } else if (dtype=='z') {
    zneg = -abs(data)
    p = -log10(pnorm(zneg))
  } else if (dtype=='p') {
    p = -log10(data)
  } else if (dtype=='PIP') {
    p = data
  } else {
    stop(paste0("Unknown dtype specification: ", dtype))
  }
  if(is.null(b)){b = rep(0,length(p))}
  if(is.null(pos)){pos = 1:length(p)}
  plot(pos,p,col="black",xlab="",ylab=ifelse(dtype=="PIP", "PIP", "-log10(p)"), pch=16, ...)
  if (!is.null(res)) {
    if (class(res) == "susie")
      res = res$alpha
    for(i in rev(1:nrow(res))){
      if (!is.null(CS$cs_index) && ! (i %in% CS$cs_index)) {
        next
      }
      if (!is.null(CS$purity) && max_cs < 1 && CS$purity[which(CS$cs_index==i),1] >= max_cs) {
        x0 = pos[CS$cs[[which(CS$cs_index==i)]]]
        y1 = p[CS$cs[[which(CS$cs_index==i)]]]
      } else if (n_in_CS(res, coverage)[i]<max_cs) {
        x0 = pos[which(in_CS(res, coverage)[i,]>0)]
        y1 = p[which(in_CS(res, coverage)[i,]>0)]
      } else {
        x0 = NULL
        y1 = NULL
      }
      if (is.null(x0)) {
        next
      }
      if (add_bar) {
        y0 = rep(0, length(x0))
        x1 = x0
        segments(x0,y0,x1,y1,lwd=1.5,col='gray')
      }
      points(x0, y1,col=i+2,cex=1.5,lwd=2.5)
    }
  }
  points(pos[b!=0],p[b!=0],col=2,pch=16)
}

#' @title Diagnostic plot for SuSiE iterations
#' @param res a susie fit, the output of `susieR::susie()`.
#' Multiple plots will be made for all iterations if `track_fit` was set to `TRUE` when running SuSiE.
#' @param L an integer, number of CS to plot
#' @param file_prefix prefix to path of output plot file
#' @param pos position of variables to display, default to all variables
#' @export
susie_iterplot = function(res, L, file_prefix, pos=NULL) {
  if(!requireNamespace("ggplot2",quietly = TRUE))
    stop("Required package ggplot2 not found")
  if(!requireNamespace("reshape",quietly = TRUE))
    stop("Required package reshape not found")
  get_layer = function(obj, k, idx, vars) {
    require(ggplot2,quietly = TRUE)
    alpha = reshape::melt(obj$alpha[1:k,vars,drop=F])
    colnames(alpha) = c('L', 'variables', 'alpha')
    alpha$L = as.factor(alpha$L)
    ggplot(alpha, aes(variables, alpha, group=L)) +
      geom_col(aes(fill=L)) +
      ggtitle(paste('Iteration', idx)) +
      theme_classic()
  }
  k = min(nrow(res$alpha), L)
  if (is.null(pos)) vars = 1:ncol(res$alpha)
  else vars = pos
  pdf(paste0(file_prefix, '.pdf'), 8, 3)
  if (is.null(res$trace)) {
    print(get_layer(res, k, res$niter, vars))
  } else {
    for (i in 2:length(res$trace)) {
      print(get_layer(res$trace[[i]], k, i-1, vars))
    }
  }
  dev.off()
  format = '.pdf'
  if (!is.null(res$trace)) {
    cmd = paste("convert -delay 30 -loop 0 -density 300 -dispose previous",
                paste0(file_prefix, '.pdf'),
                "\\( -clone 0 -set delay 300 \\) -swap 0 +delete \\( +clone -set delay 300 \\) +swap +delete -coalesce -layers optimize",
                paste0(file_prefix, '.gif'))
    cat("Creating GIF animation ...\n")
    if (file.exists(paste0(file_prefix, '.gif')))
      file.remove(paste0(file_prefix, '.gif'))
    output = try(system(cmd))
    if (class(output) == 'try-error') {
      cat("Cannot create GIF animation because `convert` command failed.\n")
    } else {
      format = '.gif'
    }
  }
  cat(paste0('Iterplot saved to ', file_prefix, format, '\n'))
}

# return residuals from Y after removing susie fit
get_R = function(X,Y,s){
  Y- X %*% coef(s)
}

#' @title Get updated prior variance
#' @param res a susie fit, the output of `susieR::susie()`.
#' @export
susie_get_prior_variance <- function(res) {
  return(res$sa2)
}

#' @title Get updated residual variance
#' @param res a susie fit, the output of `susieR::susie()`.
#' @export
susie_get_residual_variance <- function(res) {
  return(res$sigma2)
}
