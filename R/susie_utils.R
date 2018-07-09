#' @title Local false sign rate (lfsr) for susie confidence sets
#' @details This computes the average lfsr across SNPs for each l, weighted by the
#' posterior inclusion probability alpha
#' @param fitted a susie fit, the output of `susieR::susie()`
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
      value = abs(cor(X[,pos]))
    } else {
      value = abs(Xcorr[pos, pos])
    }
    c(min(value), mean(value), median(value))
  }
}

#' @title Extract confidence sets from fitted SuSiE model
#' @details It reports indices of variables in each confidence set identified,
#' as well as summaries of correlation between variables within each set.
#' @param fitted a susie fit, the output of `susieR::susie()`, or simply the posterior
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
susie_get_CS = function(fitted,
                        X = NULL, Xcorr = NULL,
                        coverage = 0.9,
                        min_abs_corr = 0.5) {
  if (class(fitted) == "susie")
    fitted = fitted$alpha
  if (!is.null(X) && !is.null(Xcorr)) {
    stop("Only one of X or Xcorr should be specified")
  }
  if (!is.null(X) && ncol(X) != length(fitted[1,])) {
    stop("Column of X matrix does not match length of fitted posterior")
  }
  if (!is.null(Xcorr) && ncol(Xcorr) != length(fitted[1,])) {
    stop("Dimension of Xcorr matrix does not match length of fitted posterior")
  }
  if (!is.null(Xcorr) && !isSymmetric(Xcorr)) {
    stop("Xcorr matrix must be symmetric")
  }
  # L by P binary matrix
  status = in_CS(fitted, coverage)
  # an L list of CS positions
  cs = lapply(1:nrow(status), function(i) which(status[i,]!=0))
  cs = cs[lapply(cs, length) > 0]
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
    } else {
      row_names = paste0("L", 1:length(cs))
    }
    names(cs) = row_names
    rownames(purity) = row_names
    return(list(cs = cs, purity = purity, cs_index = is_pure))
  }
}

#' @title Compute posterior inclusion probability (PIP) for all variables
#' @param fitted a susie fit, the output of `susieR::susie()`, or simply the posterior
#' inclusion probability matrix alpha.
#' @param include_index index of single effect models to consider when calculating PIP. Default to NULL.
#' @return a vector of posterior inclusion probability.
#' @export
susie_get_PIP = function(fitted, include_index = NULL) {
  if (class(fitted) == "susie")
    fitted = fitted$alpha
  if (!is.null(include_index)) {
    fitted = fitted[include_index,,drop=FALSE]
  }
  return(as.vector(1 - apply(1 - fitted, 2, prod)))
}

# compute standard error for regression coef
calc_stderr = function(X, residuals) {
    # S = (X'X)^-1 \Sigma
    sqrt(diag(sum(residuals^2) / (nrow(X) - 2) * chol2inv(chol(t(X) %*% X))))
  }

# univariate regression between each column of X and y
# Remove covariates if Z is not NULL
univariate_regression = function(X, y, Z=NULL, return_residue=FALSE) {
  if (!is.null(Z)) {
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
  if (return_residue) {
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
#' @param fitted a susie fit, the output of `susieR::susie()`, or simply the posterior
#' inclusion probability matrix alpha.
#' @param dtype a string indicating the input data type (choices are raw_data (for raw data input), p (for p-value input), z (for z-scores input) and PIP (for PIP input))
#' @param coverage coverage of confident sets. Default to 0.9 for 90\% confidence interval.
#' @param add_bar add horizontal bar to signals in confidence interval.
#' @param b for simulated data, specify b = true effects (highlights in red).
#' @param CSmax maximum size of CS to display.
#' @export
susie_pplot = function(data,fitted=NULL,dtype='raw_data',coverage=0.9,add_bar=FALSE,pos=NULL,b=NULL,CSmax=400,...){
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
  if (!is.null(fitted)) {
    if (class(fitted) == "susie")
      fitted = fitted$alpha
    for(i in 1:nrow(fitted)){
      if(n_in_CS(fitted, coverage)[i]<CSmax) {
        x0 = pos[which(in_CS(fitted, coverage)[i,]>0)]
        y1 = p[which(in_CS(fitted, coverage)[i,]>0)]
        if (add_bar) {
          y0 = rep(0, length(x0))
          x1 = x0
          segments(x0,y0,x1,y1,lwd=1.5,col='gray')
        }
        points(x0, y1,col=i+2,cex=1.5,lwd=2.5)
        }
    }
  }
  points(pos[b!=0],p[b!=0],col=2,pch=16)
}

# return residuals from Y after removing susie fit
get_R = function(X,Y,s){
  Y- X %*% coef(s)
}

# get number of iterations
susie_get_niter <- function(res) {
  return(res$niter)
}
