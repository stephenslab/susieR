#' @title Local false sign rate (lfsr) for susie credible sets
#' @details This computes the average lfsr across SNPs for each l, weighted by the
#' posterior inclusion probability alpha
#' @param res a susie fit, the output of `susieR::susie()`
#' @return an l vector of lfsr for credible sets
#' @importFrom stats pnorm
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
# indicating which variables are in susie credible sets
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

# subsample and compute min, mean, median and max abs corr
#
#' @importFrom stats median
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
      value = abs(muffled_corr(as.matrix(X_sub)))
    } else {
      value = abs(Xcorr[pos, pos])
    }
    c(min(value,na.rm=T), mean(value,na.rm=T), median(value,na.rm=T))
  }
}

#' @title `cor` function with specified warning muffled
#' @keywords internal
#' @importFrom stats cor
muffled_corr = function(x)
  withCallingHandlers(cor(x),
                    warning = function(w) {
                      if (grepl("the standard deviation is zero", w$message))
                        invokeRestart("muffleWarning")
                    } )

#' @title `cov2cor` function with specified warning muffled
#' @keywords internal
#' @importFrom stats cov2cor
muffled_cov2cor = function(x)
  withCallingHandlers(cov2cor(x),
                      warning = function(w) {
                        if (grepl("had 0 or NA entries; non-finite result is doubtful", w$message))
                          invokeRestart("muffleWarning")
                      } )

#' @title check for symmetric matrix
#' @keywords internal

is_symmetric_matrix = function(x) {
    res = isSymmetric(x)
    if (res == FALSE) res = isSymmetric(unname(x))
    return(res)
}

#' @title Extract credible sets from SuSiE fit
#' @details Reports indices of variables in each credible set (CS) identified,
#' as well as summaries of correlation among the variables included in each CS.
#' If desired, one can filter out CSs that do not meet a specified purity threshold (min_abs_corr);
#' to do this either `X` or `Xcorr` must be supplied.
#' @param res a susie fit, the output of `susieR::susie()`, or simply the posterior
#' inclusion probability matrix alpha.
#' @param X N by P matrix of values of the P variables (covariates) in N samples.
#' When provided, correlation between variables will be computed and used to remove
#' CSs whose minimum correlation among variables is smaller than `min_abs_corr` (see below).
#' @param Xcorr P by P matrix of correlations between variables (covariates).
#' When provided, it will be used to remove CSs
#' whose minimum correlation among variables is smaller than `min_abs_corr` (see below).
#' @param coverage a number (in [0,1]) specifying desired coverage of each CS
#' @param min_abs_corr a "purity" threshold for the CS. Any CS that contains
#' a pair of variables with correlation less than this threshold will be filtered out and not reported.
#' @param dedup If TRUE, "deduplicate" - that is remove duplicated CSs.
#' @return a list with elements:
#'
#' \item{cs}{a list, each element corresponds to a CS, and is a vector containing the indices of the variables in the CS.}
#' \item{coverage}{the nominal coverage specified for each CS.}
#' \item{purity}{(If `X` or `Xcorr` are provided), the purity of each CS.}
#' \item{cs_index}{(If `X` or `Xcorr` are provided) the index (in 1,...,L) of each reported CS in the supplied susie fit.}
#'
#' @export
susie_get_cs = function(res,
                        X = NULL, Xcorr = NULL,
                        coverage = 0.95,
                        min_abs_corr = 0.5,
                        dedup = TRUE) {
  if (class(res) == "susie") {
    null_index = res$null_index
    if (is.numeric(res$V)) include_idx = which(res$V != 0)
    else include_idx = 1:nrow(res$alpha)
    if (length(include_idx) == 0) return(list(cs = NULL,coverage=coverage))
    res = res$alpha[include_idx, , drop=F]
  } else {
    null_index = 0
  }
  if (!is.null(X) && !is.null(Xcorr)) {
    stop("Only one of X or Xcorr should be specified")
  }
  if (!is.null(Xcorr) && !is_symmetric_matrix(Xcorr)) {
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
  # drop the single effect with estimated prior zero
  # compute and filter by "purity"
  if (is.null(Xcorr) && is.null(X)) {
    return(list(cs=cs,coverage=coverage))
  } else {
    purity = data.frame(do.call(rbind, lapply(1:length(cs), function(i)
            {
                if (null_index > 0 && null_index %in% cs[[i]]) c(-9,-9,-9)
                else get_purity(cs[[i]], X, Xcorr)
            })))
    colnames(purity) = c('min.abs.corr', 'mean.abs.corr', 'median.abs.corr')
    is_pure = which(purity$min.abs.corr >= min_abs_corr)
    if (length(is_pure) > 0) {
      cs = cs[is_pure]
      purity = purity[is_pure,]
      row_names = paste0("L", is_pure)
      names(cs) = row_names
      rownames(purity) = row_names
      ## re-order CS list and purity rows based on purity
      ordering = order(purity$min.abs.corr, decreasing=T)
      return(list(cs = cs[ordering], purity = purity[ordering,], cs_index = is_pure[ordering],coverage=coverage))
    } else {
      return(list(cs = NULL,coverage=coverage))
    }
  }
}

#' @title Compute posterior inclusion probability (PIP) for all variables
#' @param res a susie fit, the output of `susieR::susie()`, or simply the posterior
#' inclusion probability matrix alpha.
#' @param prune_by_cs Whether or not to ignore single effects not in reported CS when calculating PIP. Default to FALSE.
#' @return a vector of posterior inclusion probability.
#' @export
susie_get_pip = function(res, prune_by_cs = FALSE) {
  if (class(res) == "susie") {
    # drop null weight columns
    if (res$null_index > 0) res$alpha = res$alpha[,-res$null_index]
    # drop the single effect with estimated prior zero
    if (is.numeric(res$V)) include_idx = which(res$V != 0)
    else include_idx = 1:nrow(res$alpha)
    # only consider variables in reported CS
    # this is not what we do in the SuSiE paper
    # so by default prune_by_cs = FALSE means we do not run the following code
    if (!is.null(res$sets$cs_index) && prune_by_cs)
      include_idx = intersect(include_idx, res$sets$cs_index)
    if (is.null(res$sets$cs_index) && prune_by_cs)
      include_idx = numeric(0)
    # now extract relevant rows from alpha matrix
    if (length(include_idx) > 0) {
      res = res$alpha[include_idx,,drop=FALSE]
    } else {
      res = matrix(0,1,ncol(res$alpha))
    }
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
#
#' @importFrom stats lm
#' @importFrom stats .lm.fit
#' @importFrom stats coef
#' @importFrom stats summary.lm
univariate_regression = function(X, y, Z=NULL, centered=FALSE,
                                 return_residuals=FALSE) {
  if (!centered) {
    y = y - mean(y)
    X = set_X_attributes(X, center=TRUE, scale = FALSE)
  }
  if (!is.null(Z)) {
    if (!centered) Z = set_X_attributes(Z, center=TRUE, scale=FALSE)
    y = .lm.fit(Z, y)$residuals
  }
  output = try(do.call(rbind,
                       lapply(1:ncol(X), function(i) {
                         g = .lm.fit(cbind(1, X[,i]), y)
                         return(c(coef(g)[2], calc_stderr(cbind(1, X[,i]), g$residuals)[2]))
                       })),
               silent = TRUE)
  if (class(output) == 'try-error') {
    # Exception occurs, fall back to a safer but slower calculation
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
calc_z = function(X,y,centered=FALSE){
  out = univariate_regression(X,y,centered=centered)
  return(out$betahat/out$sebetahat)
}

#' @title Plot per variable summary in SuSiE CSs
#' @param model a susie fit, the output of `susieR::susie()`.
#' It has to contain `z`, `PIP` and optionally `sets`.
#' It is also possible to take in a vector of z-score or PIP,
#' in order to plot data from other software program.
#' @param y a string indicating what to plot: z (for z-score), PIP, log10PIP
#' or a random label to plot input data as is.
#' @param add_bar add horizontal bar to signals in credible interval.
#' @param pos coordinates of variables to plot, default to all variables
#' @param b for simulated data, specify b = true effects (highlights in red).
#' @param max_cs the biggest CS to display, based on purity (set max_cs in between 0 and 1) or size (>1).
#' @param add_legend if TRUE, add a legend to annotate the size and purity of each CS discovered.
#' @importFrom stats pnorm
#' @importFrom graphics plot
#' @importFrom graphics segments
#' @importFrom graphics points
#' @importFrom graphics legend
#' @importFrom graphics par
#'
#' @export
susie_plot = function(model,y,add_bar=FALSE,pos=NULL,b=NULL,max_cs=400,add_legend=FALSE,...){
  is_susie = (class(model) == "susie")
  ylab = y
  if (y=='z') {
    if (is_susie) {
      if (is.null(model$z))
        stop('z-score not available from SuSiE fit. Please set `compute_univariate_zscore=TRUE` in `susie()` function call.')
      zneg = -abs(model$z)
    }
    else zneg = -abs(model)
    p = -log10(pnorm(zneg))
    ylab = "-log10(p)"
  } else if (y=='PIP') {
    if (is_susie) p = model$pip
    else p = model
  } else if (y=='log10PIP') {
    if (is_susie) p = log10(model$pip)
    else p = log10(model)
    ylab = "log10(PIP)"
  } else {
    if (is_susie) stop('Need to specify z or PIP or log10PIP for SuSiE fits.')
    p = model
  }
  if(is.null(b)){b = rep(0,length(p))}
  if(is.null(pos)){pos = 1:length(p)}
  legend_text = list(col = vector(), purity = vector(), size = vector())
  plot(pos,p,col="black", ylab=ylab, pch=16, ...)
  if (is_susie && !is.null(model$sets$cs)) {
    for(i in rev(1:nrow(model$alpha))){
      if (!is.null(model$sets$cs_index) && !(i %in% model$sets$cs_index)) {
        next
      }
      purity = model$sets$purity[which(model$sets$cs_index==i),1]
      if (!is.null(model$sets$purity) && max_cs < 1 && purity >= max_cs) {
        x0 = pos[model$sets$cs[[which(model$sets$cs_index==i)]]]
        y1 = p[model$sets$cs[[which(model$sets$cs_index==i)]]]
      } else if (n_in_CS(model, model$sets$coverage)[i]<max_cs) {
        x0 = pos[which(in_CS(model, model$sets$coverage)[i,]>0)]
        y1 = p[which(in_CS(model, model$sets$coverage)[i,]>0)]
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
      legend_text$col = append(i+2, legend_text$col)
      legend_text$purity = append(round(purity,4), legend_text$purity)
      legend_text$size = append(length(x0), legend_text$size)
    }
    if (length(legend_text$col) > 0 && add_legend) {
      # plot legend
      text = vector()
      for (i in 1:length(legend_text$col)) {
        if (legend_text$size[i] == 1) text[i] = paste0("L", legend_text$col[i]-2, ": C=1")
        else text[i] = paste0("L", legend_text$col[i]-2, ": C=", legend_text$size[i], "/R=", legend_text$purity[i])
      }
      legend(par("xaxp")[1], 1.1 * par("yaxp")[2], text,
        xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n", pch = 15, col = legend_text$col, cex = 0.75)
    }
  }
  points(pos[b!=0],p[b!=0],col=2,pch=16)
}

#' @title Diagnostic plot for SuSiE iterations
#' @param model a susie fit, the output of `susieR::susie()`.
#' Multiple plots will be made for all iterations if `track_fit` was set to `TRUE` when running SuSiE.
#' @param L an integer, number of CS to plot
#' @param file_prefix prefix to path of output plot file
#' @param pos position of variables to display, default to all variables
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @export
susie_plot_iteration = function(model, L, file_prefix, pos=NULL) {
  if(!requireNamespace("ggplot2",quietly = TRUE))
    stop("Required package ggplot2 not found")
  if(!requireNamespace("reshape",quietly = TRUE))
    stop("Required package reshape not found")
  get_layer = function(obj, k, idx, vars) {
    alpha = reshape::melt(obj$alpha[1:k,vars,drop=F])
    colnames(alpha) = c("L","variables","alpha")
    alpha$L = as.factor(alpha$L)
    ggplot2::ggplot(alpha,ggplot2::aes_string("variables","alpha",group="L")) +
      ggplot2::geom_col(ggplot2::aes_string(fill = "L")) +
      ggplot2::ggtitle(paste("Iteration", idx)) +
      ggplot2::theme_classic()
  }
  k = min(nrow(model$alpha), L)
  if (is.null(pos))
    vars = 1:ncol(model$alpha)
  else
    vars = pos
  pdf(paste0(file_prefix,".pdf"), 8, 3)
  if (is.null(model$trace)) {
    print(get_layer(model, k, model$niter, vars))
  } else {
    for (i in 2:length(model$trace)) {
      print(get_layer(model$trace[[i]], k, i-1, vars))
    }
  }
  dev.off()
  format = '.pdf'
  if (!is.null(model$trace)) {
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
#
#' @importFrom stats coef
get_R = function(X,Y,s){
  Y- X %*% coef(s)
}

#' @title Get updated prior variance
#' @param res a susie fit, the output of `susieR::susie()`.
#' @export
susie_get_prior_variance <- function(res) {
  return(res$V)
}

#' @title Get updated residual variance
#' @param res a susie fit, the output of `susieR::susie()`.
#' @export
susie_get_residual_variance <- function(res) {
  return(res$sigma2)
}

#' @title Get evidence lower bound (ELBO) from fitted SuSiE model
#' @param res a susie fit, the output of `susieR::susie()`.
#' @param all whether or not to get ELBO from all iterations
#' @param warning_tol warn if ELBO is non-decreasing by this tolerance level (default set to 1E-6)
#' @export
susie_get_objective <- function(res, all = FALSE, warning_tol = 1E-6) {
  if (!all(diff(res$elbo) >= (-1 * warning_tol))) {
    warning('Objective is not non-decreasing')
  }
  if (all) return(res$elbo)
  else return(res$elbo[length(res$elbo)])
}

#' @title Slim the result of fitted SuSiE model
#' @param res a susie fit, the output of `susieR::susie()`
susie_slim = function(res){
  list(alpha = res$alpha, niter = res$niter, V = res$V, sigma2 = res$sigma2)
}

#' @title Get posterior mean for coefficients from fitted SuSiE model
#' @param res a susie fit
susie_get_posterior_mean = function(res){
  colSums(res$alpha*s$mu)/res$X_column_scale_factors
}

#' @title Get posterior standard deviation for coefficients from fitted SuSiE model
#' @param res a susie fit
susie_get_posterior_sd = function(res){
  sqrt(colSums(res$alpha * res$mu2 - (res$alpha*res$mu)^2))/(res$X_column_scale_factors)
}
