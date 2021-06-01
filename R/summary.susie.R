#' @title Summarize Susie Fit.
#' 
#' @description \code{summary} method for the \dQuote{susie} class.
#' 
#' @param object A susie fit.
#'
#' @param \dots Additional arguments passed to the generic \code{summary}
#'   or \code{print.summary} method.
#' 
#' @return \code{summary.susie} returns a list containing a data frame
#'   of variables and a data frame of credible sets.
#' 
#' @method summary susie
#' 
#' @export summary.susie
#' 
#' @export
#' 
summary.susie = function (object, ...) {
  if (is.null(object$sets))
    stop("Cannot summarize SuSiE object because credible set information ",
         "is not available")
  variables = data.frame(cbind(1:length(object$pip),object$pip,-1))
  colnames(variables) = c("variable","variable_prob","cs")
  rownames(variables) = NULL
  if (object$null_index > 0)
    variables = variables[-object$null_index,]
  if (!is.null(object$sets$cs)) {
    cs = data.frame(matrix(NA,length(object$sets$cs),5))
    colnames(cs) = c("cs","cs_log10bf","cs_avg_r2","cs_min_r2","variable")
      for (i in 1:length(object$sets$cs)) {
        variables$cs[variables$variable %in% object$sets$cs[[i]]] =
          object$sets$cs_index[[i]]
        cs$cs[i] = object$sets$cs_index[[i]]
        cs$cs_log10bf[i] = log10(exp(object$lbf[cs$cs[i]]))
        cs$cs_avg_r2[i] = object$sets$purity$mean.abs.corr[i]^2
        cs$cs_min_r2[i] = object$sets$purity$min.abs.corr[i]^2
        cs$variable[i] = paste(object$sets$cs[[i]],collapse=",")
      }
      variables = variables[order(variables$variable_prob,decreasing = TRUE),]
  } else
    cs = NULL
  out = list(vars = variables,cs = cs)
  class(out) = c("summary.susie","list")
  return(out)
}

#' @rdname summary.susie
#' 
#' @param x A susie summary.
#'
#' @method print summary.susie
#' 
#' @export print.summary.susie
#' 
#' @export
#' 
print.summary.susie = function (x, ...) {
  cat("\nVariables in credible sets:\n\n")
  print.data.frame(x$vars[which(x$vars$cs > 0),],row.names = FALSE)
  cat("\nCredible sets summary:\n\n")
  print.data.frame(x$cs,row.names = FALSE)
}
