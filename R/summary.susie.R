#' @title summarize susie fit
#' @param object a susie fit
#' @return a data.frame of variables and a data.frame of credible sets
#' @export summary.susie
#' @export
summary.susie = function (object, ...) {
  if (is.null(object$sets))
    stop("Cannot summarize SuSiE object because credible set information is not available")
  variables = data.frame(cbind(1:length(object$pip), object$pip, -1))
  colnames(variables) = c('variable', 'variable_prob', 'cs')
  rownames(variables) = NULL
  if (object$null_index > 0) variables = variables[-object$null_index,]
  if (!is.null(object$sets$cs)) {
      cs = data.frame(matrix(NA, length(object$sets$cs), 5))
      colnames(cs) = c('cs', 'cs_log10bf', 'cs_avg_r2', 'cs_min_r2', 'variable')
      for (i in 1:length(object$sets$cs)) {
        variables$cs[variables$variable %in% object$sets$cs[[i]]] = object$sets$cs_index[[i]]
        cs$cs[i] = object$sets$cs_index[[i]]
        cs$cs_log10bf[i] = object$lbf[cs$cs[i]]
        cs$cs_avg_r2[i] = object$sets$purity$mean.abs.corr[i]^2
        cs$cs_min_r2[i] = object$sets$purity$min.abs.corr[i]^2
        cs$variable[i] = paste(object$sets$cs[[i]], collapse=',')
      }
      variables = variables[order(variables$variable_prob, decreasing = T),]
  } else {
      cs = NULL
  }
  return(list(vars=variables, cs=cs))
}
