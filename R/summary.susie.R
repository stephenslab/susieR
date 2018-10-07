#' @title summarize susie fit
#' @param s a susie fit
#' @return a data.frame of variables and a data.frame of credible sets
#' @method summary susie
#' @export
summary.susie = function(s) {
  if (is.null(s$sets))
    stop("Cannot summarize SuSiE object because credible set information is not available")
  cs = data.frame(matrix(NA, length(s$sets$cs), 5))
  colnames(cs) = c('cs', 'cs_log10bf', 'cs_avg_r2', 'cs_min_r2', 'variable')
  variables = data.frame(cbind(1:length(s$pip), s$pip, -1))
  colnames(variables) = c('variable', 'variable_prob', 'cs')
  for (i in 1:length(s$sets$cs)) {
    variables$cs[variables$variable %in% s$sets$cs[[i]]] = s$sets$cs_index[[i]]
    cs$cs[i] = s$sets$cs_index[[i]]
    cs$cs_log10bf[i] = s$lbf[cs$cs[i]]
    cs$cs_avg_r2[i] = s$sets$purity$mean.abs.corr[cs$cs[i]]^2
    cs$cs_min_r2[i] = s$sets$purity$min.abs.corr[cs$cs[i]]^2
    cs$variable[i] = paste(s$sets$cs[[i]], collapse=',')
  }
  variables = variables[order(variables$variable_prob, decreasing = T),]
  return(list(vars=variables, cs=cs))
}
