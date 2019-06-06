#' @title Plot changepoint data and susie fit using ggplot2
#' @details Plots original data, y, overlaid with line showing susie fitted value and shaded rectangles showing credible sets for changepoint locations
#' @param y an n vector of observations that are ordered in time or space (assumed equally-spaced)
#' @param s a susie fit obtained by applying susie_trendfilter(y,order=0) to y
#' @param line_col color for the line showing fitted values
#' @param line_size a size for the line
#' @param cs_col color for the shaded rectangles showing credible sets
#' @return a ggplot2 object for plotting
#' @examples
#' set.seed(1)
#' mu = c(rep(0,50),rep(1,50),rep(3,50),rep(-2,50),rep(0,300))
#' y = mu + rnorm(500)
#' s = susie_trendfilter(y)
#' susie_plot_changepoint(s,y) # produces ggplot with credible sets for changepoints on top of plot
#'
#' @export
#' 
susie_plot_changepoint <-
  function(s,y, line_col="blue", line_size=1.5, cs_col="red"){
  df = data.frame(x = 1:length(y),y = y, mu = predict.susie(s))
  CS = susie_get_cs(s)$cs

  p= ggplot2::ggplot(df) +
    ggplot2::geom_point(data = df, ggplot2::aes_string(x="x", y="y")) +
    ggplot2::geom_line(color=line_col,data = df,
                       ggplot2::aes_string(x = "x",y = "mu"), size=line_size)
  for(i in 1:length(CS)){
    p = p + ggplot2::annotate("rect", fill = cs_col, alpha = 0.5,
                     xmin = min(CS[[i]])-0.5, xmax = max(CS[[i]])+0.5,
                     ymin = -Inf, ymax = Inf)
  }
  p
}
