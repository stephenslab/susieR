#' @rdname susie_plots
#' 
#' @title SuSiE Plots.
#'
#' @details Plot per variable summary in SuSiE CSs.
#' 
#' @param model a susie fit, the output of `susieR::susie()`.
#' It has to contain `z`, `PIP` and optionally `sets`.
#' It is also possible to take in a vector of z-score or PIP,
#' in order to plot data from other software program.
#' @param y a string indicating what to plot: z (for z-score), PIP, log10PIP
#' or a random label to plot input data as is.
#' @param add_bar add horizontal bar to signals in credible interval.
#' @param pos can be either 1) numeric vector of indices of variables to plot, default to all variables, or
#' 2) a list of \code{list(attr = , start = , end = )} where \code{attr} is a character string of the name of 
#' index variable in \code{model} object, \code{start} and \code{end} are boundaries of indices to plot.
#' @param b for simulated data, specify b = true effects (highlights in red).
#' @param max_cs the biggest CS to display, based on purity (set max_cs in between 0 and 1) or size (>1).
#' @param add_legend if TRUE, add a legend to annotate the size and purity of each CS discovered.
#'
#' @importFrom utils head
#' @importFrom stats pnorm
#' @importFrom graphics plot
#' @importFrom graphics segments
#' @importFrom graphics points
#' @importFrom graphics legend
#' @importFrom graphics par
#'
#' @export
susie_plot = function(model,y,add_bar=FALSE,pos=NULL,b=NULL,max_cs=400,add_legend=FALSE,...){
  is_susie = inherits(model,"susie")
  ylab = y
  color = c(
    "dodgerblue2",
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )
  if (y=="z") {
    if (is_susie) {
      if (is.null(model$z))
        stop("z-score not available from SuSiE fit. Please set `compute_univariate_zscore=TRUE` in `susie()` function call.")
      zneg = -abs(model$z)
    }
    else zneg = -abs(model)
    p = -log10(pnorm(zneg))
    ylab = "-log10(p)"
  } else if (y=="PIP") {
    if (is_susie) p = model$pip
    else p = model
  } else if (y=="log10PIP") {
    if (is_susie) p = log10(model$pip)
    else p = log10(model)
    ylab = "log10(PIP)"
  } else {
    if (is_susie) stop("Need to specify z or PIP or log10PIP for SuSiE fits")
    p = model
  }
  if(is.null(b)){b = rep(0,length(p))}
  if(is.null(pos)){
    pos = 1:length(p)
  }
  if (inherits(pos, 'list')) {
    # check input
    if (is.null(pos$attr) || is.null(pos$start) || is.null(pos$end)) {
      stop("pos argument should be a list of list(attr=,start=,end=)")
    }
    if (!(pos$attr %in% names(model))) stop(paste("Cannot find attribute", pos$attr, "in input model object"))
    if (pos$start>=pos$end) stop("Position start should be smaller than end")
    start = min(min(model[[pos$attr]]), pos$start)
    end = max(max(model[[pos$attr]]), pos$end)
    # add zeros to alpha and p
    alpha = matrix(0, nrow(model$alpha), end - start + 1)
    new_p = rep(min(p), end - start + 1)
    pos_with_value = model[[pos$attr]] - start + 1
    new_p[pos_with_value] = p
    alpha[,pos_with_value] = model$alpha
    p = new_p
    model$alpha = alpha
    # adjust model$cs
    if (!is.null(model$sets$cs)) {
      for (i in 1:length(model$sets$cs)) {
        model$sets$cs[[i]] = pos_with_value[model$sets$cs[[i]]]
      }
    }
    # change pos object to be indices 
    start_adj = -min(min(model[[pos$attr]]) - pos$start, 0)
    end_adj = max(max(model[[pos$attr]]) - pos$end, 0)
    pos = (1 + start_adj):(length(p) - end_adj)
  }
  legend_text = list(col = vector(), purity = vector(), size = vector())
  options(scipen=10)
  plot(pos,p[pos],ylab=ylab, pch=16, ...)
  if (is_susie && !is.null(model$sets$cs)) {
    for(i in rev(1:nrow(model$alpha))){
      if (!is.null(model$sets$cs_index) && !(i %in% model$sets$cs_index)) {
        next
      }
      purity = model$sets$purity[which(model$sets$cs_index==i),1]
      if (!is.null(model$sets$purity) && max_cs < 1 && purity >= max_cs) {
        x0 = intersect(pos, model$sets$cs[[which(model$sets$cs_index==i)]])
        y1 = p[x0]
      } else if (n_in_CS(model, model$sets$coverage)[i]<max_cs) {
        x0 = intersect(pos, which(in_CS(model, model$sets$coverage)[i,]>0))
        y1 = p[x0]
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
        segments(x0,y0,x1,y1,lwd=1.5,col="gray")
      }
      points(x0, y1,col=head(color, 1),cex=1.5,lwd=2.5)
      legend_text$col = append(legend_text$col, head(color, 1))
      # rotate color
      color = c(color[-1], color[1])
      legend_text$purity = append(round(purity,4), legend_text$purity)
      legend_text$size = append(length(x0), legend_text$size)
    }
    if (length(legend_text$col) > 0 && add_legend) {
      # plot legend
      text = vector()
      for (i in 1:length(legend_text$col)) {
        if (legend_text$size[i] == 1) text[i] = paste0("L", i, ": C=1")
        else text[i] = paste0("L", i, ": C=", legend_text$size[i], "/R=", legend_text$purity[i])
      }
      legend(par("xaxp")[1], 1.1 * par("yaxp")[2], text,
        xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n", pch = 15, col = legend_text$col, cex = 0.75)
    }
  }
  points(pos[b!=0],p[b!=0],col=2,pch=16)
}

#' @rdname susie_plots
#' 
#' @details Diagnostic plot for SuSiE iterations
#' @param model a susie fit, the output of `susieR::susie()`.
#' Multiple plots will be made for all iterations if `track_fit` was set to `TRUE` when running SuSiE.
#' @param L an integer, number of CS to plot
#' @param file_prefix prefix to path of output plot file
#' @param pos index of variables to plot, default to all variables
#' 
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom reshape melt
#' @importFrom ggplot2 ggplot
#' 
#' @export
#' 
susie_plot_iteration = function (model, L, file_prefix, pos=NULL) {
  if(!requireNamespace("ggplot2",quietly = TRUE))
    stop("Required package ggplot2 not found")
  if(!requireNamespace("reshape",quietly = TRUE))
    stop("Required package reshape not found")
  get_layer = function(obj, k, idx, vars) {
    alpha = melt(obj$alpha[1:k,vars,drop=FALSE])
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
  format = ".pdf"
  if (!is.null(model$trace)) {
    cmd = paste("convert -delay 30 -loop 0 -density 300 -dispose previous",
                paste0(file_prefix,".pdf"),
                "\\( -clone 0 -set delay 300 \\) -swap 0 +delete \\( +clone -set delay 300 \\) +swap +delete -coalesce -layers optimize",
                paste0(file_prefix, ".gif"))
    cat("Creating GIF animation ...\n")
    if (file.exists(paste0(file_prefix,".gif")))
      file.remove(paste0(file_prefix,".gif"))
    output = try(system(cmd))
    if (inherits(output,"try-error")) {
      cat("Cannot create GIF animation because `convert` command failed.\n")
    } else {
      format = ".gif"
    }
  }
  cat(paste0("Iterplot saved to ", file_prefix, format, "\n"))
}
