coef.susie_slide <- function(object,...) {
  scale <- object$X_column_scale_factors
  ans <- cbind(additive=c(object$intercept,colSums(object$alpha*object$mu)/scale),
               heterozygote=c(0,colSums(object$alpha*object$mu*object$delta)/scale))
  snps <- colnames(object$alpha)
  if(is.null(snps)) snps <- paste0("SNP",seq_len(ncol(object$alpha)))
  rownames(ans) <- c("(Intercept)",snps)
  if(object$null_index>0) ans <- ans[-(object$null_index+1L),,drop=FALSE]
  ans
}
predict.susie_slide <- function(object,newx=NULL,type=c("response","coefficients"),...) {
  type <- match.arg(type)
  if(type=="coefficients") {
    if(!is.null(newx)) stop("Do not supply newx when requesting coefficients.")
    return(coef(object))
  }
  if(is.null(newx)) return(object$fitted)
  if((!is.matrix(newx) && !inherits(newx,"sparseMatrix")) || ncol(newx)!=object$input_p)
    stop("newx must have the same SNP columns, in the same order, as the training genotypes.")
  values <- if(inherits(newx,"sparseMatrix")) newx@x else newx
  if(any(!is.finite(values)) || any(!values %in% 0:2)) stop("newx must contain finite hard-call genotypes.")
  cf <- coef(object)
  if(!is.null(colnames(newx)) && !is.null(colnames(object$alpha)) &&
     !identical(colnames(newx),rownames(cf)[-1])) stop("newx SNP names/order differ from the training data.")
  result <- as.numeric(cf[1,1]+newx %*% cf[-1,1])
  active <- which(cf[-1,2]!=0)
  if(length(active)) {
    cap <- max(1,floor(64e6/(8*max(1,nrow(newx)))))
    for(k in .chunks(length(active),cap)) {
      j <- active[k]
      result <- result+as.numeric(((newx[,j,drop=FALSE]==1)*1) %*% cf[j+1,2])
    }
  }
  result
}
slider_cs_table <- function(fit) {
  empty <- data.frame(cs=character(),component=integer(),snp_index=integer(),
                      snp=character(),alpha=numeric(),pip=numeric(),delta=numeric(),
                      delta_forced=logical())
  if(is.null(fit$sets$cs) || !length(fit$sets$cs)) return(empty)
  snps <- colnames(fit$alpha)
  if(is.null(snps)) snps <- paste0("SNP",seq_len(ncol(fit$alpha)))
  rows <- lapply(seq_along(fit$sets$cs),function(k) {
    l <- fit$sets$cs_index[k]; j <- fit$sets$cs[[k]]
    data.frame(cs=names(fit$sets$cs)[k],component=l,snp_index=j,snp=snps[j],
      alpha=fit$alpha[l,j],pip=fit$pip[j],delta=fit$delta[l,j],
      delta_forced=unname(fit$delta_forced[j]),row.names=NULL)
  })
  do.call(rbind,rows)
}
summary.susie_slide <- function(object,...) {
  # Use the familiar summary layout; avoid the upstream null-column row drop.
  temporary <- object
  temporary$null_index <- 0
  result <- susieR::summary.susie(temporary,...)
  result$delta <- slider_cs_table(object)
  result
}
