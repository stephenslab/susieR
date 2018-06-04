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

#find how many variables in the 95% CS
# x is a probability vector
n_in_CS_x = function(x){
  sum(cumsum(sort(x,decreasing = TRUE))<0.95)+1
}

# return binary vector indicating if each point is in CS
# x is a probability vector
in_CS_x = function(x){
  n = n_in_CS_x(x)
  o = order(x,decreasing=TRUE)
  result = rep(0,length(x))
  result[o[1:n]] = 1
  return(result)
}

#' @title Variables in susie confidence sets
#' @details It returns a binary matrix indicating which variables are in CS
#' of each effect
#' @param fitted a susie fit, the output of `susieR::susie()`, or simply the posterior
#' inclusion probability matrix alpha
#' @return a l by p binary matrix
#' @export
susie_in_CS = function(res){
  if (class(res) == "susie")
    res = res$alpha
  t(apply(res,1,in_CS_x))
}

n_in_CS = function(res){
  if (class(res) == "susie")
    res = res$alpha
  apply(res,1,n_in_CS_x)
}

# computes z score for association between each
# column of X and y
calc_z = function(X,y){
  z = rep(0,ncol(X))
  for(i in 1:ncol(X)){
    z[i] = summary(lm(y ~ X[,i]))$coeff[2,3]
  }
  return(z)
}


# plot p values of data and color in the 95% CSs
# for simulated data, specify b = true effects (highlights in red)
pplot = function(X,y,res,pos=NULL,b=NULL,CSmax = 400,...){
  z = calc_z(X,y)
  zneg = -abs(z)
  logp = log10(pnorm(zneg))
  if(is.null(b)){b = rep(0,ncol(X))}
  if(is.null(pos)){pos = 1:ncol(X)}
  plot(pos,-logp,col="grey",xlab="",ylab="-log10(p)",...)
  points(pos[b!=0],-logp[b!=0],col=2,pch=16)
  for(i in 1:nrow(res$alpha)){
    if(n_in_CS(res)[i]<CSmax)
      points(pos[which(in_CS(res)[i,]>0)],-logp[which(in_CS(res)[i,]>0)],col=i+2)
  }

}


# return residuals from Y after removing susie fit
get_R = function(X,Y,s){
  Y- X %*%  coef(s)
}

# get number of iterations
susie_get_niter <- function(res) {
  return(res$niter)
}
