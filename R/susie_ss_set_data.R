#' @title Create data for susie_ss
#' @param bhat a p vector, the estimated effects from p simple linear regression: y = Xb + e
#' @param shat a p vector of corresponding standard errors
#' @param zhat a p vector of z scores. bhat/shat
#' @param XtX a p by p matrix, X'X, where columns of X are centered to have mean 0
#' @param R a p by p correlation matrix of X
#' @param var_y the (sample) variance of y
#' @param n sample size
#' @param standardize logical flag (default=TRUE) for whether to fit susie_ss with standardize columns of X
#' @param bhat.standardized whether the bhat comes from the standardized X
#'
#' @export
susie_ss_set_data = function(bhat = NULL, shat = NULL, zhat = NULL,
                             XtX = NULL, R = NULL,
                             var_y = 1, n,
                             standardize = TRUE, bhat.standardized = TRUE){
  if(is.null(bhat) && is.null(zhat)){
    stop('One of bhat or zhat must be specified.')
  }

  if(is.null(XtX) && is.null(R)){
    stop('One of XtX or R must be specified.')
  }

  if(!is.null(bhat) && !is.null(shat) && !is.null(zhat)){
    if(!all.equal(bhat/shat, zhat)){
      stop('The provided zhat is not bhat/shat.')
    }
  }

  if(!is.null(bhat) && !is.null(zhat)){
    message('Both bhat and zhat are specified. We use zhat in susie_ss.')
  }

  if(standardize){
    if(is.null(R)){
      R = cov2cor(XtX)
    }
    # z score
    if(!is.null(bhat) && !is.null(shat)){
      zhat = bhat/shat
    }
    if(!is.null(zhat)){
      R2 = zhat^2/(zhat^2 + n-2)
      sigma2 = var_y*(n-1)*(1-R2)/(n-2)
      Xty = sqrt(sigma2) * sqrt(n-1) * zhat
      return(list(XtX = (n-1)*R, Xty = Xty, var_y = var_y, n = n))
    }

    # bhat, no shat
    if(!is.null(bhat)){
      if(bhat.standardized){
        Xty = (n-1)* bhat
        return(list(XtX = (n-1)*R, Xty = Xty, var_y = var_y, n = n))
      }else{
        if(is.null(XtX)){
          stop('XtX must be specified....')
        }
        Xty = sqrt(n-1) * sqrt(diag(XtX)) * bhat
        return(list(XtX = (n-1)*R, Xty = Xty, var_y = var_y, n = n))
      }
    }

  }else{
    if(is.null(XtX)){
      if(!is.null(bhat) && !is.null(shat)){
        XtX = (n-1)*var_y/(shat^2 * (n-2) + bhat^2)
        Xty = bhat * XtX
        XtX = t(R * sqrt(XtX)) * sqrt(XtX)
        return(list(XtX = XtX, Xty = Xty, var_y = var_y, n = n))
      }else{
        stop('To have non standardized results, XtX must be specified.')
      }
    }

    # z score
    if(!is.null(bhat) && !is.null(shat)){
      zhat = bhat/shat
    }
    if(!is.null(zhat)){
      R2 = zhat^2/(zhat^2 + n-2)
      sigma2 = var_y*(n-1)*(1-R2)/(n-2)
      Xty = sigma2^(0.5) * sqrt(diag(XtX)) * zhat
      return(list(XtX = XtX, Xty = Xty, var_y = var_y, n = n))
    }

    # bhat, no shat
    if(!is.null(bhat)){
      if(bhat.standardized){
        Xty = sqrt((n-1)*diag(XtX)) * bhat
        return(list(XtX = XtX, Xty = Xty, var_y = var_y, n = n))
      }else{
        Xty = diag(XtX)* bhat
        return(list(XtX = XtX, Xty = Xty, var_y = var_y, n = n))
      }
    }
  }
}
