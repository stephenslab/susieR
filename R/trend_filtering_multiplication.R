# @title Compute unscaled X %*% b using the special structure of trend
#   filtering
# @param order is the order of trend filtering
# @param b an n=p vector
# @return an n vector
compute_tf_Xb = function(order,b) {
  for (i in 1:(order+1))
    b = rev(-1*cumsum(rev(b)))
  return(b)
}

# @title Compute unscaled t(X) %*% y using the special structure of
#   trend filtering
# @param order is the order of trend filtering
# @param y an n vector
# @return an n vector
compute_tf_Xty = function(order, y) {
  for (i in 1:(order+1))
    y = -1*cumsum(y)
  return(y)
}

# @title Compute colSums(X*X) for X under four scenarios
# @param order is the order of trend filtering
# @param n the length of y
# @param cm column means of X
# @param csd column standard deviations of X
# @param intercept a boolean denotes whether mean centering X
# @param standardize a boolean denotes whether scaling X by standard deviation
# @return an n vector
compute_tf_d = function (order, n, cm, csd, standardize = FALSE,
                         intercept = FALSE) {
  if (intercept) {
      
    # When standardize = TRUE, intercept = TRUE: by special
    # observation d = [n-1, n-1, ...]
    d = rep(n-1,n)
    if (order == 0)
      d[n] = 0

    # When standardize = FALSE, intercept = TRUE:
    # d = [n-1, n-1, ...] * (csd^2)
    if (!standardize)
      d = d*csd^2
    return(d)
  } else {
      
    # When standardize = FALSE, intercept = FALSE: d = colSums(X^2)
    base = rep(-1,n)
    if (order == 0)
      d = cumsum(base^2)
    else {
      for (i in 1:order)
        base = cumsum(base)
      d = cumsum(base^2) 
    }

    # When standardize = TRUE, intercept = TRUE:
    # d = colSums(X^2) / (csd^2)
    if (standardize)
      d = d/csd^2
    return(d)
  }
}

# @title Compute column mean of the trend filtering matrix X.
# @param order is the order of trend filtering
# @param n the length of y
# @return an n vector
compute_tf_cm = function (order, n) {
  base = rep(1,n)
  for (i in 1:(order+1))
    base = -cumsum(base)
  return(base/n)
}

# @title Compute column standard deviation of the trend filtering
#   matrix X
# @param order is the order of trend filtering
# @param n is the length of y
# @return an n vector
compute_tf_csd = function (order, n) {
  cm = compute_tf_cm(order,n)
  csd = sqrt((compute_tf_d(order,n)/n - cm^2)*n/(n-1))
  csd[which(csd == 0)] = 1
  return(csd)
}

# @title A fast way to compute colSums(X*X), where X is a
#   mean-centered and standardized trend filtering matrix.
# @param order order of trend filtering
# @param n the length of y
# @return an n vector
compute_tf_std_d = function (order, n) {
  res = rep(n-1,n)
  if (order == 0)
    res[n] = 0
  return(res)
}
