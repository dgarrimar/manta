

##' Compute asymptotic P-values
##' 
##' @param ss numerator of the pseudo-F statistic.
##' @param lambda eigenvalues.
##' @param df degrees of freedom of the numerator of the test statistic.
##' @param eps precision limit.
##' @param eps.updt factor by which precision limit is updated.
##' 
##' @author Diego Garrido-Martín
##' 
##' @keywords internal
##' 
pv.ss <- function(ss, lambda, df, eps = 1e-14, eps.updt = 2){
  
  pv  <- AS204(c = ss, lambda = lambda, mult = rep(df, length(lambda)), eps = eps)
  while (is.null(pv)) {
    eps <- eps * eps.updt
    if(eps > 1e-10){
      stop(sprintf("Precision > 1e-10"))
    }
    pv  <- AS204(c = ss, lambda = lambda, mult = rep(df, length(lambda)), eps = eps)
  }
  if(pv < eps){
    pv <- eps
  }
  return(c(pv, eps))
}

##' AS204 algorithm
##' 
##' @author Diego Garrido-Martín
##' 
##' @useDynLib mlm ruben
##' 
##' @keywords internal
##'
AS204 <- function (c, lambda, mult = rep(1, length(lambda)), delta = rep(0, length(lambda)),
                   maxit = 100000, eps = 1e-10, mode = 1) {
  
  out <- .C("ruben", lambda = as.double(lambda), mult = as.integer(mult), 
            delta = as.double(delta), n = as.integer(length(lambda)), 
            c = as.double(c), mode = as.double(mode), maxit = as.integer(maxit), 
            eps = as.double(eps), dnsty = as.double(0), ifault = as.integer(0), 
            res = as.double(0), PACKAGE = "mlm")
  
  if(out$ifault == 0){
    return(1 - out$res)
  } else if (out$ifault %in% c(4, 5, 9)){
    return(NULL)
  } else {
    stop(sprintf("Fault indicator: %s", out$ifault))
  }
}