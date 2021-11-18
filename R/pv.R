##' Asymptotic P-values
##' 
##' Computes asymptotic P-values given the numerator of the pseudo-F statistic, 
##' its degrees of freedom and the eigenvalues of the residual covariance matrix.
##' 
##' @param ss numerator of the pseudo-F statistic.
##' @param lambda eigenvalues of the residual covariance matrix.
##' @param df degrees of freedom of the numerator of the pseudo-F statistic.
##' @param eps the desired level of accuracy.
##' @param eps.updt factor by which \code{eps} is updated to retry execution of
##' algorithm AS 204 when it fails with fault indicator 4, 5 or 9.
##' @param eps.stop if \code{eps > eps.stop}, execution of algorithm AS 204 is 
##' not retried and the function raises an error. Default is \code{1e-10}.
##' 
##' @return A vector containing the P-value and the level of accuracy. 
##' 
##' @seealso \code{\link{AS204}}
##' 
##' @author Diego Garrido-Martín
##' 
##' @keywords internal
##' 
p.asympt <- function(ss, df, lambda, eps = 1e-14, eps.updt = 2, eps.stop = 1e-10){
  
  pv  <- AS204(c = ss, lambda = lambda, mult = rep(df, length(lambda)), eps = eps)
  while (is.null(pv)) {
    eps <- eps * eps.updt
    if(eps > eps.stop){
      stop(sprintf("Precision of asymptotic P-value > %.2e", eps.stop))
    }
    pv  <- AS204(c = ss, lambda = lambda, mult = rep(df, length(lambda)), eps = eps)
  }
  if(pv < eps){
    pv <- eps
  }
  return(c(pv, eps))
}

##' Algorithm AS 204
##' 
##' Distribution of a positive linear combination of \eqn{\chi^2} random variables.
##' 
##' Algorithm AS 204 evaluates the expression
##' \deqn{
##' P [X < c] = P [ \sum_{j=1}^n \lambda_j \chi^2(m_j, \delta^2_j) < c ]
##' }
##' where \eqn{\lambda_j} and \eqn{c} are positive constants and 
##' \eqn{\chi^2(m_j, \delta^2_j)} represents an independent \eqn{\chi^2} 
##' random variable with \eqn{m_j} degrees of freedom and non-centrality
##' parameter \eqn{\delta^2_j}. This can be approximated by the truncated series 
##' \deqn{
##' \sum_{k=0}^{K-1} a_k P [\chi^2(m+2k) < c/\beta]
##' } 
##' where \eqn{m = \sum_{j=1}^n m_j} and \eqn{\beta} is an arbitrary constant 
##' (as given by argument "mode"). 
##' 
##' The \code{C++} implementation of algorithm AS 204 used here is identical 
##' to the one employed by the \code{\link[CompQuadForm]{farebrother}} method
##' in the \code{CompQuadForm} package, with minor modifications.
##' 
##' @param c value point at which distribution is to be evaluated.
##' @param lambda the weights \eqn{\lambda_j}.
##' @param mult the multiplicities \eqn{m_j}.
##' @param delta the non-centrality parameters \eqn{\delta^2_j}.
##' @param maxit the maximum number of terms \eqn{K} (see Details).
##' @param eps the desired level of accuracy.
##' @param mode if "\code{mode}" > 0 then \eqn{\beta=mode\lambda_{min}},
##' otherwise \eqn{\beta=2/(1/\lambda_{min}+1/\lambda_{max})}.
##' 
##' @return The function returns the probability \eqn{P[X > c] = 1 - P[X < c]} 
##' if the AS 204 fault indicator is 0 (see Note below), and \code{NULL} if 
##' the fault indicator is 4, 5 or 9, as the corresponding faults can be
##' corrected by increasing "\code{eps}". Other faults raise an error. 
##' 
##' @note The algorithm AS 204 defines the following fault indicators:
##' \strong{-j)} one or more of the constraints \eqn{\lambda_j > 0}, 
##' \eqn{m_j > 0} and \eqn{\delta^2_j \ge 0} is not satisfied.
##' \strong{1)} non-fatal underflow of \eqn{a_0}.
##' \strong{2)} one or more of the constraints \eqn{n > 0}, 
##' \eqn{c > 0}, \eqn{maxit > 0} and \eqn{eps > 0} is not satisfied.
##' \strong{3)} the current estimate of the probability is < -1.
##' \strong{4)} the required accuracy could not be obtained in \eqn{maxit} iterations.
##' \strong{5)} the value returned by the procedure does not satisfy 
##' \eqn{0 \le P [X < c] \le 1}.
##' \strong{6)} the density of the linear form is negative.
##' \strong{9)} faults 4 and 5.
##' \strong{10)} faults 4 and 6. 
##' \strong{0)} otherwise.
##' 
##' @references 
##' P. Duchesne, P. Lafaye de Micheaux, Computing the distribution of quadratic forms: 
##' Further comparisons between the Liu-Tang-Zhang approximation and exact methods, 
##' Computational Statistics and Data Analysis, Vol. 54, (2010), 858-862 
##' 
##' Farebrother R.W., Algorithm AS 204: The distribution of a Positive Linear Combination 
##' of chi-squared random variables, Journal of the Royal Statistical Society, 
##' Series C (applied Statistics), Vol. 33, No. 3 (1984), 332-339
##' 
##' @seealso \link[CompQuadForm]{farebrother}
##' 
##' @author Diego Garrido-Martín
##' 
##' @useDynLib manta ruben
##' 
##' @keywords internal
##' 
AS204 <- function (c, lambda, mult = rep(1, length(lambda)), delta = rep(0, length(lambda)),
                   maxit = 100000, eps = 1e-14, mode = 1) {
  
  out <- .C("ruben", lambda = as.double(lambda), mult = as.integer(mult), 
            delta = as.double(delta), n = as.integer(length(lambda)), 
            c = as.double(c), mode = as.double(mode), maxit = as.integer(maxit), 
            eps = as.double(eps), dnsty = as.double(0), ifault = as.integer(0), 
            res = as.double(0), PACKAGE = "manta")
  
  if(out$ifault == 0){
    return(1 - out$res)
  } else if (out$ifault %in% c(4, 5, 9)){
    return(NULL)
  } else {
    stop(sprintf("Algorithm AS 204 failed with fault indicator: %s", out$ifault))
  }
}
