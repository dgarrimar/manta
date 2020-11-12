##' Multivariate linear models
##' 
##' Fits a multivariate linear model and computes test statistics and p-values
##' for a set of predictors in a non-parametric manner. 
##' 
##' A \code{Y} matrix is obtained after transforming and centering the original 
##' response variables. Then, the multivariate fit obtained by \code{\link{lm}} 
##' can be used to compute sums of squares (\code{I}, \code{II} or \code{III}), 
##' pseudo F statistics and asymptotic p-values for the explanatory 
##' variables in a non-parametric manner.
##' 
##' @param formula object of class "\code{\link{formula}}": a symbolic 
##' description of the model to be fitted. 
##' @param data an optional data frame containing the variables in the model. 
##' If not found in data, the variables are taken from \code{environment(formula)}, 
##' typically the environment from which \code{mlm} is called.
##' @param transform transformation of the response variables: "\code{none}", 
##' "\code{sqrt}" or "\code{log}". Default is "\code{none}".
##' @param type type of sum of squares: "\code{I}", "\code{II}" or "\code{III}". 
##' Default is "\code{II}".
##' @param contrasts an optional list. See \code{contrasts.arg} in 
##' \code{\link{model.matrix.default}}. Default is "\code{\link{contr.sum}}" 
##' for ordered factors and "\code{\link{contr.poly}}" for unordered factors. 
##' Note that this is different from the default setting in 
##' \code{\link{options}("contrasts")}.
##' @param subset subset of explanatory variables for which summary statistics will be reported.
##' @param fit logical. If \code{TRUE} the multivariate fit on transformed and centered responses
##' is returned.
##' 
##' @return \code{mlm} returns an object of \code{\link{class}} "MLM", a list containing:
##' \item{call}{the matched call.}
##' \item{aov.tab}{ANOVA table with Df, Sum Sq, Mean Sq, F values, partial R2 and P values.}
##' \item{type}{the type of sum of squares (\code{I}, \code{II} or \code{III}).}
##' \item{precision}{the precision in P value computation.}
##' \item{transform}{the transformation applied to the response variables.}
##' \item{na.omit}{incomplete cases removed (see \code{\link{na.omit}})}
##' \item{fit}{if \code{fit = TRUE} the multivariate fit done on the transformed and centered 
##' response variables is also returned.}
##' 
##' @seealso \code{\link{lm}}
##' 
##' @author Diego Garrido-Martín
##' 
##' @import stats
##' @export
##' 
mlm <- function(formula, data, transform = "none", type = "II", 
                contrasts = NULL, subset = NULL, fit = FALSE){
  
  ## Checks
  transform <- match.arg(transform, c("none", "sqrt", "log"))
  type <- match.arg(type, c("I", "II", "III"))
  
  ## Save call, build model frame, obtain responses
  cl <- match.call()
  m <- match(c("formula", "data"), names(cl), 0L)
  mf <- cl[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action = "na.omit"
    # The rows with at least one NA either in Y or X 
    # (only considering variables used in the formula) 
    # will be removed before transforming/centering
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())               
  mt <- attr(mf, "terms")
  response <- model.response(mf, "numeric")
  
  ## Checks
  if (NCOL(response) < 2){
    stop("The number of response variables should be >= 2.")
  }
  if (length(attr(mt, "term.labels")) < 1) {
    stop("The model should contain at least one predictor (excluding the intercept).")
  }
  
  ## Transform and center responses, update model frame
  if(transform == "none"){
    Y <- response
  } else if (transform == "sqrt"){
    Y <- sqrt(response)
  } else if (transform == "log"){
    Y <- log(response)
  } 
  Y <- scale(Y, center = TRUE, scale = FALSE)
  mf[[1L]] <- Y
  
  ## Define contrasts
  if(is.null(contrasts)){
    contrasts <- list(unordered = "contr.sum", ordered = "contr.poly")
    dc <- attr(mt, "dataClasses")[-1]
    contr.list <- lapply(dc, FUN = function(k){
      # No contrast for quantitaive predictors
      # Sum contrasts for unordered categorical predictors
      # Polynomial contrasts for ordered categorical predictors
      contr.type <- switch(k, "factor" = contrasts$unordered,
                           "ordered" = contrasts$ordered)
      return(contr.type)
    })
    contr.list <- contr.list[!unlist(lapply(contr.list, is.null))]
  } else {
    contr.list <- contrasts
  }

  ## Build model matrix
  X <- model.matrix(mt, mf, contr.list)
  
  ## Fit lm
  lmfit <- lm.fit(X, Y)
  class(lmfit) <- c("mlm", "lm")
  lmfit$na.action <- attr(mf, "na.action")
  lmfit$contrasts <- attr(X, "contrasts")
  lmfit$xlevels <- .getXlevels(mt, mf)
  lmfit$call <- cl[c(1L, m)]
  lmfit$call[[1L]] <- quote(lm)
  if(length(contr.list)>0) lmfit$call$contrasts <- quote(contr.list)
  lmfit$terms <- mt
  lmfit$model <- mf

  ## Get dfs, sums of squares, f tildes, partial R2s and eigenvalues 
  stats <- mlmtst(fit = lmfit, X = X, type = type, subset = subset)
  SSi <- stats$SSi
  SSe <- stats$SSe
  df.i <- stats$df.i
  df.e <- stats$df.e
  f.tilde <- stats$f.tilde
  r2 <- stats$r2
  e <- stats$e
  
  # Compute p.values
  pv.acc <- mapply(pv.ss, ss = SSi, df.i = df.i, MoreArgs = list(lambda = e))
  
  # Output
  ss <- c(SSi, Residuals = SSe)
  df <- c(df.i, Residuals = df.e)
  ms <- ss/df
  stats.l <- list(df, ss, ms, f.tilde, r2, pv.acc[1, ])
  cmat <- data.frame()
  for(i in seq(along = stats.l)) {
    for(j in names(stats.l[[i]])){
      cmat[j,i] <- stats.l[[i]][j]
    }
  }
  cmat <- as.matrix(cmat)
  colnames(cmat) <- c("Df", "Sum Sq", "Mean Sq", "F value", "R2", "Pr(>F)")

  out <- list("call" = cl,
              "aov.tab" = cmat,
              "type" = type,
              "precision" = pv.acc[2, ],
              "transform" = transform,
              "na.omit" = lmfit$na.action)
  if(fit){
    out$fit <- lmfit
  }

  ## Update class
  class(out) <- c('MLM', class(out)) # Add help for class MLM
  return(out)
}

##' Compute test statistic
##' 
##' This function computes the degrees of freedom, sum of squares, partial R2 and
##' pseudo F statistics for each explanatory variable from \code{fit}.
##' 
##' Different types of sums of squares are available.
##' 
##' @param fit multivariate fit obtained by \code{\link{lm}}.
##' @param X design matrix obtained by \code{\link{model.matrix}}.
##' @param type type of sum of squares (\code{I}, \code{II} or \code{III}). Default is \code{II}.
##' @param subset subset of explanatory variables for which summary statistics will be reported.
##' @param tol \code{e[e/sum(e) > tol]}, where \code{e} is the vector of eigenvalues
##' of the residual covariance matrix.
##' 
##' @author Diego Garrido-Martín
##' 
##' @export
##' 
mlmtst <- function(fit, X, type = "II", subset = NULL, tol = 1e-3){
  
  # Error sum-of-squares and cross-products (SSCP) matrix 
  SSCP.e <- crossprod(fit$residuals) 
  
  ## Residual sum-of-squares and df
  SSe <- sum(diag(SSCP.e))
  df.e <- fit$df.residual # df.e <- (n-1) - sum(df)
  
  ## Total sum-of-squares
  SSt <- sum(diag(crossprod(fit$model[[1L]])))
  
  ## Partial sums-of-squares
  terms <- attr(fit$terms, "term.labels") # Model terms
  n.terms <- length(terms)
  
  if(!is.null(subset)) {
    if(all(subset %in% terms)) {
      iterms <- which(terms %in% subset)
    }else {
      stop(sprintf("Unknown terms in subset: %s",
                   paste0("'", subset[which(! subset %in% terms)], "'", 
                          collapse = ", ")))
    }
  } else {
    iterms <- 1:n.terms
  }
  asgn <- fit$assign
  
  df.i <- SSi <- numeric(n.terms) # Initialize empty
  names(df.i) <- names(SSi) <- terms
  
  if (type == "I"){
    
    effects <- as.matrix(fit$effects)[seq_along(asgn), , drop = FALSE]
    
    for (i in iterms) {
      subs <- which(asgn == i) 
      SSi[i] <- sum(diag(crossprod(effects[subs, , drop = FALSE])))
      df.i[i] <- length(subs)
    }
    
  } else {
    
    sscp <- function(L, B, V){
      LB <- L %*% B
      crossprod(LB, solve(L %*% tcrossprod(V, L), LB))
    }
    
    B <- fit$coefficients     # Coefficients
    V <- solve(crossprod(X))  # V = (X'X)^{-1}
    p <- nrow(B)
    I.p <- diag(p)
    
    # In contrast to car::Anova, intercept 
    # information is not returned for
    # type III sums-of-squares
    
    if (type == "III"){
      
      for (i in iterms){
        subs <- which(asgn == i) 
        L <- I.p[subs, , drop = FALSE] # Hypothesis matrix
        SSi[i] <- sum(diag(sscp(L, B, V)))
        df.i[i] <- length(subs)
      }
      
    } else {
      
      is.relative <- function(term1, term2, factors) {
        all( !( factors[, term1] & ( !factors[, term2] ) ) )
      }
      
      fac <- attr(fit$terms, "factors") 
      for (i in iterms){
        term <- terms[i]
        subs.term <- which(asgn == i)
        if(n.terms > 1) { # Obtain relatives
          relatives <- (1:n.terms)[-i][sapply(terms[-i], 
                                              function(term2) 
                                                is.relative(term, term2, fac))]
        } else { 
          relatives <- NULL
        }
        subs.relatives <- NULL
        for (relative in relatives){
          subs.relatives <- c(subs.relatives, which(asgn == relative))
        }
        L1 <- I.p[subs.relatives, , drop = FALSE] # Hyp. matrix (relatives) 
        if (length(subs.relatives) == 0) {
          SSCP1 <- 0
        } else {
          SSCP1 <- sscp(L1, B, V)
        }
        L2 <- I.p[c(subs.relatives, subs.term), , drop = FALSE] # Hyp. matrix (relatives + term) 
        SSCP2 <- sscp(L2, B, V)
        SSi[i] <- sum(diag(SSCP2 - SSCP1))
        df.i[i] <- length(subs.term)
      }
      
    }
  }
  
  ## subset
  if(!is.null(subset)){
    SSi <- SSi[iterms]
    df.i <- df.i[iterms]
  }
  
  ## pseudo F
  f.tilde <- SSi/SSe*df.e/df.i
  
  ## r.squared
  R2 <- (SSt - SSe)/SSt
  # R2adj <- 1-( (1-R2)*(n-1) / df.e ) 
  r2 <- SSi/SSt
  # r2adj <- 1-( (1-r2)*(n-1) / df.e )
  
  # Get eigenvalues from cov(R)*(n-1)/df.e
  e <- eigen(SSCP.e/df.e, symmetric = T, only.values = T)$values
  e <- e[e/sum(e) > tol]
  
  return(list("SSi" = SSi, "SSe" = SSe, 
              "df.i" = df.i, "df.e" = df.e, 
              "f.tilde" = f.tilde, "r2" = r2,
              "e" = e))
}

##' Compute asymptotic P-values
##' 
##' Description
##' 
##' Details
##' 
##' @param ss numerator of the pseudo-F statistic.
##' @param lambda eigenvalues.
##' @param df.i degrees of freedom of the numerator of the test statistic.
##' @param eps precision limit.
##' 
##' @author Diego Garrido-Martín
##' 
##' @export
##' 
pv.ss <- function(ss, lambda, df.i, eps = 1e-14){
  
  pv.farebrother <- function(ss, lambda, df.i, maxit = 100000, eps = 1e-14){
    pv <- CompQuadForm::farebrother(ss, lambda = lambda, 
                                    h = rep(df.i, length(lambda)), 
                                    maxit = maxit, eps = eps) # Skip warning
    if(pv$ifault %in% c(0,4)){
      return(pv$Qq)
    } else if (pv$ifault %in% c(5,9)){
      return(pv)
    } else {
      print(pv)
      stop ("Unexpected fault indicator in Farebrother method.")
    }
  }
  
  eps0 <- eps
  pv  <- pv.farebrother(ss = ss, lambda = lambda, df.i = df.i, eps = eps)
  while (length(pv) > 1) {
    eps <- eps * 10
    if(eps > 1e-8){
      print(pv)
      stop("Precision > 1e-8.")
    }
    pv  <- pv.farebrother(ss = ss, lambda = lambda, df.i = df.i, eps = eps)
  }
  if (pv < eps0) {
    pv <- eps0
  }
  return(c(pv, eps))
}

##' @author Diego Garrido-Martín
##' 
##' @keywords internal
##'
##' @export
##' 
print.MLM <- function (x, digits = max(getOption("digits") - 2L, 3L), ...){
  
  ## Print Call and type of SS
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Type", x$type, "Sum of Squares\n\n")
  
  ## Print ANOVA table
  cmat <- x$aov.tab
  if (!is.null(heading <- attr(cmat, "heading"))) 
    cat(heading, sep = "\n")
  nc <- dim(cmat)[2L]
  if (is.null(cn <- colnames(cmat))) 
    stop("'anova' object must have colnames")
  has.P <- grepl("^(P|Pr)\\(", cn[nc])
  zap.i <- 1L:(if (has.P) 
    nc - 1
    else nc)
  i <- which(substr(cn, 2, 7) == " value")
  i <- c(i, which(!is.na(match(cn, "F"))))
  if (length(i)) 
    zap.i <- zap.i[!(zap.i %in% i)]
  tst.i <- i
  if (length(i <- grep("Df$", cn))) 
    zap.i <- zap.i[!(zap.i %in% i)]
  printCoefmat.mp(cmat, digits = digits, has.Pvalue = TRUE, 
                  P.values = TRUE, cs.ind = NULL, zap.ind = zap.i, 
                  tst.ind = tst.i, na.print = "", eps.Pvalue = x$precision + 1e-30, ...)
  na <- x$na.omit
  if(!is.null(na)){
    cat(sprintf("%s observation%s deleted due to missingness\n", 
                length(na), ifelse(length(na) > 1, "s", "")))
  }
  invisible(x)
}

##' Print Coefficient Matrices (multiple p-value precision limits)
##' 
##' Function \code{\link{printCoefmat}} modified to use multiple p-value 
##' precision limits in higher-level print methods.
##' 
##' @seealso \code{\link{printCoefmat}}.
##' 
##' @author Diego Garrido-Martín
##' 
##' @keywords internal
##' 
printCoefmat.mp <- function (x, digits = max(3L, getOption("digits") - 2L),
                             signif.stars = getOption("show.signif.stars"), 
                             signif.legend = signif.stars, 
                             dig.tst = max(1L,min(5L, digits - 1L)), 
                             cs.ind = 1:k, tst.ind = k + 1, zap.ind = integer(), 
                             P.values = NULL, has.Pvalue = nc >= 4 && 
                               substr(colnames(x)[nc], 1, 3) == "Pr(", 
                             eps.Pvalue = .Machine$double.eps, na.print = "NA", 
                             ...) {
  if (is.null(d <- dim(x)) || length(d) != 2L) 
    stop("'x' must be coefficient matrix/data frame")
  nc <- d[2L]
  if (is.null(P.values)) {
    scp <- getOption("show.coef.Pvalues")
    if (!is.logical(scp) || is.na(scp)) {
      warning("option \"show.coef.Pvalues\" is invalid: assuming TRUE")
      scp <- TRUE
    }
    P.values <- has.Pvalue && scp
  }
  else if (P.values && !has.Pvalue) 
    stop("'P.values' is TRUE, but 'has.Pvalue' is not")
  if (has.Pvalue && !P.values) {
    d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
    nc <- nc - 1
    has.Pvalue <- FALSE
  }
  else xm <- data.matrix(x)
  k <- nc - has.Pvalue - (if (missing(tst.ind)) 
    1
    else length(tst.ind))
  if (!missing(cs.ind) && length(cs.ind) > k) 
    stop("wrong k / cs.ind")
  Cf <- array("", dim = d, dimnames = dimnames(xm))
  ok <- !(ina <- is.na(xm))
  for (i in zap.ind) xm[, i] <- zapsmall(xm[, i], digits)
  if (length(cs.ind)) {
    acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
    if (any(ia <- is.finite(acs))) {
      digmin <- 1 + if (length(acs <- acs[ia & acs != 0])) 
        floor(log10(range(acs[acs != 0], finite = TRUE)))
      else 0
      Cf[, cs.ind] <- format(round(coef.se, max(1L, digits - 
                                                  digmin)), digits = digits)
    }
  }
  if (length(tst.ind)) 
    Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst), 
                            digits = digits)
  if (any(r.ind <- !((1L:nc) %in% c(cs.ind, tst.ind, if (has.Pvalue) nc)))) 
    for (i in which(r.ind)) Cf[, i] <- format(xm[, i], digits = digits)
  ok[, tst.ind] <- FALSE
  okP <- if (has.Pvalue) 
    ok[, -nc]
  else ok
  x1 <- Cf[okP]
  dec <- getOption("OutDec")
  if (dec != ".") 
    x1 <- chartr(dec, ".", x1)
  x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
  if (length(not.both.0 <- which(x0 & !is.na(x0)))) {
    Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(1L, 
                                                                    digits - 1L))
  }
  if (any(ina)) 
    Cf[ina] <- na.print
  if (P.values) {
    if (!is.logical(signif.stars) || is.na(signif.stars)) {
      warning("option \"show.signif.stars\" is invalid: assuming TRUE")
      signif.stars <- TRUE
    }
    if (any(okP <- ok[, nc])) {
      pv <- as.vector(xm[, nc])
      # Added ================================================================ #
      Cf[okP, nc] <- mapply(format.pval, pv = pv[okP], eps = eps.Pvalue, 
                            MoreArgs = list(digits = dig.tst))
      # Removed
      # Cf[okP, nc] <- format.pval(pv[okP], digits = dig.tst, eps = eps.Pvalue)
      # ====================================================================== #
      signif.stars <- signif.stars && any(pv[okP] < 0.1)
      if (signif.stars) {
        Signif <- symnum(pv, corr = FALSE, na = FALSE, 
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                         symbols = c("***", "**", "*", ".", " "))
        Cf <- cbind(Cf, format(Signif))
      }
    }
    else signif.stars <- FALSE
  }
  else signif.stars <- FALSE
  print.default(Cf, quote = FALSE, right = TRUE, na.print = na.print, 
                ...)
  if (signif.stars && signif.legend) {
    if ((w <- getOption("width")) < nchar(sleg <- attr(Signif, 
                                                       "legend"))) 
      sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
    cat("---\nSignif. codes:  ", sleg, sep = "", fill = w + 
          4 + max(nchar(sleg, "bytes") - nchar(sleg)))
  }
  invisible(x)
}

#' Biomarkers
#'
#' A simulated dataset containing the levels of 5 biomarkers, 
#' measured in 100 individuals, with different scales. 
#' Missing observations appear as \code{NA}.
#'
#' @format A matrix with 100 rows and 5 numerical variables:
#' \describe{
#'   \item{biomarker1}{levels of biomarker1}
#'   \item{biomarker2}{levels of biomarker2}
#'   ...
#' }
#' @author Diego Garrido-Martín
#' 
"biomarkers"

#' Patients
#'
#' A simulated dataset containing the gender, age and disease status of 100
#' individuals. Missing observations appear as \code{NA}.
#'
#' @format A matrix with 100 rows and 3 variables:
#' \describe{
#'   \item{gender}{Gender of the patient (factor with levels: \code{male} and 
#'   \code{female})}
#'   \item{age}{Age of the patient (numerical)}
#'   \item{disease}{Disease status of the patient (ordered factor with levels
#'   \code{healthy}, \code{"mild"}, \code{"severe"})}
#' }
#' @author Diego Garrido-Martín
#' 
"patients"
