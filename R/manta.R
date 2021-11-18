##' Non-parametric, Asymptotic P-values for Multivariate Linear Models
##' 
##' Fits a multivariate linear model and computes test statistics and asymptotic 
##' P-values for predictors in a non-parametric manner. 
##' 
##' A \code{Y} matrix is obtained after transforming (optionally) and centering 
##' the original response variables. Then, the multivariate fit obtained by 
##' \code{\link{lm}} can be used to compute sums of squares (type-I, type-II or 
##' type-III), pseudo-F statistics and asymptotic P-values for the terms specified
##' by the \code{formula} in a non-parametric manner. The designations "type-II" 
##' and "type-III" correspond exactly to those used in \code{\link[car]{Anova}}. 
##' "type-I" refers to sequential sums of squares.
##' 
##' @param formula object of class "\code{\link{formula}}" (or one that can be 
##' coerced to that class): a symbolic description of the model to be fitted. 
##' @param data an optional data frame, list or environment (or object coercible 
##' by \code{\link{as.data.frame}} to a data frame) containing the variables in 
##' the model. If not found in data, the variables are taken from 
##' \code{environment(formula)}, typically the environment from which \code{manta} 
##' is called.
##' @param transform transformation of the response variables: "\code{none}", 
##' "\code{sqrt}" or "\code{log}". Default is "\code{none}".
##' @param type type of sum of squares: "\code{I}", "\code{II}" or "\code{III}". 
##' Default is "\code{II}".
##' @param contrasts an optional list. See \code{contrasts.arg} in 
##' \code{\link{model.matrix.default}}. Default is "\code{\link{contr.sum}}" 
##' for ordered factors and "\code{\link{contr.poly}}" for unordered factors. 
##' Note that this is different from the default setting in \code{\link{options}("contrasts")}.
##' @param subset subset of predictors for which summary statistics will be 
##' reported. Note that this is different from the "\code{subset}" argument in \code{\link{lm}}.
##' @param fit logical. If \code{TRUE} the multivariate fit on transformed and 
##' centered responses is returned.
##' 
##' @return \code{manta} returns an object of \code{\link{class}} "manta", a list containing:
##' \item{call}{the matched call.}
##' \item{aov.tab}{ANOVA table with Df, Sum Sq, Mean Sq, F values, 
##' partial R-squared and P-values.}
##' \item{type}{the type of sum of squares (\code{"I"}, \code{"II"} or \code{"III"}).}
##' \item{precision}{the precision in P-value computation.}
##' \item{transform}{the transformation applied to the response variables.}
##' \item{na.omit}{incomplete cases removed (see \code{\link{na.omit}}).}
##' \item{fit}{if \code{fit = TRUE} the multivariate fit done on the transformed 
##' and centered response variables is also returned.}
##' 
##' @seealso \code{\link{lm}}, \code{\link[car]{Anova}}
##' 
##' @author Diego Garrido-Martín
##' 
##' @import stats
##' 
##' @export
##' 
manta <- function(formula, data, transform = "none", type = "II", 
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
    stop("The number of response variables should be >= 2")
  }
  if (length(attr(mt, "term.labels")) < 1) {
    stop("The model should contain at least one predictor (excluding the intercept)")
  }
  
  ## Transform and center responses, update model frame
  if(transform == "none"){
    Y <- response
  } else if (transform == "sqrt"){
    if (any(response < 0)) {
      stop("'sqrt' transformation requires all response values >= 0")
    }
    Y <- sqrt(response)
  } else if (transform == "log"){
    if (any(response <= 0)) {
      stop("'log' transformation requires all response values > 0")
    }
    Y <- log(response)
  }
  Y <- scale(Y, center = TRUE, scale = FALSE)
  mf[[1L]] <- Y
  
  ## Define contrasts
  if(is.null(contrasts)){
    contrasts <- list(unordered = "contr.sum", ordered = "contr.poly")
    dc <- attr(mt, "dataClasses")[-1]
    contr.list <- lapply(dc, FUN = function(k){
      # No contrast for quantitative predictors
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
  if(length(contr.list) > 0) lmfit$call$contrasts <- quote(contr.list)
  lmfit$terms <- mt
  lmfit$model <- mf

  ## Compute sums of squares, df's, pseudo-F statistics, partial R2s and eigenvalues 
  stats <- manta.ss(fit = lmfit, X = X, type = type, subset = subset)
  SS <- stats$SS
  df <- stats$df
  f.tilde <- stats$f.tilde
  r2 <- stats$r2
  e <- stats$e

  ## Compute P-values
  l <- length(df) # SS[l], df[l] correspond to Residuals 
  pv.acc <- mapply(p.asympt, ss = SS[-l], df = df[-l], MoreArgs = list(lambda = e))  
  
  ## ANOVA table
  stats.l <- list(df, SS, SS/df, f.tilde, r2, pv.acc[1, ])
  cmat <- data.frame()
  for(i in seq(along = stats.l)) {
    for(j in names(stats.l[[i]])){
      cmat[j, i] <- stats.l[[i]][j]
    }
  }
  cmat <- as.matrix(cmat)
  colnames(cmat) <- c("Df", "Sum Sq", "Mean Sq", "F value", "R2", "Pr(>F)")

  ## Output
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
  class(out) <- c('manta', class(out))
  return(out)
}

##' @author Diego Garrido-Martín
##' 
##' @keywords internal
##'
##' @export
##' 
print.manta <- function (x, digits = max(getOption("digits") - 2L, 3L), ...){ # #nocov start
  
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
    stop("'aov.tab' object must have colnames")
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
} # #nocov end

##' Print Coefficient Matrices (Multiple P-value Precision Limits)
##' 
##' Modification of \code{\link{printCoefmat}} to use multiple P-value 
##' precision limits.
##' 
##' @seealso \code{\link{printCoefmat}}
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
                             ...) { # #nocov start
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
} # #nocov end
