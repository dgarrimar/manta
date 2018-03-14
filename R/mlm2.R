##' Multivariate linear models
##' 
##' Transforms the response variables and fits a multivariate
##' linear model using \code{\link{lm}}. 
##' 
##' A "\code{Y}" matrix is obtained after transforming and centering the 
##' original response variables. Then, the multivariate fit obtained by 
##' \code{\link{lm}} can be used to obtain sums of squares (I, II, III), 
##' pseudo F statistics and asymptotic p-values for the explanatory variables in
##' a non-parametric manner.
##' 
##' @param formula object of class "\code{\link{formula}}": a symbolic 
##' description of the model to be fitted.
##' @param data an optional data frame, list or environment (or object 
##' coercible by \code{\link{as.data.frame}} to a data frame) containing the 
##' variables in the model. If not found in data, the variables are taken from 
##' \code{environment(formula)}, typically the environment from which \code{mlm}
##' is called.
##' @param distance data transformation. One of c("\code{euclidean}", "\code{hellinger}"). 
##' Default is "\code{euclidean}".
##' @param contrasts an optional list. See the \code{contrasts.arg} of 
##' \code{\link{model.matrix.default}}. Default is "\code{\link{contr.sum}}" 
##' for ordered factors and "\code{\link{contr.poly}}" for unordered factors. 
##' Note that this is different from the default setting in \code{\link{options}
##' ("contrasts")}.
##' @param ... additional arguments to be passed to \code{\link{lm}}.
##' 
##' @return \code{mlm} returns an object of \code{\link{class}} c("mlm2", "mlm", "lm"); a list containing:
##' \item{coefficients}{a named vector of coefficients.}
##' \item{residuals}{the residuals, that is response minus fitted values.}
##' \item{fitted.values}{the fitted mean values.}
##' \item{rank}{the numeric rank of the fitted linear model.}
##' \item{df.residual}{the residual degrees of freedom.}
##' \item{call}{the matched call.}
##' \item{terms}{the terms object used (transformed response).}
##' \item{contrasts}{(only where relevant) the contrasts used.}
##' \item{xlevels}{(only where relevant) a record of the levels of the factors used in fitting.}
##' \item{model}{the model frame used (transformed response).}
##' \item{na.action}{(where relevant) information returned by \code{\link{model.frame}} on the special handling of NAs.}
##' \item{transformation}{transformations of the response matrix.}
##' In addition, non-null fits will have components assign, effects and qr 
##' relating to the linear fit. See \code{\link{lm}} for further details. 
##' Note the multivariate fit is done on the transformed response \code{"Y"}.
##' 
##' @seealso \code{\link{summary.mlm2}}
##' 
##' @author Diego Garrido-Martín
##' @import stats
##' @export
mlm2 <- function(formula, data, distance = "euclidean", contrasts = NULL, ...){
  
  ## Define model frame, response and explanatory variables 
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- "na.pass"
  mf[[1L]] <- quote(stats::model.frame) # need stats:: for non-standard evaluation
  mf <- eval(mf, parent.frame())

  Y <- model.response(mf, "numeric")         # Response variables - will be transformed
  X <- mf[, colnames(mf)[-1], drop = FALSE]  # Explanatory variables - no change
  
  ## Checks
   # > 1 response variable
   # arguments
   # response = factor
   # [...]
  distance <- match.arg(distance, c("euclidean", "hellinger"))
  
  ## Transformation (distance)
  if(distance == "hellinger"){
    Y <- sqrt(Y)
  }
  
  ## Center response variables
  Y <- scale(Y, center = TRUE, scale = FALSE)

  ## Define contrasts
  if(is.null(contrasts)){
    contrasts <- list(unordered = "contr.sum", ordered = "contr.poly")
    contr.list <- lapply(1:ncol(X), FUN = function(k){
      # Default: no contrast for quantitaive predictor
      contr.type <- NULL
      # Sum contrasts for unordered categorical predictor
      if(is.factor(X[,k])){
        contr.type <- contrasts$unordered
      }
      # Polynomial contrasts for ordered categorical predictor
      if(is.ordered(X[,k])){
        contr.type <- contrasts$ordered
      }
      return(contr.type)
    })
    names(contr.list) <- colnames(X)
    contr.list <- contr.list[!unlist(lapply(contr.list, is.null))]
  } else {
    contr.list <- contrasts
  }
  
  ## Update formula
  if (!missing(data)) 
    formula <- terms(formula, data = data)
    formula <- update(formula, Y ~ .)
    ## no data? find variables in .GlobalEnv
  
  ## Fit lm 
  fit <- lm(formula, data = X, contrasts = contr.list, ...)

  ## Update object call and class
  fit$call <- cl
  fit$transformation <- 
    attributes(fit$model)$transformation <- 
    c(sprintf('distance:%s', distance), 'scaled:center')
  class(fit) <- c('mlm2', class(fit))
  
  return(fit)
}


##' @author Diego Garrido-Martín
##' @keywords internal
##' @importFrom car Anova
##' @export
print.mlm2 <- function (x, type = "II", digits = max(3L, getOption("digits") - 3L), 
                        ...){
  
  ## Checks
  type <- match.arg(type, c("I", "II", "III", 1, 2, 3))
  
  ## Print "Call" and "Terms"
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Terms (Type", type, "Sum of Squares):\n")
  
  ## Compute sums of squares 
  if (type == "I" || type == 1){
    mnv <- summary(manova(x))    # Intercept possible here, but always 0 due to centering
    SSP <- mnv$SS[-length(mnv$SS)]
    SSPE <- mnv$SS$Residuals
  } else {
    if((type == "III" || type == 3) && 
       any(unlist(x$contrasts) %in% c("contr.treatment", "contr.SAS"))){
      warning(strwrap("Type III Sum of Squares requires effect- or orthogonal
                      coding for unordered categorical variables (i.e. contr.sum,
                      contr.helmert)."))
    }
    UU <- car::Anova(x, type = type) 
    SSP <- UU$SSP
    SSPE <- UU$SSPE
  }
  SS <- lapply(SSP, function(x){diag(x)})
  SSe <- diag(SSPE)
  
  ## Degrees of freedom
  if(type == "III"){
    Df <- table(x$assign)
    names(Df) <- c("(Intercept)", attributes(x$terms)$term.labels)
  } else{
    Df <- table(x$assign)[-1]
    names(Df) <- attributes(x$terms)$term.labels
  }
  df.e <- x$df.residual # df.e <- (n-1) - sum(Df)
  
  ## Build table to print 
  ss <- do.call(cbind, c(SS, list(Residuals = SSe)))
  rse <- sapply(sqrt(ss[,"Residuals"]/df.e), format)
  ss <- apply(zapsmall(ss), 2L, format)
  rownames(ss) <- paste0("resp ", 1:nrow(ss))
  out <- rbind(ss, "Deg. of Freedom" = c(Df, df.e))
  print(out, quote = F, right = T)
  cat("\nResidual standard errors: ", rse, "\n")
  na <- attributes(x$model)$na.action
  if(!is.null(na)){
    cat(sprintf("%s observation%s deleted due to missingness\n", 
                length(na), ifelse(length(na) > 1, "s", "")))
  }
  invisible(x)
}

##' Summary method for Non-parametric Multivariate Analysis of Variance
##' 
##' A summary method for class "\code{mlm2}".
##' 
##' @param object an object of class "\code{mlm2}"
##' @param type type of sum of squares (either "I", "II", "III" or 1, 2, 3).
##' @param ... further arguments passed to or from other methods.
##' 
##' @importFrom car Anova
##' @export
summary.mlm2 <- function(object, type = "II", ...){
  
  ## Type of sum of squares
  type <- match.arg(type, c("I", "II", "III", 1, 2, 3))
  
  ## Get residuals and sample size
  R <- object$residuals
  n <- nrow(R)
  
  ## Compute sums of squares 
  if (type == "I" || type == 1){
    mnv <- summary(manova(object))    # Intercept possible here, but always 0 due to centering
    SSP <- mnv$SS[-length(mnv$SS)]
    SSPE <- mnv$SS$Residuals
  } else {
    if((type == "III" || type == 3) && any(unlist(object$contrasts) %in% c("contr.treatment", "contr.SAS"))){
      warning(strwrap("Type III Sum of Squares requires effect- or orthogonal
                      coding for unordered categorical variables (i.e. contr.sum,
                      contr.helmert).")) 
    }
    UU <- car::Anova(object, type = type) # Intercept here? McArtor
    SSP <- UU$SSP
    SSPE <- UU$SSPE
    }
  SS <- lapply(SSP, function(x){sum(diag(x))})
  SSe <- sum(diag(SSPE))
  
  ## Compute pseudo F's
  f.tilde <- unlist(lapply(SS, function(x){x/SSe}))
  
  ## Degrees of freedom
  if(type == "III"){
    Df <- table(object$assign)
    names(Df) <- c("(Intercept)", attributes(object$terms)$term.labels)
  } else{
    Df <- table(object$assign)[-1]
    names(Df) <- attributes(object$terms)$term.labels
  }
  
  df.e <- object$df.residual # df.e <- (n-1) - sum(Df)
  
  ## Get response 
  Y <- object$model[,1]
  
  ## Compute omnibus r.squared and adj.r.squared (per response variable and global)
  sscp <- crossprod(Y)
  R2 <- (sscp-SSPE)/sscp
  R2 <- c(overall = sum(diag(sscp-SSPE))/sum(diag(sscp)),diag(R2))
  R2adj <- 1-( (1-R2)*(n-1) / df.e )
  
  # Compute r.squared per explanatory variable (and response variable)
  r2 <- lapply(SSP, function(x){c(overall = sum(diag(x))/sum(diag(sscp)), diag(x/sscp))})
  r2adj <- lapply(r2, function(x){1-( (1-x)*(n-1) / df.e )})
  
  # Get eigenvalues from R
  e <- eigen(cov(R)*(n-1)/df.e, symmetric = T, only.values = T)$values
  
  # Compute p.values
  pv.acc <- mapply(pv.f, f = f.tilde, df.i = Df, MoreArgs = list(df.e = df.e, lambda = e))
  
  # Output 
  SS <- c(SSP, list(Residuals = SSPE)) 
  attributes(SS)$type <- type
  
  r.squared <- unlist(lapply(r2, function(x){as.numeric(x[1])}))
  stats.l <- list(c(Df, Residuals = df.e), f.tilde*df.e/Df, r.squared, pv.acc[1,])
  cmat <- data.frame()
  for(i in seq(along = stats.l)) {
    for(j in names(stats.l[[i]])){
      cmat[j,i] <- stats.l[[i]][j]
    }
  } 
  cmat <- as.matrix(cmat)
  colnames(cmat) <- c("Df", "F", "r2", "P(>F)")
  attributes(cmat)$precision <- pv.acc[2,]
  
  x <- list("row.names" = rownames(cmat),
            "SS" = SS,
            "Eigenvalues" = e,
            "stats" = as.matrix(cmat))
  
  class(x) <- 'summary.mlm2'
  return(x)
  }

##' Compute asymptotic P-values
##' 
##' Description
##' 
##' Details
##' 
##' @param f pseudo-F statistic.
##' @param lambda eigenvalues
##' @param df.i degrees of freedom of the variable
##' @param df.e residual degrees of freedom
##' @param acc precision limit 
##' 
##' @export
pv.f <- function(f, lambda, df.i, df.e, acc = 1e-14){
  
  pv.davies <- function(f, lambda, df.i, df.e, lim = 50000, acc = 1e-14){
    H <- c(rep(df.i, length(lambda)), rep(df.e, length(lambda)))
    pv <- CompQuadForm::davies(0, lambda = c(lambda, -f * lambda), h = H, lim = lim, acc = acc)
    if(pv$ifault != 0 || pv$Qq < 0 || pv$Qq > 1){
      return(pv)
    } else {
      return(pv$Qq)
    }
  }
  
  pv <- pv.davies(f = f, lambda = lambda, df.i = df.i, df.e = df.e, acc = acc)
  while (length(pv) > 1) {
    acc <- acc * 10
    pv  <- pv.davies(f = f, lambda = lambda, df.i = df.i, df.e = df.e, acc = acc)
  }
  if (pv < acc) {
    pv <- acc
  }
  return(c(pv, acc))
}

##' @export
print.summary.mlm2 <- function(x, digits = max(getOption("digits") - 2L, 3L), tol = 1e-30, ...){
  
  cmat <- x$stats
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
                  tst.ind = tst.i, na.print = "", eps.Pvalue = attributes(cmat)$precision + tol, ...)
  invisible(x)
}

##' Print Coefficient Matrices (multiple p-value precision limits)
##' 
##' Function \code{\link{printCoefmat}} modified to use multiple p-value 
##' precision limits in higher-level print methods, such as those for 
##' \code{\link{summary.mlm2}}.
##' 
##' @seealso \code{\link{printCoefmat}}.
##' 
##' @keywords internal
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

