## Main egpd fitting function and S3 methods
## Adapted from evgam by Ben Youngman

#' Fit Extended Generalized Pareto Distribution GAMs
#'
#' @param formula a formula or list of formulae
#' @param data a data frame
#' @param family a character string: "egpd", "degpd", or "zidegpd"
#' @param correctV logical: should variance-covariance matrix account for smoothing parameter uncertainty? Defaults to TRUE
#' @param rho0 initial log smoothing parameters
#' @param inits initial parameter values
#' @param outer outer optimization method: "bfgs" (default), "newton", "fd", or "fixed"
#' @param control a list of control parameters
#' @param removeData logical: should data be removed from the returned object? Defaults to FALSE
#' @param trace an integer controlling output verbosity
#' @param knots a list of knot values for smooth terms
#' @param maxdata maximum number of data rows
#' @param maxspline maximum number of rows for spline basis construction
#' @param compact logical: use compact representation? Defaults to FALSE
#' @param egpd.args a list of arguments for EGPD family (e.g., m=1)
#' @param degpd.args a list of arguments for DEGPD family (e.g., m=1)
#' @param zidegpd.args a list of arguments for ZIDEGPD family (e.g., m=1)
#' @param sandwich.args a list of sandwich correction arguments
#' @param custom.fns a list of custom likelihood functions
#' @param sp fixed smoothing parameters (if supplied, outer optimization is skipped)
#' @param gamma a gamma multiplier for the likelihood
#'
#' @return An object of class \code{egpd}
#'
#' @export
egpd <- function(formula, data, family="egpd", correctV=TRUE, rho0=0,
inits=NULL, outer="bfgs", control=NULL, removeData=FALSE, trace=0,
knots=NULL, maxdata=1e20, maxspline=1e20, compact=FALSE,
egpd.args=list(), degpd.args=list(), zidegpd.args=list(),
sandwich.args=list(), custom.fns=list(), sp=NULL, gamma=1) {

## setup family
family.info <- .setup.family.egpd(family, egpd.args, degpd.args, zidegpd.args, formula, custom.fns)
family <- family.info$family

## setup formulae
formula <- .setup.formulae(formula, family.info$npar, family.info$npar2, data, trace, family.info$nms)
response.name <- attr(formula, "response.name")

## setup mgcv objects and data
temp.data <- .setup.data(data, response.name, formula, family, family.info$nms,
  removeData, knots, maxdata, maxspline, compact, sandwich.args,
  tolower(outer), trace, gamma)
data <- temp.data$data

## initialise inner iteration
beta <- .setup.inner.inits(inits, temp.data$lik.data, family.info$lik.fns, family.info$npar, family)
lik.data <- .sandwich(temp.data$lik.data, beta)
if (trace > 0 & lik.data$adjust > 0) cat(paste("\n Sandwich correct lambda =", signif(lik.data$k, 3), "\n"))

## check whether any smoothing parameters need estimating
smooths <- length(temp.data$gotsmooth) > 0

if (smooths) {

## initialise outer iteration
S.data <- .joinSmooth(temp.data$gams)
nsp <- length(attr(S.data, "Sl"))
if (is.null(rho0)) {
    diagSl <- sapply(attr(S.data, "Sl"), diag)
    rho0 <- apply(diagSl, 2, function(y) uniroot(.guess, c(-1e2, 1e2), d=attr(beta, "diagH"), s=y)$root)
} else {
    if (length(rho0) == 1) rho0 <- rep(rho0, nsp)
}

## check for fixed smoothing parameters
if (!is.null(sp)) {
  rho0 <- log(sp)
  lik.data$outer <- "fixed"
}

lik.data$S <- .makeS(S.data, exp(rho0))

## perform outer iteration
fit.reml <- .outer(rho0, beta, family.info$lik.fns, lik.data, S.data, control, correctV, lik.data$outer, trace)

sp <- exp(fit.reml$par)
lik.data$S <- .makeS(S.data, sp)

} else {

S.data <- NULL
fit.reml <- .outer.nosmooth(beta, family.info$lik.fns, lik.data, control, trace)

}

## covariance matrices
VpVc <- .VpVc(fit.reml, family.info$lik.fns, lik.data, S.data, correctV=correctV, sandwich=temp.data$sandwich, smooths=smooths, trace=trace)

## effective degrees of freedom
edf <- .edf(fit.reml$beta, family.info$lik.fns, lik.data, VpVc, temp.data$sandwich)

## update mgcv objects
names(temp.data$gams) <- family.info$nms
gams <- .swap(fit.reml, temp.data$gams, lik.data, VpVc, temp.data$gotsmooth, edf, smooths)

## add extra things that make an egpd object
gams <- .finalise(gams, data, family.info$lik.fns, lik.data, S.data, fit.reml, VpVc, family, temp.data$gotsmooth, formula, response.name, removeData, edf, family.info$nms2)

return(gams)
}

#' Extract Model Fitted Values
#'
#' @param object a fitted \code{egpd} object
#' @param ... not used
#'
#' @return Fitted values extracted from the object
#'
#' @export
fitted.egpd <- function(object, ...) {
predict(object)
}

#' Simulations from a fitted \code{egpd} object
#'
#' @param object a fitted \code{egpd} object
#' @param nsim an integer giving the number of simulations
#' @param seed an integer giving the seed for simulations
#' @param newdata a data frame
#' @param type a character string: "link" or "response"
#' @param probs a scalar or vector of probabilities
#' @param threshold a scalar added to simulations
#' @param marginal logical: should uncertainty integrate out smoothing parameter uncertainty?
#' @param ... additional arguments
#'
#' @return Simulations of parameters
#'
#' @export
simulate.egpd <- function(object, nsim=1e3, seed=NULL, newdata,
  type="link", probs=NULL, threshold=0, marginal=TRUE, ...) {
if (!is.null(probs))
  type <- "quantile"
if (type %in% c("link", "response")) {
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if(is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  family <- object$family
  if (marginal) {
    V.type <- "Vc"
  } else {
    V.type <- "Vp"
  }
  B <- .pivchol_rmvn(nsim, object$coefficients, object[[V.type]])
  idpars <- object$idpars
  X <- predict.egpd(object, newdata, type="lpmatrix")
  nms <- names(X)
  B <- lapply(seq_along(X), function(i) B[idpars == i, , drop=FALSE])
  X <- lapply(seq_along(X), function(i) X[[i]] %*% B[[i]])
  names(X) <- nms
  if (type == "response") {
    unlink <- which(substr(nms, 1, 3) == "log")
    for (i in unlink) {
      X[[i]] <- exp(X[[i]])
      if (substr(nms[i], 1, 5) == "logit")
        X[[i]] <- X[[i]] / (1 + X[[i]])
    }
  }
  nms <- gsub("cloglog", "", nms)
  nms <- gsub("probit", "", nms)
  nms <- gsub("logit", "", nms)
  nms <- gsub("log", "", nms)
  names(X) <- nms
}
if (type == "quantile") {
  stop("Quantile simulations not yet implemented for egpd families.")
}
return(X)
}

#' Log-likelihood from a fitted \code{egpd} object
#'
#' @param object a fitted \code{egpd} object
#' @param ... not used
#'
#' @return A \code{logLik} object
#'
#' @export
logLik.egpd <- function(object, ...) {
if (!missing(...)) warning("extra arguments discarded")
out <- object$logLik
attr(out, "df") <- attr(object, "df")
attr(out, "nobs") <- nobs(object)
class(out) <- "logLik"
out
}

#' Plot a fitted \code{egpd} object
#'
#' @param x a fitted \code{egpd} object
#' @param onepage logical: should all plots be on one page? Defaults to TRUE
#' @param which a vector of integers identifying which smooths to plot
#' @param main a character string or vector of plot titles
#' @param ask logical: ask to show next plots?
#' @param ... extra arguments to pass to \code{mgcv::plot.gam}
#'
#' @return Plots representing smooth terms
#'
#' @export
plot.egpd <- function(x, onepage = TRUE, which = NULL, main, ask = !onepage, ...) {
x <- x[x$gotsmooth]
if (is.null(which)) {
  nplot <- sum(unlist(lapply(x, function(x) as.integer(sapply(x$smooth, function(y) y$plot.me)))))
  which <- seq_len(nplot)
} else {
  nplot <- length(which)
}
if (onepage) {
  omfrow <- par("mfrow")
  nmfrow <- rev(n2mfrow(nplot))
  par(mfrow = nmfrow)
}
if (ask) {
  oask <- par("ask")
  if (nplot > prod(par("mfrow")) && dev.interactive()) {
    par(ask = TRUE)
  } else {
    ask <- FALSE
  }
}
current <- 1
for (i in seq_along(x)) {
  for (j in seq_along(x[[i]]$smooth)) {
    if (current %in% which) {
      if (missing(main)) {
        mgcv::plot.gam(x[[i]], select = j, main = paste(names(x)[i], x[[i]]$smooth[[j]]$label, sep = ": "), ...)
      } else {
        mgcv::plot.gam(x[[i]], select = j, ...)
      }
    }
    current <- current + 1
  }
}
if (onepage)
  par(mfrow = omfrow)
if (ask)
  par(ask = oask)
}

#' Summary method for a fitted \code{egpd} object
#'
#' @param object a fitted \code{egpd} object
#' @param ... not used
#'
#' @return A \code{summary.egpd} object
#'
#' @export
summary.egpd <- function(object, ...) {
if (!missing(...)) warning("extra arguments discarded")
out <- list()
out[[1]] <- .parametric.summary.egpd(object)
out[[2]] <- .smooth.summary.egpd(object)
class(out) <- "summary.egpd"
out
}

#' @param x a \code{summary.egpd} object
#' @rdname summary.egpd
#' @export
print.summary.egpd <- function(x, ...) {
if (!missing(...)) warning("extra arguments discarded")
cat("\n")
cat("** Parametric terms **")
tab <- lapply(x[[1]], .tidyParametricTable)
cat("\n")
for (i in seq_along(tab)) {
cat("\n")
cat(names(tab)[[i]])
cat("\n")
print(tab[[i]])
}
cat("\n")
cat("** Smooth terms **")
tab <- lapply(x[[2]], .tidySmoothTable)
cat("\n")
for (i in seq_along(tab)) {
cat("\n")
cat(names(tab)[[i]])
cat("\n")
print(tab[[i]])
}
invisible(x)
}

#' Print a fitted \code{egpd} object
#'
#' @param x a fitted \code{egpd} object
#' @param ... not used
#'
#' @return The call of the object (invisibly)
#'
#' @export
print.egpd <- function(x, ...) {
if (!missing(...)) warning("extra arguments discarded")
print(x$call)
invisible(x)
}

#' Bind a list of data frames
#'
#' @param x a list of data frames
#'
#' @return A data frame
#'
#' @export
dfbind <- function(x) {
nms <- names(x[[1]])
cls <- sapply(x[[1]], class)
x <- lapply(nms, function(i) unlist(lapply(x, function(y) y[,i])))
x <- as.data.frame(x)
dt <- cls == "Date"
if (any(dt)) {
  for (i in which(dt)) x[,i] <- as.Date(x[,i], origin="1970-01-01")
}
names(x) <- nms
x
}

#' Moore-Penrose pseudo-inverse of a matrix
#'
#' @param x a matrix
#' @param tol a scalar tolerance
#'
#' @return A matrix
#'
#' @export
pinv <- function(x, tol=-1) {
armapinv(x, tol)
}

#' @rdname pinv
#' @export
ginv.egpd <- function(x, tol=sqrt(.Machine$double.eps)) {
armaginv(x, tol)
}

#' Generate a sequence between a range
#'
#' @param x a 2-vector
#' @param length an integer
#'
#' @return A vector
#'
#' @export
seq_between <- function(x, length=NULL) {
if (is.null(length)) {
    return(seq(x[1], x[2]))
    } else {
    return(seq(x[1], x[2], length=length))
    }
}

#' Number of observations
#' @param object an egpd object
#' @param ... not used
#' @return integer
#' @export
nobs.egpd <- function(object, ...) {
object$nobs
}
