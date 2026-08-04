## Main egpd fitting function and S3 methods
## Adapted from evgam by Ben Youngman

#' Fit Extended Generalized Pareto Distribution GAMs
#'
#' @param formula a formula or list of formulae
#' @param data a data frame
#' @param family a character string: "egpd", "degpd", "zidegpd",
#'   "comppareto", or "custom".
#'   Continuous zero-inflated EGPD utilities are provided elsewhere in the
#'   package, but are not fitted via \code{egpd()}.
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
#' @param degpd.args a list of arguments for the DEGPD family. \code{m} selects the
#'   G-function (model 1-6, default 1). \code{xi.max} (default \code{Inf}) caps the
#'   tail index \eqn{\xi} via a bounded shape link
#'   \eqn{\xi = \texttt{xi.max}\,(1 - \exp(-\exp(\eta)/\texttt{xi.max}))}, which reduces
#'   to the usual unbounded \eqn{\xi = \exp(\eta)} log link as \code{xi.max} tends to
#'   \code{Inf} and saturates at \code{xi.max}; this stabilises fits on short, spiky,
#'   or heavily over-dispersed series. Applies to all DEGPD models 1-6.
#' @param zidegpd.args a list of arguments for ZIDEGPD family (e.g., m=1)
#' @param gpig.args a list of arguments controlling how the GPIG/ZIGPIG pmf is
#'   evaluated, passed to \code{\link{gpig_control}}: \code{method} (one of
#'   \code{"hybrid"}, \code{"recursion"}, \code{"trunc"}, \code{"saddlepoint"},
#'   \code{"legacy"}), \code{yswitch}, \code{order}, \code{eps}, \code{sptol},
#'   \code{ymax}, \code{nquad} and \code{normalize}. Applied for the duration of the
#'   fit and then restored. Ignored, with a warning, for other families.
#' @param comppareto.args a list of arguments for the CompPareto GAM family.
#'   Currently this supports \code{spec = "lnorm"}, \code{"gamma"},
#'   \code{"weibull"}, or \code{"exp"} for the body distribution.
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
egpd.args=list(), degpd.args=list(), zidegpd.args=list(), comppareto.args=list(),
gpig.args=list(), sandwich.args=list(), custom.fns=list(), sp=NULL, gamma=1,
fixed=list(), restarts=TRUE) {

## Multi-start for DEGPD model 1. Its likelihood has a nearly flat ridge along which
## sigma -> 0 and kappa -> infinity together (they enter the upper tail only through
## kappa * sigma^(1/xi)), and from a single start the optimiser can settle on it or on
## a saddle -- which shows up as a wildly large kappa, or as a negative variance for a
## parameter that is genuinely estimated. Trying a few starts and keeping the highest
## likelihood costs little and removes both failure modes. Every candidate lies in the
## same parameter space, so this is plain multi-start MLE, not selection across models.
## Skipped when the caller supplies `inits`, or sets restarts = FALSE.
if (isTRUE(restarts) && is.null(inits) && identical(family, "degpd") &&
    (is.null(degpd.args$m) || identical(as.integer(degpd.args$m), 1L))) {
  cl <- match.call()
  ## The caller's frame must be captured HERE, not inside the lapply() below.
  ## match.call() records argument *expressions*, so cl still refers to whatever
  ## symbols the caller wrote -- and calling parent.frame() from inside the lambda
  ## resolves to lapply()'s frame, where those symbols do not exist. Every restart
  ## then failed with "object not found", tryCatch turned that into NULL, and the
  ## function fell through to a single start with no error and no warning. See #6.
  caller_env <- parent.frame()
  ybar <- tryCatch({
    fl <- if (is.list(formula)) formula[[1]] else formula
    mean(data[[all.vars(fl)[1]]], na.rm = TRUE)
  }, error = function(e) NA_real_)
  if (is.finite(ybar)) {
    idlink <- !is.null(degpd.args$link) && identical(degpd.args$link, "identity")
    ## a light and a heavier shape, alongside the package default
    starts <- c(list(NULL), lapply(c(0.05, 0.3), function(v)
                 c(log(ybar + 1), if (idlink) v else log(v), 0)))
    err <- NULL
    cand <- lapply(starts, function(s) {
      cl$inits <- s; cl$restarts <- FALSE
      tryCatch(suppressWarnings(eval(cl, caller_env)),
               error = function(e) { err <<- conditionMessage(e); NULL })
    })
    cand <- Filter(function(z) !is.null(z) && is.finite(as.numeric(z$logLik)), cand)
    if (length(cand))
      return(cand[[which.max(vapply(cand, function(z) as.numeric(z$logLik), 0))]])
    ## Never degrade quietly: a silent fall-through to one start is the failure mode
    ## that made #6 invisible, and it can cost thousands of log-likelihood units.
    warning("restarts requested but every restart failed; continuing from a single ",
            "start", if (!is.null(err)) paste0(" (", err, ")"), call. = FALSE)
  }
}

## setup family
family.info <- .setup.family.egpd(family, egpd.args, degpd.args, zidegpd.args, comppareto.args, formula, custom.fns)
family <- family.info$family

## `fixed`: hold named parameters at known values. This is sugar for the offset-only
## pinning `~ offset(rep(<link>(v), n)) - 1`, which is easy to get wrong by hand because
## the offset must be on the LINK scale (log for l-prefixed parameters, identity
## otherwise). Values here are given on the RESPONSE scale, so fixed = list(kappa = 1)
## pins kappa itself rather than log kappa.
if (length(fixed)) {
  .fx <- .setup.fixed.egpd(fixed, formula, data, family.info$nms)
  formula <- .fx$formula
  data    <- .fx$data
}

## GPIG/ZIGPIG pmf evaluation route (see gpig_control). Set for the duration of
## this fit only, so a method chosen here cannot leak into later calls.
.gpig.old <- .gpig_apply_control(family, gpig.args)
if (!is.null(.gpig.old)) on.exit(do.call(gpig_control, .gpig.old), add = TRUE)

## setup formulae
formula <- .setup.formulae(formula, family.info$npar, family.info$npar2, data, trace, family.info$nms)
response.name <- attr(formula, "response.name")

## ZIDEGPD nests DEGPD (zero-inflation pi -> 0). Warm-start the richer model from
## the DEGPD MLE so a converged ZIDEGPD can never sit below the simpler model's
## likelihood (the nesting guarantee). The last ZIDEGPD parameter is always
## logitpi, so dropping it yields the matching DEGPD formula. Seed (sigma, xi,
## kappa, ...) from DEGPD's per-parameter linear-predictor means and pi near 0.
if (family == "zidegpd" && is.null(inits) && length(response.name) == 1) {
  inits <- tryCatch({
    zi_m <- if (is.null(zidegpd.args$m)) 1 else zidegpd.args$m
    dfit <- egpd(formula[-length(formula)], data, family = "degpd",
                 correctV = FALSE, inits = NULL, outer = outer, knots = knots,
                 degpd.args = list(m = zi_m), trace = 0, gamma = gamma)
    lp <- predict(dfit, type = "link")
    p0 <- min(max(mean(data[[response.name]] == 0, na.rm = TRUE), 1e-3), 0.5)
    c(vapply(lp, mean, numeric(1)), qlogis(p0))
  }, error = function(e) {
    if (trace > 0) message("ZIDEGPD warm-start from DEGPD failed; using default inits.")
    NULL
  })
}

## setup mgcv objects and data
temp.data <- .setup.data(data, response.name, formula, family, family.info$nms,
  removeData, knots, maxdata, maxspline, compact, sandwich.args,
  tolower(outer), trace, gamma)
data <- temp.data$data

## issue #2: bounded shape link for DEGPD model 1. xi.max (via degpd.args) caps the
## tail index; the default Inf reproduces the usual unbounded exp() link. Only the
## degpd1 likelihood wrapper reads this, so it is inert for other families/models.
temp.data$lik.data$xi.max <- if (!is.null(degpd.args$xi.max)) degpd.args$xi.max else Inf
temp.data$lik.data$bounded.xi <- (family == "degpd" &&
  is.finite(temp.data$lik.data$xi.max))

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

## record any `fixed` pins so predict() can re-create their synthetic columns for
## newdata, and so the constraint is visible on the fitted object
if (length(fixed)) {
  gams$fixed.offsets <- .fx$offsets
  gams$fixed         <- .fx$values
}

## Warn on the DEGPD kappa ridge. sigma and kappa enter the upper tail only through
## kappa * sigma^(1/xi), so they are separated by the bulk alone; on a short series
## that separation is weak and the optimum drifts along a flat ridge to very large
## kappa and near-zero sigma. xi usually survives this -- it is not part of the
## confounding -- but sigma and kappa individually do not, so quantiles and any
## interval that leans on them should be treated as provisional.
if (identical(attr(family, "type"), 1) && grepl("degpd", family) && !removeData) {
  kap <- tryCatch(suppressWarnings(max(as.data.frame(
           stats::predict(gams, type = "response"))$kappa, na.rm = TRUE)),
         error = function(e) NA_real_)
  if (is.finite(kap) && kap > 1e4)
    warning(sprintf(paste0("fitted kappa reaches %.3g: the DEGPD scale/carrier pair is ",
                           "effectively unidentified (the sigma-kappa ridge). xi is ",
                           "usually still identified; sigma and kappa are not. Consider ",
                           "fixed = list(kappa = 1) for a plain discretised GPD."), kap),
            call. = FALSE)
}

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
  nms <- sub("^cloglog", "", nms)
  nms <- sub("^probit", "", nms)
  nms <- sub("^logit", "", nms)
  nms <- sub("^log", "", nms)
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
