## REML functions
## Adapted from evgam by Ben Youngman

.reml0 <- function(pars, likfns, likdata, Sdata, beta=NULL) {
if (is.null(beta))
  beta <- attr(pars, "beta")
sp <- exp(pars)
likdata$S <- .makeS(Sdata, sp)
## A refit that fails outright (e.g. non-finite derivatives at a trial smoothing
## parameter) is reported as a failed inner fit, which the outer optimisers reject,
## instead of aborting the whole fit.
fitbeta <- try({
  fb <- .newton_step(beta, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns, control=likdata$control$inner)
  ## Retry with a few short gradient steps when Newton stopped short of the gradient
  ## criterion while still predicting a real decrease. A stop at the resolution limit
  ## of the objective (a negligible Newton decrement) cannot be improved by the retry,
  ## which then only doubled the cost of the tight-tolerance finite-difference fits.
  if (!fb$gradconv && fb$convergence != 2L &&
      !isTRUE(.5 * fb$decrement <= likdata$control$inner$fntol)) {
    it0 <- likdata$control$inner$itlim
    likdata$control$inner$itlim <- 10
    fb <- .newton_step(fb$par, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns, control=likdata$control$inner, newton=FALSE, alpha0=.05)
    likdata$control$inner$itlim <- it0
    fb <- .newton_step(fb$par, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns, control=likdata$control$inner)
  }
  fb
}, silent=TRUE)
if (inherits(fitbeta, "try-error"))
  return(structure(1e20, inner.convergence=1L,
                   inner.report=paste("Inner coefficient optimization failed:",
                                      conditionMessage(attr(fitbeta, "condition")))))
logdetSdata <- .logdetS(Sdata, pars)
logdetHdata <- .d0logdetH(fitbeta)
out <- fitbeta$objective + as.numeric(fitbeta$convergence != 0) * 1e20
out <- out + .5 * logdetHdata$d0 - .5 * logdetSdata$d0
## The constants of the Laplace approximation to the negative log restricted likelihood
## (Wood, 2011; the REML score of mgcv): the part of -log|S|+/2 that .logdetS leaves
## out, and -Mp/2 log(2 pi) from the improper flat prior on the Mp unpenalised
## directions: the coefficients in the Laplace integral (those the inner fit kept, as
## in log|H|) less the rank of the total penalty. They used to be +Mp/2 log(2 pi), with
## Mp counting one per linear predictor plus each smooth's null-space dimension, and
## no eigenvalue term. They are added in one final step, so that the differences of
## the criterion -- all the outer optimisers use -- are not changed by them beyond
## the rounding of that one addition.
Mp <- sum(fitbeta$kept) - logdetSdata$rank
out <- as.vector(out - .5 * (logdetSdata$ldev + Mp * log(2 * pi)))
if (!is.finite(out))
  return(structure(1e20, inner.convergence=1L,
                   inner.report="Non-finite profile REML objective."))
attr(fitbeta$par, "dropped") <- !fitbeta$kept
attr(out, "beta") <- fitbeta$par
attr(out, "gradient") <- fitbeta$gradient
attr(out, "Hessian") <- fitbeta$Hessian
attr(out, "inner.convergence") <- fitbeta$convergence
attr(out, "inner.report") <- fitbeta$report
## a coefficient fit that stopped on a saddle point or an unbounded likelihood (e.g. an
## EGPD shape below -1); reported by the outer optimiser if its search ends here
stall <- .inner_saddle(fitbeta, .search.pen, likdata, likfns,
                       gradtol=likdata$control$inner$gradtol)
if (!is.null(stall))
  attr(out, "inner.stall") <- sprintf("the coefficient search stopped (%s) with max |gradient| %.3g; %s",
                                      fitbeta$report, max(abs(fitbeta$gradient)), stall)
return(out)
}

.reml1 <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL) {
if (!length(pars)) return(numeric(0))
if (is.null(beta)) {
  beta <- attr(pars, "beta")
} else {
  attr(pars, "beta") <- beta
}
spSl <- Map("*", attr(Sdata, "Sl"), exp(pars))
likdata$S <- .makeS(Sdata, spSl = spSl)
if (is.null(H))
  H <- .Hdata(.hess.pen(beta, likdata, likfns))
d1beta <- .d1beta(pars, beta, spSl, H)
dS <- .logdetS(Sdata, pars, deriv=1)
dH <- .d1logdetH(d1beta, likdata, likfns, spSl, H)
dbSb <- sapply(spSl, function(x) crossprod(beta, x %*% beta))
d1 <- .5 * dbSb
d1 <- d1 - .5 * dS$d1
d1 <- d1 + .5 * dH$d1
d1
}

## Differences of .reml0 need inner solutions far more accurate than the objective
## evaluations of the outer line search, so tighten the inner gradient tolerance for
## them only -- unless the user chose an inner gradtol, which is respected.
.tighten_inner <- function(likdata, tol) {
if (!isTRUE(likdata$gradtol.user) && !is.null(likdata$control$inner$gradtol))
  likdata$control$inner$gradtol <- min(likdata$control$inner$gradtol, tol)
likdata
}

.reml1.fd <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL, kept=NULL) {
if (is.null(beta)) beta <- attr(pars, "beta")
likdata <- .tighten_inner(likdata, 1e-8)
tol <- .Machine$double.eps^(1/4)
eps <- pmax(tol * abs(pars), tol)
objective <- function(rho) {
  value <- .reml0(rho, likfns, likdata, Sdata, beta=beta)
  if (!isTRUE(attr(value, "inner.convergence") == 0L))
    stop("Inner optimization failed while differentiating profile REML.",
         call.=FALSE)
  as.numeric(value)
}
f0 <- objective(pars)
f1 <- 0 * pars
for (i in seq_along(pars)) {
  parsi <- pars
  parsi[i] <- parsi[i] + eps[i]
  f1[i] <- objective(parsi)
}
(f1 - f0) / eps
}

.reml12 <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL) {
if (is.null(beta)) beta <- attr(pars, "beta")
likdata <- .tighten_inner(likdata, 1e-8)
gradient <- function(rho) {
  value <- .reml0(rho, likfns, likdata, Sdata, beta=beta)
  bhat <- attr(value, "beta")
  if (is.null(bhat) || !isTRUE(attr(value, "inner.convergence") == 0L))
    stop("Inner optimization failed while differentiating profile REML.",
         call.=FALSE)
  .reml1(rho, likfns, likdata, Sdata, beta=bhat)
}
d1 <- gradient(pars)
tol <- .Machine$double.eps^(1/4)
eps <- pmax(tol * abs(pars), tol)
d2 <- matrix(NA, length(pars), length(pars))
for (i in seq_along(pars)) {
  plus <- replace(pars, i, pars[i] + eps[i])
  minus <- replace(pars, i, pars[i] - eps[i])
  d2[, i] <- (gradient(plus) - gradient(minus)) / (2 * eps[i])
}
d2 <- .5 * (d2 + t(d2))
list(d1 = d1, d2 = d2)
}

.reml2.fd <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL, kept=NULL) {
.reml12(pars, likfns, likdata, Sdata, beta=beta)[[2]]
}

## Hessian of the profile REML criterion by second differences of .reml0, for the
## families without analytic third derivatives. .reml0 is only smooth to about 1e-8
## (its 0.5 log|H| term is first order in the error of the inner solution), so the
## old steps of ~2e-4 put that noise into denominators of ~1e-8 and could return a
## curvature of the wrong sign (-0.196 against 2.20 on one ZIDEGPD fit). Second
## differences D(h), D(h/2) on the log smoothing parameter scale combined by
## Richardson extrapolation, R(h) = (4 D(h/2) - D(h)) / 3, cancel the O(h^2)
## truncation error while keeping the noise contribution near 1e-5; the inner fits use
## a tight gradient tolerance. A third step gives R(h/2), and R(h/2) - R(h) estimates
## the error of the returned R(h) entry by entry (attribute "error"): its O(h^4)
## truncation error plus noise, or a jump where the inner fit is degenerate.
.reml2.fdfd <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL, kept=NULL, h=0.2) {
if (is.null(beta)) beta <- attr(pars, "beta")
likdata <- .tighten_inner(likdata, 1e-10)
objective <- function(rho) {
  value <- .reml0(rho, likfns, likdata, Sdata, beta=beta)
  if (is.null(attr(value, "beta")) ||
      !isTRUE(attr(value, "inner.convergence") == 0L))
    stop("Inner optimization failed while differentiating profile REML.",
         call.=FALSE)
  as.numeric(value)
}
np <- length(pars)
f0 <- objective(pars)
second <- function(hh) {
  out <- matrix(0, np, np)
  for (i in seq_len(np)) for (j in seq_len(i)) {
    ei <- replace(numeric(np), i, hh)
    ej <- replace(numeric(np), j, hh)
    if (i == j) {
      out[i, i] <- (objective(pars + ei) - 2 * f0 + objective(pars - ei)) / hh^2
    } else {
      out[i, j] <- out[j, i] <-
        (objective(pars + ei + ej) - objective(pars + ei - ej) -
         objective(pars - ei + ej) + objective(pars - ei - ej)) / (4 * hh^2)
    }
  }
  out
}
d1 <- second(h)
d2 <- second(h / 2)
d3 <- second(h / 4)
out <- (4 * d2 - d1) / 3
attr(out, "error") <- (4 * d3 - d2) / 3 - out
out
}

## Covariance of the log smoothing parameters from the REML Hessian. The Hessian is
## symmetrised and each eigen-direction is judged on its own: its curvature counts
## as measured if it exceeds both twice the estimated error of the Hessian in that
## direction, |u' E u| (attribute "error" of .reml2.fdfd), and the noise floor
## `noise`; without an error estimate (finite differences of the analytic gradient)
## the floor is `tol`. Measured positive curvature is inverted, however small.
## Directions that are not measured get zero variance, as mgcv gives non-positive
## curvature: silently when their error is at noise level (a flat criterion, e.g. a
## smoothing parameter on its way to infinity), with a warning when the error is
## larger than `tol` (the criterion is not smooth there). Measured negative
## curvature means the outer optimisation did not stop at a minimum and is warned
## about. A single error bound for the whole matrix used to be dominated by the most
## curved direction and discarded accurately measured small curvatures elsewhere.
.reml_vcov <- function(H, tol=1e-4, noise=1e-5) {
E <- attr(H, "error")
H <- .5 * (H + t(H))
if (any(!is.finite(H)))
  stop("non-finite REML Hessian", call.=FALSE)
ev <- eigen(H, symmetric=TRUE)
d <- ev$values
U <- ev$vectors
if (is.matrix(E) && all(dim(E) == dim(H)) && all(is.finite(E))) {
  E <- .5 * (E + t(E))
  err <- abs(colSums(U * (E %*% U)))
  floor <- noise
} else {
  err <- rep(0, length(d))
  floor <- tol
}
measured <- abs(d) > pmax(2 * err, floor)
if (any(measured & d < 0)) {
  warning(sprintf(paste0("the REML Hessian with respect to the log smoothing parameters ",
                         "is not positive definite (smallest eigenvalue %.3g): the outer ",
                         "optimisation did not stop at a minimum, or the criterion is not ",
                         "smooth there; the smoothing-parameter correction to Vc ignores ",
                         "those directions"),
                  min(d)), call.=FALSE)
}
rough <- !measured & err > tol
if (any(rough)) {
  warning(sprintf(paste0("the REML criterion is not smooth enough near the estimated ",
                         "smoothing parameters to resolve their curvature in %d ",
                         "direction(s) (curvature %s, estimated error %s); the ",
                         "smoothing-parameter correction to Vc ignores them"),
                  sum(rough), paste(signif(d[rough], 3), collapse = ", "),
                  paste(signif(err[rough], 3), collapse = ", ")), call.=FALSE)
}
pos <- measured & d > 0
V <- matrix(0, nrow(H), ncol(H))
if (any(pos)) {
  Up <- U[, pos, drop=FALSE]
  V <- Up %*% (t(Up) / d[pos])
}
attr(V, "eigen") <- d
attr(V, "error") <- err
V
}

.search.reml <- function(pars, likfns, likdata, Sdata, H=NULL, kept=NULL) {
gH <- .reml12(pars, likfns, likdata, Sdata, H=H)
.search.dir(gH[[1]], gH[[2]], !logical(length(gH[[1]])))
}

## fixed smoothing parameters

.reml0_fixed <- function(pars, likfns, likdata, Sdata, beta=NULL) {
if (is.null(beta))
  beta <- attr(pars, "beta")
sp <- exp(pars)
likdata$S <- .makeS(Sdata, sp)
fitbeta <- .fit_inner(beta, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns)
fitbeta$beta <- fitbeta$par
fitbeta$par <- pars
fitbeta
}
