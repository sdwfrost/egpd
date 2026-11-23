## Discrete EGPD negative log-likelihood functions

## issue #2: bounded shape link for DEGPD models 2-6. The C++ d12 is fed
## lxi_eff = log(xi_bounded) via its forward map, so it returns per-observation
## derivatives w.r.t. lxi_eff. This converts the shape (parameter 2) gradient and
## Hessian columns to the raw predictor eta with the closed-form chain rule
## (m'(eta), m''(eta)); no-op when xi.max = Inf. Column packing matches the C++
## d12 layout: gradient cols 1..npar, then upper-triangular Hessian row-major.
## (Model 1 applies this same chain rule inside C++.)
.bounded_xi_chain <- function(out, pars, likdata) {
  a <- likdata$xi.max
  if (is.null(a) || !is.finite(a)) return(out)
  npar <- length(likdata$X)
  xcol <- 2L                                   # shape is parameter 2 (1-based)
  beta <- split(pars, factor(likdata$idpars, levels = seq_along(likdata$X)))[[xcol]]
  eta <- as.numeric(likdata$X[[xcol]] %*% beta)
  if (isTRUE(likdata$duplicate == 1)) eta <- eta[likdata$dupid + 1]
  if (length(likdata$offsets[[xcol]]) > 0) eta <- eta + likdata$offsets[[xcol]]
  u <- exp(eta); tt <- u / a
  small <- !is.na(tt) & (tt < 1e-4)
  mp <- ifelse(small, 1 - 0.5 * tt + tt^2 / 12, (u * exp(-tt)) / (a * -expm1(-tt)))
  mpp <- mp * (1 - tt - mp)
  sat <- is.na(mp) | !is.finite(mp) | (tt > 700); mp[sat] <- 0; mpp[sat] <- 0
  ## packed upper-tri column (1-based) for parameter pair (i,j), i<=j, 0-based
  col_ij <- function(i, j) {
    off <- 0L; if (i > 0) for (r in 0:(i - 1)) off <- off + (npar - r)
    npar + off + (j - i) + 1L
  }
  xi0 <- 1L                                    # 0-based index of the shape parameter
  g <- out[, xcol]; cdd <- col_ij(xi0, xi0); hdd <- out[, cdd]
  out[, xcol] <- g * mp
  for (k in 0:(npar - 1)) if (k != xi0) {
    cc <- col_ij(min(k, xi0), max(k, xi0)); out[, cc] <- out[, cc] * mp
  }
  out[, cdd] <- hdd * mp^2 + g * mpp
  out
}

## model 1 ##

.degpd1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  xm <- if (is.null(likdata$xi.max)) Inf else likdata$xi.max
  degpd1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, xm)
}

.degpd1.d12 <- function(pars, likdata) {
  xm <- if (is.null(likdata$xi.max)) Inf else likdata$xi.max
  degpd1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, xm)
}

## The generated high-order expressions lose precision near xi = 0, and in cells too
## narrow for the generated forms (below). Signal these numerical limitations so
## .outer() can retry with finite differences.
.degpd1.d34.unavailable <- function(message) {
  stop(structure(list(message = message, call = NULL),
                 class = c("egpd_d34_unavailable", "error", "condition")))
}

## Narrow cells: the generated d34 evaluate t = 1 + xi y/sigma at both ends of the cell
## [y, y + 1], whose width |xi|/sigma is then resolved only to eps t sigma/|xi| =
## eps (sigma + xi y)/|xi| (eps * (distance to the endpoint) for xi < 0). The third
## derivatives lose about that much, times up to 3e4 (both links), and the REML
## gradient error grows with it: below eps (sigma + xi y)/|xi| = 1e-5 it stayed below
## 2e-5 in designed fits (under the finite-difference gradient's own error, 6e-5 to
## 3e-4); above, it reached 1.5e-4 at 1e10 counts and 5e-2 (sp 5% off) at 1e13.
## Ordinary counts and scales are many orders of magnitude below the limit.
.degpd1.d34.cell.limit <- 1e-5

.degpd1.check.d34 <- function(pars, likdata, identity = FALSE) {
  b <- split(pars, factor(likdata$idpars, levels = seq_along(likdata$X)))
  predictor <- function(k) {
    eta <- as.numeric(likdata$X[[k]] %*% b[[k]])
    if (isTRUE(likdata$duplicate == 1)) eta <- eta[likdata$dupid + 1L]
    if (length(likdata$offsets[[k]])) eta <- eta + likdata$offsets[[k]]
    eta
  }
  eta <- predictor(2)
  small <- if (identity) abs(eta) < 0.05 else eta < log(1e-3)
  if (any(!is.finite(eta) | small))
    .degpd1.d34.unavailable("analytic DEGPD d34 is unstable near zero shape")
  xi <- if (identity) eta else exp(eta)
  cells <- (exp(predictor(1)) + xi * likdata$y[, 1]) / abs(xi)
  if (any(.Machine$double.eps * cells > .degpd1.d34.cell.limit, na.rm = TRUE))
    .degpd1.d34.unavailable(sprintf(paste(
      "analytic DEGPD (%s shape link) third/fourth derivatives lose accuracy at very",
      "large counts or scales: (sigma + xi*y)/|xi| reaches %.3g, above %.3g, so the",
      "cells [y, y + 1] are resolved to worse than %.0e"),
      if (identity) "identity" else "log", max(cells, na.rm = TRUE),
      .degpd1.d34.cell.limit / .Machine$double.eps, .degpd1.d34.cell.limit))
}

## These derivatives assume an unbounded log link, not the xi.max link.
.degpd1.d34 <- function(pars, likdata) {
  xm <- if (is.null(likdata$xi.max)) Inf else likdata$xi.max
  if (is.finite(xm))
    .degpd1.d34.unavailable("analytic DEGPD d34 is unavailable for finite xi.max")
  .degpd1.check.d34(pars, likdata)
  out <- degpd1d34(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, xm)
  if (any(!is.finite(out)))
    .degpd1.d34.unavailable("analytic DEGPD d34 is non-finite")
  out
}

.iG1_degpd <- function(v, kappa) v^(1/kappa)

.degpd1fns <- list(d0=.degpd1.d0, d120=.degpd1.d12, d340=.degpd1.d34, m=1, iG=.iG1_degpd)

## model 1, identity link on xi (unconstrained shape; xi may be < 0, bounded tail).
## Parameter 2 IS xi (not log xi), so NO .bounded_xi_chain conversion is applied.
.degpd1id.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  degpd1id0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, Inf)
}
.degpd1id.d12 <- function(pars, likdata) {
  degpd1id12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, Inf)
}
.degpd1id.d34 <- function(pars, likdata) {
  .degpd1.check.d34(pars, likdata, identity = TRUE)
  out <- degpd1id34(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, Inf)
  if (any(!is.finite(out)))
    .degpd1.d34.unavailable("analytic DEGPD d34 is non-finite")
  out
}
.degpd1idfns <- list(d0=.degpd1id.d0, d120=.degpd1id.d12, d340=.degpd1id.d34, m=1, iG=.iG1_degpd, identity.xi=TRUE)

## model 2 ##

.G2_degpd <- function(v, kappa1, kappa2, p) {
  p * v^kappa1 + (1 - p) * v^kappa2
}

.iG2i_degpd <- function(v, kappa1, kappa2, p) {
  q.G(v, type = 6, kappa = kappa1, delta = kappa2, prob = p)
}

.iG2_degpd <- function(v, kappa1, kappa2, p) {
  mapply(.iG2i_degpd, v, kappa1, kappa2, p, USE.NAMES = FALSE)
}

.degpd2.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  xm <- if (is.null(likdata$xi.max)) Inf else likdata$xi.max
  degpd2d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, xm)
}

.degpd2.d12 <- function(pars, likdata) {
  xm <- if (is.null(likdata$xi.max)) Inf else likdata$xi.max
  out <- degpd2d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, xm)
  .bounded_xi_chain(out, pars, likdata)
}

.degpd2fns <- list(d0=.degpd2.d0, d120=.degpd2.d12, d340=NULL, m=2, iG=.iG2_degpd)

## model 3 ##

.degpd3.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  xm <- if (is.null(likdata$xi.max)) Inf else likdata$xi.max
  degpd3d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, xm)
}

.degpd3.d12 <- function(pars, likdata) {
  xm <- if (is.null(likdata$xi.max)) Inf else likdata$xi.max
  out <- degpd3d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, xm)
  .bounded_xi_chain(out, pars, likdata)
}

.iG3_degpd <- function(v, delta) q.G(v, type = 4L, delta = delta)

.degpd3fns <- list(d0=.degpd3.d0, d120=.degpd3.d12, d340=NULL, m=3, iG=.iG3_degpd)

## model 4 ##

.degpd4.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  xm <- if (is.null(likdata$xi.max)) Inf else likdata$xi.max
  degpd4d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, xm)
}

.degpd4.d12 <- function(pars, likdata) {
  xm <- if (is.null(likdata$xi.max)) Inf else likdata$xi.max
  out <- degpd4d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, xm)
  .bounded_xi_chain(out, pars, likdata)
}

.iG4_degpd <- function(v, delta, kappa) q.G(v, type = 5L, delta = delta, kappa = kappa)

.degpd4fns <- list(d0=.degpd4.d0, d120=.degpd4.d12, d340=NULL, m=4, iG=.iG4_degpd)

## model 5 — truncated normal ##

.degpd5.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  xm <- if (is.null(likdata$xi.max)) Inf else likdata$xi.max
  degpd5d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, xm)
}

.degpd5.d12 <- function(pars, likdata) {
  xm <- if (is.null(likdata$xi.max)) Inf else likdata$xi.max
  out <- degpd5d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, xm)
  .bounded_xi_chain(out, pars, likdata)
}

.iG5_degpd <- function(v, kappa) q.G(v, type = 2, kappa = kappa)

.degpd5fns <- list(d0=.degpd5.d0, d120=.degpd5.d12, d340=NULL, m=5, iG=.iG5_degpd)

## model 6 — truncated beta ##

.degpd6.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  xm <- if (is.null(likdata$xi.max)) Inf else likdata$xi.max
  degpd6d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, xm)
}

.degpd6.d12 <- function(pars, likdata) {
  xm <- if (is.null(likdata$xi.max)) Inf else likdata$xi.max
  out <- degpd6d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets, xm)
  .bounded_xi_chain(out, pars, likdata)
}

.iG6_degpd <- function(v, kappa) q.G(v, type = 3, kappa = kappa)

.degpd6fns <- list(d0=.degpd6.d0, d120=.degpd6.d12, d340=NULL, m=6, iG=.iG6_degpd)

## ---------------------------------------------------------------------------
## IDENTITY LINK ON xi, CARRIERS 2-6.
##
## xi_identity = TRUE dispatches to direct predictor derivatives in
## src/degpd1_identity.cpp. The former inverse log-xi chain divided cancellation
## errors by xi^2 and could not represent the true xi=0 likelihood/derivatives.
## All carriers now share log1p(z)/z with an analytic Taylor continuation.
## Carriers 5/6 integrate their density over the survival interval, checking
## quadrature error in both the objective and its first/second derivatives.

## Build the identity-link wrappers mechanically: same C++ entry points, xi.max
## forced to Inf (a bounding link on xi is meaningless when xi may be negative),
## xi_identity = TRUE. Native results already use the identity predictors.
.make_degpd_id_fns <- function(m, d0fun, d12fun, nX, iG) {
  cargs <- function(likdata, pars)
    c(list(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X)))),
      likdata$X[seq_len(nX)],
      list(likdata$y[, 1], likdata$dupid, likdata$duplicate, likdata$offsets,
           Inf, TRUE))
  d0 <- function(pars, likdata) {
    if (likdata$censored)
      stop("Censored likelihoods not currently available for extended GPDs.")
    do.call(d0fun, cargs(likdata, pars))
  }
  d12 <- function(pars, likdata)
    do.call(d12fun, cargs(likdata, pars))
  list(d0 = d0, d120 = d12, d340 = NULL, m = m, iG = iG, identity.xi = TRUE)
}

.degpd2idfns <- .make_degpd_id_fns(2, degpd2d0, degpd2d12, 5, .iG2_degpd)
.degpd3idfns <- .make_degpd_id_fns(3, degpd3d0, degpd3d12, 3, .iG3_degpd)
.degpd4idfns <- .make_degpd_id_fns(4, degpd4d0, degpd4d12, 4, .iG4_degpd)
.degpd5idfns <- .make_degpd_id_fns(5, degpd5d0, degpd5d12, 3, .iG5_degpd)
.degpd6idfns <- .make_degpd_id_fns(6, degpd6d0, degpd6d12, 3, .iG6_degpd)
