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
  mp <- (u * exp(-tt)) / (a * -expm1(-tt)); mpp <- mp * (1 - tt - mp)
  sat <- !is.finite(mp) | (tt > 700); mp[sat] <- 0; mpp[sat] <- 0
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

.iG1_degpd <- function(v, kappa) v^(1/kappa)

.degpd1fns <- list(d0=.degpd1.d0, d120=.degpd1.d12, d340=NULL, m=1, iG=.iG1_degpd)

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
.degpd1idfns <- list(d0=.degpd1id.d0, d120=.degpd1id.d12, d340=NULL, m=1, iG=.iG1_degpd, identity.xi=TRUE)

## model 2 ##

.G2_degpd <- function(v, kappa1, kappa2, p) {
  p * v^kappa1 + (1 - p) * v^kappa2
}

.iG2i_degpd <- function(v, kappa1, kappa2, p) {
  if (v <= 0)
    return(0)
  if (v >= 1)
    return(1)
  if (abs(kappa1 - kappa2) < sqrt(.Machine$double.eps))
    return(v^(1 / kappa1))
  vv <- range(c(v^(1 / kappa1), v^(1 / kappa2)))
  d <- max(diff(vv), .Machine$double.eps)
  lo <- vv[1]
  while(.G2_degpd(lo, kappa1, kappa2, p) - v > 0 && lo > 0) lo <- max(0, lo - d)
  hi <- vv[2]
  while(.G2_degpd(hi, kappa1, kappa2, p) - v < 0 && hi < 1) hi <- min(1, hi + d)
  uniroot(function(x) .G2_degpd(x, kappa1, kappa2, p) - v, c(lo, hi))$root
}

.iG2_degpd <- function(v, kappa1, kappa2, p) {
  n <- length(v)
  vapply(seq_len(n), function(i) .iG2i_degpd(v[i], kappa1[i], kappa2[i], p[i]), double(1))
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

.iG3_degpd <- function(v, delta) 1 - qbeta(1 - v, 1/delta, 2)^(1/delta)

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

.iG4_degpd <- function(v, delta, kappa) 1 - qbeta(1 - v^(2/kappa), 1/delta, 2)^(1/delta)

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
