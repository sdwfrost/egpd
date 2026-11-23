## Continuous EGPD negative log-likelihood functions
## Adapted from evgam by Ben Youngman
## Native continuous objectives use the exponential limit at xi=0. Where the
## generated closed-form derivatives are unreliable (xi near 0, y/sigma near 0,
## delta * A near 0 for carriers 3/4) or overflow (e.g. next to the xi < 0
## endpoint), d12 uses exact forward-mode derivatives of that same objective,
## and d34 (where available) uses predictor-scale stencils kept inside the support.
## Carrier-3 d34 rows next to the endpoint or far in the tail instead combine the
## carrier-1 forms (the GPD part) with exact derivatives of the carrier's delta term.

## model 1 ##

.egpd1.d0 <- function(pars, likdata) {
if (likdata$censored)
  stop("Censored likelihoods not currently available for extended GPDs.")
egpd1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.egpd1.d12 <- function(pars, likdata) {
egpd1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.egpd1.d34 <- function(pars, likdata) {
egpd1d34(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.iG1_egpd <- function(v, kappa) v^(1/kappa)

.egpd1fns <- list(d0=.egpd1.d0, d120=.egpd1.d12, d340=.egpd1.d34, m=1, iG=.iG1_egpd)

## model 2 ##

.G2_egpd <- function(v, kappa1, kappa2, p) {
p * v^kappa1 + (1 - p) * v^kappa2
}

.iG2i_egpd <- function(v, kappa1, kappa2, p) {
 q.G(v, type = 6, kappa = kappa1, delta = kappa2, prob = p)
}

.iG2_egpd <- function(v, kappa1, kappa2, p) {
mapply(.iG2i_egpd, v, kappa1, kappa2, p, USE.NAMES = FALSE)
}

.egpd2.d0 <- function(pars, likdata) {
if (likdata$censored)
  stop("Censored likelihoods not currently available for extended GPDs.")
egpd2d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.egpd2.d12 <- function(pars, likdata) {
egpd2d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.egpd2.d34 <- function(pars, likdata) {
egpd2d34(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.egpd2fns <- list(d0=.egpd2.d0, d120=.egpd2.d12, d340=NULL, m=2, iG=.iG2_egpd)

## model 3 ##

.egpd3.d0 <- function(pars, likdata) {
if (likdata$censored)
  stop("Censored likelihoods not currently available for extended GPDs.")
egpd3d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.egpd3.d12 <- function(pars, likdata) {
egpd3d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.egpd3.d34 <- function(pars, likdata) {
egpd3d34(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.iG3_egpd <- function(v, delta) q.G(v, type = 4L, delta = delta)

.egpd3fns <- list(d0=.egpd3.d0, d120=.egpd3.d12, d340=.egpd3.d34, m=3, iG=.iG3_egpd)

## model 4 ##

.egpd4.d0 <- function(pars, likdata) {
if (likdata$censored)
  stop("Censored likelihoods not currently available for extended GPDs.")
egpd4d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.egpd4.d12 <- function(pars, likdata) {
egpd4d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.egpd4.d34 <- function(pars, likdata) {
egpd4d34(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.iG4_egpd <- function(v, delta, kappa) q.G(v, type = 5L, delta = delta, kappa = kappa)

.egpd4fns <- list(d0=.egpd4.d0, d120=.egpd4.d12, d340=NULL, m=4, iG=.iG4_egpd)

## model 5 — truncated normal ##

.egpd5.d0 <- function(pars, likdata) {
if (likdata$censored)
  stop("Censored likelihoods not currently available for extended GPDs.")
egpd5d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.egpd5.d12 <- function(pars, likdata) {
egpd5d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.iG5_egpd <- function(v, kappa) q.G(v, type = 2, kappa = kappa)

.egpd5fns <- list(d0=.egpd5.d0, d120=.egpd5.d12, d340=NULL, m=5, iG=.iG5_egpd)

## model 6 — truncated beta ##

.egpd6.d0 <- function(pars, likdata) {
if (likdata$censored)
  stop("Censored likelihoods not currently available for extended GPDs.")
egpd6d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.egpd6.d12 <- function(pars, likdata) {
egpd6d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.iG6_egpd <- function(v, kappa) q.G(v, type = 3, kappa = kappa)

.egpd6fns <- list(d0=.egpd6.d0, d120=.egpd6.d12, d340=NULL, m=6, iG=.iG6_egpd)
