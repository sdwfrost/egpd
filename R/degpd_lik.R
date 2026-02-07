## Discrete EGPD negative log-likelihood functions

## model 1 ##

.degpd1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  degpd1d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.degpd1.d12 <- function(pars, likdata) {
  degpd1d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.iG1_degpd <- function(v, kappa) v^(1/kappa)

.degpd1fns <- list(d0=.degpd1.d0, d120=.degpd1.d12, d340=NULL, m=1, iG=.iG1_degpd)

## model 2 ##

.G2_degpd <- function(v, kappa1, kappa2, p) {
  p * v^kappa1 + (1 - p) * v^kappa2
}

.iG2i_degpd <- function(v, kappa1, kappa2, p) {
  vv <- range(c(v^1/kappa1, v^1/kappa2))
  d <- diff(vv)
  lo <- vv[1]
  while(.G2_degpd(lo, kappa1, kappa2, p) - v > 0) lo <- max(0, lo - d)
  hi <- vv[2]
  while(.G2_degpd(hi, kappa1, kappa2, p) - v > 0) hi <- min(1, hi + d)
  uniroot(function(x) .G2_degpd(x, kappa1, kappa2, p) - v, c(lo, hi))$root
}

.iG2_degpd <- function(v, kappa1, kappa2, p) {
  n <- length(v)
  vapply(seq_len(n), function(i) .iG2i_degpd(v[i], kappa1[i], kappa2[i], p[i]), double(1))
}

.degpd2.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  degpd2d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.degpd2.d12 <- function(pars, likdata) {
  degpd2d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.degpd2fns <- list(d0=.degpd2.d0, d120=.degpd2.d12, d340=NULL, m=2, iG=.iG2_degpd)

## model 3 ##

.degpd3.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  degpd3d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.degpd3.d12 <- function(pars, likdata) {
  degpd3d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.iG3_degpd <- function(v, delta) 1 - qbeta(1 - v, 1/delta, 2)^(1/delta)

.degpd3fns <- list(d0=.degpd3.d0, d120=.degpd3.d12, d340=NULL, m=3, iG=.iG3_degpd)

## model 4 ##

.degpd4.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  degpd4d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.degpd4.d12 <- function(pars, likdata) {
  degpd4d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.iG4_degpd <- function(v, delta, kappa) 1 - qbeta(1 - v^(2/kappa), 1/delta, 2)^(1/delta)

.degpd4fns <- list(d0=.degpd4.d0, d120=.degpd4.d12, d340=NULL, m=4, iG=.iG4_degpd)
