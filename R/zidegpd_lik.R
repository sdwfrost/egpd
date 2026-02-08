## Zero-Inflated Discrete EGPD negative log-likelihood functions

## model 1 ##

.zidegpd1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  zidegpd1d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.zidegpd1.d12 <- function(pars, likdata) {
  zidegpd1d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.iG1_zidegpd <- function(v, kappa) v^(1/kappa)

.zidegpd1fns <- list(d0=.zidegpd1.d0, d120=.zidegpd1.d12, d340=NULL, m=1, iG=.iG1_zidegpd)

## model 2 ##

.G2_zidegpd <- function(v, kappa1, kappa2, p) {
  p * v^kappa1 + (1 - p) * v^kappa2
}

.iG2i_zidegpd <- function(v, kappa1, kappa2, p) {
  vv <- range(c(v^1/kappa1, v^1/kappa2))
  d <- diff(vv)
  lo <- vv[1]
  while(.G2_zidegpd(lo, kappa1, kappa2, p) - v > 0) lo <- max(0, lo - d)
  hi <- vv[2]
  while(.G2_zidegpd(hi, kappa1, kappa2, p) - v > 0) hi <- min(1, hi + d)
  uniroot(function(x) .G2_zidegpd(x, kappa1, kappa2, p) - v, c(lo, hi))$root
}

.iG2_zidegpd <- function(v, kappa1, kappa2, p) {
  n <- length(v)
  vapply(seq_len(n), function(i) .iG2i_zidegpd(v[i], kappa1[i], kappa2[i], p[i]), double(1))
}

.zidegpd2.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  zidegpd2d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$X[[6]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.zidegpd2.d12 <- function(pars, likdata) {
  zidegpd2d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$X[[6]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.zidegpd2fns <- list(d0=.zidegpd2.d0, d120=.zidegpd2.d12, d340=NULL, m=2, iG=.iG2_zidegpd)

## model 3 ##

.zidegpd3.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  zidegpd3d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.zidegpd3.d12 <- function(pars, likdata) {
  zidegpd3d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.iG3_zidegpd <- function(v, delta) 1 - qbeta(1 - v, 1/delta, 2)^(1/delta)

.zidegpd3fns <- list(d0=.zidegpd3.d0, d120=.zidegpd3.d12, d340=NULL, m=3, iG=.iG3_zidegpd)

## model 4 ##

.zidegpd4.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  zidegpd4d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.zidegpd4.d12 <- function(pars, likdata) {
  zidegpd4d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.iG4_zidegpd <- function(v, delta, kappa) 1 - qbeta(1 - v^(2/kappa), 1/delta, 2)^(1/delta)

.zidegpd4fns <- list(d0=.zidegpd4.d0, d120=.zidegpd4.d12, d340=NULL, m=4, iG=.iG4_zidegpd)

## model 5 — truncated normal ##

.zidegpd5.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  zidegpd5d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.zidegpd5.d12 <- function(pars, likdata) {
  zidegpd5d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.iG5_zidegpd <- function(v, kappa) q.G(v, type = 2, kappa = kappa)

.zidegpd5fns <- list(d0=.zidegpd5.d0, d120=.zidegpd5.d12, d340=NULL, m=5, iG=.iG5_zidegpd)

## model 6 — truncated beta ##

.zidegpd6.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  zidegpd6d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.zidegpd6.d12 <- function(pars, likdata) {
  zidegpd6d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$offsets)
}

.iG6_zidegpd <- function(v, kappa) q.G(v, type = 3, kappa = kappa)

.zidegpd6fns <- list(d0=.zidegpd6.d0, d120=.zidegpd6.d12, d340=NULL, m=6, iG=.iG6_zidegpd)
