## Generalised Poisson-inverse Gaussian (GPIG; Zhu & Joe 2009) and zero-inflated
## GPIG (ZIGPIG) likelihood bundles for the native GAM fitter. These mirror the
## PIG wrappers in R/pig_lik.R: each R wrapper splits the flat coefficient vector
## by parameter and calls the corresponding Rcpp routine in src/gpig.cpp.
##
## GPIG is a three-parameter heavy-tailed count family (a in (0,1] tail exponent,
## b > 0 level, c in [0,1] down-weighting). Two parameterisations are offered:
##   "gpig"    mean-based  (mu = mean, a, c); b implied by b = mu(1-c)^(1-a)/(ac)
##   "gpignat" native      (a, b, c) exactly as in the paper
## and likewise "zigpig" / "zigpignat" add a logit-link zero-inflation pi.
##
## The score is exact (paper eqs 3.4-3.5, with the two misprinted p0 base cases
## corrected); the Hessian is the empirical-Fisher / BHHH approximation -- see
## the header of src/gpig.cpp.

## ---- inverse-G (quantile) helpers, used for prediction/simulation ----
## qgpig() (native a,b,c) is defined in R/gpig_dist.R.
.iG_gpig    <- function(v, mu, a, c) {              # mean parameterisation
  b <- mu * (1 - c)^(1 - a) / (a * c)
  qgpig(v, a = a, b = b, c = c)
}
.iG_gpignat <- function(v, a, b, c) qgpig(v, a = a, b = b, c = c)

## ---- GPIG, mean parameterisation (3 par: log mu, logit a, logit c) ----

.gpig1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for GPIG.")
  gpig1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
          likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[, 1],
          likdata$dupid, likdata$duplicate, likdata$offsets)
}
.gpig1.d12 <- function(pars, likdata) {
  gpig1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
           likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[, 1],
           likdata$dupid, likdata$duplicate, likdata$offsets)
}
.gpig1fns <- list(d0 = .gpig1.d0, d120 = .gpig1.d12, d340 = NULL, m = 1, iG = .iG_gpig)

## ---- GPIG, native parameterisation (3 par: logit a, log b, logit c) ----

.gpignat1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for GPIG.")
  gpignat1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
             likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[, 1],
             likdata$dupid, likdata$duplicate, likdata$offsets)
}
.gpignat1.d12 <- function(pars, likdata) {
  gpignat1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
              likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[, 1],
              likdata$dupid, likdata$duplicate, likdata$offsets)
}
.gpignat1fns <- list(d0 = .gpignat1.d0, d120 = .gpignat1.d12, d340 = NULL, m = 1,
                     iG = .iG_gpignat)

## ---- ZIGPIG, mean parameterisation (4 par: log mu, logit a, logit c, logit pi) ----

.zigpig1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for ZIGPIG.")
  zigpig1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
            likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]],
            likdata$y[, 1], likdata$dupid, likdata$duplicate, likdata$offsets)
}
.zigpig1.d12 <- function(pars, likdata) {
  zigpig1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
             likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]],
             likdata$y[, 1], likdata$dupid, likdata$duplicate, likdata$offsets)
}
.zigpig1fns <- list(d0 = .zigpig1.d0, d120 = .zigpig1.d12, d340 = NULL, m = 1,
                    iG = .iG_gpig)

## ---- ZIGPIG, native parameterisation (4 par: logit a, log b, logit c, logit pi) ----

.zigpignat1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for ZIGPIG.")
  zigpignat1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
               likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]],
               likdata$y[, 1], likdata$dupid, likdata$duplicate, likdata$offsets)
}
.zigpignat1.d12 <- function(pars, likdata) {
  zigpignat1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
                likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]],
                likdata$y[, 1], likdata$dupid, likdata$duplicate, likdata$offsets)
}
.zigpignat1fns <- list(d0 = .zigpignat1.d0, d120 = .zigpignat1.d12, d340 = NULL,
                       m = 1, iG = .iG_gpignat)
