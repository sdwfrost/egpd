## Poisson-inverse Gaussian (PIG) and zero-inflated PIG (ZIPIG) likelihood
## bundles for the native GAM fitter. These mirror the DEGPD/ZIDEGPD wrappers
## in R/degpd_lik.R: each R wrapper splits the flat coefficient vector by
## parameter and calls the corresponding Rcpp routine in src/pig.cpp.
##
## PIG is a mixed-Poisson count model (not an EGPD): mu > 0 is the mean and
## sigma > 0 the dispersion, with Var = mu + sigma * mu^2. ZIPIG adds a
## logit-link zero-inflation probability. The score is exact; the Hessian is
## the empirical-Fisher / BHHH approximation (Fisher scoring), matching
## gamlss.dist::PIG -- see the header of src/pig.cpp.

## ---- PIG (2 parameters: log mu, log sigma) ----

.pig1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for PIG.")
  pig1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
         likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
         likdata$dupid, likdata$duplicate, likdata$offsets)
}

.pig1.d12 <- function(pars, likdata) {
  pig1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
          likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
          likdata$dupid, likdata$duplicate, likdata$offsets)
}

## quantile-type inverse transform (mean-scale); used for prediction/simulation.
.iG_pig <- function(v, mu, sigma) gamlss.dist::qPIG(v, mu = mu, sigma = sigma)

.pig1fns <- list(d0 = .pig1.d0, d120 = .pig1.d12, d340 = NULL, m = 1, iG = .iG_pig)

## ---- ZIPIG (3 parameters: log mu, log sigma, logit pi) ----

.zipig1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for ZIPIG.")
  zipig1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
           likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[, 1],
           likdata$dupid, likdata$duplicate, likdata$offsets)
}

.zipig1.d12 <- function(pars, likdata) {
  zipig1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
            likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[, 1],
            likdata$dupid, likdata$duplicate, likdata$offsets)
}

.zipig1fns <- list(d0 = .zipig1.d0, d120 = .zipig1.d12, d340 = NULL, m = 1, iG = .iG_pig)
