## Bell (Castellares, Ferrari & Lemonte 2018) and zero-inflated Bell (ZIBell)
## likelihood bundles for the native GAM fitter. These mirror the PIG/GPIG
## wrappers: each R wrapper splits the flat coefficient vector by parameter and
## calls the corresponding Rcpp routine in src/bell.cpp.
##
## Bell is a one-parameter overdispersed count family. Two parameterisations are
## offered:
##   "bell"    mean-based  (mu = mean, log link); internally theta = W0(mu).
##   "bellnat" native      (theta > 0, log link) exactly as in the paper.
## and likewise "zibell" / "zibellnat" add a logit-link zero-inflation pi.
##
## Because the Bell log-likelihood is closed form, both the score and the
## observed-information Hessian returned by d12 are exact (see src/bell.cpp).

## ---- inverse-G (quantile) helpers, used for prediction/simulation ----
## qbell() (native theta) is defined in R/bell_dist.R; .bell_theta_from_mean()
## (theta = W0(mu)) in R/bell_gamlss.R.
.iG_bell    <- function(v, mu)    qbell(v, theta = .bell_theta_from_mean(mu))  # mean
.iG_bellnat <- function(v, theta) qbell(v, theta = theta)                     # native

## ---- Bell, mean parameterisation (1 par: log mu) ----

.bell1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for Bell.")
  bell1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
          likdata$X[[1]], likdata$y[, 1],
          likdata$dupid, likdata$duplicate, likdata$offsets)
}
.bell1.d12 <- function(pars, likdata) {
  bell1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
           likdata$X[[1]], likdata$y[, 1],
           likdata$dupid, likdata$duplicate, likdata$offsets)
}
.bell1fns <- list(d0 = .bell1.d0, d120 = .bell1.d12, d340 = NULL, m = 1, iG = .iG_bell)

## ---- Bell, native parameterisation (1 par: log theta) ----

.bellnat1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for Bell.")
  bellnat1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
             likdata$X[[1]], likdata$y[, 1],
             likdata$dupid, likdata$duplicate, likdata$offsets)
}
.bellnat1.d12 <- function(pars, likdata) {
  bellnat1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
              likdata$X[[1]], likdata$y[, 1],
              likdata$dupid, likdata$duplicate, likdata$offsets)
}
.bellnat1fns <- list(d0 = .bellnat1.d0, d120 = .bellnat1.d12, d340 = NULL, m = 1,
                     iG = .iG_bellnat)

## ---- ZIBell, mean parameterisation (2 par: log mu, logit pi) ----

.zibell1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for ZIBell.")
  zibell1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
            likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
            likdata$dupid, likdata$duplicate, likdata$offsets)
}
.zibell1.d12 <- function(pars, likdata) {
  zibell1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
             likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
             likdata$dupid, likdata$duplicate, likdata$offsets)
}
.zibell1fns <- list(d0 = .zibell1.d0, d120 = .zibell1.d12, d340 = NULL, m = 1,
                    iG = .iG_bell)

## ---- ZIBell, native parameterisation (2 par: log theta, logit pi) ----

.zibellnat1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for ZIBell.")
  zibellnat1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
               likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
               likdata$dupid, likdata$duplicate, likdata$offsets)
}
.zibellnat1.d12 <- function(pars, likdata) {
  zibellnat1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
                likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
                likdata$dupid, likdata$duplicate, likdata$offsets)
}
.zibellnat1fns <- list(d0 = .zibellnat1.d0, d120 = .zibellnat1.d12, d340 = NULL,
                       m = 1, iG = .iG_bellnat)
