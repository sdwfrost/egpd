## Generalized Poisson, generalized Waring and Poisson-lognormal: d/p/q/r functions and
## the likelihood bundles for the native GAM fitter. The log pmfs live in
## src/countfams.cpp; see its header for the parameterisations, links and tail behaviour.
##
## The R names avoid masking stats::plnorm (the lognormal cdf): the Poisson-lognormal
## functions are d/p/q/rpoislnorm, while the egpd family string remains "plnorm".

## ---------------------------------------------------------------- generalized Poisson

#' Generalized Poisson distribution
#'
#' Consul & Jain's generalized Poisson in the mean parameterisation:
#' \code{theta = mu (1 - lambda)}, so \code{E[Y] = mu} and
#' \code{Var[Y] = mu / (1 - lambda)^2}. The tail decays geometrically, so this family
#' captures overdispersion but not a heavy tail.
#'
#' @param x,q vector of counts
#' @param p vector of probabilities
#' @param n number of draws
#' @param mu mean, \code{mu > 0}
#' @param lambda dispersion, \code{lambda} in \code{[0, 1)}; \code{lambda = 0} is Poisson
#' @param log,log.p logical; return the log density / probability
#' @param lower.tail logical; if TRUE probabilities are \code{P[X <= x]}
#' @param max.value largest count searched when inverting the cdf
#' @return a numeric vector
#' @name gpois
NULL

#' @rdname gpois
#' @export
dgpois <- function(x, mu = 1, lambda = 0.5, log = FALSE) {
  stopifnot(mu > 0, lambda >= 0, lambda < 1)
  out <- gpois_logpmf_cpp(as.numeric(x), mu, max(lambda, 1e-10))
  out[x < 0 | x != floor(x)] <- -Inf
  if (log) out else exp(out)
}
#' @rdname gpois
#' @export
pgpois <- function(q, mu = 1, lambda = 0.5, lower.tail = TRUE, log.p = FALSE) {
  hi <- max(0, floor(max(q)))
  cdf <- cumsum(dgpois(0:hi, mu, lambda))
  out <- ifelse(q < 0, 0, pmin(cdf[pmax(floor(q), 0) + 1L], 1))
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}
#' @rdname gpois
#' @export
qgpois <- function(p, mu = 1, lambda = 0.5, lower.tail = TRUE, log.p = FALSE,
                   max.value = 1e6) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  .qsearch(p, function(n) cumsum(dgpois(0:n, mu, lambda)), mu / (1 - lambda), max.value)
}
#' @rdname gpois
#' @export
rgpois <- function(n, mu = 1, lambda = 0.5) qgpois(stats::runif(n), mu, lambda)

## ---------------------------------------------------------------- generalized Waring

#' Generalized Waring distribution
#'
#' The three-parameter generalized Waring distribution in the mean parameterisation, with
#' \code{a = mu (rho - 1) / k}, so \code{E[Y] = mu}. Its tail is regularly varying,
#' \code{P(Y = y) ~ y^-(rho + 1)}: a power law with extreme-value index \code{xi = 1/rho},
#' placing it in the Frechet domain, unlike the negative binomial, PIG, generalized Poisson
#' or Poisson-lognormal. The mean exists for \code{rho > 1} and the variance for
#' \code{rho > 2}. Setting \code{k = 1} recovers the Waring distribution.
#'
#' @param x,q vector of counts
#' @param p vector of probabilities
#' @param n number of draws
#' @param mu mean, \code{mu > 0}
#' @param k shape, \code{k > 0}
#' @param rho tail parameter, \code{rho > 1}
#' @param log,log.p logical; return the log density / probability
#' @param lower.tail logical; if TRUE probabilities are \code{P[X <= x]}
#' @param max.value largest count searched when inverting the cdf
#' @return a numeric vector
#' @name gwaring
NULL

#' @rdname gwaring
#' @export
dgwaring <- function(x, mu = 1, k = 1, rho = 2.5, log = FALSE) {
  stopifnot(mu > 0, k > 0, rho > 1)
  out <- gwaring_logpmf_cpp(as.numeric(x), mu, k, rho)
  out[x < 0 | x != floor(x)] <- -Inf
  if (log) out else exp(out)
}
#' @rdname gwaring
#' @export
pgwaring <- function(q, mu = 1, k = 1, rho = 2.5, lower.tail = TRUE, log.p = FALSE) {
  hi <- max(0, floor(max(q)))
  cdf <- cumsum(dgwaring(0:hi, mu, k, rho))
  out <- ifelse(q < 0, 0, pmin(cdf[pmax(floor(q), 0) + 1L], 1))
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}
#' @rdname gwaring
#' @export
qgwaring <- function(p, mu = 1, k = 1, rho = 2.5, lower.tail = TRUE, log.p = FALSE,
                     max.value = 1e7) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  .qsearch(p, function(n) cumsum(dgwaring(0:n, mu, k, rho)), mu, max.value)
}
#' @rdname gwaring
#' @export
rgwaring <- function(n, mu = 1, k = 1, rho = 2.5) qgwaring(stats::runif(n), mu, k, rho)

## ---------------------------------------------------------------- Poisson-lognormal

#' Poisson-lognormal distribution
#'
#' \code{Y | Z ~ Poisson(exp(Z))} with \code{Z ~ N(muz, sigma^2)} and
#' \code{muz = log(mu) - sigma^2 / 2}, so \code{E[Y] = mu}. The pmf has no closed form and
#' is evaluated by adaptive Gauss-Hermite quadrature centred on the mode of the integrand.
#' The tail is lognormal-like: heavier than the negative binomial's, but still in the
#' Gumbel domain (\code{xi = 0}), so this family models a heavy body rather than a heavy
#' tail.
#'
#' @param x,q vector of counts
#' @param p vector of probabilities
#' @param n number of draws
#' @param mu mean, \code{mu > 0}
#' @param sigma standard deviation of the log-intensity, \code{sigma > 0}
#' @param log,log.p logical; return the log density / probability
#' @param lower.tail logical; if TRUE probabilities are \code{P[X <= x]}
#' @param max.value largest count searched when inverting the cdf
#' @return a numeric vector
#' @name poislnorm
NULL

#' @rdname poislnorm
#' @export
dpoislnorm <- function(x, mu = 1, sigma = 1, log = FALSE) {
  stopifnot(mu > 0, sigma > 0)
  out <- plnorm_logpmf_cpp(as.numeric(x), mu, sigma)
  out[x < 0 | x != floor(x)] <- -Inf
  if (log) out else exp(out)
}
#' @rdname poislnorm
#' @export
ppoislnorm <- function(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  hi <- max(0, floor(max(q)))
  cdf <- cumsum(dpoislnorm(0:hi, mu, sigma))
  out <- ifelse(q < 0, 0, pmin(cdf[pmax(floor(q), 0) + 1L], 1))
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}
#' @rdname poislnorm
#' @export
qpoislnorm <- function(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE,
                       max.value = 1e6) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  .qsearch(p, function(n) cumsum(dpoislnorm(0:n, mu, sigma)), mu, max.value)
}
#' @rdname poislnorm
#' @export
rpoislnorm <- function(n, mu = 1, sigma = 1) {
  z <- stats::rnorm(n, log(mu) - sigma^2 / 2, sigma)
  stats::rpois(n, exp(z))
}

## Shared cdf inversion: grow the grid geometrically from a scale hint until it covers
## max(p), then look up. Returns NA above max.value rather than looping forever.
.qsearch <- function(p, cdf_upto, scale, max.value) {
  n <- max(32, ceiling(4 * scale))
  repeat {
    cdf <- cdf_upto(n)
    if (max(cdf) >= max(p, na.rm = TRUE) - 1e-12 || n >= max.value) break
    n <- min(max.value, n * 4L)
  }
  vapply(p, function(pp) {
    if (is.na(pp)) return(NA_real_)
    i <- which(cdf >= pp - 1e-12)[1]
    if (is.na(i)) NA_real_ else i - 1
  }, 0)
}

## ---------------------------------------------------------------- likelihood bundles
## Each wrapper splits the flat coefficient vector by parameter and calls the Rcpp routine,
## mirroring R/pig_lik.R and R/gpig_lik.R.

.iG_gpois   <- function(v, mu, lambda) qgpois(v, mu = mu, lambda = lambda)
.iG_gwaring <- function(v, mu, k, rho)  qgwaring(v, mu = mu, k = k, rho = rho)
.iG_plnorm  <- function(v, mu, sigma)   qpoislnorm(v, mu = mu, sigma = sigma)

.gpois1.d0 <- function(pars, likdata) {
  if (likdata$censored) stop("Censored likelihoods not available for gpois.")
  gpois1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
           likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
           likdata$dupid, likdata$duplicate, likdata$offsets)
}
.gpois1.d12 <- function(pars, likdata) {
  gpois1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
            likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
            likdata$dupid, likdata$duplicate, likdata$offsets)
}
.gpois1fns <- list(d0 = .gpois1.d0, d120 = .gpois1.d12, d340 = NULL, m = 1, iG = .iG_gpois)

.gwaring1.d0 <- function(pars, likdata) {
  if (likdata$censored) stop("Censored likelihoods not available for gwaring.")
  gwaring1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
             likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[, 1],
             likdata$dupid, likdata$duplicate, likdata$offsets)
}
.gwaring1.d12 <- function(pars, likdata) {
  gwaring1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
              likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[, 1],
              likdata$dupid, likdata$duplicate, likdata$offsets)
}
.gwaring1fns <- list(d0 = .gwaring1.d0, d120 = .gwaring1.d12, d340 = NULL, m = 1,
                     iG = .iG_gwaring)

.plnorm1.d0 <- function(pars, likdata) {
  if (likdata$censored) stop("Censored likelihoods not available for plnorm.")
  plnorm1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
            likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
            likdata$dupid, likdata$duplicate, likdata$offsets)
}
.plnorm1.d12 <- function(pars, likdata) {
  plnorm1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
             likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
             likdata$dupid, likdata$duplicate, likdata$offsets)
}
.plnorm1fns <- list(d0 = .plnorm1.d0, d120 = .plnorm1.d12, d340 = NULL, m = 1,
                    iG = .iG_plnorm)
