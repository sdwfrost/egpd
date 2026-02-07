## gamlss family constructors for EGPD models
##
## These allow fitting EGPD, DEGPD, ZIEGPD, and ZIDEGPD models using the
## gamlss package. Each family maps gamlss parameters (mu, sigma, nu, tau)
## to the underlying egpd distribution parameters.
##
## Parameter mapping (egpd -> gamlss):
##   m=1: sigma->mu, xi->sigma, kappa->nu
##   m=3: sigma->mu, xi->sigma, delta->nu
##   m=4: sigma->mu, xi->sigma, delta->nu, kappa->tau
##   ZI adds: pi->tau (m=1,3 only; m=4 would need 5 params)

## ---------------------------------------------------------------
## Internal helpers
## ---------------------------------------------------------------

## Numerical first derivative of log-density via central differences
.gamlss_nd <- function(logd_fn, y, par, par_idx, h_scale = 1e-4) {
  h <- pmax(abs(par) * h_scale, 1e-8)
  par_plus <- par_minus <- par
  par_plus <- par + h
  par_minus <- par - h
  (logd_fn(y, par_plus, par_idx) - logd_fn(y, par_minus, par_idx)) / (2 * h)
}

## ---------------------------------------------------------------
## d/p/q/r wrappers: Continuous EGPD
## ---------------------------------------------------------------

#' @title gamlss Distribution Functions for Continuous EGPD Model 1
#' @name EGPD1
#' @description Density, distribution function, quantile function, and random
#'   generation for the continuous EGPD with G-transformation type 1
#'   (power: \eqn{G(u) = u^\kappa}), parameterised for use with \code{gamlss}.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi), positive
#' @param nu G-transformation parameter (kappa), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @return Numeric vector
#' @export
dEGPD1 <- function(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE) {
  degpd_density(x, sigma = mu, xi = sigma, kappa = nu, type = 1L, log = log)
}

#' @rdname EGPD1
#' @export
pEGPD1 <- function(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  p <- pegpd(q, sigma = mu, xi = sigma, kappa = nu, type = 1L)
  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname EGPD1
#' @export
qEGPD1 <- function(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  qegpd(p, sigma = mu, xi = sigma, kappa = nu, type = 1L)
}

#' @rdname EGPD1
#' @export
rEGPD1 <- function(n, mu = 1, sigma = 0.5, nu = 1) {
  regpd(n, sigma = mu, xi = sigma, kappa = nu, type = 1L)
}

#' @title gamlss Distribution Functions for Continuous EGPD Model 3
#' @name EGPD3
#' @description Density, distribution function, quantile function, and random
#'   generation for the continuous EGPD with G-transformation type 3
#'   (incomplete beta), parameterised for use with \code{gamlss}.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi), positive
#' @param nu G-transformation parameter (delta), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @return Numeric vector
#' @export
dEGPD3 <- function(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE) {
  degpd_density(x, sigma = mu, xi = sigma, delta = nu, type = 4L, log = log)
}

#' @rdname EGPD3
#' @export
pEGPD3 <- function(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  p <- pegpd(q, sigma = mu, xi = sigma, delta = nu, type = 4L)
  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname EGPD3
#' @export
qEGPD3 <- function(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  qegpd(p, sigma = mu, xi = sigma, delta = nu, type = 4L)
}

#' @rdname EGPD3
#' @export
rEGPD3 <- function(n, mu = 1, sigma = 0.5, nu = 1) {
  regpd(n, sigma = mu, xi = sigma, delta = nu, type = 4L)
}

#' @title gamlss Distribution Functions for Continuous EGPD Model 4
#' @name EGPD4
#' @description Density, distribution function, quantile function, and random
#'   generation for the continuous EGPD with G-transformation type 4
#'   (power-beta), parameterised for use with \code{gamlss}.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi), positive
#' @param nu G-transformation parameter (delta), positive
#' @param tau G-transformation parameter (kappa), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @return Numeric vector
#' @export
dEGPD4 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 1, log = FALSE) {
  degpd_density(x, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L, log = log)
}

#' @rdname EGPD4
#' @export
pEGPD4 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE) {
  p <- pegpd(q, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L)
  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname EGPD4
#' @export
qEGPD4 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  qegpd(p, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L)
}

#' @rdname EGPD4
#' @export
rEGPD4 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 1) {
  regpd(n, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L)
}

## ---------------------------------------------------------------
## d/p/q/r wrappers: Discrete EGPD
## ---------------------------------------------------------------

#' @title gamlss Distribution Functions for Discrete EGPD Model 1
#' @name DEGPD1
#' @description PMF, distribution function, quantile function, and random
#'   generation for the discrete EGPD with G-transformation type 1
#'   (power), parameterised for use with \code{gamlss}.
#' @param x,q vector of quantiles (non-negative integers)
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi), positive
#' @param nu G-transformation parameter (kappa), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @return Numeric vector
#' @export
dDEGPD1 <- function(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE) {
  d <- ddiscegpd(x, sigma = mu, xi = sigma, kappa = nu, type = 1L)
  if (log) d <- ifelse(d > 0, log(d), -Inf)
  d
}

#' @rdname DEGPD1
#' @export
pDEGPD1 <- function(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  p <- pdiscegpd(q, sigma = mu, xi = sigma, kappa = nu, type = 1L)
  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname DEGPD1
#' @export
qDEGPD1 <- function(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  qdiscegpd(p, sigma = mu, xi = sigma, kappa = nu, type = 1L)
}

#' @rdname DEGPD1
#' @export
rDEGPD1 <- function(n, mu = 1, sigma = 0.5, nu = 1) {
  rdiscegpd(n, sigma = mu, xi = sigma, kappa = nu, type = 1L)
}

#' @title gamlss Distribution Functions for Discrete EGPD Model 3
#' @name DEGPD3
#' @description PMF, distribution function, quantile function, and random
#'   generation for the discrete EGPD with G-transformation type 3
#'   (incomplete beta), parameterised for use with \code{gamlss}.
#' @param x,q vector of quantiles (non-negative integers)
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi), positive
#' @param nu G-transformation parameter (delta), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @return Numeric vector
#' @export
dDEGPD3 <- function(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE) {
  d <- ddiscegpd(x, sigma = mu, xi = sigma, delta = nu, type = 4L)
  if (log) d <- ifelse(d > 0, log(d), -Inf)
  d
}

#' @rdname DEGPD3
#' @export
pDEGPD3 <- function(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  p <- pdiscegpd(q, sigma = mu, xi = sigma, delta = nu, type = 4L)
  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname DEGPD3
#' @export
qDEGPD3 <- function(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  qdiscegpd(p, sigma = mu, xi = sigma, delta = nu, type = 4L)
}

#' @rdname DEGPD3
#' @export
rDEGPD3 <- function(n, mu = 1, sigma = 0.5, nu = 1) {
  rdiscegpd(n, sigma = mu, xi = sigma, delta = nu, type = 4L)
}

#' @title gamlss Distribution Functions for Discrete EGPD Model 4
#' @name DEGPD4
#' @description PMF, distribution function, quantile function, and random
#'   generation for the discrete EGPD with G-transformation type 4
#'   (power-beta), parameterised for use with \code{gamlss}.
#' @param x,q vector of quantiles (non-negative integers)
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi), positive
#' @param nu G-transformation parameter (delta), positive
#' @param tau G-transformation parameter (kappa), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @return Numeric vector
#' @export
dDEGPD4 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 1, log = FALSE) {
  d <- ddiscegpd(x, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L)
  if (log) d <- ifelse(d > 0, log(d), -Inf)
  d
}

#' @rdname DEGPD4
#' @export
pDEGPD4 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE) {
  p <- pdiscegpd(q, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L)
  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname DEGPD4
#' @export
qDEGPD4 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  qdiscegpd(p, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L)
}

#' @rdname DEGPD4
#' @export
rDEGPD4 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 1) {
  rdiscegpd(n, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L)
}

## ---------------------------------------------------------------
## d/p/q/r wrappers: Zero-Inflated Continuous EGPD
## ---------------------------------------------------------------

#' @title gamlss Distribution Functions for Zero-Inflated EGPD Model 1
#' @name ZIEGPD1
#' @description Density, distribution function, quantile function, and random
#'   generation for the zero-inflated continuous EGPD with G-transformation
#'   type 1 (power), parameterised for use with \code{gamlss}.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi), positive
#' @param nu G-transformation parameter (kappa), positive
#' @param tau zero-inflation probability (pi), in (0,1)
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @return Numeric vector
#' @export
dZIEGPD1 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE) {
  dziegpd(x, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L, log = log)
}

#' @rdname ZIEGPD1
#' @export
pZIEGPD1 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  p <- pziegpd(q, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L)
  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname ZIEGPD1
#' @export
qZIEGPD1 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  qziegpd(p, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L)
}

#' @rdname ZIEGPD1
#' @export
rZIEGPD1 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1) {
  rziegpd(n, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L)
}

#' @title gamlss Distribution Functions for Zero-Inflated EGPD Model 3
#' @name ZIEGPD3
#' @description Density, distribution function, quantile function, and random
#'   generation for the zero-inflated continuous EGPD with G-transformation
#'   type 3 (incomplete beta), parameterised for use with \code{gamlss}.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi), positive
#' @param nu G-transformation parameter (delta), positive
#' @param tau zero-inflation probability (pi), in (0,1)
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @return Numeric vector
#' @export
dZIEGPD3 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE) {
  dziegpd(x, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L, log = log)
}

#' @rdname ZIEGPD3
#' @export
pZIEGPD3 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  p <- pziegpd(q, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L)
  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname ZIEGPD3
#' @export
qZIEGPD3 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  qziegpd(p, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L)
}

#' @rdname ZIEGPD3
#' @export
rZIEGPD3 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1) {
  rziegpd(n, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L)
}

## ---------------------------------------------------------------
## d/p/q/r wrappers: Zero-Inflated Discrete EGPD
## ---------------------------------------------------------------

#' @title gamlss Distribution Functions for Zero-Inflated Discrete EGPD Model 1
#' @name ZIDEGPD1
#' @description PMF, distribution function, quantile function, and random
#'   generation for the zero-inflated discrete EGPD with G-transformation
#'   type 1 (power), parameterised for use with \code{gamlss}.
#' @param x,q vector of quantiles (non-negative integers)
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi), positive
#' @param nu G-transformation parameter (kappa), positive
#' @param tau zero-inflation probability (pi), in (0,1)
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @return Numeric vector
#' @export
dZIDEGPD1 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE) {
  d <- dzidiscegpd(x, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L)
  if (log) d <- ifelse(d > 0, log(d), -Inf)
  d
}

#' @rdname ZIDEGPD1
#' @export
pZIDEGPD1 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  p <- pzidiscegpd(q, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L)
  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname ZIDEGPD1
#' @export
qZIDEGPD1 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  qzidiscegpd(p, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L)
}

#' @rdname ZIDEGPD1
#' @export
rZIDEGPD1 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1) {
  rzidiscegpd(n, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L)
}

#' @title gamlss Distribution Functions for Zero-Inflated Discrete EGPD Model 3
#' @name ZIDEGPD3
#' @description PMF, distribution function, quantile function, and random
#'   generation for the zero-inflated discrete EGPD with G-transformation
#'   type 3 (incomplete beta), parameterised for use with \code{gamlss}.
#' @param x,q vector of quantiles (non-negative integers)
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi), positive
#' @param nu G-transformation parameter (delta), positive
#' @param tau zero-inflation probability (pi), in (0,1)
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @return Numeric vector
#' @export
dZIDEGPD3 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE) {
  d <- dzidiscegpd(x, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L)
  if (log) d <- ifelse(d > 0, log(d), -Inf)
  d
}

#' @rdname ZIDEGPD3
#' @export
pZIDEGPD3 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  p <- pzidiscegpd(q, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L)
  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname ZIDEGPD3
#' @export
qZIDEGPD3 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  qzidiscegpd(p, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L)
}

#' @rdname ZIDEGPD3
#' @export
rZIDEGPD3 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1) {
  rzidiscegpd(n, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L)
}


## ===============================================================
## gamlss Family Constructors
## ===============================================================

## Numerical derivative of log-density via central differences.
## dfn is the actual density function (resolved at family construction time).
## These are package-level helpers called from derivative functions.
.egpd_nd3 <- function(y, mu, sigma, nu, dfn, idx, h_scale = 1e-4) {
  params <- list(mu, sigma, nu)
  par_val <- params[[idx]]
  h <- pmax(abs(par_val) * h_scale, 1e-8)
  params_p <- params_m <- params
  params_p[[idx]] <- par_val + h
  params_m[[idx]] <- par_val - h
  ld_p <- dfn(y, mu = params_p[[1]], sigma = params_p[[2]], nu = params_p[[3]], log = TRUE)
  ld_m <- dfn(y, mu = params_m[[1]], sigma = params_m[[2]], nu = params_m[[3]], log = TRUE)
  dl <- (ld_p - ld_m) / (2 * h)
  dl[!is.finite(dl)] <- 0
  dl
}

.egpd_nd4 <- function(y, mu, sigma, nu, tau, dfn, idx, h_scale = 1e-4) {
  params <- list(mu, sigma, nu, tau)
  par_val <- params[[idx]]
  h <- pmax(abs(par_val) * h_scale, 1e-8)
  params_p <- params_m <- params
  params_p[[idx]] <- par_val + h
  params_m[[idx]] <- par_val - h
  ld_p <- dfn(y, mu = params_p[[1]], sigma = params_p[[2]], nu = params_p[[3]], tau = params_p[[4]], log = TRUE)
  ld_m <- dfn(y, mu = params_m[[1]], sigma = params_m[[2]], nu = params_m[[3]], tau = params_m[[4]], log = TRUE)
  dl <- (ld_p - ld_m) / (2 * h)
  dl[!is.finite(dl)] <- 0
  dl
}

## Internal factory for building gamlss.family objects.
##
## gamlss extracts function bodies via body() and evaluates them as expressions
## in its own environment. This means all function bodies must contain literal
## function names (symbols), not variable references. We use bquote() to
## substitute the actual function symbols into the function bodies at
## construction time.

.make_gamlss_family <- function(
    dfun, pfun,
    family_name, dist_type,
    npar,
    default_links,
    mu.link = NULL, sigma.link = NULL, nu.link = NULL, tau.link = NULL,
    mass.p = NULL
) {
  ## Resolve links
  mu_link    <- if (!is.null(mu.link))    mu.link    else default_links[["mu"]]
  sigma_link <- if (!is.null(sigma.link)) sigma.link else default_links[["sigma"]]
  nu_link    <- if (!is.null(nu.link))    nu.link    else default_links[["nu"]]

  ## Validate links
  gamlss.dist::checklink("mu.link",    family_name, substitute(mu_link),    c("log", "identity", "inverse"))
  gamlss.dist::checklink("sigma.link", family_name, substitute(sigma_link), c("log", "identity", "inverse"))
  gamlss.dist::checklink("nu.link",    family_name, substitute(nu_link),    c("log", "logit", "identity"))

  if (npar == 4) {
    tau_link <- if (!is.null(tau.link)) tau.link else default_links[["tau"]]
    gamlss.dist::checklink("tau.link", family_name, substitute(tau_link), c("log", "logit", "identity"))
  }

  ## Convert string names to symbols for bquote substitution
  dfun_sym <- as.name(dfun)
  pfun_str <- pfun

  if (npar == 3) {
    ## ---- 3-parameter family ----

    ## Build functions with literal density function names baked in
    G.dev.incr_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      -2 * .(dfun_sym)(y, mu = mu, sigma = sigma, nu = nu, log = TRUE)
    }))

    dldm_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 1L)
    }))
    dldd_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 2L)
    }))
    dldv_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 3L)
    }))

    d2ldm2_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      dl <- egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 1L)
      d2 <- -dl * dl
      ifelse(d2 < -1e-15, d2, -1e-15)
    }))
    d2ldd2_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      dl <- egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 2L)
      d2 <- -dl * dl
      ifelse(d2 < -1e-15, d2, -1e-15)
    }))
    d2ldv2_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      dl <- egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 3L)
      d2 <- -dl * dl
      ifelse(d2 < -1e-15, d2, -1e-15)
    }))

    d2ldmdd_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      -egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 1L) *
        egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 2L)
    }))
    d2ldmdv_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      -egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 1L) *
        egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 3L)
    }))
    d2ldddv_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      -egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 2L) *
        egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 3L)
    }))

    if (dist_type == "Discrete") {
      rqres_expr <- eval(bquote(expression(
        rqres(pfun = .(pfun_str), type = "Discrete", ymin = 0, y = y,
              mu = mu, sigma = sigma, nu = nu)
      )))
    } else {
      rqres_expr <- eval(bquote(expression(
        rqres(pfun = .(pfun_str), type = "Continuous", y = y,
              mu = mu, sigma = sigma, nu = nu)
      )))
    }

    fam <- list(
      family = c(family_name, "gamlss.family"),
      parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE),
      nopar = 3L,
      type = dist_type,

      mu.link    = mu_link,
      sigma.link = sigma_link,
      nu.link    = nu_link,

      mu.linkfun    = make.link(mu_link)$linkfun,
      mu.linkinv    = make.link(mu_link)$linkinv,
      mu.dr         = make.link(mu_link)$mu.eta,
      sigma.linkfun = make.link(sigma_link)$linkfun,
      sigma.linkinv = make.link(sigma_link)$linkinv,
      sigma.dr      = make.link(sigma_link)$mu.eta,
      nu.linkfun    = make.link(nu_link)$linkfun,
      nu.linkinv    = make.link(nu_link)$linkinv,
      nu.dr         = make.link(nu_link)$mu.eta,

      mu.initial    = expression(mu    <- rep(max(sd(y), 0.1), length(y))),
      sigma.initial = expression(sigma <- rep(0.3, length(y))),
      nu.initial    = expression(nu    <- rep(1, length(y))),

      G.dev.incr = G.dev.incr_fn,
      rqres      = rqres_expr,

      mu.valid    = function(mu) all(mu > 0),
      sigma.valid = function(sigma) all(sigma > 0),
      nu.valid    = function(nu) all(nu > 0),
      y.valid     = function(y) all(y >= 0),

      dldm = dldm_fn,
      dldd = dldd_fn,
      dldv = dldv_fn,
      d2ldm2  = d2ldm2_fn,
      d2ldd2  = d2ldd2_fn,
      d2ldv2  = d2ldv2_fn,
      d2ldmdd = d2ldmdd_fn,
      d2ldmdv = d2ldmdv_fn,
      d2ldddv = d2ldddv_fn
    )

  } else {
    ## ---- 4-parameter family ----

    G.dev.incr_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      -2 * .(dfun_sym)(y, mu = mu, sigma = sigma, nu = nu, tau = tau, log = TRUE)
    }))

    dldm_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 1L)
    }))
    dldd_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 2L)
    }))
    dldv_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 3L)
    }))
    dldt_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 4L)
    }))

    d2ldm2_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      dl <- egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 1L)
      d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
    }))
    d2ldd2_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      dl <- egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 2L)
      d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
    }))
    d2ldv2_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      dl <- egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 3L)
      d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
    }))
    d2ldt2_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      dl <- egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 4L)
      d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
    }))

    d2ldmdd_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      -egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 1L) *
        egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 2L)
    }))
    d2ldmdv_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      -egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 1L) *
        egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 3L)
    }))
    d2ldmdt_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      -egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 1L) *
        egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 4L)
    }))
    d2ldddv_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      -egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 2L) *
        egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 3L)
    }))
    d2ldddt_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      -egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 2L) *
        egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 4L)
    }))
    d2ldvdt_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      -egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 3L) *
        egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = 4L)
    }))

    if (dist_type == "Discrete") {
      rqres_expr <- eval(bquote(expression(
        rqres(pfun = .(pfun_str), type = "Discrete", ymin = 0, y = y,
              mu = mu, sigma = sigma, nu = nu, tau = tau)
      )))
    } else if (dist_type == "Mixed") {
      ## For Mixed (zero-inflated continuous): mass at 0 with prob tau
      rqres_expr <- eval(bquote(expression(
        rqres(pfun = .(pfun_str), type = "Mixed",
              mass.p = 0, prob.mp = cbind(tau),
              y = y, mu = mu, sigma = sigma, nu = nu, tau = tau)
      )))
    } else {
      rqres_expr <- eval(bquote(expression(
        rqres(pfun = .(pfun_str), type = "Continuous", y = y,
              mu = mu, sigma = sigma, nu = nu, tau = tau)
      )))
    }

    fam <- list(
      family = c(family_name, "gamlss.family"),
      parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE, tau = TRUE),
      nopar = 4L,
      type = dist_type,

      mu.link    = mu_link,
      sigma.link = sigma_link,
      nu.link    = nu_link,
      tau.link   = tau_link,

      mu.linkfun    = make.link(mu_link)$linkfun,
      mu.linkinv    = make.link(mu_link)$linkinv,
      mu.dr         = make.link(mu_link)$mu.eta,
      sigma.linkfun = make.link(sigma_link)$linkfun,
      sigma.linkinv = make.link(sigma_link)$linkinv,
      sigma.dr      = make.link(sigma_link)$mu.eta,
      nu.linkfun    = make.link(nu_link)$linkfun,
      nu.linkinv    = make.link(nu_link)$linkinv,
      nu.dr         = make.link(nu_link)$mu.eta,
      tau.linkfun   = make.link(tau_link)$linkfun,
      tau.linkinv   = make.link(tau_link)$linkinv,
      tau.dr        = make.link(tau_link)$mu.eta,

      mu.initial    = expression(mu    <- rep(max(sd(y), 0.1), length(y))),
      sigma.initial = expression(sigma <- rep(0.3, length(y))),
      nu.initial    = expression(nu    <- rep(1, length(y))),
      tau.initial   = expression(tau   <- rep(0.1, length(y))),

      G.dev.incr = G.dev.incr_fn,
      rqres      = rqres_expr,

      mu.valid    = function(mu) all(mu > 0),
      sigma.valid = function(sigma) all(sigma > 0),
      nu.valid    = function(nu) all(nu > 0),
      tau.valid   = function(tau) all(tau > 0),
      y.valid     = function(y) all(y >= 0),

      dldm = dldm_fn,
      dldd = dldd_fn,
      dldv = dldv_fn,
      dldt = dldt_fn,
      d2ldm2  = d2ldm2_fn,
      d2ldd2  = d2ldd2_fn,
      d2ldv2  = d2ldv2_fn,
      d2ldt2  = d2ldt2_fn,
      d2ldmdd = d2ldmdd_fn,
      d2ldmdv = d2ldmdv_fn,
      d2ldmdt = d2ldmdt_fn,
      d2ldddv = d2ldddv_fn,
      d2ldddt = d2ldddt_fn,
      d2ldvdt = d2ldvdt_fn
    )
  }

  class(fam) <- "gamlss.family"
  fam
}


## ---------------------------------------------------------------
## Continuous EGPD families
## ---------------------------------------------------------------

#' gamlss Family for Continuous EGPD Model 1
#'
#' Creates a \code{gamlss.family} object for fitting continuous
#' Extended Generalized Pareto Distribution models with model 1
#' (power G-transformation: \eqn{G(u) = u^\kappa}).
#'
#' The gamlss parameters map to EGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation).
#'
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"})
#' @param nu.link link function for nu (default \code{"log"})
#'
#' @return An object of class \code{gamlss.family}
#' @seealso \code{\link{dEGPD1}}, \code{\link{egpd_bamlss}}
#' @export
EGPD1 <- function(mu.link = "log", sigma.link = "log", nu.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dEGPD1", pfun = "pEGPD1",
    family_name = "EGPD1", dist_type = "Continuous",
    npar = 3,
    default_links = c(mu = "log", sigma = "log", nu = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link
  )
  fam$mu.initial    <- expression(mu    <- rep(max(sd(y), 0.1), length(y)))
  fam$sigma.initial <- expression(sigma <- rep(0.3, length(y)))
  fam$nu.initial    <- expression(nu    <- rep(max(mean(y), 0.1), length(y)))
  fam
}

#' gamlss Family for Continuous EGPD Model 3
#'
#' Creates a \code{gamlss.family} object for fitting continuous
#' Extended Generalized Pareto Distribution models with model 3
#' (incomplete beta G-transformation).
#'
#' The gamlss parameters map to EGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = delta (G-transformation).
#'
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"})
#' @param nu.link link function for nu (default \code{"log"})
#'
#' @return An object of class \code{gamlss.family}
#' @seealso \code{\link{dEGPD3}}, \code{\link{egpd_bamlss}}
#' @export
EGPD3 <- function(mu.link = "log", sigma.link = "log", nu.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dEGPD3", pfun = "pEGPD3",
    family_name = "EGPD3", dist_type = "Continuous",
    npar = 3,
    default_links = c(mu = "log", sigma = "log", nu = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link
  )
  fam$mu.initial    <- expression(mu    <- rep(max(sd(y), 0.1), length(y)))
  fam$sigma.initial <- expression(sigma <- rep(0.3, length(y)))
  fam$nu.initial    <- expression(nu    <- rep(max(0.5 * mean(y), 0.1), length(y)))
  fam
}

#' gamlss Family for Continuous EGPD Model 4
#'
#' Creates a \code{gamlss.family} object for fitting continuous
#' Extended Generalized Pareto Distribution models with model 4
#' (power-beta G-transformation).
#'
#' The gamlss parameters map to EGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = delta (G-transformation), \code{tau} = kappa (G-transformation).
#'
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"})
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"log"})
#'
#' @return An object of class \code{gamlss.family}
#' @seealso \code{\link{dEGPD4}}, \code{\link{egpd_bamlss}}
#' @export
EGPD4 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dEGPD4", pfun = "pEGPD4",
    family_name = "EGPD4", dist_type = "Continuous",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam$mu.initial    <- expression(mu    <- rep(max(sd(y), 0.1), length(y)))
  fam$sigma.initial <- expression(sigma <- rep(0.3, length(y)))
  fam$nu.initial    <- expression(nu    <- rep(max(0.5 * mean(y), 0.1), length(y)))
  fam$tau.initial   <- expression(tau   <- rep(max(mean(y), 0.1), length(y)))
  fam$tau.valid     <- function(tau) all(tau > 0)
  fam
}

## ---------------------------------------------------------------
## Discrete EGPD families
## ---------------------------------------------------------------

#' gamlss Family for Discrete EGPD Model 1
#'
#' Creates a \code{gamlss.family} object for fitting discrete
#' Extended Generalized Pareto Distribution models with model 1
#' (power G-transformation).
#'
#' The gamlss parameters map to DEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation).
#'
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"})
#' @param nu.link link function for nu (default \code{"log"})
#'
#' @return An object of class \code{gamlss.family}
#' @seealso \code{\link{dDEGPD1}}, \code{\link{degpd_bamlss}}
#' @export
DEGPD1 <- function(mu.link = "log", sigma.link = "log", nu.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dDEGPD1", pfun = "pDEGPD1",
    family_name = "DEGPD1", dist_type = "Discrete",
    npar = 3,
    default_links = c(mu = "log", sigma = "log", nu = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link
  )
  fam$mu.initial    <- expression(mu    <- rep(max(sd(y), 0.1), length(y)))
  fam$sigma.initial <- expression(sigma <- rep(0.3, length(y)))
  fam$nu.initial    <- expression(nu    <- rep(max(mean(y), 0.1), length(y)))
  fam
}

#' gamlss Family for Discrete EGPD Model 3
#'
#' Creates a \code{gamlss.family} object for fitting discrete
#' Extended Generalized Pareto Distribution models with model 3
#' (incomplete beta G-transformation).
#'
#' The gamlss parameters map to DEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = delta (G-transformation).
#'
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"})
#' @param nu.link link function for nu (default \code{"log"})
#'
#' @return An object of class \code{gamlss.family}
#' @seealso \code{\link{dDEGPD3}}, \code{\link{degpd_bamlss}}
#' @export
DEGPD3 <- function(mu.link = "log", sigma.link = "log", nu.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dDEGPD3", pfun = "pDEGPD3",
    family_name = "DEGPD3", dist_type = "Discrete",
    npar = 3,
    default_links = c(mu = "log", sigma = "log", nu = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link
  )
  fam$mu.initial    <- expression(mu    <- rep(max(sd(y), 0.1), length(y)))
  fam$sigma.initial <- expression(sigma <- rep(0.3, length(y)))
  fam$nu.initial    <- expression(nu    <- rep(max(0.5 * mean(y), 0.1), length(y)))
  fam
}

#' gamlss Family for Discrete EGPD Model 4
#'
#' Creates a \code{gamlss.family} object for fitting discrete
#' Extended Generalized Pareto Distribution models with model 4
#' (power-beta G-transformation).
#'
#' The gamlss parameters map to DEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = delta (G-transformation), \code{tau} = kappa (G-transformation).
#'
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"})
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"log"})
#'
#' @return An object of class \code{gamlss.family}
#' @seealso \code{\link{dDEGPD4}}, \code{\link{degpd_bamlss}}
#' @export
DEGPD4 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dDEGPD4", pfun = "pDEGPD4",
    family_name = "DEGPD4", dist_type = "Discrete",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam$mu.initial    <- expression(mu    <- rep(max(sd(y), 0.1), length(y)))
  fam$sigma.initial <- expression(sigma <- rep(0.3, length(y)))
  fam$nu.initial    <- expression(nu    <- rep(max(0.5 * mean(y), 0.1), length(y)))
  fam$tau.initial   <- expression(tau   <- rep(max(mean(y), 0.1), length(y)))
  fam$tau.valid     <- function(tau) all(tau > 0)
  fam
}

## ---------------------------------------------------------------
## Zero-Inflated Continuous EGPD families
## ---------------------------------------------------------------

#' gamlss Family for Zero-Inflated Continuous EGPD Model 1
#'
#' Creates a \code{gamlss.family} object for fitting zero-inflated
#' continuous Extended Generalized Pareto Distribution models with model 1
#' (power G-transformation).
#'
#' The gamlss parameters map to ZIEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation), \code{tau} = pi (zero-inflation).
#'
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"})
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"logit"})
#'
#' @return An object of class \code{gamlss.family}
#' @seealso \code{\link{dZIEGPD1}}, \code{\link{ziegpd_bamlss}}
#' @export
ZIEGPD1 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "logit") {
  fam <- .make_gamlss_family(
    dfun = "dZIEGPD1", pfun = "pZIEGPD1",
    family_name = "ZIEGPD1", dist_type = "Mixed",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "logit"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam$mu.initial    <- expression(mu    <- rep(max(sd(y[y > 0]), 0.1), length(y)))
  fam$sigma.initial <- expression(sigma <- rep(0.3, length(y)))
  fam$nu.initial    <- expression(nu    <- rep(max(mean(y), 0.1), length(y)))
  fam$tau.initial   <- expression(tau   <- rep(max(mean(y == 0), 0.01), length(y)))
  fam$tau.valid     <- function(tau) all(tau > 0 & tau < 1)
  fam
}

#' gamlss Family for Zero-Inflated Continuous EGPD Model 3
#'
#' Creates a \code{gamlss.family} object for fitting zero-inflated
#' continuous Extended Generalized Pareto Distribution models with model 3
#' (incomplete beta G-transformation).
#'
#' The gamlss parameters map to ZIEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = delta (G-transformation), \code{tau} = pi (zero-inflation).
#'
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"})
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"logit"})
#'
#' @return An object of class \code{gamlss.family}
#' @seealso \code{\link{dZIEGPD3}}, \code{\link{ziegpd_bamlss}}
#' @export
ZIEGPD3 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "logit") {
  fam <- .make_gamlss_family(
    dfun = "dZIEGPD3", pfun = "pZIEGPD3",
    family_name = "ZIEGPD3", dist_type = "Mixed",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "logit"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam$mu.initial    <- expression(mu    <- rep(max(sd(y[y > 0]), 0.1), length(y)))
  fam$sigma.initial <- expression(sigma <- rep(0.3, length(y)))
  fam$nu.initial    <- expression(nu    <- rep(max(0.5 * mean(y), 0.1), length(y)))
  fam$tau.initial   <- expression(tau   <- rep(max(mean(y == 0), 0.01), length(y)))
  fam$tau.valid     <- function(tau) all(tau > 0 & tau < 1)
  fam
}

## ---------------------------------------------------------------
## Zero-Inflated Discrete EGPD families
## ---------------------------------------------------------------

#' gamlss Family for Zero-Inflated Discrete EGPD Model 1
#'
#' Creates a \code{gamlss.family} object for fitting zero-inflated
#' discrete Extended Generalized Pareto Distribution models with model 1
#' (power G-transformation).
#'
#' The gamlss parameters map to ZIDEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation), \code{tau} = pi (zero-inflation).
#'
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"})
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"logit"})
#'
#' @return An object of class \code{gamlss.family}
#' @seealso \code{\link{dZIDEGPD1}}, \code{\link{zidegpd_bamlss}}
#' @export
ZIDEGPD1 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "logit") {
  fam <- .make_gamlss_family(
    dfun = "dZIDEGPD1", pfun = "pZIDEGPD1",
    family_name = "ZIDEGPD1", dist_type = "Discrete",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "logit"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam$mu.initial    <- expression(mu    <- rep(max(sd(y[y > 0]), 0.1), length(y)))
  fam$sigma.initial <- expression(sigma <- rep(0.3, length(y)))
  fam$nu.initial    <- expression(nu    <- rep(max(mean(y), 0.1), length(y)))
  fam$tau.initial   <- expression(tau   <- rep(max(mean(y == 0), 0.01), length(y)))
  fam$tau.valid     <- function(tau) all(tau > 0 & tau < 1)
  fam
}

#' gamlss Family for Zero-Inflated Discrete EGPD Model 3
#'
#' Creates a \code{gamlss.family} object for fitting zero-inflated
#' discrete Extended Generalized Pareto Distribution models with model 3
#' (incomplete beta G-transformation).
#'
#' The gamlss parameters map to ZIDEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = delta (G-transformation), \code{tau} = pi (zero-inflation).
#'
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"})
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"logit"})
#'
#' @return An object of class \code{gamlss.family}
#' @seealso \code{\link{dZIDEGPD3}}, \code{\link{zidegpd_bamlss}}
#' @export
ZIDEGPD3 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "logit") {
  fam <- .make_gamlss_family(
    dfun = "dZIDEGPD3", pfun = "pZIDEGPD3",
    family_name = "ZIDEGPD3", dist_type = "Discrete",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "logit"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam$mu.initial    <- expression(mu    <- rep(max(sd(y[y > 0]), 0.1), length(y)))
  fam$sigma.initial <- expression(sigma <- rep(0.3, length(y)))
  fam$nu.initial    <- expression(nu    <- rep(max(0.5 * mean(y), 0.1), length(y)))
  fam$tau.initial   <- expression(tau   <- rep(max(mean(y == 0), 0.01), length(y)))
  fam$tau.valid     <- function(tau) all(tau > 0 & tau < 1)
  fam
}
