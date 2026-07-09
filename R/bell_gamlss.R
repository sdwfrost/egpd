## gamlss family for the Bell distribution (Castellares, Ferrari & Lemonte 2018)
## and zero-inflated Bell (ZIBell). Because gamlss.dist has no Bell, these build
## on the package functions (dbell/pbell, R/bell_dist.R).
##
## Mean parameterisation (consistent with the EGPD/GPIG gamlss families):
##   mu    = mean            (log link),   internally theta = W0(mu)
##   sigma = pi, zero-inflation (ZIBell only, logit link)
## Bell is a one-parameter family, so the plain BELL() gamlss family has only mu;
## ZIBELL() adds sigma as the zero-inflation probability. Derivatives of the
## log-likelihood are obtained by central differences (.egpd_nd1 / .egpd_nd2),
## matching the approach used for the EGPD gamlss families.

## theta = W0(mu) (Lambert W); vectorised via the C++ helper in src/bell.cpp.
.bell_theta_from_mean <- function(mu) bell_W0_cpp(as.numeric(mu))

## ---- gamlss-convention d/p/q/r (mu = mean [, sigma = pi]) ----

#' @rdname BELL
#' @export
dBELL <- function(x, mu = 1, log = FALSE) {
  dbell(x, theta = .bell_theta_from_mean(mu), log = log)
}
#' @rdname BELL
#' @export
pBELL <- function(q, mu = 1, lower.tail = TRUE, log.p = FALSE) {
  pbell(q, theta = .bell_theta_from_mean(mu), lower.tail = lower.tail, log.p = log.p)
}
#' @rdname BELL
#' @export
qBELL <- function(p, mu = 1, lower.tail = TRUE, log.p = FALSE, max.value = 1e5) {
  qbell(p, theta = .bell_theta_from_mean(mu), lower.tail = lower.tail,
        log.p = log.p, max.value = max.value)
}
#' @rdname BELL
#' @export
rBELL <- function(n, mu = 1) {
  rbell(n, theta = .bell_theta_from_mean(mu))
}

#' @rdname ZIBELL
#' @export
dZIBELL <- function(x, mu = 1, sigma = 0.1, log = FALSE) {
  dzibell(x, theta = .bell_theta_from_mean(mu), pi = sigma, log = log)
}
#' @rdname ZIBELL
#' @export
pZIBELL <- function(q, mu = 1, sigma = 0.1, lower.tail = TRUE, log.p = FALSE) {
  pzibell(q, theta = .bell_theta_from_mean(mu), pi = sigma,
          lower.tail = lower.tail, log.p = log.p)
}
#' @rdname ZIBELL
#' @export
qZIBELL <- function(p, mu = 1, sigma = 0.1, lower.tail = TRUE, log.p = FALSE,
                    max.value = 1e5) {
  qzibell(p, theta = .bell_theta_from_mean(mu), pi = sigma,
          lower.tail = lower.tail, log.p = log.p, max.value = max.value)
}
#' @rdname ZIBELL
#' @export
rZIBELL <- function(n, mu = 1, sigma = 0.1) {
  rzibell(n, theta = .bell_theta_from_mean(mu), pi = sigma)
}

## ---- family builders ----

#' gamlss family for the Bell distribution
#'
#' \code{BELL()} builds a \code{gamlss.family} object for fitting the Bell
#' distribution of Castellares, Ferrari and Lemonte (2018) with \code{gamlss()},
#' in the mean parameterisation \code{mu} = mean (log link; internally
#' \eqn{\theta = W_0(\mu)}). Bell is a one-parameter family, so \code{mu} is the
#' only distribution parameter. \code{dBELL}, \code{pBELL}, \code{qBELL} and
#' \code{rBELL} are the corresponding density, distribution, quantile and
#' simulation functions in the same parameterisation.
#'
#' @param mu.link link function for the mean parameter.
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu mean parameter (\eqn{>0}).
#' @param log,log.p,lower.tail,max.value as in \code{\link{bell}}.
#' @return \code{BELL} returns a \code{gamlss.family} object.
#' @references Castellares, F., Ferrari, S. L. P. and Lemonte, A. J. (2018). On
#'   the Bell distribution and its associated regression model for count data.
#'   \emph{Applied Mathematical Modelling} 56, 172-185.
#' @seealso \code{\link{bell}} for the native \code{theta} functions.
#' @name BELL
#' @export
BELL <- function(mu.link = "log") {
  .make_bell_gamlss("dBELL", "pBELL", "BELL", 1L, mu.link, NULL)
}

#' gamlss family for the zero-inflated Bell distribution
#'
#' \code{ZIBELL()} adds a logit-link zero-inflation probability (\code{sigma}) to
#' \code{\link{BELL}}. \code{dZIBELL}, \code{pZIBELL}, \code{qZIBELL} and
#' \code{rZIBELL} are the corresponding d/p/q/r functions, with \code{sigma} the
#' zero-inflation probability.
#'
#' @inheritParams BELL
#' @param sigma.link link function for the zero-inflation parameter.
#' @param sigma zero-inflation probability \eqn{\in[0,1)}.
#' @return \code{ZIBELL} returns a \code{gamlss.family} object.
#' @seealso \code{\link{BELL}}; \code{\link{zibell}}.
#' @name ZIBELL
#' @export
ZIBELL <- function(mu.link = "log", sigma.link = "logit") {
  .make_bell_gamlss("dZIBELL", "pZIBELL", "ZIBELL", 2L, mu.link, sigma.link)
}

.make_bell_gamlss <- function(dfun, pfun, family_name, npar, mu.link, sigma.link) {
  ## Force the argument promises into locals so substitute() resolves to the
  ## link string (not an unbound symbol from the caller's frame).
  mu_link <- mu.link
  gamlss.dist::checklink("mu.link", family_name, substitute(mu_link),
                         c("log", "identity", "sqrt"))
  dfun_sym <- as.name(dfun); pfun_str <- pfun

  if (npar == 1L) {
    G.dev.incr_fn <- eval(bquote(function(y, mu, ...) {
      -2 * .(dfun_sym)(y, mu = mu, log = TRUE) }))
    dldm_fn <- eval(bquote(function(y, mu, ...)
      egpd::.egpd_nd1(y, mu, .(dfun_sym))))
    d2ldm2_fn <- eval(bquote(function(y, mu, ...) {
      dl <- egpd::.egpd_nd1(y, mu, .(dfun_sym))
      ifelse(-dl * dl < -1e-15, -dl * dl, -1e-15) }))
    rqres_expr <- eval(bquote(expression(
      rqres(pfun = .(pfun_str), type = "Discrete", ymin = 0, y = y, mu = mu))))

    fam <- list(
      family = c(family_name, "gamlss.family"),
      parameters = list(mu = TRUE), nopar = 1L, type = "Discrete",
      mu.link = mu_link,
      mu.linkfun = make.link(mu_link)$linkfun, mu.linkinv = make.link(mu_link)$linkinv,
      mu.dr = make.link(mu_link)$mu.eta,
      mu.initial = expression(mu <- rep(max(mean(y), 0.1), length(y))),
      G.dev.incr = G.dev.incr_fn, rqres = rqres_expr,
      mu.valid = function(mu) all(mu > 0),
      y.valid  = function(y) all(y >= 0),
      dldm = dldm_fn, d2ldm2 = d2ldm2_fn)

  } else {
    sigma_link <- sigma.link
    gamlss.dist::checklink("sigma.link", family_name, substitute(sigma_link),
                           c("logit", "log", "identity"))
    G.dev.incr_fn <- eval(bquote(function(y, mu, sigma, ...) {
      -2 * .(dfun_sym)(y, mu = mu, sigma = sigma, log = TRUE) }))
    dldm_fn <- eval(bquote(function(y, mu, sigma, ...)
      egpd::.egpd_nd2(y, mu, sigma, .(dfun_sym), idx = 1L)))
    dldd_fn <- eval(bquote(function(y, mu, sigma, ...)
      egpd::.egpd_nd2(y, mu, sigma, .(dfun_sym), idx = 2L)))
    d2ldm2_fn <- eval(bquote(function(y, mu, sigma, ...) {
      dl <- egpd::.egpd_nd2(y, mu, sigma, .(dfun_sym), idx = 1L)
      ifelse(-dl * dl < -1e-15, -dl * dl, -1e-15) }))
    d2ldd2_fn <- eval(bquote(function(y, mu, sigma, ...) {
      dl <- egpd::.egpd_nd2(y, mu, sigma, .(dfun_sym), idx = 2L)
      ifelse(-dl * dl < -1e-15, -dl * dl, -1e-15) }))
    d2ldmdd_fn <- eval(bquote(function(y, mu, sigma, ...)
      -egpd::.egpd_nd2(y, mu, sigma, .(dfun_sym), 1L) *
        egpd::.egpd_nd2(y, mu, sigma, .(dfun_sym), 2L)))
    rqres_expr <- eval(bquote(expression(
      rqres(pfun = .(pfun_str), type = "Discrete", ymin = 0, y = y,
            mu = mu, sigma = sigma))))

    fam <- list(
      family = c(family_name, "gamlss.family"),
      parameters = list(mu = TRUE, sigma = TRUE), nopar = 2L, type = "Discrete",
      mu.link = mu_link, sigma.link = sigma_link,
      mu.linkfun = make.link(mu_link)$linkfun, mu.linkinv = make.link(mu_link)$linkinv,
      mu.dr = make.link(mu_link)$mu.eta,
      sigma.linkfun = make.link(sigma_link)$linkfun, sigma.linkinv = make.link(sigma_link)$linkinv,
      sigma.dr = make.link(sigma_link)$mu.eta,
      mu.initial    = expression(mu    <- rep(max(mean(y), 0.1), length(y))),
      sigma.initial = expression(sigma <- rep(min(max(mean(y == 0), 0.01), 0.5), length(y))),
      G.dev.incr = G.dev.incr_fn, rqres = rqres_expr,
      mu.valid    = function(mu) all(mu > 0),
      sigma.valid = function(sigma) all(sigma >= 0 & sigma < 1),
      y.valid     = function(y) all(y >= 0),
      dldm = dldm_fn, dldd = dldd_fn,
      d2ldm2 = d2ldm2_fn, d2ldd2 = d2ldd2_fn, d2ldmdd = d2ldmdd_fn)
  }
  class(fam) <- "gamlss.family"
  fam
}
