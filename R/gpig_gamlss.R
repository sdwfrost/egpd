## gamlss family for the generalised Poisson-inverse Gaussian (GPIG; Zhu & Joe
## 2009) and zero-inflated GPIG (ZIGPIG). Because gamlss.dist has no GPIG, these
## build on the package recursion (dgpig/pgpig, R/gpig_dist.R).
##
## Mean parameterisation (consistent with the EGPD gamlss families):
##   mu    = mean            (log link)
##   sigma = a, tail exponent in (0,1)   (logit link)
##   nu    = c, down-weight  in (0,1)    (logit link)
##   tau   = pi, zero-inflation (ZIGPIG only, logit link)
## internally b = mu (1-c)^(1-a) / (a c). Derivatives of the log-likelihood are
## obtained by central differences (reusing .egpd_nd3 / .egpd_nd4), matching the
## approach used for the EGPD gamlss families.

## ---- gamlss-convention d/p/q/r (mu = mean, sigma = a, nu = c[, tau = pi]) ----

.gpig_b_from_mean <- function(mu, a, c) mu * (1 - c)^(1 - a) / (a * c)

#' @rdname GPIG
#' @export
dGPIG <- function(x, mu = 1, sigma = 0.5, nu = 0.5, log = FALSE) {
  b <- .gpig_b_from_mean(mu, sigma, nu)
  dgpig(x, a = sigma, b = b, c = nu, log = log)
}
#' @rdname GPIG
#' @export
pGPIG <- function(q, mu = 1, sigma = 0.5, nu = 0.5, lower.tail = TRUE, log.p = FALSE) {
  b <- .gpig_b_from_mean(mu, sigma, nu)
  pgpig(q, a = sigma, b = b, c = nu, lower.tail = lower.tail, log.p = log.p)
}
#' @rdname GPIG
#' @export
qGPIG <- function(p, mu = 1, sigma = 0.5, nu = 0.5, lower.tail = TRUE, log.p = FALSE,
                  max.value = 1e5) {
  b <- .gpig_b_from_mean(mu, sigma, nu)
  qgpig(p, a = sigma, b = b, c = nu, lower.tail = lower.tail, log.p = log.p,
        max.value = max.value)
}
#' @rdname GPIG
#' @export
rGPIG <- function(n, mu = 1, sigma = 0.5, nu = 0.5) {
  b <- .gpig_b_from_mean(mu, sigma, nu)
  rgpig(n, a = sigma, b = b, c = nu)
}

#' @rdname ZIGPIG
#' @export
dZIGPIG <- function(x, mu = 1, sigma = 0.5, nu = 0.5, tau = 0.1, log = FALSE) {
  b <- .gpig_b_from_mean(mu, sigma, nu)
  dzigpig(x, a = sigma, b = b, c = nu, pi = tau, log = log)
}
#' @rdname ZIGPIG
#' @export
pZIGPIG <- function(q, mu = 1, sigma = 0.5, nu = 0.5, tau = 0.1,
                    lower.tail = TRUE, log.p = FALSE) {
  b <- .gpig_b_from_mean(mu, sigma, nu)
  pzigpig(q, a = sigma, b = b, c = nu, pi = tau, lower.tail = lower.tail, log.p = log.p)
}
#' @rdname ZIGPIG
#' @export
qZIGPIG <- function(p, mu = 1, sigma = 0.5, nu = 0.5, tau = 0.1,
                    lower.tail = TRUE, log.p = FALSE, max.value = 1e5) {
  b <- .gpig_b_from_mean(mu, sigma, nu)
  qzigpig(p, a = sigma, b = b, c = nu, pi = tau, lower.tail = lower.tail,
          log.p = log.p, max.value = max.value)
}
#' @rdname ZIGPIG
#' @export
rZIGPIG <- function(n, mu = 1, sigma = 0.5, nu = 0.5, tau = 0.1) {
  b <- .gpig_b_from_mean(mu, sigma, nu)
  rzigpig(n, a = sigma, b = b, c = nu, pi = tau)
}

## ---- family builders ----
## Mirror .make_gamlss_family but with logit links / (0,1) shape parameters and
## GPIG-appropriate initial values (c must stay < 1).

#' gamlss family for the generalised Poisson-inverse Gaussian distribution
#'
#' \code{GPIG()} builds a \code{gamlss.family} object for fitting the
#' generalised Poisson-inverse Gaussian distribution of Zhu and Joe (2009) with
#' \code{gamlss()}, in the mean parameterisation \code{mu} = mean (log link),
#' \code{sigma} = tail exponent \eqn{a\in(0,1)} (logit link) and \code{nu} =
#' down-weighting \eqn{c\in(0,1)} (logit link). \code{dGPIG}, \code{pGPIG},
#' \code{qGPIG} and \code{rGPIG} are the corresponding density, distribution,
#' quantile and simulation functions in the same parameterisation.
#'
#' @param mu.link,sigma.link,nu.link link functions for the parameters.
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu mean parameter (\eqn{>0}).
#' @param sigma tail exponent \eqn{a\in(0,1)}.
#' @param nu down-weighting parameter \eqn{c\in(0,1)}.
#' @param log,log.p,lower.tail,max.value as in \code{\link{gpig}}.
#' @return \code{GPIG} returns a \code{gamlss.family} object.
#' @references Zhu, R. and Joe, H. (2009). Modelling heavy-tailed count data
#'   using a generalised Poisson-inverse Gaussian family. \emph{Statistics and
#'   Probability Letters} 79, 1695-1703.
#' @seealso \code{\link{gpig}} for the native \code{(a,b,c)} functions.
#' @name GPIG
#' @export
GPIG <- function(mu.link = "log", sigma.link = "logit", nu.link = "logit") {
  .make_gpig_gamlss("dGPIG", "pGPIG", "GPIG", 3L,
                    mu.link, sigma.link, nu.link, NULL)
}

#' gamlss family for the zero-inflated generalised Poisson-inverse Gaussian
#'
#' \code{ZIGPIG()} adds a logit-link zero-inflation probability \code{tau} to
#' \code{\link{GPIG}}. \code{dZIGPIG}, \code{pZIGPIG}, \code{qZIGPIG} and
#' \code{rZIGPIG} are the corresponding d/p/q/r functions.
#'
#' @inheritParams GPIG
#' @param tau.link link function for the zero-inflation parameter.
#' @param tau zero-inflation probability \eqn{\in[0,1)}.
#' @return \code{ZIGPIG} returns a \code{gamlss.family} object.
#' @seealso \code{\link{GPIG}}; \code{\link{zigpig}}.
#' @name ZIGPIG
#' @export
ZIGPIG <- function(mu.link = "log", sigma.link = "logit", nu.link = "logit",
                   tau.link = "logit") {
  .make_gpig_gamlss("dZIGPIG", "pZIGPIG", "ZIGPIG", 4L,
                    mu.link, sigma.link, nu.link, tau.link)
}

.make_gpig_gamlss <- function(dfun, pfun, family_name, npar,
                              mu.link, sigma.link, nu.link, tau.link) {
  ## Force the argument promises into locals so substitute() resolves to the
  ## link string (not an unbound symbol from the caller's frame).
  mu_link <- mu.link; sigma_link <- sigma.link; nu_link <- nu.link
  gamlss.dist::checklink("mu.link",    family_name, substitute(mu_link),    c("log", "identity"))
  gamlss.dist::checklink("sigma.link", family_name, substitute(sigma_link), c("logit", "log", "identity"))
  gamlss.dist::checklink("nu.link",    family_name, substitute(nu_link),    c("logit", "log", "identity"))
  dfun_sym <- as.name(dfun); pfun_str <- pfun

  if (npar == 3L) {
    G.dev.incr_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      -2 * .(dfun_sym)(y, mu = mu, sigma = sigma, nu = nu, log = TRUE)
    }))
    dldm_fn <- eval(bquote(function(y, mu, sigma, nu, ...)
      egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 1L)))
    dldd_fn <- eval(bquote(function(y, mu, sigma, nu, ...)
      egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 2L)))
    dldv_fn <- eval(bquote(function(y, mu, sigma, nu, ...)
      egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 3L)))
    d2ldm2_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      dl <- egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 1L)
      ifelse(-dl * dl < -1e-15, -dl * dl, -1e-15) }))
    d2ldd2_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      dl <- egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 2L)
      ifelse(-dl * dl < -1e-15, -dl * dl, -1e-15) }))
    d2ldv2_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      dl <- egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), idx = 3L)
      ifelse(-dl * dl < -1e-15, -dl * dl, -1e-15) }))
    d2ldmdd_fn <- eval(bquote(function(y, mu, sigma, nu, ...)
      -egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), 1L) *
        egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), 2L)))
    d2ldmdv_fn <- eval(bquote(function(y, mu, sigma, nu, ...)
      -egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), 1L) *
        egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), 3L)))
    d2ldddv_fn <- eval(bquote(function(y, mu, sigma, nu, ...)
      -egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), 2L) *
        egpd::.egpd_nd3(y, mu, sigma, nu, .(dfun_sym), 3L)))
    rqres_expr <- eval(bquote(expression(
      rqres(pfun = .(pfun_str), type = "Discrete", ymin = 0, y = y,
            mu = mu, sigma = sigma, nu = nu))))

    fam <- list(
      family = c(family_name, "gamlss.family"),
      parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE),
      nopar = 3L, type = "Discrete",
      mu.link = mu_link, sigma.link = sigma_link, nu.link = nu_link,
      mu.linkfun = make.link(mu_link)$linkfun, mu.linkinv = make.link(mu_link)$linkinv,
      mu.dr = make.link(mu_link)$mu.eta,
      sigma.linkfun = make.link(sigma_link)$linkfun, sigma.linkinv = make.link(sigma_link)$linkinv,
      sigma.dr = make.link(sigma_link)$mu.eta,
      nu.linkfun = make.link(nu_link)$linkfun, nu.linkinv = make.link(nu_link)$linkinv,
      nu.dr = make.link(nu_link)$mu.eta,
      mu.initial    = expression(mu    <- rep(max(mean(y), 0.1), length(y))),
      sigma.initial = expression(sigma <- rep(0.5, length(y))),
      nu.initial    = expression(nu    <- rep(0.7, length(y))),
      G.dev.incr = G.dev.incr_fn, rqres = rqres_expr,
      mu.valid    = function(mu) all(mu > 0),
      sigma.valid = function(sigma) all(sigma > 0 & sigma < 1),
      nu.valid    = function(nu) all(nu > 0 & nu < 1),
      y.valid     = function(y) all(y >= 0),
      dldm = dldm_fn, dldd = dldd_fn, dldv = dldv_fn,
      d2ldm2 = d2ldm2_fn, d2ldd2 = d2ldd2_fn, d2ldv2 = d2ldv2_fn,
      d2ldmdd = d2ldmdd_fn, d2ldmdv = d2ldmdv_fn, d2ldddv = d2ldddv_fn)

  } else {
    tau_link <- tau.link
    gamlss.dist::checklink("tau.link", family_name, substitute(tau_link), c("logit", "log", "identity"))
    G.dev.incr_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      -2 * .(dfun_sym)(y, mu = mu, sigma = sigma, nu = nu, tau = tau, log = TRUE) }))
    nd <- function(idx) eval(bquote(function(y, mu, sigma, nu, tau, ...)
      egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = .(idx))))
    dldm_fn <- nd(1L); dldd_fn <- nd(2L); dldv_fn <- nd(3L); dldt_fn <- nd(4L)
    d2self <- function(idx) eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      dl <- egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), idx = .(idx))
      ifelse(-dl * dl < -1e-15, -dl * dl, -1e-15) }))
    d2cross <- function(i, j) eval(bquote(function(y, mu, sigma, nu, tau, ...)
      -egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), .(i)) *
        egpd::.egpd_nd4(y, mu, sigma, nu, tau, .(dfun_sym), .(j))))
    rqres_expr <- eval(bquote(expression(
      rqres(pfun = .(pfun_str), type = "Discrete", ymin = 0, y = y,
            mu = mu, sigma = sigma, nu = nu, tau = tau))))

    fam <- list(
      family = c(family_name, "gamlss.family"),
      parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE, tau = TRUE),
      nopar = 4L, type = "Discrete",
      mu.link = mu_link, sigma.link = sigma_link, nu.link = nu_link, tau.link = tau_link,
      mu.linkfun = make.link(mu_link)$linkfun, mu.linkinv = make.link(mu_link)$linkinv,
      mu.dr = make.link(mu_link)$mu.eta,
      sigma.linkfun = make.link(sigma_link)$linkfun, sigma.linkinv = make.link(sigma_link)$linkinv,
      sigma.dr = make.link(sigma_link)$mu.eta,
      nu.linkfun = make.link(nu_link)$linkfun, nu.linkinv = make.link(nu_link)$linkinv,
      nu.dr = make.link(nu_link)$mu.eta,
      tau.linkfun = make.link(tau_link)$linkfun, tau.linkinv = make.link(tau_link)$linkinv,
      tau.dr = make.link(tau_link)$mu.eta,
      mu.initial    = expression(mu    <- rep(max(mean(y), 0.1), length(y))),
      sigma.initial = expression(sigma <- rep(0.5, length(y))),
      nu.initial    = expression(nu    <- rep(0.7, length(y))),
      tau.initial   = expression(tau   <- rep(min(max(mean(y == 0), 0.01), 0.5), length(y))),
      G.dev.incr = G.dev.incr_fn, rqres = rqres_expr,
      mu.valid    = function(mu) all(mu > 0),
      sigma.valid = function(sigma) all(sigma > 0 & sigma < 1),
      nu.valid    = function(nu) all(nu > 0 & nu < 1),
      tau.valid   = function(tau) all(tau >= 0 & tau < 1),
      y.valid     = function(y) all(y >= 0),
      dldm = dldm_fn, dldd = dldd_fn, dldv = dldv_fn, dldt = dldt_fn,
      d2ldm2 = d2self(1L), d2ldd2 = d2self(2L), d2ldv2 = d2self(3L), d2ldt2 = d2self(4L),
      d2ldmdd = d2cross(1L, 2L), d2ldmdv = d2cross(1L, 3L), d2ldmdt = d2cross(1L, 4L),
      d2ldddv = d2cross(2L, 3L), d2ldddt = d2cross(2L, 4L), d2ldvdt = d2cross(3L, 4L))
  }
  class(fam) <- "gamlss.family"
  fam
}
