## gamlss family for the Bell distribution (Castellares, Ferrari & Lemonte 2018)
## and zero-inflated Bell (ZIBell). Because gamlss.dist has no Bell, these build
## on the package functions (dbell/pbell, R/bell_dist.R).
##
## Mean parameterisation (consistent with the EGPD/GPIG gamlss families):
##   mu    = mean            (log link),   internally theta = W0(mu)
##   sigma = pi, zero-inflation (ZIBell only, logit link)
## Bell is a one-parameter family, so the plain BELL() gamlss family has only mu;
## ZIBELL() adds sigma as the zero-inflation probability. Both use analytic
## natural-scale scores and the per-observation BHHH curvature.

## theta = W0(mu) (Lambert W); vectorised via the C++ helper in src/bell.cpp.
.bell_theta_from_mean <- function(mu) bell_W0_cpp(as.numeric(mu))

## d log P(y) / d mu = (y - mu) / Var(Y). For a mixture zero, weight
## this by the posterior probability of the Bell component. The right
## derivative in sigma at zero is included (no division by sigma).
.bell_gamlss_score <- function(y, mu, sigma = 0, inflation = FALSE) {
  p <- .gpig_recycle(y = y, mu = mu, sigma = sigma)
  y <- p$y; mu <- p$mu; sigma <- p$sigma
  theta <- .bell_theta_from_mean(mu)
  valid <- is.finite(theta) & theta > 0 & is.finite(sigma) & sigma >= 0 & sigma < 1 &
    is.finite(y) & y >= 0 & y == floor(y)
  score <- if (inflation) -1 / (1 - sigma) else ((y - mu) / (1 + theta)) / mu
  zero <- which(valid & y == 0 & (inflation | sigma > 0))
  if (length(zero)) {
    lp0 <- -expm1(theta[zero])
    pi <- sigma[zero]
    score[zero] <- if (inflation) {
      -expm1(lp0) / (pi + (1 - pi) * exp(lp0))
    } else {
      -plogis(lp0 - qlogis(pi)) / (1 + theta[zero])
    }
  }
  score[!valid] <- NaN
  missing <- Reduce(`|`, lapply(p, function(x) is.na(x) & !is.nan(x)))
  score[missing] <- NA_real_
  if (any(is.nan(score)))
    warning("Bell score undefined for invalid parameters or observations outside the support.",
            call. = FALSE)
  score
}

## ---- gamlss-convention d/p/q/r (mu = mean [, sigma = pi]) ----

#' @rdname bell
#' @export
dBELL <- function(x, mu = 1, log = FALSE) {
  dbell(x, theta = .bell_theta_from_mean(mu), log = log)
}
#' @rdname bell
#' @export
pBELL <- function(q, mu = 1, lower.tail = TRUE, log.p = FALSE) {
  pbell(q, theta = .bell_theta_from_mean(mu), lower.tail = lower.tail, log.p = log.p)
}
#' @rdname bell
#' @export
qBELL <- function(p, mu = 1, lower.tail = TRUE, log.p = FALSE, max.value = 1e5) {
  qbell(p, theta = .bell_theta_from_mean(mu), lower.tail = lower.tail,
        log.p = log.p, max.value = max.value)
}
#' @rdname bell
#' @export
rBELL <- function(n, mu = 1) {
  rbell(n, theta = .bell_theta_from_mean(mu))
}

#' @rdname zibell
#' @export
dZIBELL <- function(x, mu = 1, sigma = 0.1, log = FALSE) {
  dzibell(x, theta = .bell_theta_from_mean(mu), pi = sigma, log = log)
}
#' @rdname zibell
#' @export
pZIBELL <- function(q, mu = 1, sigma = 0.1, lower.tail = TRUE, log.p = FALSE) {
  pzibell(q, theta = .bell_theta_from_mean(mu), pi = sigma,
          lower.tail = lower.tail, log.p = log.p)
}
#' @rdname zibell
#' @export
qZIBELL <- function(p, mu = 1, sigma = 0.1, lower.tail = TRUE, log.p = FALSE,
                    max.value = 1e5) {
  qzibell(p, theta = .bell_theta_from_mean(mu), pi = sigma,
          lower.tail = lower.tail, log.p = log.p, max.value = max.value)
}
#' @rdname zibell
#' @export
rZIBELL <- function(n, mu = 1, sigma = 0.1) {
  rzibell(n, theta = .bell_theta_from_mean(mu), pi = sigma)
}

## ---- family builders ----

#' @details
#' \code{BELL()} builds a \code{gamlss.family} object for fitting the Bell
#' distribution with \code{gamlss()}, in the mean parameterisation \code{mu} =
#' mean (log link; internally \eqn{\theta = W_0(\mu)}). \code{dBELL},
#' \code{pBELL}, \code{qBELL} and \code{rBELL} are the corresponding density,
#' distribution, quantile and simulation functions in the same parameterisation.
#' The natural-scale mean score is \eqn{(y-\mu)/[\mu(1+\theta)]}.
#' Curvature uses the negative outer product of scores.
#' @inheritSection EGPD1 GAMLSS namespace use
#'
#' @param mu.link link function for the mean parameter.
#' @param mu mean parameter (\eqn{>0}).
#' @rdname bell
#' @export
BELL <- function(mu.link = "log") {
  .make_bell_gamlss("dBELL", "pBELL", "BELL", 1L, mu.link, NULL)
}

#' @details
#' \code{ZIBELL()} adds a logit-link zero-inflation probability (\code{sigma}) to
#' \code{\link{BELL}}. \code{dZIBELL}, \code{pZIBELL}, \code{qZIBELL} and
#' \code{rZIBELL} are the corresponding mean-parameterised d/p/q/r functions.
#' For a zero count the Bell mean score is weighted by the conditional
#' probability that the zero arose from the Bell component.
#'
#' @param mu.link link function for the mean parameter.
#' @param sigma.link link function for the zero-inflation parameter.
#' @param mu mean parameter (\eqn{>0}).
#' @param sigma zero-inflation probability \eqn{\in[0,1)}.
#' @rdname zibell
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
  dfun_sym <- call("::", quote(egpd), as.name(dfun))
  rqres_expr <- .gamlss_rqres_expr(pfun, "Discrete", npar)
  ## gamlss rebinds derivative bodies outside our namespace.
  score <- quote(get(".bell_gamlss_score", envir = asNamespace("egpd")))

  if (npar == 1L) {
    G.dev.incr_fn <- eval(bquote(function(y, mu, ...) {
      -2 * .(dfun_sym)(y, mu = mu, log = TRUE) }))
    sm <- bquote(.(score)(y, mu))
    dldm_fn <- eval(bquote(function(y, mu, ...) .(sm)))
    d2ldm2_fn <- eval(bquote(function(y, mu, ...) pmin(-.(sm)^2, -1e-15)))

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
    sm <- bquote(.(score)(y, mu, sigma))
    sd <- bquote(.(score)(y, mu, sigma, inflation = TRUE))
    dldm_fn <- eval(bquote(function(y, mu, sigma, ...) .(sm)))
    dldd_fn <- eval(bquote(function(y, mu, sigma, ...) .(sd)))
    d2ldm2_fn <- eval(bquote(function(y, mu, sigma, ...) pmin(-.(sm)^2, -1e-15)))
    d2ldd2_fn <- eval(bquote(function(y, mu, sigma, ...) pmin(-.(sd)^2, -1e-15)))
    d2ldmdd_fn <- eval(bquote(function(y, mu, sigma, ...) -.(sm) * .(sd)))

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
