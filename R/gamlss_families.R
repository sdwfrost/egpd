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

#' @rdname EGPD1
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (kappa), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dEGPD1 <- function(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE) {
  degpd_density(x, sigma = mu, xi = sigma, kappa = nu, type = 1L, log = log)
}

#' @rdname EGPD1
#' @export
pEGPD1 <- function(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  pegpd(q, sigma = mu, xi = sigma, kappa = nu, type = 1L,
        lower.tail = lower.tail, log.p = log.p)
}

#' @rdname EGPD1
#' @export
qEGPD1 <- function(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  qegpd(p, sigma = mu, xi = sigma, kappa = nu, type = 1L,
        lower.tail = lower.tail, log.p = log.p)
}

#' @rdname EGPD1
#' @export
rEGPD1 <- function(n, mu = 1, sigma = 0.5, nu = 1) {
  regpd(n, sigma = mu, xi = sigma, kappa = nu, type = 1L)
}

#' @rdname EGPD3
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (delta), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dEGPD3 <- function(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE) {
  degpd_density(x, sigma = mu, xi = sigma, delta = nu, type = 4L, log = log)
}

#' @rdname EGPD3
#' @export
pEGPD3 <- function(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  pegpd(q, sigma = mu, xi = sigma, delta = nu, type = 4L,
        lower.tail = lower.tail, log.p = log.p)
}

#' @rdname EGPD3
#' @export
qEGPD3 <- function(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  qegpd(p, sigma = mu, xi = sigma, delta = nu, type = 4L,
        lower.tail = lower.tail, log.p = log.p)
}

#' @rdname EGPD3
#' @export
rEGPD3 <- function(n, mu = 1, sigma = 0.5, nu = 1) {
  regpd(n, sigma = mu, xi = sigma, delta = nu, type = 4L)
}

#' @rdname EGPD4
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (delta), positive
#' @param tau G-transformation parameter (kappa), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dEGPD4 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 1, log = FALSE) {
  degpd_density(x, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L, log = log)
}

#' @rdname EGPD4
#' @export
pEGPD4 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE) {
  pegpd(q, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L,
        lower.tail = lower.tail, log.p = log.p)
}

#' @rdname EGPD4
#' @export
qEGPD4 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE) {
  qegpd(p, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L,
        lower.tail = lower.tail, log.p = log.p)
}

#' @rdname EGPD4
#' @export
rEGPD4 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 1) {
  regpd(n, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L)
}

## ---------------------------------------------------------------
## d/p/q/r wrappers: Discrete EGPD
## ---------------------------------------------------------------

#' @rdname DEGPD1
#' @param x,q vector of quantiles (non-negative integers)
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (kappa), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dDEGPD1 <- function(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE) {
  ddiscegpd(x, sigma = mu, xi = sigma, kappa = nu, type = 1L, log = log)
}

#' @rdname DEGPD1
#' @export
pDEGPD1 <- function(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  pdiscegpd(q, sigma = mu, xi = sigma, kappa = nu, type = 1L,
            lower.tail = lower.tail, log.p = log.p)
}

#' @rdname DEGPD1
#' @export
qDEGPD1 <- function(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  qdiscegpd(p, sigma = mu, xi = sigma, kappa = nu, type = 1L,
            lower.tail = lower.tail, log.p = log.p)
}

#' @rdname DEGPD1
#' @export
rDEGPD1 <- function(n, mu = 1, sigma = 0.5, nu = 1) {
  rdiscegpd(n, sigma = mu, xi = sigma, kappa = nu, type = 1L)
}

#' @rdname DEGPD3
#' @param x,q vector of quantiles (non-negative integers)
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (delta), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dDEGPD3 <- function(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE) {
  ddiscegpd(x, sigma = mu, xi = sigma, delta = nu, type = 4L, log = log)
}

#' @rdname DEGPD3
#' @export
pDEGPD3 <- function(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  pdiscegpd(q, sigma = mu, xi = sigma, delta = nu, type = 4L,
            lower.tail = lower.tail, log.p = log.p)
}

#' @rdname DEGPD3
#' @export
qDEGPD3 <- function(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  qdiscegpd(p, sigma = mu, xi = sigma, delta = nu, type = 4L,
            lower.tail = lower.tail, log.p = log.p)
}

#' @rdname DEGPD3
#' @export
rDEGPD3 <- function(n, mu = 1, sigma = 0.5, nu = 1) {
  rdiscegpd(n, sigma = mu, xi = sigma, delta = nu, type = 4L)
}

#' @rdname DEGPD4
#' @param x,q vector of quantiles (non-negative integers)
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (delta), positive
#' @param tau G-transformation parameter (kappa), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dDEGPD4 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 1, log = FALSE) {
  ddiscegpd(x, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L, log = log)
}

#' @rdname DEGPD4
#' @export
pDEGPD4 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE) {
  pdiscegpd(q, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L,
            lower.tail = lower.tail, log.p = log.p)
}

#' @rdname DEGPD4
#' @export
qDEGPD4 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE) {
  qdiscegpd(p, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L,
            lower.tail = lower.tail, log.p = log.p)
}

#' @rdname DEGPD4
#' @export
rDEGPD4 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 1) {
  rdiscegpd(n, sigma = mu, xi = sigma, delta = nu, kappa = tau, type = 5L)
}

## ---------------------------------------------------------------
## d/p/q/r wrappers: Zero-Inflated Continuous EGPD
## ---------------------------------------------------------------

#' @rdname ZIEGPD1
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (kappa), positive
#' @param tau zero-inflation probability (pi), in (0,1)
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dZIEGPD1 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE) {
  dziegpd(x, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L, log = log)
}

#' @rdname ZIEGPD1
#' @export
pZIEGPD1 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  pziegpd(q, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L,
          lower.tail = lower.tail, log.p = log.p)
}

#' @rdname ZIEGPD1
#' @export
qZIEGPD1 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  qziegpd(p, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L,
          lower.tail = lower.tail, log.p = log.p)
}

#' @rdname ZIEGPD1
#' @export
rZIEGPD1 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1) {
  rziegpd(n, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L)
}

#' @rdname ZIEGPD3
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (delta), positive
#' @param tau zero-inflation probability (pi), in (0,1)
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dZIEGPD3 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE) {
  dziegpd(x, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L, log = log)
}

#' @rdname ZIEGPD3
#' @export
pZIEGPD3 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  pziegpd(q, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L,
          lower.tail = lower.tail, log.p = log.p)
}

#' @rdname ZIEGPD3
#' @export
qZIEGPD3 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  qziegpd(p, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L,
          lower.tail = lower.tail, log.p = log.p)
}

#' @rdname ZIEGPD3
#' @export
rZIEGPD3 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1) {
  rziegpd(n, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L)
}

## ---------------------------------------------------------------
## d/p/q/r wrappers: Zero-Inflated Discrete EGPD
## ---------------------------------------------------------------

#' @rdname ZIDEGPD1
#' @param x,q vector of quantiles (non-negative integers)
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (kappa), positive
#' @param tau zero-inflation probability (pi), in (0,1)
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dZIDEGPD1 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE) {
  dzidiscegpd(x, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L, log = log)
}

#' @rdname ZIDEGPD1
#' @export
pZIDEGPD1 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  pzidiscegpd(q, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L,
              lower.tail = lower.tail, log.p = log.p)
}

#' @rdname ZIDEGPD1
#' @export
qZIDEGPD1 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  qzidiscegpd(p, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L,
              lower.tail = lower.tail, log.p = log.p)
}

#' @rdname ZIDEGPD1
#' @export
rZIDEGPD1 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1) {
  rzidiscegpd(n, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 1L)
}

#' @rdname ZIDEGPD3
#' @param x,q vector of quantiles (non-negative integers)
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (delta), positive
#' @param tau zero-inflation probability (pi), in (0,1)
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dZIDEGPD3 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE) {
  dzidiscegpd(x, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L, log = log)
}

#' @rdname ZIDEGPD3
#' @export
pZIDEGPD3 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  pzidiscegpd(q, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L,
              lower.tail = lower.tail, log.p = log.p)
}

#' @rdname ZIDEGPD3
#' @export
qZIDEGPD3 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  qzidiscegpd(p, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L,
              lower.tail = lower.tail, log.p = log.p)
}

#' @rdname ZIDEGPD3
#' @export
rZIDEGPD3 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1) {
  rzidiscegpd(n, pi = tau, sigma = mu, xi = sigma, delta = nu, type = 4L)
}


## ===============================================================
## gamlss Family Constructors
## ===============================================================

## Numerical derivatives of a log-density with respect to one gamlss
## parameter.  dfn is the density function (resolved at family construction
## time); these package-level helpers are called from the family derivative
## functions.
##
## Central differences are used where both probes give a finite log-density.
## When a probe leaves the parameter domain (for example pi - h < 0 for a
## zero-inflation probability, or a probe that moves y outside the support)
## a second-order one-sided difference is used instead, and if neither probe
## is usable the step is reduced.  Observations whose own log-density is not
## finite (zero likelihood or invalid parameters) get NaN with a warning; NA
## inputs give NA.  No non-finite derivative is replaced by zero.
.egpd_nd_core <- function(ld, par, h) {
  n <- length(par)
  rows <- seq_len(n)
  ld0 <- ld(par, rows)
  dl <- rep(NaN, n)
  fin0 <- which(is.finite(ld0))
  failed <- 0L
  if (length(fin0)) {
    hh <- h[fin0]
    ldp <- suppressWarnings(ld(par[fin0] + hh, fin0))
    ldm <- suppressWarnings(ld(par[fin0] - hh, fin0))
    cen <- is.finite(ldp) & is.finite(ldm)
    dl[fin0[cen]] <- (ldp[cen] - ldm[cen]) / (2 * hh[cen])
    for (i in fin0[!cen]) {
      dl[i] <- .egpd_nd_boundary(ld, par[i], h[i], ld0[i], i)
      if (!is.finite(dl[i])) failed <- failed + 1L
    }
  }
  na_only <- is.na(ld0) & !is.nan(ld0)
  undefined <- sum(!is.finite(ld0) & !na_only) + failed
  if (undefined) {
    warning(sprintf(paste0(
      "derivative undefined for %d observation(s) whose log-density is not ",
      "finite near the current parameters; NaN returned"), undefined),
      call. = FALSE)
  }
  dl[which(na_only)] <- NA_real_
  dl[!is.finite(dl) & !na_only] <- NaN
  dl
}

## Derivative at a point whose central stencil leaves the domain (a parameter
## boundary or a support edge lies within h).  Steps h, h/10, ... are tried,
## using central differences when both probes are usable and second-order
## one-sided differences otherwise; the estimate at which successive steps
## agree best is returned.
.egpd_nd_boundary <- function(ld, p, h, f0, i) {
  est <- numeric(0)
  for (k in 0:7) {
    hk <- h * 10^(-k)
    fp <- suppressWarnings(ld(p + hk, i))
    fm <- suppressWarnings(ld(p - hk, i))
    d <- if (is.finite(fp) && is.finite(fm)) {
      (fp - fm) / (2 * hk)
    } else if (is.finite(fp)) {
      fpp <- suppressWarnings(ld(p + 2 * hk, i))
      if (is.finite(fpp)) (-3 * f0 + 4 * fp - fpp) / (2 * hk) else (fp - f0) / hk
    } else if (is.finite(fm)) {
      fmm <- suppressWarnings(ld(p - 2 * hk, i))
      if (is.finite(fmm)) (3 * f0 - 4 * fm + fmm) / (2 * hk) else (f0 - fm) / hk
    } else {
      NA_real_
    }
    if (!is.finite(d)) next
    est <- c(est, d)
    m <- length(est)
    if (m >= 2L && abs(est[m] - est[m - 1L]) <= 1e-8 * max(1, abs(est[m]))) return(est[m])
  }
  m <- length(est)
  if (!m) return(NaN)
  if (m == 1L) return(est)
  est[which.min(abs(diff(est))) + 1L]
}

.egpd_nd_step <- function(par_val, idx, h_scale) {
  h <- pmax(abs(par_val) * h_scale, 1e-8)
  if (idx == 2L) {
    ## sigma is the GPD shape xi in the EGPD families: crossing zero is valid
    h <- ifelse(abs(par_val) < 1e-4, 1e-5, h)
  }
  h[!is.finite(h)] <- 1e-8
  h
}

.egpd_nd_generic <- function(y, params, dfn, idx, h_scale, xi_index = 2L) {
  lens <- c(length(y), lengths(params))
  if (min(lens) == 0L) return(numeric(0))
  n <- max(lens)
  y <- rep_len(y, n)
  params <- lapply(params, rep_len, length.out = n)
  par_val <- as.numeric(params[[idx]])
  h <- .egpd_nd_step(par_val, if (idx == xi_index) 2L else 0L, h_scale)
  ld <- function(p, i) {
    args <- lapply(params, function(v) v[i])
    args[[idx]] <- p
    do.call(dfn, c(list(y[i]), args, list(log = TRUE)))
  }
  .egpd_nd_core(ld, par_val, h)
}

.egpd_nd3 <- function(y, mu, sigma, nu, dfn, idx, h_scale = 1e-4) {
  .egpd_nd_generic(y, list(mu = mu, sigma = sigma, nu = nu), dfn, idx, h_scale)
}

.egpd_nd4 <- function(y, mu, sigma, nu, tau, dfn, idx, h_scale = 1e-4) {
  .egpd_nd_generic(y, list(mu = mu, sigma = sigma, nu = nu, tau = tau), dfn, idx, h_scale)
}

## One- and two-parameter variants for densities without a shape xi.
.egpd_nd1 <- function(y, mu, dfn, h_scale = 1e-4) {
  .egpd_nd_generic(y, list(mu = mu), dfn, 1L, h_scale, xi_index = 0L)
}

.egpd_nd2 <- function(y, mu, sigma, dfn, idx, h_scale = 1e-4) {
  .egpd_nd_generic(y, list(mu = mu, sigma = sigma), dfn, idx, h_scale, xi_index = 0L)
}

## gamlss() and get.rqres() eval the stored expression in different frames.
## Neither has access to the constructor's bindings, nor need egpd be attached.
.gamlss_rqres_expr <- function(pfun, type, npar) {
  nms <- c("mu", "sigma", "nu", "tau")[seq_len(npar)]
  params <- stats::setNames(lapply(nms, as.name), nms)
  par <- as.call(c(list(quote(list)), params))
  support <- switch(type, Continuous = "continuous", Discrete = "count", Mixed = "mixed")
  cdf <- call("::", quote(egpd), as.name(pfun))
  as.expression(list(bquote(get(".gamlss_rqres", envir = asNamespace("egpd"))(
    y, .(par), .(cdf), .(support)))))
}

.gamlss_rqres <- function(y, par, pfun, support) {
  par <- lapply(par, rep_len, length(y))
  cdf <- function(q, i, lower) {
    do.call(pfun, c(list(q = q, lower.tail = lower, log.p = TRUE), lapply(par, `[`, i)))
  }
  ## GAMLSS randomizes a mixed zero using the uniform for that observation's
  ## row. The promise is not forced for continuous samples without jumps.
  .rqres_logtails(y, cdf, support, stats::runif(length(y)))
}

## gamlss also evaluates body(G.dev.incr) and rebinds derivative bodies in
## its own frames: bake qualified density references into those bodies.
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

  dfun_sym <- call("::", quote(egpd), as.name(dfun))
  rqres_expr <- .gamlss_rqres_expr(pfun, dist_type, npar)
  sigma_valid_fn <- if (identical(sigma_link, "identity")) {
    function(sigma) all(is.finite(sigma))
  } else {
    function(sigma) all(sigma > 0)
  }

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

      mu.initial    = expression({ yk <- y[!is.na(y)]; s0 <- if (length(yk) >= 2L) stats::sd(yk) else 1; mu <- rep(max(s0, 0.1), length(y)) }),
      sigma.initial = expression(sigma <- rep(0.3, length(y))),
      nu.initial    = expression(nu    <- rep(1, length(y))),

      G.dev.incr = G.dev.incr_fn,
      rqres      = rqres_expr,

      mu.valid    = function(mu) all(mu > 0),
      sigma.valid = sigma_valid_fn,
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

      mu.initial    = expression({ yk <- y[!is.na(y)]; s0 <- if (length(yk) >= 2L) stats::sd(yk) else 1; mu <- rep(max(s0, 0.1), length(y)) }),
      sigma.initial = expression(sigma <- rep(0.3, length(y))),
      nu.initial    = expression(nu    <- rep(1, length(y))),
      tau.initial   = expression(tau   <- rep(0.1, length(y))),

      G.dev.incr = G.dev.incr_fn,
      rqres      = rqres_expr,

      mu.valid    = function(mu) all(mu > 0),
      sigma.valid = sigma_valid_fn,
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


## Starting values for the EGPD gamlss families.
##
## The (sigma, xi, kappa) likelihood surface of the EGPD families has long,
## flat ridges, along which the coordinate-wise RS algorithm moves slowly and
## may stop early with converged = TRUE.  The families therefore start from
## the maximum likelihood estimate of the intercept-only model (fitted to at
## most `max_n` evenly spaced observations plus the smallest and largest),
## falling back to simple moment starts when that optimisation fails or the
## estimate gives a non-finite log-likelihood on the full data.  Returns
## c(mu, sigma, nu[, tau]) on the gamlss parameter scale.
##
## gamlssML() evaluates each *.initial expression in a fresh environment, so
## the most recent result is cached at package level, keyed on the inputs and
## the full response (summary statistics of different responses can coincide).
.egpd_start_cache <- new.env(parent = emptyenv())

.egpd_gamlss_start <- function(y, family, sigma_link = "log", max_n = 2000L,
                               nu_link = "log", tau_link = "log") {
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  ckey <- list(family, sigma_link, max_n)
  hit <- .egpd_start_cache$last
  if (!is.null(hit) && identical(hit$key, ckey) && identical(hit$y, y)) {
    ## reuse only a start that gives every observation positive likelihood
    dfun <- get(paste0("d", family), envir = asNamespace("egpd"), mode = "function")
    ld <- suppressWarnings(do.call(dfun, c(list(y), as.list(hit$est), list(log = TRUE))))
    if (all(is.finite(ld))) return(.egpd_start_in_link(hit$est, nu_link, tau_link))
  }
  ycache <- y
  key <- tolower(family)
  info <- .egpd_score_info(key)
  if (!length(y)) y <- 1
  ypos <- y[y > 0]
  m0 <- mean(y)
  s0 <- if (info$zi) {
    if (length(ypos) >= 2L) stats::sd(ypos) else 1
  } else {
    if (length(y) >= 2L) stats::sd(y) else 1
  }
  half_nu <- info$type %in% c(4L, 5L)
  default <- c(mu = max(s0, 0.1), sigma = 0.3,
               nu = if (half_nu) max(0.5 * m0, 0.1) else max(m0, 0.1))
  if (info$npar == 4L) {
    default <- c(default, tau = if (info$zi) min(max(mean(y == 0), 0.01), 0.99) else max(m0, 0.1))
  }
  y_full <- y
  if (length(y) > max_n) {
    y <- y[sort(unique(c(round(seq(1, length(y), length.out = max_n)),
                         which.min(y), which.max(y))))]
  }
  dfun <- get(paste0("d", family), envir = asNamespace("egpd"), mode = "function")
  xi_free <- identical(sigma_link, "identity")
  to_par <- function(th) {
    p <- exp(th)
    if (xi_free) p[2] <- th[2]
    if (info$zi) p[4] <- stats::plogis(th[4])
    p
  }
  th0 <- log(default)
  if (xi_free) th0[2] <- default[["sigma"]]
  if (info$zi) th0[4] <- stats::qlogis(default[["tau"]])
  pnames <- names(default)
  nll_par <- function(p, yy) {
    if (any(!is.finite(p))) return(1e100)
    args <- stats::setNames(as.list(p), pnames)
    ld <- suppressWarnings(do.call(dfun, c(list(yy), args, list(log = TRUE))))
    if (any(!is.finite(ld))) return(1e100)
    -sum(ld)
  }
  nll <- function(th) nll_par(to_par(th), y)
  ## with xi free (identity link) the support can end just beyond max(y),
  ## where finite-difference gradients straddle the boundary: use the analytic
  ## gradient from the score engine there
  gr_par <- function(th, yy) {
    p <- to_par(th)
    tau <- if (info$npar == 4L) p[[4]] else NULL
    g <- vapply(seq_along(p), function(i) {
      sum(.egpd_score(yy, p[[1]], p[[2]], p[[3]], tau, key, i))
    }, numeric(1))
    jac <- p
    if (xi_free) jac[2] <- 1
    if (info$zi) jac[4] <- p[[4]] * (1 - p[[4]])
    out <- -g * jac
    if (!all(is.finite(out))) stop("non-finite gradient")
    out
  }
  fit_from <- function(th, yy, maxit) {
    fn <- function(t) nll_par(to_par(t), yy)
    ctrl <- list(maxit = maxit, reltol = 1e-10)
    opt <- NULL
    if (xi_free) {
      opt <- tryCatch(
        suppressWarnings(stats::optim(th, fn, function(t) gr_par(t, yy), method = "BFGS", control = ctrl)),
        error = function(e) NULL
      )
    }
    if (is.null(opt)) {
      opt <- tryCatch(suppressWarnings(stats::optim(th, fn, method = "BFGS", control = ctrl)),
                      error = function(e) NULL)
    }
    if (!is.null(opt) && all(is.finite(opt$par)) && is.finite(opt$value) && opt$value < 1e100) opt else NULL
  }
  ## the carrier-shape profile can be multimodal: try a few shape starts
  starts <- list(th0)
  for (sh in list(c(1, 1), c(3, 2.5))) {
    th <- th0
    th[3] <- log(sh[1])
    if (info$type == 5L) th[4] <- log(sh[2])
    starts[[length(starts) + 1L]] <- th
  }
  best <- NULL
  for (th in starts) {
    opt <- fit_from(th, y, 500L)
    if (!is.null(opt) && (is.null(best) || opt$value < best$value)) best <- opt
  }
  est <- default
  if (!is.null(best) && best$value <= nll(th0)) {
    est <- stats::setNames(to_par(best$par), pnames)
    ## a fitted bounded tail puts the endpoint close to max(y), where a start
    ## that is not stationary for the full data can send the first gamlss step
    ## past max(y): polish such estimates on all observations
    if (xi_free && est[["sigma"]] < 0 && length(y_full) > length(y)) {
      pol <- fit_from(best$par, y_full, 200L)
      if (!is.null(pol) && pol$value < nll_par(est, y_full)) {
        est <- stats::setNames(to_par(pol$par), pnames)
      }
    }
  }
  ## the GPD scale is clamped to a range relative to the data's scale (a fixed
  ## cap would move the endpoint of a bounded tail inside data on a large scale)
  sc <- max(s0, m0, 1)
  est[["mu"]] <- min(max(est[["mu"]], 1e-8 * sc), 1e8 * sc)
  est[["sigma"]] <- if (identical(sigma_link, "log")) min(max(est[["sigma"]], 1e-4), 1e3) else est[["sigma"]]
  est[["nu"]] <- min(max(est[["nu"]], 1e-3), 1e3)
  if (info$npar == 4L) {
    est[["tau"]] <- if (info$zi) min(max(est[["tau"]], 1e-4), 1 - 1e-4) else min(max(est[["tau"]], 1e-3), 1e3)
  }
  ## the start must give every observation positive likelihood (e.g. lie
  ## beyond max(y) for a bounded tail), or gamlss() fails at iteration 1
  if (nll_par(est, y_full) >= 1e100) est <- default
  .egpd_start_cache$last <- list(key = ckey, y = ycache, est = est)
  .egpd_start_in_link(est, nu_link, tau_link)
}

## Keep a start inside the range of its link.  The carrier shapes are clamped
## to [1e-3, 1e3] above, which a "logit" nu or tau link (range (0, 1)) does
## not contain: gamlss() would stop with "Value ... out of range".
.egpd_start_in_link <- function(est, nu_link, tau_link) {
  in01 <- function(v) min(max(v, 1e-4), 1 - 1e-4)
  if (identical(nu_link, "logit")) est[["nu"]] <- in01(est[["nu"]])
  if (identical(tau_link, "logit") && "tau" %in% names(est)) est[["tau"]] <- in01(est[["tau"]])
  est
}

## Install the initial-value expressions (evaluated by gamlss() in its own
## frame, or by gamlssML() in a list environment) on an EGPD gamlss family.
.egpd_set_initial <- function(fam) {
  fname <- fam$family[1]
  slink <- fam$sigma.link
  nlink <- fam$nu.link
  tlink <- if (is.null(fam$tau.link)) "log" else fam$tau.link
  start_call <- bquote(egpd::.egpd_gamlss_start(y, .(fname), sigma_link = .(slink),
                                                nu_link = .(nlink), tau_link = .(tlink)))
  get_start <- bquote(if (!exists(".egpd_start", inherits = FALSE)) .egpd_start <- .(start_call))
  pars <- names(fam$parameters)
  for (p in pars) {
    fam[[paste0(p, ".initial")]] <- as.expression(bquote({
      .(get_start)
      .(as.name(p)) <- rep(.egpd_start[[.(p)]], length(y))
    }))
  }
  fam
}

## ---------------------------------------------------------------
## Continuous EGPD families
## ---------------------------------------------------------------

#' gamlss Family for Continuous EGPD Model 1
#'
#' The \code{EGPD1()} function defines a \code{gamlss.family} distribution for
#' fitting a continuous Extended Generalized Pareto Distribution with model 1
#' (power G-transformation: \eqn{G(u) = u^\kappa}). The functions
#' \code{dEGPD1}, \code{pEGPD1}, \code{qEGPD1} and \code{rEGPD1} give the
#' density, distribution function, quantile function and random generation for
#' this distribution, using the gamlss parameter names.
#'
#' The gamlss parameters map to EGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @section GAMLSS namespace use:
#' The families can be passed directly to \code{gamlss::gamlss()}, e.g.
#' \code{family = egpd::DEGPD1()}, without attaching \pkg{egpd}.
#' For families defined by \pkg{egpd}, fitted quantile residuals and
#' \code{gamlss::get.rqres()} use log CDFs and, near one, direct log survivors.
#' Counts randomize uniformly between the two jump endpoints; mixed families
#' randomize only at the zero atom. Interpolation and normal quantiles use
#' log probabilities, without clipping small tails. Precision is limited by
#' the family's log tails; a true zero tail probability gives an infinite
#' normal quantile.
#' Some separate GAMLSS helpers, including \code{gamlssML()},
#' \code{devianceIncr()}, \code{predictAll()} at new data and the likelihood
#' used by \code{vcov()}, look up the family or density by bare name instead
#' of using the family callbacks.
#' Attaching \pkg{egpd} can resolve those lookups but does not correct a
#' GAMLSS version that reconstructs fixed parameters with default links in
#' \code{vcov()}. The \pkg{egpd} \code{\link{predictive_coverage}} method
#' instead computes its own held-out conditional information. It does not
#' change upstream GAMLSS helpers, the re-exported \pkg{gamlss.dist} families
#' \code{PIG}/\code{ZIPIG}, or residual expressions in previously saved fits.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#'
#' @return \code{EGPD1()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dEGPD1} gives the density, \code{pEGPD1} gives the distribution
#'   function, \code{qEGPD1} gives the quantile function, and \code{rEGPD1}
#'   generates random deviates.
#' @references
#' Abbas, A., Ahmad, T. and Ahmad, I. (2025).
#' Modeling zero-inflated precipitation extremes.
#' \emph{Communications in Statistics - Simulation and Computation}, 1--17.
#' \doi{10.1080/03610918.2025.2585398}
#'
#' Naveau, P., Huser, R., Ribereau, P. and Hannart, A. (2016).
#' Modeling jointly low, moderate, and heavy rainfall intensities without a
#' threshold selection. \emph{Water Resources Research}, 52(4), 2753--2769.
#' \doi{10.1002/2015WR018552}
#' @seealso \code{\link{egpd_bamlss}}, \code{\link{degpd_density}},
#'   \code{\link{pegpd}}
#' @examples
#' \dontrun{
#'   library(gamlss)
#'   set.seed(1)
#'   y <- rEGPD1(500, mu = 2, sigma = 0.2, nu = 1.5)
#'   fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
#'                 family = EGPD1())
#' }
#' @order 1
#' @export
EGPD1 <- function(mu.link = "log", sigma.link = "log", nu.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dEGPD1", pfun = "pEGPD1",
    family_name = "EGPD1", dist_type = "Continuous",
    npar = 3,
    default_links = c(mu = "log", sigma = "log", nu = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link
  )
  fam <- .override_derivs3(fam, "egpd1")
  fam <- .egpd_set_initial(fam)
  fam
}

#' gamlss Family for Continuous EGPD Model 3
#'
#' The \code{EGPD3()} function defines a \code{gamlss.family} distribution for
#' fitting a continuous Extended Generalized Pareto Distribution with model 3
#' (incomplete beta G-transformation, \code{type = 4}). The functions
#' \code{dEGPD3}, \code{pEGPD3}, \code{qEGPD3} and \code{rEGPD3} give the
#' density, distribution function, quantile function and random generation for
#' this distribution, using the gamlss parameter names.
#'
#' The gamlss parameters map to EGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = delta (G-transformation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#'
#' @return \code{EGPD3()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dEGPD3} gives the density, \code{pEGPD3} gives the distribution
#'   function, \code{qEGPD3} gives the quantile function, and \code{rEGPD3}
#'   generates random deviates.
#' @references
#' Abbas, A., Ahmad, T. and Ahmad, I. (2025).
#' Modeling zero-inflated precipitation extremes.
#' \emph{Communications in Statistics - Simulation and Computation}, 1--17.
#' \doi{10.1080/03610918.2025.2585398}
#'
#' Naveau, P., Huser, R., Ribereau, P. and Hannart, A. (2016).
#' Modeling jointly low, moderate, and heavy rainfall intensities without a
#' threshold selection. \emph{Water Resources Research}, 52(4), 2753--2769.
#' \doi{10.1002/2015WR018552}
#' @seealso \code{\link{EGPD1}}, \code{\link{EGPD4}},
#'   \code{\link{egpd_bamlss}}
#' @order 1
#' @export
EGPD3 <- function(mu.link = "log", sigma.link = "log", nu.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dEGPD3", pfun = "pEGPD3",
    family_name = "EGPD3", dist_type = "Continuous",
    npar = 3,
    default_links = c(mu = "log", sigma = "log", nu = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link
  )
  fam <- .override_derivs3(fam, "egpd3")
  fam <- .egpd_set_initial(fam)
  fam
}

#' gamlss Family for Continuous EGPD Model 4
#'
#' The \code{EGPD4()} function defines a \code{gamlss.family} distribution for
#' fitting a continuous Extended Generalized Pareto Distribution with model 4
#' (power-beta G-transformation, \code{type = 5}). This is a 4-parameter
#' family. The functions \code{dEGPD4}, \code{pEGPD4}, \code{qEGPD4} and
#' \code{rEGPD4} give the density, distribution function, quantile function
#' and random generation for this distribution, using the gamlss parameter
#' names.
#'
#' The gamlss parameters map to EGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = delta (G-transformation), \code{tau} = kappa (G-transformation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"log"})
#'
#' @return \code{EGPD4()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dEGPD4} gives the density, \code{pEGPD4} gives the distribution
#'   function, \code{qEGPD4} gives the quantile function, and \code{rEGPD4}
#'   generates random deviates.
#' @references
#' Abbas, A., Ahmad, T. and Ahmad, I. (2025).
#' Modeling zero-inflated precipitation extremes.
#' \emph{Communications in Statistics - Simulation and Computation}, 1--17.
#' \doi{10.1080/03610918.2025.2585398}
#'
#' Naveau, P., Huser, R., Ribereau, P. and Hannart, A. (2016).
#' Modeling jointly low, moderate, and heavy rainfall intensities without a
#' threshold selection. \emph{Water Resources Research}, 52(4), 2753--2769.
#' \doi{10.1002/2015WR018552}
#' @seealso \code{\link{EGPD1}}, \code{\link{EGPD3}},
#'   \code{\link{egpd_bamlss}}
#' @order 1
#' @export
EGPD4 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dEGPD4", pfun = "pEGPD4",
    family_name = "EGPD4", dist_type = "Continuous",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam <- .override_derivs4(fam, "egpd4")
  fam <- .egpd_set_initial(fam)
  fam$tau.valid     <- function(tau) all(tau > 0)
  fam
}

## ---------------------------------------------------------------
## Discrete EGPD families
## ---------------------------------------------------------------

#' gamlss Family for Discrete EGPD Model 1
#'
#' The \code{DEGPD1()} function defines a \code{gamlss.family} distribution
#' for fitting a discrete Extended Generalized Pareto Distribution with model
#' 1 (power G-transformation: \eqn{G(u) = u^\kappa}). The functions
#' \code{dDEGPD1}, \code{pDEGPD1}, \code{qDEGPD1} and \code{rDEGPD1} give the
#' PMF, distribution function, quantile function and random generation for
#' this distribution, using the gamlss parameter names.
#'
#' The gamlss parameters map to DEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#'
#' @return \code{DEGPD1()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dDEGPD1} gives the PMF, \code{pDEGPD1} gives the distribution
#'   function, \code{qDEGPD1} gives the quantile function, and \code{rDEGPD1}
#'   generates random deviates.
#' @references
#' Ahmad, T. and Arshad, I. A. (2026).
#' New flexible versions of the extended generalized Pareto distribution
#' for count data. \emph{Journal of Applied Statistics}, 1--21.
#' \doi{10.1080/02664763.2026.2634785} (preprint: arXiv:2409.18719)
#'
#' Ahmad, T., Gaetan, C. and Naveau, P. (2025).
#' An extended generalized Pareto regression model for count data.
#' \emph{Statistical Modelling}, 25(5), 416--431.
#' \doi{10.1177/1471082X241266729}
#'
#' Naveau, P., Huser, R., Ribereau, P. and Hannart, A. (2016).
#' Modeling jointly low, moderate, and heavy rainfall intensities without a
#' threshold selection. \emph{Water Resources Research}, 52(4), 2753--2769.
#' \doi{10.1002/2015WR018552}
#' @seealso \code{\link{DEGPD3}}, \code{\link{DEGPD4}},
#'   \code{\link{degpd_bamlss}}, \code{\link{ddiscegpd}}
#' @examples
#' \dontrun{
#'   library(gamlss)
#'   set.seed(1)
#'   y <- rDEGPD1(500, mu = 3, sigma = 0.15, nu = 2)
#'   fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
#'                 family = DEGPD1())
#' }
#' @order 1
#' @export
DEGPD1 <- function(mu.link = "log", sigma.link = "log", nu.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dDEGPD1", pfun = "pDEGPD1",
    family_name = "DEGPD1", dist_type = "Discrete",
    npar = 3,
    default_links = c(mu = "log", sigma = "log", nu = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link
  )
  fam <- .override_derivs3(fam, "degpd1")
  fam <- .egpd_set_initial(fam)
  fam
}

#' gamlss Family for Discrete EGPD Model 3
#'
#' The \code{DEGPD3()} function defines a \code{gamlss.family} distribution
#' for fitting a discrete Extended Generalized Pareto Distribution with model
#' 3 (incomplete beta G-transformation, \code{type = 4}). The functions
#' \code{dDEGPD3}, \code{pDEGPD3}, \code{qDEGPD3} and \code{rDEGPD3} give the
#' PMF, distribution function, quantile function and random generation for
#' this distribution, using the gamlss parameter names.
#'
#' The gamlss parameters map to DEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = delta (G-transformation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#'
#' @return \code{DEGPD3()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dDEGPD3} gives the PMF, \code{pDEGPD3} gives the distribution
#'   function, \code{qDEGPD3} gives the quantile function, and \code{rDEGPD3}
#'   generates random deviates.
#' @references
#' Ahmad, T. and Arshad, I. A. (2026).
#' New flexible versions of the extended generalized Pareto distribution
#' for count data. \emph{Journal of Applied Statistics}, 1--21.
#' \doi{10.1080/02664763.2026.2634785} (preprint: arXiv:2409.18719)
#'
#' Ahmad, T., Gaetan, C. and Naveau, P. (2025).
#' An extended generalized Pareto regression model for count data.
#' \emph{Statistical Modelling}, 25(5), 416--431.
#' \doi{10.1177/1471082X241266729}
#'
#' Naveau, P., Huser, R., Ribereau, P. and Hannart, A. (2016).
#' Modeling jointly low, moderate, and heavy rainfall intensities without a
#' threshold selection. \emph{Water Resources Research}, 52(4), 2753--2769.
#' \doi{10.1002/2015WR018552}
#' @seealso \code{\link{DEGPD1}}, \code{\link{DEGPD4}},
#'   \code{\link{degpd_bamlss}}
#' @order 1
#' @export
DEGPD3 <- function(mu.link = "log", sigma.link = "log", nu.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dDEGPD3", pfun = "pDEGPD3",
    family_name = "DEGPD3", dist_type = "Discrete",
    npar = 3,
    default_links = c(mu = "log", sigma = "log", nu = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link
  )
  fam <- .override_derivs3(fam, "degpd3")
  fam <- .egpd_set_initial(fam)
  fam
}

#' gamlss Family for Discrete EGPD Model 4
#'
#' The \code{DEGPD4()} function defines a \code{gamlss.family} distribution
#' for fitting a discrete Extended Generalized Pareto Distribution with model
#' 4 (power-beta G-transformation, \code{type = 5}). This is a 4-parameter
#' family. The functions \code{dDEGPD4}, \code{pDEGPD4}, \code{qDEGPD4} and
#' \code{rDEGPD4} give the PMF, distribution function, quantile function and
#' random generation for this distribution, using the gamlss parameter names.
#'
#' The gamlss parameters map to DEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = delta (G-transformation), \code{tau} = kappa (G-transformation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"log"})
#'
#' @return \code{DEGPD4()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dDEGPD4} gives the PMF, \code{pDEGPD4} gives the distribution
#'   function, \code{qDEGPD4} gives the quantile function, and \code{rDEGPD4}
#'   generates random deviates.
#' @references
#' Ahmad, T. and Arshad, I. A. (2026).
#' New flexible versions of the extended generalized Pareto distribution
#' for count data. \emph{Journal of Applied Statistics}, 1--21.
#' \doi{10.1080/02664763.2026.2634785} (preprint: arXiv:2409.18719)
#'
#' Ahmad, T., Gaetan, C. and Naveau, P. (2025).
#' An extended generalized Pareto regression model for count data.
#' \emph{Statistical Modelling}, 25(5), 416--431.
#' \doi{10.1177/1471082X241266729}
#'
#' Naveau, P., Huser, R., Ribereau, P. and Hannart, A. (2016).
#' Modeling jointly low, moderate, and heavy rainfall intensities without a
#' threshold selection. \emph{Water Resources Research}, 52(4), 2753--2769.
#' \doi{10.1002/2015WR018552}
#' @seealso \code{\link{DEGPD1}}, \code{\link{DEGPD3}},
#'   \code{\link{degpd_bamlss}}
#' @order 1
#' @export
DEGPD4 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dDEGPD4", pfun = "pDEGPD4",
    family_name = "DEGPD4", dist_type = "Discrete",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam <- .override_derivs4(fam, "degpd4")
  fam <- .egpd_set_initial(fam)
  fam$tau.valid     <- function(tau) all(tau > 0)
  fam
}

## ---------------------------------------------------------------
## Zero-Inflated Continuous EGPD families
## ---------------------------------------------------------------

#' gamlss Family for Zero-Inflated Continuous EGPD Model 1
#'
#' The \code{ZIEGPD1()} function defines a \code{gamlss.family} distribution
#' for fitting a zero-inflated continuous Extended Generalized Pareto
#' Distribution with model 1 (power G-transformation: \eqn{G(u) = u^\kappa}).
#' The functions \code{dZIEGPD1}, \code{pZIEGPD1}, \code{qZIEGPD1} and
#' \code{rZIEGPD1} give the density, distribution function, quantile function
#' and random generation for this distribution, using the gamlss parameter
#' names.
#'
#' The gamlss parameters map to ZIEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation), \code{tau} = pi (zero-inflation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"logit"})
#'
#' @return \code{ZIEGPD1()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dZIEGPD1} gives the density, \code{pZIEGPD1} gives the distribution
#'   function, \code{qZIEGPD1} gives the quantile function, and
#'   \code{rZIEGPD1} generates random deviates.
#' @references
#' Abbas, A., Ahmad, T. and Ahmad, I. (2025).
#' Modeling zero-inflated precipitation extremes.
#' \emph{Communications in Statistics - Simulation and Computation}, 1--17.
#' \doi{10.1080/03610918.2025.2585398}
#'
#' Naveau, P., Huser, R., Ribereau, P. and Hannart, A. (2016).
#' Modeling jointly low, moderate, and heavy rainfall intensities without a
#' threshold selection. \emph{Water Resources Research}, 52(4), 2753--2769.
#' \doi{10.1002/2015WR018552}
#' @seealso \code{\link{ZIEGPD3}}, \code{\link{ziegpd_bamlss}},
#'   \code{\link{dziegpd}}
#' @order 1
#' @export
ZIEGPD1 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "logit") {
  fam <- .make_gamlss_family(
    dfun = "dZIEGPD1", pfun = "pZIEGPD1",
    family_name = "ZIEGPD1", dist_type = "Mixed",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "logit"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam <- .override_derivs4(fam, "ziegpd1")
  fam <- .egpd_set_initial(fam)
  fam$tau.valid     <- function(tau) all(tau > 0 & tau < 1)
  fam
}

#' gamlss Family for Zero-Inflated Continuous EGPD Model 3
#'
#' The \code{ZIEGPD3()} function defines a \code{gamlss.family} distribution
#' for fitting a zero-inflated continuous Extended Generalized Pareto
#' Distribution with model 3 (incomplete beta G-transformation,
#' \code{type = 4}). The functions \code{dZIEGPD3}, \code{pZIEGPD3},
#' \code{qZIEGPD3} and \code{rZIEGPD3} give the density, distribution
#' function, quantile function and random generation for this distribution,
#' using the gamlss parameter names.
#'
#' The gamlss parameters map to ZIEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = delta (G-transformation), \code{tau} = pi (zero-inflation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"logit"})
#'
#' @return \code{ZIEGPD3()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dZIEGPD3} gives the density, \code{pZIEGPD3} gives the distribution
#'   function, \code{qZIEGPD3} gives the quantile function, and
#'   \code{rZIEGPD3} generates random deviates.
#' @references
#' Abbas, A., Ahmad, T. and Ahmad, I. (2025).
#' Modeling zero-inflated precipitation extremes.
#' \emph{Communications in Statistics - Simulation and Computation}, 1--17.
#' \doi{10.1080/03610918.2025.2585398}
#'
#' Naveau, P., Huser, R., Ribereau, P. and Hannart, A. (2016).
#' Modeling jointly low, moderate, and heavy rainfall intensities without a
#' threshold selection. \emph{Water Resources Research}, 52(4), 2753--2769.
#' \doi{10.1002/2015WR018552}
#' @seealso \code{\link{ZIEGPD1}}, \code{\link{ziegpd_bamlss}},
#'   \code{\link{dziegpd}}
#' @order 1
#' @export
ZIEGPD3 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "logit") {
  fam <- .make_gamlss_family(
    dfun = "dZIEGPD3", pfun = "pZIEGPD3",
    family_name = "ZIEGPD3", dist_type = "Mixed",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "logit"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam <- .override_derivs4(fam, "ziegpd3")
  fam <- .egpd_set_initial(fam)
  fam$tau.valid     <- function(tau) all(tau > 0 & tau < 1)
  fam
}

## ---------------------------------------------------------------
## Zero-Inflated Discrete EGPD families
## ---------------------------------------------------------------

#' gamlss Family for Zero-Inflated Discrete EGPD Model 1
#'
#' The \code{ZIDEGPD1()} function defines a \code{gamlss.family} distribution
#' for fitting a zero-inflated discrete Extended Generalized Pareto
#' Distribution with model 1 (power G-transformation: \eqn{G(u) = u^\kappa}).
#' The functions \code{dZIDEGPD1}, \code{pZIDEGPD1}, \code{qZIDEGPD1} and
#' \code{rZIDEGPD1} give the PMF, distribution function, quantile function and
#' random generation for this distribution, using the gamlss parameter names.
#'
#' The gamlss parameters map to ZIDEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation), \code{tau} = pi (zero-inflation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"logit"})
#'
#' @return \code{ZIDEGPD1()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dZIDEGPD1} gives the PMF, \code{pZIDEGPD1} gives the distribution
#'   function, \code{qZIDEGPD1} gives the quantile function, and
#'   \code{rZIDEGPD1} generates random deviates.
#' @references
#' Ahmad, T. and Arshad, I. A. (2026).
#' New flexible versions of the extended generalized Pareto distribution
#' for count data. \emph{Journal of Applied Statistics}, 1--21.
#' \doi{10.1080/02664763.2026.2634785} (preprint: arXiv:2409.18719)
#'
#' Ahmad, T., Gaetan, C. and Naveau, P. (2025).
#' An extended generalized Pareto regression model for count data.
#' \emph{Statistical Modelling}, 25(5), 416--431.
#' \doi{10.1177/1471082X241266729}
#'
#' Ahmad, T. and Hussain, A. (2025).
#' Flexible model for varying levels of zeros and outliers in count data.
#' \emph{arXiv preprint} arXiv:2510.27365.
#' \url{https://arxiv.org/abs/2510.27365}
#'
#' Naveau, P., Huser, R., Ribereau, P. and Hannart, A. (2016).
#' Modeling jointly low, moderate, and heavy rainfall intensities without a
#' threshold selection. \emph{Water Resources Research}, 52(4), 2753--2769.
#' \doi{10.1002/2015WR018552}
#' @seealso \code{\link{ZIDEGPD3}}, \code{\link{zidegpd_bamlss}},
#'   \code{\link{dzidiscegpd}}
#' @examples
#' \dontrun{
#'   library(gamlss)
#'   set.seed(1)
#'   y <- rZIDEGPD1(500, mu = 3, sigma = 0.15, nu = 2, tau = 0.25)
#'   fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
#'                 tau.formula = ~ 1, family = ZIDEGPD1())
#' }
#' @order 1
#' @export
ZIDEGPD1 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "logit") {
  fam <- .make_gamlss_family(
    dfun = "dZIDEGPD1", pfun = "pZIDEGPD1",
    family_name = "ZIDEGPD1", dist_type = "Discrete",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "logit"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam <- .override_derivs4(fam, "zidegpd1")
  fam <- .egpd_set_initial(fam)
  fam$tau.valid     <- function(tau) all(tau > 0 & tau < 1)
  fam
}

#' gamlss Family for Zero-Inflated Discrete EGPD Model 3
#'
#' The \code{ZIDEGPD3()} function defines a \code{gamlss.family} distribution
#' for fitting a zero-inflated discrete Extended Generalized Pareto
#' Distribution with model 3 (incomplete beta G-transformation,
#' \code{type = 4}). The functions \code{dZIDEGPD3}, \code{pZIDEGPD3},
#' \code{qZIDEGPD3} and \code{rZIDEGPD3} give the PMF, distribution function,
#' quantile function and random generation for this distribution, using the
#' gamlss parameter names.
#'
#' The gamlss parameters map to ZIDEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = delta (G-transformation), \code{tau} = pi (zero-inflation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"logit"})
#'
#' @return \code{ZIDEGPD3()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dZIDEGPD3} gives the PMF, \code{pZIDEGPD3} gives the distribution
#'   function, \code{qZIDEGPD3} gives the quantile function, and
#'   \code{rZIDEGPD3} generates random deviates.
#' @references
#' Ahmad, T. and Arshad, I. A. (2026).
#' New flexible versions of the extended generalized Pareto distribution
#' for count data. \emph{Journal of Applied Statistics}, 1--21.
#' \doi{10.1080/02664763.2026.2634785} (preprint: arXiv:2409.18719)
#'
#' Ahmad, T., Gaetan, C. and Naveau, P. (2025).
#' An extended generalized Pareto regression model for count data.
#' \emph{Statistical Modelling}, 25(5), 416--431.
#' \doi{10.1177/1471082X241266729}
#'
#' Ahmad, T. and Hussain, A. (2025).
#' Flexible model for varying levels of zeros and outliers in count data.
#' \emph{arXiv preprint} arXiv:2510.27365.
#' \url{https://arxiv.org/abs/2510.27365}
#'
#' Naveau, P., Huser, R., Ribereau, P. and Hannart, A. (2016).
#' Modeling jointly low, moderate, and heavy rainfall intensities without a
#' threshold selection. \emph{Water Resources Research}, 52(4), 2753--2769.
#' \doi{10.1002/2015WR018552}
#' @seealso \code{\link{ZIDEGPD1}}, \code{\link{zidegpd_bamlss}},
#'   \code{\link{dzidiscegpd}}
#' @order 1
#' @export
ZIDEGPD3 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "logit") {
  fam <- .make_gamlss_family(
    dfun = "dZIDEGPD3", pfun = "pZIDEGPD3",
    family_name = "ZIDEGPD3", dist_type = "Discrete",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "logit"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam <- .override_derivs4(fam, "zidegpd3")
  fam <- .egpd_set_initial(fam)
  fam$tau.valid     <- function(tau) all(tau > 0 & tau < 1)
  fam
}

## ---------------------------------------------------------------
## d/p/q/r wrappers: Continuous EGPD Model 5 (truncated normal)
## ---------------------------------------------------------------

#' @rdname EGPD5
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (kappa), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dEGPD5 <- function(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE) {
  degpd_density(x, sigma = mu, xi = sigma, kappa = nu, type = 2L, log = log)
}

#' @rdname EGPD5
#' @export
pEGPD5 <- function(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  pegpd(q, sigma = mu, xi = sigma, kappa = nu, type = 2L,
        lower.tail = lower.tail, log.p = log.p)
}

#' @rdname EGPD5
#' @export
qEGPD5 <- function(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  qegpd(p, sigma = mu, xi = sigma, kappa = nu, type = 2L,
        lower.tail = lower.tail, log.p = log.p)
}

#' @rdname EGPD5
#' @export
rEGPD5 <- function(n, mu = 1, sigma = 0.5, nu = 1) {
  regpd(n, sigma = mu, xi = sigma, kappa = nu, type = 2L)
}

#' @rdname EGPD6
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (kappa), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dEGPD6 <- function(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE) {
  degpd_density(x, sigma = mu, xi = sigma, kappa = nu, type = 3L, log = log)
}

#' @rdname EGPD6
#' @export
pEGPD6 <- function(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  pegpd(q, sigma = mu, xi = sigma, kappa = nu, type = 3L,
        lower.tail = lower.tail, log.p = log.p)
}

#' @rdname EGPD6
#' @export
qEGPD6 <- function(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  qegpd(p, sigma = mu, xi = sigma, kappa = nu, type = 3L,
        lower.tail = lower.tail, log.p = log.p)
}

#' @rdname EGPD6
#' @export
rEGPD6 <- function(n, mu = 1, sigma = 0.5, nu = 1) {
  regpd(n, sigma = mu, xi = sigma, kappa = nu, type = 3L)
}

## ---------------------------------------------------------------
## d/p/q/r wrappers: Discrete EGPD Models 5 and 6
## ---------------------------------------------------------------

#' @rdname DEGPD5
#' @param x,q vector of quantiles (non-negative integers)
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (kappa), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dDEGPD5 <- function(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE) {
  ddiscegpd(x, sigma = mu, xi = sigma, kappa = nu, type = 2L, log = log)
}

#' @rdname DEGPD5
#' @export
pDEGPD5 <- function(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  pdiscegpd(q, sigma = mu, xi = sigma, kappa = nu, type = 2L,
            lower.tail = lower.tail, log.p = log.p)
}

#' @rdname DEGPD5
#' @export
qDEGPD5 <- function(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  qdiscegpd(p, sigma = mu, xi = sigma, kappa = nu, type = 2L,
            lower.tail = lower.tail, log.p = log.p)
}

#' @rdname DEGPD5
#' @export
rDEGPD5 <- function(n, mu = 1, sigma = 0.5, nu = 1) {
  rdiscegpd(n, sigma = mu, xi = sigma, kappa = nu, type = 2L)
}

#' @rdname DEGPD6
#' @param x,q vector of quantiles (non-negative integers)
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (kappa), positive
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dDEGPD6 <- function(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE) {
  ddiscegpd(x, sigma = mu, xi = sigma, kappa = nu, type = 3L, log = log)
}

#' @rdname DEGPD6
#' @export
pDEGPD6 <- function(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  pdiscegpd(q, sigma = mu, xi = sigma, kappa = nu, type = 3L,
            lower.tail = lower.tail, log.p = log.p)
}

#' @rdname DEGPD6
#' @export
qDEGPD6 <- function(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  qdiscegpd(p, sigma = mu, xi = sigma, kappa = nu, type = 3L,
            lower.tail = lower.tail, log.p = log.p)
}

#' @rdname DEGPD6
#' @export
rDEGPD6 <- function(n, mu = 1, sigma = 0.5, nu = 1) {
  rdiscegpd(n, sigma = mu, xi = sigma, kappa = nu, type = 3L)
}

## ---------------------------------------------------------------
## Zero-Inflated EGPD Models 5 and 6
## ---------------------------------------------------------------

#' @rdname ZIEGPD5
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (kappa), positive
#' @param tau zero-inflation probability (pi), in (0,1)
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dZIEGPD5 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE) {
  dziegpd(x, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 2L, log = log)
}
#' @rdname ZIEGPD5
#' @export
pZIEGPD5 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  pziegpd(q, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 2L,
          lower.tail = lower.tail, log.p = log.p)
}
#' @rdname ZIEGPD5
#' @export
qZIEGPD5 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  qziegpd(p, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 2L,
          lower.tail = lower.tail, log.p = log.p)
}
#' @rdname ZIEGPD5
#' @export
rZIEGPD5 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1) {
  rziegpd(n, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 2L)
}

#' @rdname ZIEGPD6
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (kappa), positive
#' @param tau zero-inflation probability (pi), in (0,1)
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dZIEGPD6 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE) {
  dziegpd(x, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 3L, log = log)
}
#' @rdname ZIEGPD6
#' @export
pZIEGPD6 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  pziegpd(q, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 3L,
          lower.tail = lower.tail, log.p = log.p)
}
#' @rdname ZIEGPD6
#' @export
qZIEGPD6 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  qziegpd(p, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 3L,
          lower.tail = lower.tail, log.p = log.p)
}
#' @rdname ZIEGPD6
#' @export
rZIEGPD6 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1) {
  rziegpd(n, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 3L)
}

## ---------------------------------------------------------------
## Zero-Inflated Discrete EGPD Models 5 and 6
## ---------------------------------------------------------------

#' @rdname ZIDEGPD5
#' @param x,q vector of quantiles (non-negative integers)
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (kappa), positive
#' @param tau zero-inflation probability (pi), in (0,1)
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dZIDEGPD5 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE) {
  dzidiscegpd(x, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 2L, log = log)
}
#' @rdname ZIDEGPD5
#' @export
pZIDEGPD5 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  pzidiscegpd(q, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 2L,
              lower.tail = lower.tail, log.p = log.p)
}
#' @rdname ZIDEGPD5
#' @export
qZIDEGPD5 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  qzidiscegpd(p, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 2L,
              lower.tail = lower.tail, log.p = log.p)
}
#' @rdname ZIDEGPD5
#' @export
rZIDEGPD5 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1) {
  rzidiscegpd(n, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 2L)
}

#' @rdname ZIDEGPD6
#' @param x,q vector of quantiles (non-negative integers)
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu GPD scale parameter (sigma), positive
#' @param sigma GPD shape parameter (xi); a negative value gives a bounded
#'   upper tail
#' @param nu G-transformation parameter (kappa), positive
#' @param tau zero-inflation probability (pi), in (0,1)
#' @param log,log.p logical; if TRUE, probabilities/densities are given as log
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x)
#' @export
dZIDEGPD6 <- function(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE) {
  dzidiscegpd(x, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 3L, log = log)
}
#' @rdname ZIDEGPD6
#' @export
pZIDEGPD6 <- function(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  pzidiscegpd(q, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 3L,
              lower.tail = lower.tail, log.p = log.p)
}
#' @rdname ZIDEGPD6
#' @export
qZIDEGPD6 <- function(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE, log.p = FALSE) {
  qzidiscegpd(p, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 3L,
              lower.tail = lower.tail, log.p = log.p)
}
#' @rdname ZIDEGPD6
#' @export
rZIDEGPD6 <- function(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1) {
  rzidiscegpd(n, pi = tau, sigma = mu, xi = sigma, kappa = nu, type = 3L)
}

## ---------------------------------------------------------------
## gamlss Family Constructors for new models
## ---------------------------------------------------------------

#' gamlss Family for Continuous EGPD Model 5 (truncated normal)
#'
#' The \code{EGPD5()} function defines a \code{gamlss.family} distribution for
#' fitting a continuous Extended Generalized Pareto Distribution with model 5
#' (truncated normal G-transformation, \code{type = 2}:
#' \eqn{G(u) = (\Phi(\sqrt{\kappa}(u-1)) - \Phi(-\sqrt{\kappa})) /
#' (0.5 - \Phi(-\sqrt{\kappa}))}). The functions \code{dEGPD5}, \code{pEGPD5},
#' \code{qEGPD5} and \code{rEGPD5} give the density, distribution function,
#' quantile function and random generation for this distribution, using the
#' gamlss parameter names.
#'
#' The gamlss parameters map to EGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#'
#' @return \code{EGPD5()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dEGPD5} gives the density, \code{pEGPD5} gives the distribution
#'   function, \code{qEGPD5} gives the quantile function, and \code{rEGPD5}
#'   generates random deviates.
#' @references
#' Ahmad, T. and Arshad, I. A. (2026).
#' New flexible versions of the extended generalized Pareto distribution
#' for count data. \emph{Journal of Applied Statistics}, 1--21.
#' \doi{10.1080/02664763.2026.2634785} (preprint: arXiv:2409.18719)
#'
#' Gamet, P. and Jalbert, J. (2022).
#' A flexible extended generalized Pareto distribution for tail estimation.
#' \emph{Environmetrics}, 33(6), e2744. \doi{10.1002/env.2744}
#' @seealso \code{\link{EGPD1}}, \code{\link{EGPD6}},
#'   \code{\link{egpd_bamlss}}, \code{\link{pegpd}}
#' @examples
#' \dontrun{
#'   library(gamlss)
#'   set.seed(1)
#'   y <- rEGPD5(500, mu = 2, sigma = 0.2, nu = 1.5)
#'   fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
#'                 family = EGPD5())
#' }
#' @order 1
#' @export
EGPD5 <- function(mu.link = "log", sigma.link = "log", nu.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dEGPD5", pfun = "pEGPD5",
    family_name = "EGPD5", dist_type = "Continuous",
    npar = 3,
    default_links = c(mu = "log", sigma = "log", nu = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link
  )
  fam <- .override_derivs3(fam, "egpd5")
  fam <- .egpd_set_initial(fam)
  fam
}

#' gamlss Family for Continuous EGPD Model 6 (truncated beta)
#'
#' The \code{EGPD6()} function defines a \code{gamlss.family} distribution for
#' fitting a continuous Extended Generalized Pareto Distribution with model 6
#' (truncated beta G-transformation, \code{type = 3}). The functions
#' \code{dEGPD6}, \code{pEGPD6}, \code{qEGPD6} and \code{rEGPD6} give the
#' density, distribution function, quantile function and random generation for
#' this distribution, using the gamlss parameter names.
#'
#' The gamlss parameters map to EGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#'
#' @return \code{EGPD6()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dEGPD6} gives the density, \code{pEGPD6} gives the distribution
#'   function, \code{qEGPD6} gives the quantile function, and \code{rEGPD6}
#'   generates random deviates.
#' @references
#' Ahmad, T. and Arshad, I. A. (2026).
#' New flexible versions of the extended generalized Pareto distribution
#' for count data. \emph{Journal of Applied Statistics}, 1--21.
#' \doi{10.1080/02664763.2026.2634785} (preprint: arXiv:2409.18719)
#'
#' Gamet, P. and Jalbert, J. (2022).
#' A flexible extended generalized Pareto distribution for tail estimation.
#' \emph{Environmetrics}, 33(6), e2744. \doi{10.1002/env.2744}
#' @seealso \code{\link{EGPD1}}, \code{\link{EGPD5}},
#'   \code{\link{egpd_bamlss}}, \code{\link{pegpd}}
#' @examples
#' \dontrun{
#'   library(gamlss)
#'   set.seed(1)
#'   y <- rEGPD6(500, mu = 2, sigma = 0.2, nu = 1.5)
#'   fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
#'                 family = EGPD6())
#' }
#' @order 1
#' @export
EGPD6 <- function(mu.link = "log", sigma.link = "log", nu.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dEGPD6", pfun = "pEGPD6",
    family_name = "EGPD6", dist_type = "Continuous",
    npar = 3,
    default_links = c(mu = "log", sigma = "log", nu = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link
  )
  fam <- .override_derivs3(fam, "egpd6")
  fam <- .egpd_set_initial(fam)
  fam
}

#' gamlss Family for Discrete EGPD Model 5 (truncated normal)
#'
#' The \code{DEGPD5()} function defines a \code{gamlss.family} distribution
#' for fitting a discrete Extended Generalized Pareto Distribution with model
#' 5 (truncated normal G-transformation, \code{type = 2}:
#' \eqn{G(u) = (\Phi(\sqrt{\kappa}(u-1)) - \Phi(-\sqrt{\kappa})) /
#' (0.5 - \Phi(-\sqrt{\kappa}))}). The functions \code{dDEGPD5},
#' \code{pDEGPD5}, \code{qDEGPD5} and \code{rDEGPD5} give the PMF,
#' distribution function, quantile function and random generation for this
#' distribution, using the gamlss parameter names.
#'
#' The gamlss parameters map to DEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#'
#' @return \code{DEGPD5()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dDEGPD5} gives the PMF, \code{pDEGPD5} gives the distribution
#'   function, \code{qDEGPD5} gives the quantile function, and \code{rDEGPD5}
#'   generates random deviates.
#' @references
#' Ahmad, T. and Arshad, I. A. (2026).
#' New flexible versions of the extended generalized Pareto distribution
#' for count data. \emph{Journal of Applied Statistics}, 1--21.
#' \doi{10.1080/02664763.2026.2634785} (preprint: arXiv:2409.18719)
#'
#' Ahmad, T., Gaetan, C. and Naveau, P. (2025).
#' An extended generalized Pareto regression model for count data.
#' \emph{Statistical Modelling}, 25(5), 416--431.
#' \doi{10.1177/1471082X241266729}
#'
#' Gamet, P. and Jalbert, J. (2022).
#' A flexible extended generalized Pareto distribution for tail estimation.
#' \emph{Environmetrics}, 33(6), e2744. \doi{10.1002/env.2744}
#' @seealso \code{\link{DEGPD1}}, \code{\link{DEGPD6}},
#'   \code{\link{degpd_bamlss}}, \code{\link{ddiscegpd}}
#' @examples
#' \dontrun{
#'   library(gamlss)
#'   set.seed(1)
#'   y <- rDEGPD5(500, mu = 3, sigma = 0.15, nu = 2)
#'   fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
#'                 family = DEGPD5())
#' }
#' @order 1
#' @export
DEGPD5 <- function(mu.link = "log", sigma.link = "log", nu.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dDEGPD5", pfun = "pDEGPD5",
    family_name = "DEGPD5", dist_type = "Discrete",
    npar = 3,
    default_links = c(mu = "log", sigma = "log", nu = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link
  )
  fam <- .override_derivs3(fam, "degpd5")
  fam <- .egpd_set_initial(fam)
  fam
}

#' gamlss Family for Discrete EGPD Model 6 (truncated beta)
#'
#' The \code{DEGPD6()} function defines a \code{gamlss.family} distribution
#' for fitting a discrete Extended Generalized Pareto Distribution with model
#' 6 (truncated beta G-transformation, \code{type = 3}). The functions
#' \code{dDEGPD6}, \code{pDEGPD6}, \code{qDEGPD6} and \code{rDEGPD6} give the
#' PMF, distribution function, quantile function and random generation for
#' this distribution, using the gamlss parameter names.
#'
#' The gamlss parameters map to DEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#'
#' @return \code{DEGPD6()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dDEGPD6} gives the PMF, \code{pDEGPD6} gives the distribution
#'   function, \code{qDEGPD6} gives the quantile function, and \code{rDEGPD6}
#'   generates random deviates.
#' @references
#' Ahmad, T. and Arshad, I. A. (2026).
#' New flexible versions of the extended generalized Pareto distribution
#' for count data. \emph{Journal of Applied Statistics}, 1--21.
#' \doi{10.1080/02664763.2026.2634785} (preprint: arXiv:2409.18719)
#'
#' Ahmad, T., Gaetan, C. and Naveau, P. (2025).
#' An extended generalized Pareto regression model for count data.
#' \emph{Statistical Modelling}, 25(5), 416--431.
#' \doi{10.1177/1471082X241266729}
#'
#' Gamet, P. and Jalbert, J. (2022).
#' A flexible extended generalized Pareto distribution for tail estimation.
#' \emph{Environmetrics}, 33(6), e2744. \doi{10.1002/env.2744}
#' @seealso \code{\link{DEGPD1}}, \code{\link{DEGPD5}},
#'   \code{\link{degpd_bamlss}}, \code{\link{ddiscegpd}}
#' @examples
#' \dontrun{
#'   library(gamlss)
#'   set.seed(1)
#'   y <- rDEGPD6(500, mu = 3, sigma = 0.15, nu = 2)
#'   fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
#'                 family = DEGPD6())
#' }
#' @order 1
#' @export
DEGPD6 <- function(mu.link = "log", sigma.link = "log", nu.link = "log") {
  fam <- .make_gamlss_family(
    dfun = "dDEGPD6", pfun = "pDEGPD6",
    family_name = "DEGPD6", dist_type = "Discrete",
    npar = 3,
    default_links = c(mu = "log", sigma = "log", nu = "log"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link
  )
  fam <- .override_derivs3(fam, "degpd6")
  fam <- .egpd_set_initial(fam)
  fam
}

#' gamlss Family for Zero-Inflated Continuous EGPD Model 5
#'
#' The \code{ZIEGPD5()} function defines a \code{gamlss.family} distribution
#' for fitting a zero-inflated continuous Extended Generalized Pareto
#' Distribution with model 5 (truncated normal G-transformation,
#' \code{type = 2}:
#' \eqn{G(u) = (\Phi(\sqrt{\kappa}(u-1)) - \Phi(-\sqrt{\kappa})) /
#' (0.5 - \Phi(-\sqrt{\kappa}))}). The functions \code{dZIEGPD5},
#' \code{pZIEGPD5}, \code{qZIEGPD5} and \code{rZIEGPD5} give the density,
#' distribution function, quantile function and random generation for this
#' distribution, using the gamlss parameter names.
#'
#' The gamlss parameters map to ZIEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation), \code{tau} = pi (zero-inflation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"logit"})
#'
#' @return \code{ZIEGPD5()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dZIEGPD5} gives the density, \code{pZIEGPD5} gives the distribution
#'   function, \code{qZIEGPD5} gives the quantile function, and
#'   \code{rZIEGPD5} generates random deviates.
#' @references
#' Abbas, A., Ahmad, T. and Ahmad, I. (2025).
#' Modeling zero-inflated precipitation extremes.
#' \emph{Communications in Statistics - Simulation and Computation}, 1--17.
#' \doi{10.1080/03610918.2025.2585398}
#'
#' Ahmad, T. and Arshad, I. A. (2026).
#' New flexible versions of the extended generalized Pareto distribution
#' for count data. \emph{Journal of Applied Statistics}, 1--21.
#' \doi{10.1080/02664763.2026.2634785} (preprint: arXiv:2409.18719)
#'
#' Gamet, P. and Jalbert, J. (2022).
#' A flexible extended generalized Pareto distribution for tail estimation.
#' \emph{Environmetrics}, 33(6), e2744. \doi{10.1002/env.2744}
#' @seealso \code{\link{ZIEGPD1}}, \code{\link{ZIEGPD6}},
#'   \code{\link{ziegpd_bamlss}}, \code{\link{dziegpd}}
#' @order 1
#' @export
ZIEGPD5 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "logit") {
  fam <- .make_gamlss_family(
    dfun = "dZIEGPD5", pfun = "pZIEGPD5",
    family_name = "ZIEGPD5", dist_type = "Mixed",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "logit"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam <- .override_derivs4(fam, "ziegpd5")
  fam <- .egpd_set_initial(fam)
  fam$tau.valid     <- function(tau) all(tau > 0 & tau < 1)
  fam
}

#' gamlss Family for Zero-Inflated Continuous EGPD Model 6
#'
#' The \code{ZIEGPD6()} function defines a \code{gamlss.family} distribution
#' for fitting a zero-inflated continuous Extended Generalized Pareto
#' Distribution with model 6 (truncated beta G-transformation,
#' \code{type = 3}). The functions \code{dZIEGPD6}, \code{pZIEGPD6},
#' \code{qZIEGPD6} and \code{rZIEGPD6} give the density, distribution
#' function, quantile function and random generation for this distribution,
#' using the gamlss parameter names.
#'
#' The gamlss parameters map to ZIEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation), \code{tau} = pi (zero-inflation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"logit"})
#'
#' @return \code{ZIEGPD6()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dZIEGPD6} gives the density, \code{pZIEGPD6} gives the distribution
#'   function, \code{qZIEGPD6} gives the quantile function, and
#'   \code{rZIEGPD6} generates random deviates.
#' @references
#' Abbas, A., Ahmad, T. and Ahmad, I. (2025).
#' Modeling zero-inflated precipitation extremes.
#' \emph{Communications in Statistics - Simulation and Computation}, 1--17.
#' \doi{10.1080/03610918.2025.2585398}
#'
#' Ahmad, T. and Arshad, I. A. (2026).
#' New flexible versions of the extended generalized Pareto distribution
#' for count data. \emph{Journal of Applied Statistics}, 1--21.
#' \doi{10.1080/02664763.2026.2634785} (preprint: arXiv:2409.18719)
#'
#' Gamet, P. and Jalbert, J. (2022).
#' A flexible extended generalized Pareto distribution for tail estimation.
#' \emph{Environmetrics}, 33(6), e2744. \doi{10.1002/env.2744}
#' @seealso \code{\link{ZIEGPD1}}, \code{\link{ZIEGPD5}},
#'   \code{\link{ziegpd_bamlss}}, \code{\link{dziegpd}}
#' @order 1
#' @export
ZIEGPD6 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "logit") {
  fam <- .make_gamlss_family(
    dfun = "dZIEGPD6", pfun = "pZIEGPD6",
    family_name = "ZIEGPD6", dist_type = "Mixed",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "logit"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam <- .override_derivs4(fam, "ziegpd6")
  fam <- .egpd_set_initial(fam)
  fam$tau.valid     <- function(tau) all(tau > 0 & tau < 1)
  fam
}

#' gamlss Family for Zero-Inflated Discrete EGPD Model 5
#'
#' The \code{ZIDEGPD5()} function defines a \code{gamlss.family} distribution
#' for fitting a zero-inflated discrete Extended Generalized Pareto
#' Distribution with model 5 (truncated normal G-transformation,
#' \code{type = 2}:
#' \eqn{G(u) = (\Phi(\sqrt{\kappa}(u-1)) - \Phi(-\sqrt{\kappa})) /
#' (0.5 - \Phi(-\sqrt{\kappa}))}). The functions \code{dZIDEGPD5},
#' \code{pZIDEGPD5}, \code{qZIDEGPD5} and \code{rZIDEGPD5} give the PMF,
#' distribution function, quantile function and random generation for this
#' distribution, using the gamlss parameter names.
#'
#' The gamlss parameters map to ZIDEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation), \code{tau} = pi (zero-inflation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"logit"})
#'
#' @return \code{ZIDEGPD5()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dZIDEGPD5} gives the PMF, \code{pZIDEGPD5} gives the distribution
#'   function, \code{qZIDEGPD5} gives the quantile function, and
#'   \code{rZIDEGPD5} generates random deviates.
#' @references
#' Ahmad, T. and Arshad, I. A. (2026).
#' New flexible versions of the extended generalized Pareto distribution
#' for count data. \emph{Journal of Applied Statistics}, 1--21.
#' \doi{10.1080/02664763.2026.2634785} (preprint: arXiv:2409.18719)
#'
#' Ahmad, T., Gaetan, C. and Naveau, P. (2025).
#' An extended generalized Pareto regression model for count data.
#' \emph{Statistical Modelling}, 25(5), 416--431.
#' \doi{10.1177/1471082X241266729}
#'
#' Ahmad, T. and Hussain, A. (2025).
#' Flexible model for varying levels of zeros and outliers in count data.
#' \emph{arXiv preprint} arXiv:2510.27365.
#' \url{https://arxiv.org/abs/2510.27365}
#'
#' Gamet, P. and Jalbert, J. (2022).
#' A flexible extended generalized Pareto distribution for tail estimation.
#' \emph{Environmetrics}, 33(6), e2744. \doi{10.1002/env.2744}
#' @seealso \code{\link{ZIDEGPD1}}, \code{\link{ZIDEGPD6}},
#'   \code{\link{zidegpd_bamlss}}, \code{\link{dzidiscegpd}}
#' @examples
#' \dontrun{
#'   library(gamlss)
#'   set.seed(1)
#'   y <- rZIDEGPD5(500, mu = 3, sigma = 0.15, nu = 2, tau = 0.25)
#'   fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
#'                 tau.formula = ~ 1, family = ZIDEGPD5())
#' }
#' @order 1
#' @export
ZIDEGPD5 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "logit") {
  fam <- .make_gamlss_family(
    dfun = "dZIDEGPD5", pfun = "pZIDEGPD5",
    family_name = "ZIDEGPD5", dist_type = "Discrete",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "logit"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam <- .override_derivs4(fam, "zidegpd5")
  fam <- .egpd_set_initial(fam)
  fam$tau.valid     <- function(tau) all(tau > 0 & tau < 1)
  fam
}

#' gamlss Family for Zero-Inflated Discrete EGPD Model 6
#'
#' The \code{ZIDEGPD6()} function defines a \code{gamlss.family} distribution
#' for fitting a zero-inflated discrete Extended Generalized Pareto
#' Distribution with model 6 (truncated beta G-transformation,
#' \code{type = 3}). The functions \code{dZIDEGPD6}, \code{pZIDEGPD6},
#' \code{qZIDEGPD6} and \code{rZIDEGPD6} give the PMF, distribution function,
#' quantile function and random generation for this distribution, using the
#' gamlss parameter names.
#'
#' The gamlss parameters map to ZIDEGPD parameters as follows:
#' \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
#' \code{nu} = kappa (G-transformation), \code{tau} = pi (zero-inflation).
#'
#' @details The scores used by \code{gamlss()} are analytic and numerically
#'   stable, including zero counts, the far tail and shapes near zero; second
#'   derivatives use the per-observation approximation \eqn{-(\partial l /
#'   \partial \theta)^2}. An observation with zero likelihood at the current
#'   parameters (e.g. a count beyond a finite upper endpoint when
#'   \code{sigma.link = "identity"}) has no score: \code{NaN} is returned with
#'   a warning and \code{gamlss()} stops rather than silently ignoring it.
#'   Starting values are the maximum likelihood estimates of the
#'   intercept-only model (on at most 2000 evenly spaced observations plus the
#'   smallest and largest), since the RS algorithm moves slowly along the flat
#'   scale/shape/carrier ridges of these families. With
#'   \code{sigma.link = "identity"} and a fitted bounded tail (\eqn{\xi < 0})
#'   they are refined on all observations, so that the start gives every
#'   observation positive likelihood; a start that does not falls back to
#'   moment-based values. They are computed once per response and reused by
#'   the other \code{*.initial} expressions (also under \code{gamlssML()});
#'   supply \code{mu.start} etc. to override them. With a \code{"logit"}
#'   link for \code{nu} or \code{tau}, a start outside (0, 1) is moved inside it.
#' @param mu.link link function for mu (default \code{"log"})
#' @param sigma.link link function for sigma (default \code{"log"}); use
#'   \code{"identity"} to allow a negative shape (bounded upper tail)
#' @param nu.link link function for nu (default \code{"log"})
#' @param tau.link link function for tau (default \code{"logit"})
#'
#' @return \code{ZIDEGPD6()} returns an object of class \code{"gamlss.family"}.
#'
#'   \code{dZIDEGPD6} gives the PMF, \code{pZIDEGPD6} gives the distribution
#'   function, \code{qZIDEGPD6} gives the quantile function, and
#'   \code{rZIDEGPD6} generates random deviates.
#' @references
#' Ahmad, T. and Arshad, I. A. (2026).
#' New flexible versions of the extended generalized Pareto distribution
#' for count data. \emph{Journal of Applied Statistics}, 1--21.
#' \doi{10.1080/02664763.2026.2634785} (preprint: arXiv:2409.18719)
#'
#' Ahmad, T., Gaetan, C. and Naveau, P. (2025).
#' An extended generalized Pareto regression model for count data.
#' \emph{Statistical Modelling}, 25(5), 416--431.
#' \doi{10.1177/1471082X241266729}
#'
#' Ahmad, T. and Hussain, A. (2025).
#' Flexible model for varying levels of zeros and outliers in count data.
#' \emph{arXiv preprint} arXiv:2510.27365.
#' \url{https://arxiv.org/abs/2510.27365}
#'
#' Gamet, P. and Jalbert, J. (2022).
#' A flexible extended generalized Pareto distribution for tail estimation.
#' \emph{Environmetrics}, 33(6), e2744. \doi{10.1002/env.2744}
#' @seealso \code{\link{ZIDEGPD1}}, \code{\link{ZIDEGPD5}},
#'   \code{\link{zidegpd_bamlss}}, \code{\link{dzidiscegpd}}
#' @examples
#' \dontrun{
#'   library(gamlss)
#'   set.seed(1)
#'   y <- rZIDEGPD6(500, mu = 3, sigma = 0.15, nu = 2, tau = 0.25)
#'   fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
#'                 tau.formula = ~ 1, family = ZIDEGPD6())
#' }
#' @order 1
#' @export
ZIDEGPD6 <- function(mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "logit") {
  fam <- .make_gamlss_family(
    dfun = "dZIDEGPD6", pfun = "pZIDEGPD6",
    family_name = "ZIDEGPD6", dist_type = "Discrete",
    npar = 4,
    default_links = c(mu = "log", sigma = "log", nu = "log", tau = "logit"),
    mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link, tau.link = tau.link
  )
  fam <- .override_derivs4(fam, "zidegpd6")
  fam <- .egpd_set_initial(fam)
  fam$tau.valid     <- function(tau) all(tau > 0 & tau < 1)
  fam
}
