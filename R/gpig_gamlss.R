## gamlss family for the generalised Poisson-inverse Gaussian (GPIG; Zhu & Joe
## 2009) and zero-inflated GPIG (ZIGPIG). Because gamlss.dist has no GPIG, these
## build on the package recursion (dgpig/pgpig, R/gpig_dist.R).
##
## Mean parameterisation (consistent with the EGPD gamlss families):
##   mu    = mean            (log link)
##   sigma = a, tail exponent in (0,1)   (logit link)
##   nu    = c, down-weight  in (0,1)    (logit link)
##   tau   = pi, zero-inflation (ZIGPIG only, logit link)
## internally b = mu (1-c)^(1-a) / (a c).

## ---- gamlss-convention d/p/q/r (mu = mean, sigma = a, nu = c[, tau = pi]) ----

.gpig_b_from_mean <- function(mu, a, c, n = max(length(mu), length(a), length(c))) {
  if (!length(mu) || !length(a) || !length(c)) return(numeric(0))
  ## Transform the jointly recycled rows. The derived b need not have the
  ## same repetition period as any one of its three inputs.
  mu <- rep_len(mu, n); a <- rep_len(a, n); c <- rep_len(c, n)
  na_par <- (is.na(mu) & !is.nan(mu)) | (is.na(a) & !is.nan(a)) | (is.na(c) & !is.nan(c))
  ok <- is.finite(mu) & mu > 0 & is.finite(a) & a > 0 & a < 1 &
    is.finite(c) & c > 0 & c < 1
  out <- rep(NaN, n)
  out[na_par] <- NA_real_
  out[ok] <- mu[ok] * (1 - c[ok])^(1 - a[ok]) / (a[ok] * c[ok])
  out
}

#' @rdname gpig
#' @export
dGPIG <- function(x, mu = 1, sigma = 0.5, nu = 0.5, log = FALSE) {
  b <- .gpig_b_from_mean(mu, sigma, nu, max(lengths(list(x, mu, sigma, nu))))
  dgpig(x, a = sigma, b = b, c = nu, log = log)
}
#' @rdname gpig
#' @export
pGPIG <- function(q, mu = 1, sigma = 0.5, nu = 0.5, lower.tail = TRUE, log.p = FALSE) {
  b <- .gpig_b_from_mean(mu, sigma, nu, max(lengths(list(q, mu, sigma, nu))))
  pgpig(q, a = sigma, b = b, c = nu, lower.tail = lower.tail, log.p = log.p)
}
#' @rdname gpig
#' @export
qGPIG <- function(p, mu = 1, sigma = 0.5, nu = 0.5, lower.tail = TRUE, log.p = FALSE,
                  max.value = 1e5) {
  b <- .gpig_b_from_mean(mu, sigma, nu, max(lengths(list(p, mu, sigma, nu))))
  qgpig(p, a = sigma, b = b, c = nu, lower.tail = lower.tail, log.p = log.p,
        max.value = max.value)
}
#' @rdname gpig
#' @export
rGPIG <- function(n, mu = 1, sigma = 0.5, nu = 0.5, max.value = 1e5) {
  n <- .count_rn(n)
  b <- .gpig_b_from_mean(mu, sigma, nu, n)
  rgpig(n, a = sigma, b = b, c = nu, max.value = max.value)
}

#' @rdname zigpig
#' @export
dZIGPIG <- function(x, mu = 1, sigma = 0.5, nu = 0.5, tau = 0.1, log = FALSE) {
  b <- .gpig_b_from_mean(mu, sigma, nu, max(lengths(list(x, mu, sigma, nu, tau))))
  dzigpig(x, a = sigma, b = b, c = nu, pi = tau, log = log)
}
#' @rdname zigpig
#' @export
pZIGPIG <- function(q, mu = 1, sigma = 0.5, nu = 0.5, tau = 0.1,
                    lower.tail = TRUE, log.p = FALSE) {
  b <- .gpig_b_from_mean(mu, sigma, nu, max(lengths(list(q, mu, sigma, nu, tau))))
  pzigpig(q, a = sigma, b = b, c = nu, pi = tau, lower.tail = lower.tail, log.p = log.p)
}
#' @rdname zigpig
#' @export
qZIGPIG <- function(p, mu = 1, sigma = 0.5, nu = 0.5, tau = 0.1,
                    lower.tail = TRUE, log.p = FALSE, max.value = 1e5) {
  b <- .gpig_b_from_mean(mu, sigma, nu, max(lengths(list(p, mu, sigma, nu, tau))))
  qzigpig(p, a = sigma, b = b, c = nu, pi = tau, lower.tail = lower.tail,
          log.p = log.p, max.value = max.value)
}
#' @rdname zigpig
#' @export
rZIGPIG <- function(n, mu = 1, sigma = 0.5, nu = 0.5, tau = 0.1, max.value = 1e5) {
  n <- .count_rn(n)
  b <- .gpig_b_from_mean(mu, sigma, nu, n)
  rzigpig(n, a = sigma, b = b, c = nu, pi = tau, max.value = max.value)
}

## ---- family builders ----
## Mirror .make_gamlss_family but with logit links / (0,1) shape parameters and
## GPIG-appropriate initial values (c must stay < 1).

## log P(0) = -mu*K. Differentiate K = [(1-c)^(1-a) - (1-c)]/(a*c)
## without differencing nearly equal powers at small a or c.
.gpig_zero_score <- function(mu, a, c, idx) {
  L <- -log1p(-c)
  x <- a * L
  K <- (L / c) * exp(-(1 - a) * L) * .expm1_over_x(x)
  if (idx == 1L) return(-K)
  if (idx == 4L) return(-mu * K)  # log P(0) itself for mixture weights
  if (idx == 2L) {
    dlogK <- L * if (x < 1e-3) 1/2 + x/12 - x^3/720 else 1/(-expm1(-x)) - 1/x
    return(-mu * K * dlogK)
  }
  ## -K' = ((1-a*c)/(1-c)^a - 1)/(a*c^2).
  if (c < 1e-3)
    return(mu * (1-a) * (1/2 + c*(1+a)*(1/3 + c*(2+a)*(1/8 + c*(3+a)/30))))
  t <- if (a < 0.5) a*L + log1p(-a*c) else log1p((1-a)*c/(1-c)) - (1-a)*L
  mu * (expm1(t) / a) / c^2
}

## The log PGF is log P(0) + sum_j lambda_j z^j, with lambda_1 =
## mu*(1-c)^(1-a), lambda_{j+1}/lambda_j = c*(j-a)/(j+1).
## Differentiating the PGF gives the positive-count correction to the zero
## score: sum_j w_j d log(lambda_j), w_j = lambda_j P(y-j)/P(y).
## Compute the weights from the *same exact ladder as dGPIG*, normalising
## by sum(j*w_j) = y. No fitting approximation or gpig_control setting enters.
.gpig_positive_score <- function(y, mu, a, c, b, idx) {
  n <- max(y)
  j <- seq_len(n)
  k <- j[-n]
  logrel <- c(0, cumsum(log(c) + log1p(-a/k) - log1p(1/k)))
  lp <- .gpig_fam$lpmf(n, list(a = a, b = b, c = c))
  if (idx == 2L) da <- -log1p(-c) - c(0, cumsum(1/(k-a)))
  vapply(y, function(i) {
    h <- seq_len(i)
    prev <- lp[i + 1 - h]
    lw <- logrel[h] + (prev - max(prev))
    w <- exp(lw - max(lw))
    w <- w * (i / sum(h*w))
    switch(idx, sum(w)/mu, sum(w*da[h]), sum((h-1)*w)/c - sum(w)*(1-a)/(1-c))
  }, numeric(1))
}

.gpig_gamlss_score <- function(y, mu, sigma, nu, tau = 0, idx) {
  z <- .gpig_recycle(y = y, mu = mu, a = sigma, c = nu, pi = tau)
  y <- z$y; pi <- z$pi
  b <- .gpig_b_from_mean(z$mu, z$a, z$c)
  good <- which(is.finite(b) & b > 0 & is.finite(pi) & pi >= 0 & pi < 1 &
                  is.finite(y) & y >= 0 & y == floor(y))
  out <- rep(NaN, length(y))
  for (group in .count_groups(lapply(z[c("mu", "a", "c")], `[`, good))) {
    i <- good[group]; h <- i[1L]
    mu <- z$mu[h]; a <- z$a[h]; c <- z$c[h]
    iz <- i[y[i] == 0]; ip <- i[y[i] > 0]
    if (idx == 4L) {
      lp0 <- .gpig_zero_score(mu, a, c, 4L)
      out[iz] <- -expm1(lp0) / (pi[iz] + (1 - pi[iz]) * exp(lp0))
      out[ip] <- -1 / (1 - pi[ip])
    } else {
      s0 <- .gpig_zero_score(mu, a, c, idx)
      out[i] <- s0
      mix <- iz[pi[iz] > 0]
      if (length(mix)) {
        lp0 <- .gpig_zero_score(mu, a, c, 4L)
        out[mix] <- s0 * stats::plogis(lp0 - stats::qlogis(pi[mix]))
      }
      if (length(ip)) {
        u <- unique(y[ip])
        s <- if (max(u) == 1) switch(idx, 1/mu, -log1p(-c), -(1-a)/(1-c)) else
          .gpig_positive_score(u, mu, a, c, b[h], idx)
        out[ip] <- s0 + s[match(y[ip], u)]
      }
    }
  }
  missing <- Reduce(`|`, lapply(z, function(x) is.na(x) & !is.nan(x)))
  out[missing] <- NA_real_
  if (any(is.nan(out)))
    warning("GPIG score undefined for invalid parameters or observations outside the support.",
            call. = FALSE)
  out
}

#' @details
#' \code{GPIG()} builds a \code{gamlss.family} object for fitting the
#' generalised Poisson-inverse Gaussian distribution with \code{gamlss()}, in
#' the mean parameterisation \code{mu} = mean (log link), \code{sigma} = tail
#' exponent \eqn{a\in(0,1)} (logit link) and \code{nu} = down-weighting
#' \eqn{c\in(0,1)} (logit link). \code{dGPIG}, \code{pGPIG}, \code{qGPIG} and
#' \code{rGPIG} are the corresponding density, distribution, quantile and
#' simulation functions in the same parameterisation.
#' Scores differentiate the exact probability-generating function, including
#' at small means and near-unit down-weighting. Curvature uses the negative
#' outer product of scores.
#'
#' The starting values of \code{GPIG()} and \code{ZIGPIG()} are the maximum likelihood
#' estimates of the intercept-only model (from the fixed values \code{mean(y)},
#' \code{0.5}, \code{0.7} the RS algorithm can drift to \eqn{a \to 0} on heavy-tailed
#' data). They maximise the native GPIG likelihood under the default
#' \code{\link{gpig_control}()} settings, fall back to the fixed values if that fails,
#' and are computed once per response (also under \code{gamlssML()}); supply
#' \code{mu.start} etc. to override them.
#' @inheritSection EGPD1 GAMLSS namespace use
#'
#' @param mu.link,sigma.link,nu.link link functions for the parameters.
#' @param mu mean parameter (\eqn{>0}).
#' @param sigma tail exponent \eqn{a\in(0,1)}.
#' @param nu down-weighting parameter \eqn{c\in(0,1)}.
#' @rdname gpig
#' @export
GPIG <- function(mu.link = "log", sigma.link = "logit", nu.link = "logit") {
  .make_gpig_gamlss("dGPIG", "pGPIG", "GPIG", 3L,
                    mu.link, sigma.link, nu.link, NULL)
}

#' @details
#' \code{ZIGPIG()} adds a logit-link zero-inflation probability \code{tau} to
#' \code{\link{GPIG}}. \code{dZIGPIG}, \code{pZIGPIG}, \code{qZIGPIG} and
#' \code{rZIGPIG} are the corresponding mean-parameterised d/p/q/r functions.
#'
#' @param mu.link,sigma.link,nu.link link functions for the parameters.
#' @param tau.link link function for the zero-inflation parameter.
#' @param mu mean parameter (\eqn{>0}).
#' @param sigma tail exponent \eqn{a\in(0,1)}.
#' @param nu down-weighting parameter \eqn{c\in(0,1)}.
#' @param tau zero-inflation probability \eqn{\in[0,1)}.
#' @rdname zigpig
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
  dfun_sym <- call("::", quote(egpd), as.name(dfun))
  rqres_expr <- .gamlss_rqres_expr(pfun, "Discrete", npar)
  ## gamlss rebinds derivative bodies outside our namespace.
  score <- function(idx) bquote(get(".gpig_gamlss_score", envir = asNamespace("egpd"))(
    y, mu, sigma, nu, .(if (npar == 4L) quote(tau) else 0), idx = .(idx)))
  callback <- function(expr) {
    f <- if (npar == 3L) function(y, mu, sigma, nu, ...) NULL else
      function(y, mu, sigma, nu, tau, ...) NULL
    body(f) <- expr
    f
  }
  d1 <- function(i) callback(score(i))
  d2self <- function(i) callback(bquote(pmin(-.(score(i))^2, -1e-15)))
  d2cross <- function(i, j) callback(bquote(-.(score(i)) * .(score(j))))

  if (npar == 3L) {
    G.dev.incr_fn <- eval(bquote(function(y, mu, sigma, nu, ...) {
      -2 * .(dfun_sym)(y, mu = mu, sigma = sigma, nu = nu, log = TRUE)
    }))

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
      mu.initial    = .gpig_initial("mu", family_name),
      sigma.initial = .gpig_initial("sigma", family_name),
      nu.initial    = .gpig_initial("nu", family_name),
      G.dev.incr = G.dev.incr_fn, rqres = rqres_expr,
      mu.valid    = function(mu) all(mu > 0),
      sigma.valid = function(sigma) all(sigma > 0 & sigma < 1),
      nu.valid    = function(nu) all(nu > 0 & nu < 1),
      y.valid     = function(y) all(y >= 0),
      dldm = d1(1L), dldd = d1(2L), dldv = d1(3L),
      d2ldm2 = d2self(1L), d2ldd2 = d2self(2L), d2ldv2 = d2self(3L),
      d2ldmdd = d2cross(1L, 2L), d2ldmdv = d2cross(1L, 3L), d2ldddv = d2cross(2L, 3L))

  } else {
    tau_link <- tau.link
    gamlss.dist::checklink("tau.link", family_name, substitute(tau_link), c("logit", "log", "identity"))
    G.dev.incr_fn <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
      -2 * .(dfun_sym)(y, mu = mu, sigma = sigma, nu = nu, tau = tau, log = TRUE) }))

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
      mu.initial    = .gpig_initial("mu", family_name),
      sigma.initial = .gpig_initial("sigma", family_name),
      nu.initial    = .gpig_initial("nu", family_name),
      tau.initial   = .gpig_initial("tau", family_name),
      G.dev.incr = G.dev.incr_fn, rqres = rqres_expr,
      mu.valid    = function(mu) all(mu > 0),
      sigma.valid = function(sigma) all(sigma > 0 & sigma < 1),
      nu.valid    = function(nu) all(nu > 0 & nu < 1),
      tau.valid   = function(tau) all(tau >= 0 & tau < 1),
      y.valid     = function(y) all(y >= 0),
      dldm = d1(1L), dldd = d1(2L), dldv = d1(3L), dldt = d1(4L),
      d2ldm2 = d2self(1L), d2ldd2 = d2self(2L), d2ldv2 = d2self(3L), d2ldt2 = d2self(4L),
      d2ldmdd = d2cross(1L, 2L), d2ldmdv = d2cross(1L, 3L), d2ldmdt = d2cross(1L, 4L),
      d2ldddv = d2cross(2L, 3L), d2ldddt = d2cross(2L, 4L), d2ldvdt = d2cross(3L, 4L))
  }
  class(fam) <- "gamlss.family"
  fam
}

## Starting values for GPIG() and ZIGPIG(): the maximum likelihood estimate of the
## intercept-only model, as for the EGPD families (.egpd_gamlss_start). From the fixed
## starts (mean(y), a = 0.5, c = 0.7) the RS algorithm can drift to a -> 0 on heavy-tailed
## data and stop far below the maximum. The estimate maximises the native GPIG likelihood
## (src/gpig.cpp) with its exact score, under the default gpig_control() settings, and the
## fixed starts are kept if that fails or does no better. Returns c(mu, sigma, nu[, tau])
## on the gamlss parameter scale. gamlssML() evaluates each *.initial expression in a fresh
## environment, so the latest result is cached, keyed on the family and the full response.
.gpig_start_cache <- new.env(parent = emptyenv())

.gpig_gamlss_start <- function(y, family) {
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  n <- length(y)
  zi <- identical(family, "ZIGPIG")
  ckey <- list(family)
  hit <- .gpig_start_cache$last
  if (!is.null(hit) && identical(hit$key, ckey) && identical(hit$y, y)) return(hit$est)
  default <- c(mu = max(if (n) mean(y) else 0, 0.1), sigma = 0.5, nu = 0.7)
  if (zi) default <- c(default, tau = min(max(if (n) mean(y == 0) else 0, 0.01), 0.5))
  est <- default
  if (n >= 2L && all(y >= 0 & y == floor(y)) && any(y > 0)) {
    old <- gpig_control()
    on.exit(do.call(gpig_control, old), add = TRUE)
    np <- length(default)
    lik <- c(list(y, integer(0), 0L, rep(list(numeric(0)), np)))
    Xs <- rep(list(matrix(1, n, 1L)), np)
    d0 <- if (zi) zigpig1d0 else gpig1d0
    d12 <- if (zi) zigpig1d12 else gpig1d12
    nll <- function(th) {
      v <- tryCatch(do.call(d0, c(list(as.list(th)), Xs, lik)), error = function(e) NA_real_)
      if (length(v) == 1L && is.finite(v)) v else 1e100
    }
    gr <- function(th) {
      g <- tryCatch(colSums(do.call(d12, c(list(as.list(th)), Xs, lik))[, seq_len(np),
                                                                         drop = FALSE]),
                    error = function(e) rep(NA_real_, np))
      if (all(is.finite(g))) g else stop("non-finite GPIG score")
    }
    th0 <- c(log(default[["mu"]]), stats::qlogis(default[-1L]))
    opt <- tryCatch(
      suppressWarnings(stats::optim(th0, nll, gr, method = "BFGS",
                                    control = list(maxit = 500L, reltol = 1e-10))),
      error = function(e) NULL)
    if (!is.null(opt) && all(is.finite(opt$par)) && opt$value < nll(th0)) {
      p <- c(exp(opt$par[1L]), stats::plogis(opt$par[-1L]))
      p[1L] <- max(p[1L], 1e-8)
      p[2L] <- min(max(p[2L], 1e-6), 1 - 1e-6)
      p[3L] <- min(max(p[3L], 1e-6), 1 - 1e-8)
      if (zi) p[4L] <- min(max(p[4L], 1e-4), 1 - 1e-4)
      if (all(is.finite(p))) est <- stats::setNames(p, names(default))
    }
  }
  .gpig_start_cache$last <- list(key = ckey, y = y, est = est)
  est
}

## Initial-value expression for one parameter, evaluated by gamlss() in its own frame (the
## first one computes the start, the others reuse it) or by gamlssML() in a fresh one.
.gpig_initial <- function(par, family_name) {
  as.expression(bquote({
    if (!exists(".gpig_start", inherits = FALSE))
      .gpig_start <- get(".gpig_gamlss_start", envir = asNamespace("egpd"),
                         mode = "function")(y, .(family_name))
    .(as.name(par)) <- rep(.gpig_start[[.(par)]], length(y))
  }))
}
