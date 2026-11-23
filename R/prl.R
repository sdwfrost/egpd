## Poisson Ramos-Louzada (PRL), Alkhairy (2023). One parameter, tau >= 2. See src/prl.cpp
## for the pmf and the two structural facts that make it worth fitting: the mean cannot go
## below 4, and the dispersion index is pinned to the mean (DI ~ mu).

#' Poisson Ramos-Louzada distribution
#'
#' The one-parameter mixed Poisson of Alkhairy (2023): \code{X | theta ~ Poisson(theta)}
#' with \code{theta} following a Ramos-Louzada distribution with scale \code{tau >= 2}.
#' Its mean is \code{tau^2/(tau-1)} -- minimised at \code{4}, so no smaller mean is
#' representable -- and its dispersion index is \code{(tau^2+tau-3)/(tau-1)}, approximately
#' the mean. It is therefore geometric-like: strictly stronger than the Bell distribution
#' (whose dispersion index is about 1), but with no free dispersion parameter.
#'
#' @param x,q vector of counts
#' @param p vector of probabilities
#' @param n number of draws; if \code{length(n) != 1}, the length is taken to be
#'   the number required (\code{numeric(0)} gives none), as in \code{\link[stats]{rpois}}
#' @param tau scale parameter, \code{tau >= 2}
#' @param log,log.p logical; return the log density / probability
#' @param lower.tail logical; if TRUE probabilities are \code{P[X <= x]}
#' @param max.value largest count searched when inverting the cdf; a nonattaining
#'   cap returns \code{NA} with a warning.
#' @details The survivor has the closed form of their Eq. (3)/(4), so both tails keep full
#'   relative accuracy. Invalid input follows \code{\link[stats]{dpois}}: an \code{NA}
#'   parameter gives \code{NA}; an invalid one (\code{tau < 2}, \code{NaN}, infinite)
#'   gives \code{NaN} with a warning; \code{x} that is negative or infinite has mass zero,
#'   and a non-integer \code{x} has mass zero with a warning; probabilities outside
#'   \eqn{[0, 1]} give \code{NaN} with a warning. Quantile endpoints are 0 and \code{Inf}.
#' @return a numeric vector
#' @references Alkhairy, I. (2023) Classical and Bayesian inference for the discrete
#'   Poisson Ramos-Louzada distribution with application to COVID-19 data.
#'   \emph{Mathematical Biosciences and Engineering} 20(8), 14061-14080.
#' @name prl
NULL

#' @rdname prl
#' @export
dprl <- function(x, tau = 2.5, log = FALSE) {
  .count_d(.prl_fam, x, list(tau = tau), log = log)
}
#' @rdname prl
#' @export
pprl <- function(q, tau = 2.5, lower.tail = TRUE, log.p = FALSE) {
  .count_p(.prl_fam, q, list(tau = tau), lower.tail, log.p)
}
#' @rdname prl
#' @export
qprl <- function(p, tau = 2.5, lower.tail = TRUE, log.p = FALSE, max.value = 1e7) {
  .count_q(.prl_fam, p, list(tau = tau), lower.tail, log.p, max.value)
}
#' @rdname prl
#' @export
rprl <- function(n, tau = 2.5) {
  .count_r(.prl_fam, n, list(tau = tau), max.value = 1e7)
}

## Their Eq. (3)/(4): S(x) = P(X > x) = (1+1/tau)^-x tau (x+tau^2) / ((tau-1)(1+tau)^2), the
## STRICT survivor (P(X = x) = S(x-1) - S(x)), used at every x.
.prl_lsurv <- function(x, tau)
  pmin(-x * log1p(1 / tau) + log(tau) + log(x + tau^2) - log(tau - 1) - 2 * log1p(tau), 0)

.prl_fam <- list(
  pars = "tau",
  valid = function(p) is.finite(p$tau) & p$tau >= 2,
  lpmf = function(n, p) prl_logpmf_cpp(0:n, p$tau),
  lpmf_at = function(x, p) prl_logpmf_cpp(as.numeric(x), p$tau),
  lsurv = function(k, p) .prl_lsurv(k, p$tau),
  surv_all = TRUE,
  init = function(p) max(32, ceiling(4 * p$tau^2 / (p$tau - 1))))

## mean <-> tau, used for starting values and by the fitter's clamp
.prl_tau_from_mu <- function(mu) {
  mu <- pmax(mu, 4 + 1e-4)
  (mu + sqrt(mu^2 - 4 * mu)) / 2
}

.iG_prl <- function(v, mu) qprl(v, tau = .prl_tau_from_mu(mu))

.prl1.d0 <- function(pars, likdata) {
  if (likdata$censored) stop("Censored likelihoods not available for prl.")
  prl1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
         likdata$X[[1]], likdata$y[, 1], likdata$dupid, likdata$duplicate, likdata$offsets)
}
.prl1.d12 <- function(pars, likdata) {
  prl1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
          likdata$X[[1]], likdata$y[, 1], likdata$dupid, likdata$duplicate, likdata$offsets)
}
.prl1fns <- list(d0 = .prl1.d0, d120 = .prl1.d12, d340 = NULL, m = 1, iG = .iG_prl)
