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
#' @param n number of draws
#' @param tau scale parameter, \code{tau >= 2}
#' @param log,log.p logical; return the log density / probability
#' @param lower.tail logical; if TRUE probabilities are \code{P[X <= x]}
#' @param max.value largest count searched when inverting the cdf
#' @return a numeric vector
#' @references Alkhairy, I. (2023) Classical and Bayesian inference for the discrete
#'   Poisson Ramos-Louzada distribution with application to COVID-19 data.
#'   \emph{Mathematical Biosciences and Engineering} 20(8), 14061-14080.
#' @name prl
NULL

#' @rdname prl
#' @export
dprl <- function(x, tau = 2.5, log = FALSE) {
  stopifnot(tau >= 2)
  out <- prl_logpmf_cpp(as.numeric(x), tau)
  out[x < 0 | x != floor(x)] <- -Inf
  if (log) out else exp(out)
}
#' @rdname prl
#' @export
pprl <- function(q, tau = 2.5, lower.tail = TRUE, log.p = FALSE) {
  ## Their Eq. (3)/(4): S(x) = P(X > x) = (1+1/tau)^-x tau(x+tau^2) / ((tau-1)(1+tau)^2),
  ## so F(x) = P(X <= x) = 1 - S(x). Note S is the STRICT survivor: P(X=x) = S(x-1) - S(x).
  x <- floor(q)
  s <- (1 + 1 / tau)^(-x) * tau * (x + tau^2) / ((tau - 1) * (1 + tau)^2)
  out <- ifelse(q < 0, 0, pmin(pmax(1 - s, 0), 1))
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}
#' @rdname prl
#' @export
qprl <- function(p, tau = 2.5, lower.tail = TRUE, log.p = FALSE, max.value = 1e7) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  mu <- tau^2 / (tau - 1)
  n <- max(32, ceiling(4 * mu))
  repeat {
    cdf <- cumsum(dprl(0:n, tau))
    if (max(cdf) >= max(p, na.rm = TRUE) - 1e-12 || n >= max.value) break
    n <- min(max.value, n * 4L)
  }
  vapply(p, function(pp) {
    if (is.na(pp)) return(NA_real_)
    i <- which(cdf >= pp - 1e-12)[1]
    if (is.na(i)) NA_real_ else i - 1
  }, 0)
}
#' @rdname prl
#' @export
rprl <- function(n, tau = 2.5) qprl(stats::runif(n), tau)

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
