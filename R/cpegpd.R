## Compound Poisson-EGPD distribution
## Panjer recursion machinery and d/p/q/r functions

#' Panjer recursion for Poisson case (a=0, b=lambda)
#' @param lambda Poisson rate parameter
#' @param f discretized severity PMF vector (f[1] = P(X in [0, h)), etc.)
#' @return g vector of compound distribution PMF
#' @noRd
.panjer_poisson <- function(lambda, f) {
  K <- length(f)
  g <- numeric(K)
  g[1] <- exp(-lambda * (1 - f[1]))

  for (k in 2:K) {
    idx <- seq_len(k - 1)
    s <- sum(idx * f[idx + 1] * g[k - idx])
    g[k] <- (lambda / (k - 1)) * s
  }

  g
}


#' Discretize EGPD severity on grid {0, h, 2h, ..., K*h}
#'
#' Uses the mean-preserving rounding (midpoint) method:
#'   f_0 = F(h/2)
#'   f_k = F((k+0.5)*h) - F((k-0.5)*h) for k >= 1
#' This preserves the mean of the continuous distribution and gives
#' much better Panjer recursion accuracy than the lower discretization.
#'
#' @param sigma GPD scale
#' @param xi GPD shape
#' @param kappa G-transformation parameter
#' @param delta G-transformation parameter
#' @param prob G-transformation parameter
#' @param type integer 1-6
#' @param h bin width
#' @param K number of bins
#' @return numeric vector of length K+1 with bin probabilities
#' @noRd
.discretize_egpd <- function(sigma, xi, kappa = NA, delta = NA, prob = NA,
                              type, h, K) {
  ## Midpoint boundaries: h/2, 3h/2, 5h/2, ..., (K+0.5)*h
  midbreaks <- (seq_len(K + 1) - 0.5) * h
  cdf_mid <- pegpd(midbreaks, prob = prob, kappa = kappa, delta = delta,
                    sigma = sigma, xi = xi, type = type)
  f <- numeric(K + 1)
  f[1] <- cdf_mid[1]                      # F(h/2) - F(0) = F(h/2)
  f[2:(K + 1)] <- diff(cdf_mid)           # F((k+0.5)h) - F((k-0.5)h)
  f <- pmax(f, 0)
  f
}


#' Density (PMF) of Compound Poisson-EGPD
#'
#' Computes the probability mass function of the Compound Poisson-EGPD
#' distribution using Panjer recursion on a discretized EGPD severity.
#'
#' @param x numeric vector of quantiles
#' @param lambda Poisson rate parameter (> 0)
#' @param prob G-transformation parameter (type 6 only)
#' @param kappa G-transformation parameter
#' @param delta G-transformation parameter
#' @param sigma GPD scale parameter (> 0)
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G-transformation type
#' @param h bin width for discretization (default 0.2)
#' @param K number of bins (default: ceiling(max(x)/h) + 50)
#' @param log logical; if TRUE, return log-probabilities
#'
#' @return numeric vector of (log-)probabilities
#' @export
dcpegpd <- function(x, lambda, prob = NA, kappa = NA, delta = NA,
                     sigma, xi, type = 1, h = 0.2, K = NULL,
                     log = FALSE) {
  if (is.null(K)) {
    K <- max(ceiling(max(c(x, 1)) / h) + 50, 100)
  }

  f <- .discretize_egpd(sigma = sigma, xi = xi, kappa = kappa,
                         delta = delta, prob = prob, type = type,
                         h = h, K = K)
  g <- .panjer_poisson(lambda, f)

  bins <- round(x / h) + 1L
  bins <- pmax(bins, 1L)
  bins <- pmin(bins, length(g))

  out <- g[bins]
  out <- pmax(out, .Machine$double.xmin)

  if (log) return(log(out))
  out
}


#' Distribution function of Compound Poisson-EGPD
#'
#' @inheritParams dcpegpd
#' @param q numeric vector of quantiles
#'
#' @return numeric vector of cumulative probabilities
#' @export
pcpegpd <- function(q, lambda, prob = NA, kappa = NA, delta = NA,
                     sigma, xi, type = 1, h = 0.2, K = NULL) {
  if (is.null(K)) {
    K <- max(ceiling(max(c(q, 1)) / h) + 50, 100)
  }

  f <- .discretize_egpd(sigma = sigma, xi = xi, kappa = kappa,
                         delta = delta, prob = prob, type = type,
                         h = h, K = K)
  g <- .panjer_poisson(lambda, f)
  cg <- cumsum(g)

  bins <- floor(q / h) + 1L
  bins <- pmax(bins, 1L)
  bins <- pmin(bins, length(cg))

  cg[bins]
}


#' Quantile function of Compound Poisson-EGPD
#'
#' @inheritParams dcpegpd
#' @param p numeric vector of probabilities
#'
#' @return numeric vector of quantiles
#' @export
qcpegpd <- function(p, lambda, prob = NA, kappa = NA, delta = NA,
                     sigma, xi, type = 1, h = 0.2, K = NULL) {
  if (is.null(K)) K <- 5000

  f <- .discretize_egpd(sigma = sigma, xi = xi, kappa = kappa,
                         delta = delta, prob = prob, type = type,
                         h = h, K = K)
  g <- .panjer_poisson(lambda, f)
  cg <- cumsum(g)

  vapply(p, function(pi) {
    idx <- which(cg >= pi)[1]
    if (is.na(idx)) return(K * h)
    (idx - 1) * h
  }, numeric(1))
}


#' Random generation from Compound Poisson-EGPD
#'
#' Direct simulation: draw N ~ Poisson(lambda), then sum N iid EGPD draws.
#'
#' @inheritParams dcpegpd
#' @param n number of observations
#'
#' @return numeric vector of length n
#' @export
rcpegpd <- function(n, lambda, prob = NA, kappa = NA, delta = NA,
                     sigma, xi, type = 1) {
  N <- rpois(n, lambda)
  vapply(N, function(ni) {
    if (ni == 0) return(0)
    sum(regpd(ni, prob = prob, kappa = kappa, delta = delta,
              sigma = sigma, xi = xi, type = type))
  }, numeric(1))
}


############################################
## Compound Poisson-Discrete EGPD        ##
## (exact integer-valued compound sum)    ##
############################################

#' Density (PMF) of Compound Poisson-Discrete EGPD
#'
#' Computes the probability mass function of the Compound Poisson-Discrete
#' EGPD distribution using Panjer recursion on the exact discrete EGPD
#' severity PMF (no discretization step).
#'
#' The model is \eqn{S = X_1 + \cdots + X_N} where
#' \eqn{N \sim \mathrm{Poisson}(\lambda)} and
#' \eqn{X_i \sim \mathrm{Discrete\mbox{-}EGPD}(\sigma, \xi, \kappa, \ldots)}.
#' Since the severity is already integer-valued, the Panjer recursion
#' produces an exact compound distribution with no discretization error.
#'
#' @param x numeric vector of non-negative integer quantiles
#' @param lambda Poisson rate parameter (> 0)
#' @param prob G-transformation parameter (type 6 only)
#' @param kappa G-transformation parameter
#' @param delta G-transformation parameter
#' @param sigma GPD scale parameter (> 0)
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G-transformation type
#' @param K grid size for Panjer recursion (default: max(max(x), 100) + 50)
#' @param log logical; if TRUE, return log-probabilities
#'
#' @return numeric vector of (log-)probabilities
#' @export
#' @seealso \code{\link{pcpdegpd}}, \code{\link{qcpdegpd}},
#'   \code{\link{rcpdegpd}}, \code{\link{fitegpd}}
dcpdegpd <- function(x, lambda, prob = NA, kappa = NA, delta = NA,
                      sigma, xi, type = 1, K = NULL, log = FALSE) {
  if (is.null(K)) {
    K <- max(max(c(x, 1)), 100) + 50
  }

  ## Exact integer PMF from discrete EGPD
  f <- ddiscegpd(0:K, prob = prob, kappa = kappa, delta = delta,
                  sigma = sigma, xi = xi, type = type)
  f <- pmax(f, 0)

  ## Panjer recursion
  g <- .panjer_poisson(lambda, f)

  ## Map x values to g (R is 1-indexed: g[1] = P(S=0))
  idx <- as.integer(floor(x)) + 1L
  idx <- pmax(idx, 1L)
  idx <- pmin(idx, length(g))

  out <- g[idx]
  out <- pmax(out, .Machine$double.xmin)

  if (log) return(log(out))
  out
}


#' Distribution function of Compound Poisson-Discrete EGPD
#'
#' @inheritParams dcpdegpd
#' @param q numeric vector of quantiles
#'
#' @return numeric vector of cumulative probabilities
#' @export
#' @seealso \code{\link{dcpdegpd}}, \code{\link{qcpdegpd}},
#'   \code{\link{rcpdegpd}}
pcpdegpd <- function(q, lambda, prob = NA, kappa = NA, delta = NA,
                      sigma, xi, type = 1, K = NULL) {
  if (is.null(K)) {
    K <- max(max(c(q, 1)), 100) + 50
  }

  f <- ddiscegpd(0:K, prob = prob, kappa = kappa, delta = delta,
                  sigma = sigma, xi = xi, type = type)
  f <- pmax(f, 0)
  g <- .panjer_poisson(lambda, f)
  cg <- cumsum(g)

  idx <- as.integer(floor(q)) + 1L
  idx <- pmax(idx, 1L)
  idx <- pmin(idx, length(cg))

  cg[idx]
}


#' Quantile function of Compound Poisson-Discrete EGPD
#'
#' @inheritParams dcpdegpd
#' @param p numeric vector of probabilities
#'
#' @return numeric vector of quantiles (non-negative integers)
#' @export
#' @seealso \code{\link{dcpdegpd}}, \code{\link{pcpdegpd}},
#'   \code{\link{rcpdegpd}}
qcpdegpd <- function(p, lambda, prob = NA, kappa = NA, delta = NA,
                      sigma, xi, type = 1, K = NULL) {
  if (is.null(K)) K <- 5000

  f <- ddiscegpd(0:K, prob = prob, kappa = kappa, delta = delta,
                  sigma = sigma, xi = xi, type = type)
  f <- pmax(f, 0)
  g <- .panjer_poisson(lambda, f)
  cg <- cumsum(g)

  vapply(p, function(pi) {
    idx <- which(cg >= pi)[1]
    if (is.na(idx)) return(K)
    idx - 1L
  }, numeric(1))
}


#' Random generation from Compound Poisson-Discrete EGPD
#'
#' Direct simulation: draw \eqn{N \sim \mathrm{Poisson}(\lambda)}, then
#' sum \eqn{N} i.i.d. Discrete-EGPD draws.
#'
#' @inheritParams dcpdegpd
#' @param n number of observations
#'
#' @return integer vector of length \code{n}
#' @export
#' @seealso \code{\link{dcpdegpd}}, \code{\link{pcpdegpd}},
#'   \code{\link{qcpdegpd}}
rcpdegpd <- function(n, lambda, prob = NA, kappa = NA, delta = NA,
                      sigma, xi, type = 1) {
  N <- rpois(n, lambda)
  vapply(N, function(ni) {
    if (ni == 0) return(0L)
    sum(rdiscegpd(ni, prob = prob, kappa = kappa, delta = delta,
                   sigma = sigma, xi = xi, type = type))
  }, numeric(1))
}
