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
