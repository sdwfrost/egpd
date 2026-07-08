## Distribution functions for the generalised Poisson-inverse Gaussian (GPIG;
## Zhu & Joe 2009, Stat. Prob. Lett. 79:1695-1703) and zero-inflated GPIG
## (ZIGPIG) count distributions.
##
## Unlike PIG, gamlss.dist has no GPIG, so the pmf is computed here from the
## paper's forward recursion (src/gpig.cpp, gpig_pmf_cpp). The natural (a,b,c)
## parameterisation of the paper is used: 0 < a <= 1 (tail exponent), b > 0
## (level), 0 <= c <= 1 (geometric down-weighting). Special cases: Poisson
## (a=1, mean bc), PIG (a=1/2), discrete stable (c=1).

## Recycle args to a common length.
.gpig_recycle <- function(...) {
  a <- list(...)
  n <- max(lengths(a))
  lapply(a, rep_len, n)
}

## pmf at (integer) x-values, scalar or vector (a,b,c).
.gpig_pmf_at <- function(x, a, b, c) {
  z <- .gpig_recycle(x = x, a = a, b = b, c = c)
  x <- z$x; a <- z$a; b <- z$b; c <- z$c
  n <- length(x)
  out <- numeric(n)
  ok <- is.finite(x) & x >= 0 & x == floor(x)
  scalarpar <- (length(unique(a[ok])) <= 1 && length(unique(b[ok])) <= 1 &&
                length(unique(c[ok])) <= 1)
  if (any(ok) && scalarpar) {
    i0 <- which(ok)[1]
    xm <- max(x[ok])
    pmf <- gpig_pmf_cpp(as.integer(xm), a[i0], b[i0], c[i0])
    out[ok] <- pmf[x[ok] + 1L]
  } else if (any(ok)) {
    for (i in which(ok))
      out[i] <- gpig_pmf_cpp(as.integer(x[i]), a[i], b[i], c[i])[x[i] + 1L]
  }
  out
}

## CDF at (integer) q-values.
.gpig_cdf_at <- function(q, a, b, c) {
  z <- .gpig_recycle(q = q, a = a, b = b, c = c)
  q <- z$q; a <- z$a; b <- z$b; c <- z$c
  n <- length(q); out <- numeric(n)
  qf <- floor(q)
  scalarpar <- (length(unique(a)) <= 1 && length(unique(b)) <= 1 &&
                length(unique(c)) <= 1)
  if (scalarpar && n > 0) {
    hi <- max(qf[is.finite(qf)], 0)
    cdf <- cumsum(gpig_pmf_cpp(as.integer(max(hi, 0)), a[1], b[1], c[1]))
    for (i in seq_len(n)) {
      if (!is.finite(qf[i])) { out[i] <- NA; next }
      out[i] <- if (qf[i] < 0) 0 else min(cdf[qf[i] + 1L], 1)
    }
  } else {
    for (i in seq_len(n)) {
      if (!is.finite(qf[i])) { out[i] <- NA; next }
      if (qf[i] < 0) { out[i] <- 0; next }
      out[i] <- min(sum(gpig_pmf_cpp(as.integer(qf[i]), a[i], b[i], c[i])), 1)
    }
  }
  out
}

## scalar quantile: smallest k with F(k) >= p. A small tolerance keeps discrete
## round-trips q(p(k)) == k robust to floating-point error (e.g. the ZI
## probability rescaling (p - pi)/(1 - pi)).
.gpig_q1 <- function(p, a, b, c, max.value, tol = 1e-9) {
  if (!is.finite(p) || p <= 0) return(0)
  if (p >= 1) return(max.value)
  N <- 64L
  repeat {
    cdf <- cumsum(gpig_pmf_cpp(as.integer(N), a, b, c))
    if (cdf[length(cdf)] >= p - tol || N >= max.value) break
    N <- min(N * 2L, as.integer(max.value))
  }
  k <- which(cdf >= p - tol)[1] - 1L
  if (is.na(k)) as.integer(N) else as.integer(k)
}

#' The generalised Poisson-inverse Gaussian distribution
#'
#' Density (PMF), distribution function, quantile function and random
#' generation for the generalised Poisson-inverse Gaussian (GPIG) distribution
#' of Zhu and Joe (2009), a three-parameter heavy-tailed count family with
#' probability generating function
#' \eqn{G(s) = \exp\{b[(1-c)^a - (1-cs)^a]\}}. It nests the Poisson
#' (\eqn{a=1}, mean \eqn{bc}), the Poisson-inverse Gaussian (\eqn{a=1/2}) and
#' the discrete-stable (\eqn{c=1}) distributions. The PMF has no closed form and
#' is evaluated by the paper's forward recursion.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param a vector of tail exponents in \eqn{(0, 1]}.
#' @param b vector of positive level parameters.
#' @param c vector of down-weighting parameters in \eqn{[0, 1]}.
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \eqn{p} are
#'   given as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param max.value integer cap on the support searched by \code{qgpig}
#'   (default \code{1e5}); relevant for very heavy tails.
#' @return \code{dgpig} gives the PMF, \code{pgpig} the CDF, \code{qgpig} the
#'   quantile function, and \code{rgpig} random deviates.
#' @references Zhu, R. and Joe, H. (2009). Modelling heavy-tailed count data
#'   using a generalised Poisson-inverse Gaussian family. \emph{Statistics and
#'   Probability Letters} 79, 1695-1703.
#' @seealso \code{\link{dpig}} for the (a=1/2) Poisson-inverse Gaussian special
#'   case; \code{\link{dzigpig}} for the zero-inflated version.
#' @examples
#' dgpig(0:5, a = 0.4, b = 2.4, c = 0.9)
#' qgpig(pgpig(3, a = 0.4, b = 2.4, c = 0.9), a = 0.4, b = 2.4, c = 0.9)  # == 3
#' @name gpig
NULL

#' @rdname gpig
#' @export
dgpig <- function(x, a = 0.5, b = 1, c = 0.5, log = FALSE) {
  out <- .gpig_pmf_at(x, a, b, c)
  if (log) log(out) else out
}

#' @rdname gpig
#' @export
pgpig <- function(q, a = 0.5, b = 1, c = 0.5, lower.tail = TRUE, log.p = FALSE) {
  out <- .gpig_cdf_at(q, a, b, c)
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}

#' @rdname gpig
#' @export
qgpig <- function(p, a = 0.5, b = 1, c = 0.5, lower.tail = TRUE, log.p = FALSE,
                  max.value = 1e5) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  z <- .gpig_recycle(p = p, a = a, b = b, c = c)
  mapply(.gpig_q1, z$p, z$a, z$b, z$c, MoreArgs = list(max.value = max.value))
}

#' @rdname gpig
#' @export
rgpig <- function(n, a = 0.5, b = 1, c = 0.5) {
  if (length(n) > 1L) n <- length(n)
  qgpig(stats::runif(n), a = a, b = b, c = c)
}


#' The zero-inflated generalised Poisson-inverse Gaussian distribution
#'
#' Density (PMF), distribution function, quantile function and random
#' generation for the zero-inflated generalised Poisson-inverse Gaussian
#' (ZIGPIG) distribution: a \code{\link{gpig}} distribution mixed with a point
#' mass at zero of probability \code{pi}.
#'
#' @inheritParams gpig
#' @param pi vector of zero-inflation probabilities in \eqn{[0, 1)}.
#' @return \code{dzigpig} gives the PMF, \code{pzigpig} the CDF, \code{qzigpig}
#'   the quantile function, and \code{rzigpig} random deviates.
#' @seealso \code{\link{gpig}} for the non-inflated version.
#' @name zigpig
NULL

#' @rdname zigpig
#' @export
dzigpig <- function(x, a = 0.5, b = 1, c = 0.5, pi = 0.1, log = FALSE) {
  z <- .gpig_recycle(x = x, a = a, b = b, c = c, pi = pi)
  base <- .gpig_pmf_at(z$x, z$a, z$b, z$c)
  is0 <- is.finite(z$x) & z$x == 0
  out <- (1 - z$pi) * base
  out[is0] <- z$pi[is0] + (1 - z$pi[is0]) * base[is0]
  if (log) log(out) else out
}

#' @rdname zigpig
#' @export
pzigpig <- function(q, a = 0.5, b = 1, c = 0.5, pi = 0.1,
                    lower.tail = TRUE, log.p = FALSE) {
  z <- .gpig_recycle(q = q, a = a, b = b, c = c, pi = pi)
  base <- .gpig_cdf_at(z$q, z$a, z$b, z$c)
  out <- z$pi + (1 - z$pi) * base
  out[is.finite(z$q) & floor(z$q) < 0] <- 0
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}

#' @rdname zigpig
#' @export
qzigpig <- function(p, a = 0.5, b = 1, c = 0.5, pi = 0.1,
                    lower.tail = TRUE, log.p = FALSE, max.value = 1e5) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  z <- .gpig_recycle(p = p, a = a, b = b, c = c, pi = pi)
  ## invert the mixture: for p <= pi the quantile is 0, else invert the GPIG
  ## part at the rescaled probability (p - pi)/(1 - pi).
  padj <- ifelse(z$p <= z$pi, 0, (z$p - z$pi) / (1 - z$pi))
  mapply(function(pp, aa, bb, cc) if (pp <= 0) 0L else .gpig_q1(pp, aa, bb, cc, max.value),
         padj, z$a, z$b, z$c)
}

#' @rdname zigpig
#' @export
rzigpig <- function(n, a = 0.5, b = 1, c = 0.5, pi = 0.1) {
  if (length(n) > 1L) n <- length(n)
  z0 <- stats::runif(n) < pi
  out <- integer(n)
  if (any(!z0))
    out[!z0] <- qgpig(stats::runif(sum(!z0)),
                      a = rep_len(a, n)[!z0], b = rep_len(b, n)[!z0],
                      c = rep_len(c, n)[!z0])
  out
}
