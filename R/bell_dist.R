## Distribution functions for the Bell (Castellares, Ferrari & Lemonte 2018,
## Appl. Math. Modelling 56:172-185) and zero-inflated Bell (ZIBell) count
## distributions.
##
## The pmf is computed from the closed form via src/bell.cpp (bell_pmf_cpp),
## using log(B_k/k!) log-space ratios for numerical stability. The natural
## theta parameterisation of the paper (Definition 1) is used: theta > 0, with
##   Pr(Y = y) = theta^y exp(1 - e^theta) B_y / y!.
## Mean = theta e^theta, variance = theta(1+theta) e^theta.

## Recycle args to a common length.
.bell_recycle <- function(...) {
  a <- list(...)
  n <- max(lengths(a))
  lapply(a, rep_len, n)
}

## pmf at (integer) x-values, scalar or vector theta.
.bell_pmf_at <- function(x, theta) {
  z <- .bell_recycle(x = x, theta = theta)
  x <- z$x; theta <- z$theta
  n <- length(x)
  out <- numeric(n)
  ok <- is.finite(x) & x >= 0 & x == floor(x) & is.finite(theta) & theta > 0
  scalarpar <- length(unique(theta[ok])) <= 1
  if (any(ok) && scalarpar) {
    i0 <- which(ok)[1]
    xm <- max(x[ok])
    pmf <- bell_pmf_cpp(as.integer(xm), theta[i0])
    out[ok] <- pmf[x[ok] + 1L]
  } else if (any(ok)) {
    for (i in which(ok))
      out[i] <- bell_pmf_cpp(as.integer(x[i]), theta[i])[x[i] + 1L]
  }
  out
}

## CDF at (integer) q-values.
.bell_cdf_at <- function(q, theta) {
  z <- .bell_recycle(q = q, theta = theta)
  q <- z$q; theta <- z$theta
  n <- length(q); out <- numeric(n)
  qf <- floor(q)
  scalarpar <- length(unique(theta)) <= 1
  if (scalarpar && n > 0) {
    hi <- max(qf[is.finite(qf)], 0)
    cdf <- cumsum(bell_pmf_cpp(as.integer(max(hi, 0)), theta[1]))
    for (i in seq_len(n)) {
      if (!is.finite(qf[i])) { out[i] <- NA; next }
      out[i] <- if (qf[i] < 0) 0 else min(cdf[qf[i] + 1L], 1)
    }
  } else {
    for (i in seq_len(n)) {
      if (!is.finite(qf[i])) { out[i] <- NA; next }
      if (qf[i] < 0) { out[i] <- 0; next }
      out[i] <- min(sum(bell_pmf_cpp(as.integer(qf[i]), theta[i])), 1)
    }
  }
  out
}

## scalar quantile: smallest k with F(k) >= p. A small tolerance keeps discrete
## round-trips q(p(k)) == k robust to floating-point error (e.g. the ZI
## probability rescaling (p - pi)/(1 - pi)).
.bell_q1 <- function(p, theta, max.value, tol = 1e-12) {
  if (!is.finite(p) || p <= 0) return(0)
  if (p >= 1) return(max.value)
  N <- 64L
  repeat {
    cdf <- cumsum(bell_pmf_cpp(as.integer(N), theta))
    if (cdf[length(cdf)] >= p - tol || N >= max.value) break
    N <- min(N * 2L, as.integer(max.value))
  }
  k <- which(cdf >= p - tol)[1] - 1L
  if (is.na(k)) as.integer(N) else as.integer(k)
}

#' The Bell distribution
#'
#' Density (PMF), distribution function, quantile function and random generation
#' for the Bell distribution of Castellares, Ferrari and Lemonte (2018), a
#' one-parameter overdispersed count family with
#' \eqn{\Pr(Y = y) = \theta^y e^{1 - e^\theta} B_y / y!}, \eqn{y = 0, 1, 2,
#' \ldots}, where \eqn{B_y} are the Bell numbers. The mean is
#' \eqn{\theta e^\theta} and the variance \eqn{\theta(1+\theta)e^\theta}, so the
#' index of dispersion is \eqn{1 + \theta > 1} (always overdispersed).
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param theta vector of positive shape parameters.
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \eqn{p} are
#'   given as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param max.value integer cap on the support searched by \code{qbell}
#'   (default \code{1e5}).
#' @return \code{dbell} gives the PMF, \code{pbell} the CDF, \code{qbell} the
#'   quantile function, and \code{rbell} random deviates.
#' @references Castellares, F., Ferrari, S. L. P. and Lemonte, A. J. (2018). On
#'   the Bell distribution and its associated regression model for count data.
#'   \emph{Applied Mathematical Modelling} 56, 172-185.
#' @seealso \code{\link{zibell}} for the zero-inflated version; \code{\link{BELL}}
#'   for the mean-parameterised \code{gamlss} family.
#' @examples
#' dbell(0:5, theta = 0.7)
#' qbell(pbell(3, theta = 0.7), theta = 0.7)  # == 3
#' @name bell
NULL

#' @rdname bell
#' @export
dbell <- function(x, theta = 1, log = FALSE) {
  out <- .bell_pmf_at(x, theta)
  if (log) log(out) else out
}

#' @rdname bell
#' @export
pbell <- function(q, theta = 1, lower.tail = TRUE, log.p = FALSE) {
  out <- .bell_cdf_at(q, theta)
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}

#' @rdname bell
#' @export
qbell <- function(p, theta = 1, lower.tail = TRUE, log.p = FALSE, max.value = 1e5) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  z <- .bell_recycle(p = p, theta = theta)
  mapply(.bell_q1, z$p, z$theta, MoreArgs = list(max.value = max.value))
}

#' @rdname bell
#' @export
rbell <- function(n, theta = 1) {
  if (length(n) > 1L) n <- length(n)
  qbell(stats::runif(n), theta = theta)
}


#' The zero-inflated Bell distribution
#'
#' Density (PMF), distribution function, quantile function and random generation
#' for the zero-inflated Bell (ZIBell) distribution: a \code{\link{bell}}
#' distribution mixed with a point mass at zero of probability \code{pi}.
#'
#' @inheritParams bell
#' @param pi vector of zero-inflation probabilities in \eqn{[0, 1)}.
#' @return \code{dzibell} gives the PMF, \code{pzibell} the CDF, \code{qzibell}
#'   the quantile function, and \code{rzibell} random deviates.
#' @seealso \code{\link{bell}} for the non-inflated version.
#' @name zibell
NULL

#' @rdname zibell
#' @export
dzibell <- function(x, theta = 1, pi = 0.1, log = FALSE) {
  z <- .bell_recycle(x = x, theta = theta, pi = pi)
  base <- .bell_pmf_at(z$x, z$theta)
  is0 <- is.finite(z$x) & z$x == 0
  out <- (1 - z$pi) * base
  out[is0] <- z$pi[is0] + (1 - z$pi[is0]) * base[is0]
  if (log) log(out) else out
}

#' @rdname zibell
#' @export
pzibell <- function(q, theta = 1, pi = 0.1, lower.tail = TRUE, log.p = FALSE) {
  z <- .bell_recycle(q = q, theta = theta, pi = pi)
  base <- .bell_cdf_at(z$q, z$theta)
  out <- z$pi + (1 - z$pi) * base
  out[is.finite(z$q) & floor(z$q) < 0] <- 0
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}

#' @rdname zibell
#' @export
qzibell <- function(p, theta = 1, pi = 0.1, lower.tail = TRUE, log.p = FALSE,
                    max.value = 1e5) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  z <- .bell_recycle(p = p, theta = theta, pi = pi)
  ## invert the mixture: for p <= pi the quantile is 0, else invert the Bell
  ## part at the rescaled probability (p - pi)/(1 - pi).
  padj <- ifelse(z$p <= z$pi, 0, (z$p - z$pi) / (1 - z$pi))
  mapply(function(pp, th) if (pp <= 0) 0L else .bell_q1(pp, th, max.value),
         padj, z$theta)
}

#' @rdname zibell
#' @export
rzibell <- function(n, theta = 1, pi = 0.1) {
  if (length(n) > 1L) n <- length(n)
  z0 <- stats::runif(n) < pi
  out <- integer(n)
  if (any(!z0))
    out[!z0] <- qbell(stats::runif(sum(!z0)), theta = rep_len(theta, n)[!z0])
  out
}
