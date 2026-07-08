## Distribution functions for the Poisson-inverse Gaussian (PIG) and
## zero-inflated PIG (ZIPIG) count distributions.
##
## These are thin wrappers around the corresponding gamlss.dist functions,
## provided for naming consistency with the other egpd distributions (e.g.
## ddiscegpd/pdiscegpd/qdiscegpd). The mixed-Poisson PMF/CDF/quantile maths is
## deliberately NOT reimplemented: gamlss.dist (already a hard dependency)
## supplies mature, numerically careful versions. The native GAM likelihood in
## src/pig.cpp has its own (independently tested) implementation.

#' The Poisson-inverse Gaussian distribution
#'
#' Density (PMF), distribution function, quantile function and random
#' generation for the Poisson-inverse Gaussian (PIG) distribution: a
#' mixed-Poisson model for overdispersed counts with mean \code{mu} and
#' dispersion \code{sigma} (variance \eqn{\mu + \sigma\mu^2}). These are thin
#' wrappers around \code{\link[gamlss.dist]{dPIG}} etc.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu vector of positive means.
#' @param sigma vector of positive dispersion parameters.
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \eqn{p} are
#'   given as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param ... further arguments passed to the underlying \code{gamlss.dist}
#'   function (e.g. \code{max.value} for \code{qpig}).
#' @return \code{dpig} gives the PMF, \code{ppig} the CDF, \code{qpig} the
#'   quantile function, and \code{rpig} random deviates.
#' @seealso \code{\link[gamlss.dist]{PIG}} for the gamlss family;
#'   \code{\link{dzipig}} for the zero-inflated version.
#' @examples
#' dpig(0:5, mu = 2, sigma = 0.8)
#' qpig(ppig(3, mu = 2, sigma = 0.8), mu = 2, sigma = 0.8)  # == 3
#' @name pig
NULL

#' @rdname pig
#' @export
dpig <- function(x, mu = 1, sigma = 1, log = FALSE)
  gamlss.dist::dPIG(x, mu = mu, sigma = sigma, log = log)

#' @rdname pig
#' @export
ppig <- function(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
  gamlss.dist::pPIG(q, mu = mu, sigma = sigma, lower.tail = lower.tail, log.p = log.p)

#' @rdname pig
#' @export
qpig <- function(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE, ...)
  gamlss.dist::qPIG(p, mu = mu, sigma = sigma, lower.tail = lower.tail, log.p = log.p, ...)

#' @rdname pig
#' @export
rpig <- function(n, mu = 1, sigma = 1)
  gamlss.dist::rPIG(n, mu = mu, sigma = sigma)


#' The zero-inflated Poisson-inverse Gaussian distribution
#'
#' Density (PMF), distribution function, quantile function and random
#' generation for the zero-inflated Poisson-inverse Gaussian (ZIPIG)
#' distribution: a PIG (mean \code{mu}, dispersion \code{sigma}) mixed with a
#' point mass at zero of probability \code{pi}. These are thin wrappers around
#' \code{\link[gamlss.dist]{dZIPIG}} etc. (whose zero-inflation argument is
#' named \code{nu}).
#'
#' @inheritParams pig
#' @param pi vector of zero-inflation probabilities in \eqn{(0, 1)} (passed to
#'   \code{gamlss.dist} as \code{nu}).
#' @return \code{dzipig} gives the PMF, \code{pzipig} the CDF, \code{qzipig} the
#'   quantile function, and \code{rzipig} random deviates.
#' @seealso \code{\link[gamlss.dist]{ZIPIG}} for the gamlss family;
#'   \code{\link{dpig}} for the non-inflated version.
#' @name zipig
NULL

#' @rdname zipig
#' @export
dzipig <- function(x, mu = 1, sigma = 1, pi = 0.1, log = FALSE)
  gamlss.dist::dZIPIG(x, mu = mu, sigma = sigma, nu = pi, log = log)

#' @rdname zipig
#' @export
pzipig <- function(q, mu = 1, sigma = 1, pi = 0.1, lower.tail = TRUE, log.p = FALSE)
  gamlss.dist::pZIPIG(q, mu = mu, sigma = sigma, nu = pi,
                      lower.tail = lower.tail, log.p = log.p)

#' @rdname zipig
#' @export
qzipig <- function(p, mu = 1, sigma = 1, pi = 0.1, lower.tail = TRUE, log.p = FALSE, ...)
  gamlss.dist::qZIPIG(p, mu = mu, sigma = sigma, nu = pi,
                      lower.tail = lower.tail, log.p = log.p, ...)

#' @rdname zipig
#' @export
rzipig <- function(n, mu = 1, sigma = 1, pi = 0.1)
  gamlss.dist::rZIPIG(n, mu = mu, sigma = sigma, nu = pi)
