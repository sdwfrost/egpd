## Distribution functions for the Poisson-inverse Gaussian (PIG) and
## zero-inflated PIG (ZIPIG) count distributions, in the parameterisation of
## gamlss.dist::PIG (mean mu, dispersion sigma, Var = mu + sigma mu^2).
##
## The pmf is computed natively, by the same cancellation-free recursion that the GAM
## likelihood uses (src/pig.cpp, pig_logpmf_ladder_cpp), for every sigma: it is exact to
## rounding from the Poisson limit to heavy overdispersion. The tail decays geometrically
## with limiting ratio 2 sigma mu / (1 + 2 sigma mu), which bounds the remainder of the
## survivor sums (R/countfams.R); it is 1 where 2 sigma mu overflows.
.pig_fam <- list(
  pars = c("mu", "sigma"),
  valid = function(p) is.finite(p$mu) & p$mu > 0 & is.finite(p$sigma) & p$sigma > 0,
  lpmf = function(n, p) pig_logpmf_ladder_cpp(as.integer(n), p$mu, p$sigma),
  rlim = function(p) { s <- 2 * p$sigma * p$mu; if (is.finite(s)) s / (1 + s) else 1 },
  init = function(p) max(32, ceiling(4 * p$mu)))

#' The Poisson-inverse Gaussian distribution
#'
#' Density (PMF), distribution function, quantile function and random
#' generation for the Poisson-inverse Gaussian (PIG) distribution: a
#' mixed-Poisson model for overdispersed counts with mean \code{mu} and
#' dispersion \code{sigma} (variance \eqn{\mu + \sigma\mu^2}), in the
#' parameterisation of \code{\link[gamlss.dist]{PIG}}. The PMF is computed natively
#' by a cancellation-free forward recursion that is exact to rounding for every
#' \code{sigma}, including the Poisson limit \code{sigma -> 0}.
#'
#' @details Both tails are computed from sums of the PMF taken in the direction that
#'   keeps full relative accuracy (the upper tail from the survivor itself), so
#'   \code{lower.tail = FALSE} and \code{log.p = TRUE} stay accurate far into the tail.
#'   The quantile function inverts these probabilities up to rounding tolerance;
#'   near one it compares on the log scale to resolve the small upper tail. Its
#'   endpoints are 0 and \code{Inf}, and a nonattaining \code{max.value} gives \code{NA}
#'   with a warning. Invalid input follows \code{\link[stats]{dpois}}: an \code{NA}
#'   parameter gives \code{NA}; an invalid one (non-positive, \code{NaN} or infinite)
#'   gives \code{NaN} with a warning, except that an \code{NA} or \code{NaN} \code{x},
#'   \code{q} or \code{p} is returned silently; \code{x} that is negative or infinite has mass
#'   zero, and a non-integer \code{x} has mass zero with a warning; probabilities
#'   outside \eqn{[0, 1]} give \code{NaN} with a warning.
#'
#'   Random generation is exact and has no cap: \code{rpig} draws the mixing variable
#'   \eqn{Z}, inverse Gaussian with mean 1 and variance \code{sigma}, by the method of
#'   Michael, Schucany and Haas (1976), and then a Poisson count with mean
#'   \eqn{\mu Z}.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return; if \code{length(n) != 1}, the length is
#'   taken to be the number required (\code{numeric(0)} gives none), as in
#'   \code{\link[stats]{rpois}}.
#' @param mu vector of positive means.
#' @param sigma vector of positive dispersion parameters.
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \eqn{p} are
#'   given as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param max.value integer cap on the support searched by \code{qpig}
#'   (default \code{1e5}). It does not limit \code{rpig}, which samples exactly.
#' @return \code{dpig} gives the PMF, \code{ppig} the CDF, \code{qpig} the
#'   quantile function, and \code{rpig} random deviates.
#' @references Michael, J. R., Schucany, W. R. and Haas, R. W. (1976). Generating
#'   random variates using transformations with multiple roots. \emph{The American
#'   Statistician} 30, 88-90.
#' @seealso \code{\link[gamlss.dist]{PIG}} for the gamlss family;
#'   \code{\link{dzipig}} for the zero-inflated version.
#' @examples
#' dpig(0:5, mu = 2, sigma = 0.8)
#' qpig(ppig(3, mu = 2, sigma = 0.8), mu = 2, sigma = 0.8)  # == 3
#' @name pig
NULL

#' @rdname pig
#' @export
dpig <- function(x, mu = 1, sigma = 1, log = FALSE) {
  .count_d(.pig_fam, x, list(mu = mu, sigma = sigma), log = log)
}

#' @rdname pig
#' @export
ppig <- function(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .count_p(.pig_fam, q, list(mu = mu, sigma = sigma), lower.tail, log.p)
}

#' @rdname pig
#' @export
qpig <- function(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE, max.value = 1e5) {
  .count_q(.pig_fam, p, list(mu = mu, sigma = sigma), lower.tail, log.p, max.value)
}

#' @rdname pig
#' @export
rpig <- function(n, mu = 1, sigma = 1) {
  .count_r_exact(.pig_fam, n, list(mu = mu, sigma = sigma), .rpig_exact)
}

## Exact PIG draws from the mixture Y | Z ~ Poisson(mu Z), where Z is inverse Gaussian with
## mean 1 and variance sigma (shape 1 / sigma), so that Var Y = mu + sigma mu^2. Z is drawn
## by Michael, Schucany and Haas (1976): with w = sigma chi^2_1, the smaller root of their
## quadratic, 1 + w/2 - sqrt(w + w^2/4), is formed as its reciprocal-conjugate (no
## cancellation), and is kept with probability 1 / (1 + x), otherwise replaced by 1 / x.
.rpig_exact <- function(p) {
  m <- length(p$mu)
  w <- p$sigma * stats::rnorm(m)^2
  x <- 1 / (1 + w / 2 + sqrt(w) * sqrt(1 + w / 4))
  z <- ifelse(stats::runif(m) <= 1 / (1 + x), x, 1 / x)
  stats::rpois(m, p$mu * z)
}


#' The zero-inflated Poisson-inverse Gaussian distribution
#'
#' Density (PMF), distribution function, quantile function and random
#' generation for the zero-inflated Poisson-inverse Gaussian (ZIPIG)
#' distribution: a PIG (mean \code{mu}, dispersion \code{sigma}) mixed with a
#' point mass at zero of probability \code{pi} (the \code{nu} of
#' \code{\link[gamlss.dist]{ZIPIG}}). All four functions are computed natively from
#' the PIG recursion; tails, quantiles and invalid input are handled as for
#' \code{\link{pig}}, and \code{pi} outside \eqn{[0, 1)} is invalid. Random generation
#' draws the structural zero with probability \code{pi} and otherwise an exact PIG
#' draw (see \code{\link{rpig}}).
#'
#' @inheritParams pig
#' @param pi vector of zero-inflation probabilities in \eqn{[0, 1)}.
#' @return \code{dzipig} gives the PMF, \code{pzipig} the CDF, \code{qzipig} the
#'   quantile function, and \code{rzipig} random deviates.
#' @seealso \code{\link[gamlss.dist]{ZIPIG}} for the gamlss family;
#'   \code{\link{dpig}} for the non-inflated version.
#' @name zipig
NULL

#' @rdname zipig
#' @export
dzipig <- function(x, mu = 1, sigma = 1, pi = 0.1, log = FALSE) {
  .count_d(.pig_fam, x, list(mu = mu, sigma = sigma), log = log, pi = pi)
}

#' @rdname zipig
#' @export
pzipig <- function(q, mu = 1, sigma = 1, pi = 0.1, lower.tail = TRUE, log.p = FALSE) {
  .count_p(.pig_fam, q, list(mu = mu, sigma = sigma), lower.tail, log.p, pi = pi)
}

#' @rdname zipig
#' @export
qzipig <- function(p, mu = 1, sigma = 1, pi = 0.1, lower.tail = TRUE, log.p = FALSE,
                   max.value = 1e5) {
  .count_q(.pig_fam, p, list(mu = mu, sigma = sigma), lower.tail, log.p, max.value,
           pi = pi)
}

#' @rdname zipig
#' @export
rzipig <- function(n, mu = 1, sigma = 1, pi = 0.1) {
  .count_r_exact(.pig_fam, n, list(mu = mu, sigma = sigma), .rpig_exact, pi = pi)
}
