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
  if (any(lengths(a) == 0L)) return(lapply(a, function(x) x[FALSE]))
  n <- max(lengths(a))
  lapply(a, rep_len, n)
}

## The log pmf ladder comes from the exact recursion (gpig_pmf_cpp). Its cost is
## quadratic in the general case, so the survivor sum is capped and formed only
## when needed. The a = 1/2 subfamily has the exact linear-time PIG recurrence.
.gpig_fam <- list(
  pars = c("a", "b", "c"),
  valid = function(p) is.finite(p$a) & p$a > 0 & p$a <= 1 & is.finite(p$b) & p$b > 0 &
    is.finite(p$c) & p$c >= 0 & p$c <= 1,
  lpmf = function(n, p) {
    if (p$a == 0.5 && p$c > 0 && p$c < 1) {
      ## Match the PGFs: 1/sigma = b*s and 2*sigma*mu = c/(1-c).
      s <- sqrt(1 - p$c)
      mu <- p$b * p$c / (2 * s)
      sigma <- (1 / s) / p$b
      if (is.finite(mu) && mu >= .Machine$double.xmin &&
          is.finite(sigma) && sigma >= .Machine$double.xmin)
        return(pig_logpmf_ladder_cpp(as.integer(n), mu, sigma))
    }
    attr(gpig_pmf_cpp(as.integer(n), p$a, p$b, p$c), "logpmf")
  },
  rlim = function(p) if (p$a >= 1) 0 else p$c,
  maxN = function(p, need) max(2 * need, need + 20000),
  slow_ladder = TRUE,
  point0 = function(p) p$c == 0,
  init = function(p) 16)

#' The generalised Poisson-inverse Gaussian distribution
#'
#' Density (PMF), distribution function, quantile function and random
#' generation for the generalised Poisson-inverse Gaussian (GPIG) distribution
#' of Zhu and Joe (2009), a three-parameter count family with
#' probability generating function
#' \eqn{G(s) = \exp\{b[(1-c)^a - (1-cs)^a]\}}. It nests the Poisson
#' (\eqn{a=1}, mean \eqn{bc}), the Poisson-inverse Gaussian (\eqn{a=1/2}) and
#' the discrete-stable (\eqn{c=1}) distributions. The general PMF is evaluated
#' by the paper's forward recursion.
#' A power-law PMF, proportional asymptotically to \eqn{y^{-1-a}}, occurs only
#' at \eqn{c=1, 0<a<1}; \eqn{c<1} exponentially tempers this tail.
#' The \eqn{a=1} case is Poisson, not power-law.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return; if \code{length(n) != 1}, the length is
#'   taken to be the number required (\code{numeric(0)} gives none), as in
#'   \code{\link[stats]{rpois}}.
#' @param a vector of tail exponents in \eqn{(0, 1]}.
#' @param b vector of positive level parameters.
#' @param c vector of down-weighting parameters in \eqn{[0, 1]}.
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \eqn{p} are
#'   given as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param max.value integer cap on the support searched by \code{qgpig}, and by
#'   \code{rgpig}, which draws by inverting the distribution function (default
#'   \code{1e5}); relevant for very heavy tails. A quantile or draw above the cap is
#'   \code{NA} with a warning, not a fabricated value. Raising the cap costs time: the
#'   general recursion is quadratic in the largest count it reaches. At
#'   \eqn{a=1/2}, the linear-time PIG recurrence is used when the finite-mean
#'   PIG parameters can be represented.
#' @details Both tails are computed from sums of the PMF taken in the direction that
#'   keeps full relative accuracy (the upper tail from the survivor itself, summed until
#'   the geometric decay \eqn{c^k} bounds the remainder), so \code{lower.tail = FALSE}
#'   and \code{log.p = TRUE} stay accurate far into the tail. When \code{c} is so close
#'   to 1 that the survivor cannot be resolved within the summation limit, it is taken
#'   from \eqn{1 - F} while that is accurate and is \code{NaN} with a warning beyond.
#'   Quantiles invert the distribution, not a table of rounded lower-tail probabilities:
#'   near one, a reported lower CDF can round above the actual jump. The quantile
#'   endpoints are 0 and \code{Inf} (both zero for the point mass at zero when
#'   \code{c = 0}). Invalid input follows \code{\link[stats]{dpois}}: an \code{NA}
#'   parameter gives \code{NA}; an invalid one gives \code{NaN} with a warning, except
#'   that an \code{NA} or \code{NaN} \code{x}, \code{q} or \code{p} is returned silently;
#'   \code{x} that is negative or infinite has mass zero, and a non-integer \code{x} has
#'   mass zero with a warning; probabilities outside \eqn{[0, 1]} give \code{NaN} with
#'   a warning.
#' @return \code{dgpig} gives the PMF, \code{pgpig} the CDF, \code{qgpig} the
#'   quantile function, and \code{rgpig} random deviates; \code{dGPIG},
#'   \code{pGPIG}, \code{qGPIG} and \code{rGPIG} give the same in the mean
#'   parameterisation, and \code{GPIG()} returns a \code{gamlss.family} object.
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
  .count_d(.gpig_fam, x, list(a = a, b = b, c = c), log = log)
}

#' @rdname gpig
#' @export
pgpig <- function(q, a = 0.5, b = 1, c = 0.5, lower.tail = TRUE, log.p = FALSE) {
  .count_p(.gpig_fam, q, list(a = a, b = b, c = c), lower.tail, log.p)
}

#' @rdname gpig
#' @export
qgpig <- function(p, a = 0.5, b = 1, c = 0.5, lower.tail = TRUE, log.p = FALSE,
                  max.value = 1e5) {
  .count_q(.gpig_fam, p, list(a = a, b = b, c = c), lower.tail, log.p, max.value)
}

#' @rdname gpig
#' @export
rgpig <- function(n, a = 0.5, b = 1, c = 0.5, max.value = 1e5) {
  .count_r(.gpig_fam, n, list(a = a, b = b, c = c), max.value = max.value)
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
#' @details Tails, quantiles and invalid input are handled as for \code{\link{gpig}};
#'   \code{pi} outside \eqn{[0, 1)} is invalid.
#' @return \code{dzigpig} gives the PMF, \code{pzigpig} the CDF, \code{qzigpig}
#'   the quantile function, and \code{rzigpig} random deviates; \code{dZIGPIG},
#'   \code{pZIGPIG}, \code{qZIGPIG} and \code{rZIGPIG} give the same in the mean
#'   parameterisation, and \code{ZIGPIG()} returns a \code{gamlss.family} object.
#' @seealso \code{\link{gpig}} for the non-inflated version.
#' @name zigpig
NULL

#' @rdname zigpig
#' @export
dzigpig <- function(x, a = 0.5, b = 1, c = 0.5, pi = 0.1, log = FALSE) {
  .count_d(.gpig_fam, x, list(a = a, b = b, c = c), log = log, pi = pi)
}

#' @rdname zigpig
#' @export
pzigpig <- function(q, a = 0.5, b = 1, c = 0.5, pi = 0.1,
                    lower.tail = TRUE, log.p = FALSE) {
  .count_p(.gpig_fam, q, list(a = a, b = b, c = c), lower.tail, log.p, pi = pi)
}

#' @rdname zigpig
#' @export
qzigpig <- function(p, a = 0.5, b = 1, c = 0.5, pi = 0.1,
                    lower.tail = TRUE, log.p = FALSE, max.value = 1e5) {
  .count_q(.gpig_fam, p, list(a = a, b = b, c = c), lower.tail, log.p, max.value,
           pi = pi)
}

#' @rdname zigpig
#' @export
rzigpig <- function(n, a = 0.5, b = 1, c = 0.5, pi = 0.1, max.value = 1e5) {
  .count_r(.gpig_fam, n, list(a = a, b = b, c = c), max.value = max.value, pi = pi)
}
