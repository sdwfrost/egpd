## Conway-Maxwell-Poisson (CMP), mean-parameterised. See src/cmp.cpp for why Z is computed
## exactly (Gaunt's asymptotic is invalid at the nu these data force) and why that is still
## affordable (nu is intercept-only, so the exact sums are gridded once per likelihood call;
## the grid is validated against exact solves, which replace it where it is inaccurate).

#' Mean-parameterised Conway-Maxwell-Poisson distribution
#'
#' Density, distribution function, quantile function and random generation for the
#' Conway-Maxwell-Poisson (CMP) distribution in the mean parameterisation:
#' \code{P(Y = y) = lambda^y / ((y!)^nu Z(lambda, nu))}, so the mass is proportional to
#' \eqn{\lambda^y/(y!)^\nu}, with \code{lambda} chosen so that \code{E[Y] = mu}.
#' \code{nu < 1} is overdispersed, \code{nu = 1} is Poisson, \code{nu > 1} underdispersed.
#' The tail decays like \code{exp(-nu y log y)}, lighter than geometric for any
#' \code{nu > 0}, so the extreme-value index is \code{0}: CMP models dispersion, not a heavy
#' tail.
#'
#' The normalising constant and the mean inversion are computed exactly. Failure to invert
#' the mean or to resolve the normalising sum within the native computational limit is an
#' error.
#'
#' @details CDFs at minus/plus infinity are zero/one, and \code{q < 0} gives zero; a
#'   non-integer \code{q} is truncated. The upper tail (\code{lower.tail = FALSE}) and
#'   \code{log.p} are computed from the survivor itself, so they keep full relative
#'   accuracy far into the tail. Quantile endpoints are zero and \code{Inf}; the quantile
#'   is the generalised inverse of the reported distribution function (of the survivor for
#'   \code{lower.tail = FALSE}).
#'
#'   Invalid input follows \code{\link[stats]{dpois}}: \code{NA} parameters give
#'   \code{NA}; invalid parameters (\code{mu <= 0}, \code{nu <= 0}, \code{NaN},
#'   infinite) give \code{NaN} with a warning, except that an \code{NA} or \code{NaN}
#'   \code{x}, \code{q} or \code{p} is returned silently; \code{x} that is negative or
#'   infinite has mass zero, and a non-integer \code{x} has mass zero with a warning;
#'   probabilities outside \eqn{[0, 1]} give \code{NaN} with a warning.
#'
#' @param x,q vector of counts.
#' @param p vector of probabilities.
#' @param n number of draws; if \code{length(n) != 1}, the length is taken to be
#'   the number required (\code{numeric(0)} gives none), as in \code{\link[stats]{rpois}}.
#' @param mu mean, \code{mu > 0}.
#' @param nu dispersion parameter, \code{nu > 0}.
#' @param log,log.p logical; return the log density/probability, or supply log
#'   probabilities to \code{qcmp}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le q]},
#'   otherwise \eqn{P[X > q]}.
#' @param max.value largest count searched when inverting the cdf. A nonattaining cap
#'   returns \code{NA} with a warning, not a capped quantile.
#' @return A numeric vector of densities, probabilities, quantiles or draws.
#' @name cmp
#' @rdname cmp-count
NULL

#' @rdname cmp-count
#' @export
dcmp <- function(x, mu = 1, nu = 1, log = FALSE) {
  .count_d(.cmp_fam, x, list(mu = mu, nu = nu), log = log)
}
#' @rdname cmp-count
#' @export
pcmp <- function(q, mu = 1, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  .count_p(.cmp_fam, q, list(mu = mu, nu = nu), lower.tail, log.p)
}
#' @rdname cmp-count
#' @export
qcmp <- function(p, mu = 1, nu = 1, lower.tail = TRUE, log.p = FALSE, max.value = 1e6) {
  .count_q(.cmp_fam, p, list(mu = mu, nu = nu), lower.tail, log.p, max.value)
}
#' @rdname cmp-count
#' @export
rcmp <- function(n, mu = 1, nu = 1) {
  .count_r(.cmp_fam, n, list(mu = mu, nu = nu), max.value = 1e6)
}

## The tail is sub-geometric, so the survivor is summed until a ratio bound certifies the
## remainder.
.cmp_fam <- list(
  pars = c("mu", "nu"),
  valid = function(p) is.finite(p$mu) & p$mu > 0 & is.finite(p$nu) & p$nu > 0,
  lpmf = function(n, p) cmp_logpmf_cpp(0:n, p$mu, p$nu),
  lpmf_at = function(x, p) cmp_logpmf_cpp(as.numeric(x), p$mu, p$nu),
  init = function(p) max(32, ceiling(4 * p$mu)))

.iG_cmp <- function(v, mu, nu) qcmp(v, mu = mu, nu = nu)

.cmp1.d0 <- function(pars, likdata) {
  if (likdata$censored) stop("Censored likelihoods not available for cmp.")
  cmp1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
         likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
         likdata$dupid, likdata$duplicate, likdata$offsets)
}
.cmp1.d12 <- function(pars, likdata) {
  cmp1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
          likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
          likdata$dupid, likdata$duplicate, likdata$offsets)
}
.cmp1fns <- list(d0 = .cmp1.d0, d120 = .cmp1.d12, d340 = NULL, m = 1, iG = .iG_cmp)
