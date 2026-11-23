## Distribution functions for the Bell (Castellares, Ferrari & Lemonte 2018,
## Appl. Math. Modelling 56:172-185) and zero-inflated Bell (ZIBell) count
## distributions.
##
## The pmf is computed from the closed form via src/bell.cpp (bell_pmf_cpp),
## using log(B_k/k!) log-space ratios for numerical stability. The natural
## theta parameterisation of the paper (Definition 1) is used: theta > 0, with
##   Pr(Y = y) = theta^y exp(1 - e^theta) B_y / y!.
## Mean = theta e^theta, variance = theta(1+theta) e^theta.

## The log pmf ladder log p_0..log p_n comes from bell_pmf_cpp; its theta-free table
## log(B_k/k!) is cached for the session. The tail is super-geometric, so survivors are
## summed directly until a ratio bound certifies the remainder (R/countfams.R).
.bell_fam <- list(
  pars = "theta",
  valid = function(p) is.finite(p$theta) & p$theta > 0,
  lpmf = function(n, p) attr(bell_pmf_cpp(as.integer(n), p$theta), "logpmf"),
  init = function(p) max(32, ceiling(4 * p$theta * exp(p$theta))))

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
#' @param n number of random values to return; if \code{length(n) != 1}, the length is
#'   taken to be the number required (\code{numeric(0)} gives none), as in
#'   \code{\link[stats]{rpois}}.
#' @param theta vector of positive shape parameters.
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \eqn{p} are
#'   given as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param max.value integer cap on the support searched by \code{qbell}
#'   (default \code{1e5}). A nonattaining cap returns \code{NA} with a warning. It does
#'   not limit \code{rbell}, which samples exactly.
#' @details Both tails are computed from sums of the PMF taken in the direction that
#'   keeps full relative accuracy (the upper tail from the survivor itself), so
#'   \code{lower.tail = FALSE} and \code{log.p = TRUE} stay accurate far into the tail.
#'   Quantiles invert the distribution, not a table of rounded lower-tail
#'   probabilities. Near one, distinct counts can have the same reported CDF,
#'   and \code{qbell(pbell(k))} need not return \code{k}.
#'   Invalid input follows \code{\link[stats]{dpois}}: an \code{NA} parameter gives
#'   \code{NA}; an invalid one gives \code{NaN} with a warning, except that an \code{NA}
#'   or \code{NaN} \code{x}, \code{q} or \code{p} is returned silently; \code{x} that is
#'   negative or infinite has mass zero, and a non-integer \code{x} has mass zero with a
#'   warning; probabilities outside \eqn{[0, 1]} give \code{NaN} with a warning. Quantile
#'   endpoints are 0 and \code{Inf}.
#'
#'   Random generation is exact and has no cap. The probability generating function
#'   \eqn{\exp(e^{\theta s} - e^\theta)} is that of a Neyman type A law, so
#'   \code{rbell} draws \eqn{K \sim} Poisson\eqn{(e^\theta)} and then a Poisson count
#'   with mean \eqn{\theta K}; equivalently, \eqn{Y} is a Poisson\eqn{(e^\theta - 1)}
#'   sum of zero-truncated Poisson\eqn{(\theta)} counts.
#' @return \code{dbell} gives the PMF, \code{pbell} the CDF, \code{qbell} the
#'   quantile function, and \code{rbell} random deviates; \code{dBELL}, \code{pBELL},
#'   \code{qBELL} and \code{rBELL} give the same in the mean parameterisation, and
#'   \code{BELL()} returns a \code{gamlss.family} object.
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
  .count_d(.bell_fam, x, list(theta = theta), log = log)
}

#' @rdname bell
#' @export
pbell <- function(q, theta = 1, lower.tail = TRUE, log.p = FALSE) {
  .count_p(.bell_fam, q, list(theta = theta), lower.tail, log.p)
}

#' @rdname bell
#' @export
qbell <- function(p, theta = 1, lower.tail = TRUE, log.p = FALSE, max.value = 1e5) {
  .count_q(.bell_fam, p, list(theta = theta), lower.tail, log.p, max.value)
}

#' @rdname bell
#' @export
rbell <- function(n, theta = 1) {
  .count_r_exact(.bell_fam, n, list(theta = theta), .rbell_exact)
}

## Exact Bell draws. The pgf exp(e^(theta s) - e^theta) is that of a Neyman type A law with
## K ~ Poisson(e^theta) clusters of Poisson(theta) counts: E[exp(theta K (s - 1))] =
## exp(e^theta (e^(theta (s - 1)) - 1)). (Dropping the empty clusters gives the compound
## Poisson form, Poisson(e^theta - 1) clusters of zero-truncated Poisson(theta) counts.)
## Given K the total is Poisson(theta K), so two Poisson draws give each count, at any theta.
.rbell_exact <- function(p) {
  m <- length(p$theta)
  stats::rpois(m, p$theta * stats::rpois(m, exp(p$theta)))
}


#' The zero-inflated Bell distribution
#'
#' Density (PMF), distribution function, quantile function and random generation
#' for the zero-inflated Bell (ZIBell) distribution: a \code{\link{bell}}
#' distribution mixed with a point mass at zero of probability \code{pi}.
#'
#' @inheritParams bell
#' @param pi vector of zero-inflation probabilities in \eqn{[0, 1)}.
#' @details Tails, quantiles and invalid input are handled as for \code{\link{bell}};
#'   \code{pi} outside \eqn{[0, 1)} is invalid. Random generation draws the structural
#'   zero with probability \code{pi} and otherwise an exact Bell draw (see
#'   \code{\link{rbell}}).
#' @return \code{dzibell} gives the PMF, \code{pzibell} the CDF, \code{qzibell}
#'   the quantile function, and \code{rzibell} random deviates; \code{dZIBELL},
#'   \code{pZIBELL}, \code{qZIBELL} and \code{rZIBELL} give the same in the mean
#'   parameterisation, and \code{ZIBELL()} returns a \code{gamlss.family} object.
#' @seealso \code{\link{bell}} for the non-inflated version.
#' @name zibell
NULL

#' @rdname zibell
#' @export
dzibell <- function(x, theta = 1, pi = 0.1, log = FALSE) {
  .count_d(.bell_fam, x, list(theta = theta), log = log, pi = pi)
}

#' @rdname zibell
#' @export
pzibell <- function(q, theta = 1, pi = 0.1, lower.tail = TRUE, log.p = FALSE) {
  .count_p(.bell_fam, q, list(theta = theta), lower.tail, log.p, pi = pi)
}

#' @rdname zibell
#' @export
qzibell <- function(p, theta = 1, pi = 0.1, lower.tail = TRUE, log.p = FALSE,
                    max.value = 1e5) {
  .count_q(.bell_fam, p, list(theta = theta), lower.tail, log.p, max.value, pi = pi)
}

#' @rdname zibell
#' @export
rzibell <- function(n, theta = 1, pi = 0.1) {
  .count_r_exact(.bell_fam, n, list(theta = theta), .rbell_exact, pi = pi)
}
