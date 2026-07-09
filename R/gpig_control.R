## Evaluation-method options for the GPIG / ZIGPIG likelihoods (src/gpig.cpp).
##
## Zhu & Joe's forward recursion reaches p_y by convolution, costing O(y^2) per
## observation, and carries the ladder in raw probability space seeded by
## p_0 = exp(b[(1-c)^a - 1]). For large counts that is both slow (a likelihood
## evaluation costs sum_i y_i^2) and numerically fatal (p_0 underflows to zero once
## the exponent drops below -745, i.e. a mean of order 1e3, after which the pmf is
## identically zero). gpig_control() selects among faster / safer routes.

GPIG_METHODS <- c(legacy = 0L, recursion = 1L, trunc = 2L, saddlepoint = 3L, hybrid = 4L)

#' Evaluation method for the GPIG and ZIGPIG likelihoods
#'
#' Selects how the GPIG probability mass function and its score are evaluated when
#' fitting \code{family = "gpig"}, \code{"gpignat"}, \code{"zigpig"} or
#' \code{"zigpignat"} with \code{\link{egpd}}. The setting is global and persists
#' until changed; \code{egpd()} sets and restores it around a fit when passed
#' \code{gpig.args}.
#'
#' @param method one of
#'   \describe{
#'     \item{\code{"hybrid"}}{(default) exact recursion for \code{y <= yswitch};
#'       above it the saddlepoint, but only where its own asymptotic series has
#'       converged (see \code{sptol}), otherwise the exact recursion again.}
#'     \item{\code{"recursion"}}{the exact \eqn{O(y^2)} recursion, carried in log
#'       space with dynamic rescaling. Immune to the \eqn{p_0} underflow.}
#'     \item{\code{"trunc"}}{as \code{"recursion"}, but the severity convolution is
#'       accumulated from its largest terms downwards and stopped once a term falls
#'       below \code{eps} times the partial sum. The criterion is deliberately
#'       relative: \eqn{p_y} can be of order \eqn{e^{-300}}, so an absolute bound on
#'       the discarded mass says nothing about the accuracy of \eqn{\log p_y}.}
#'     \item{\code{"saddlepoint"}}{Daniels' saddlepoint applied to the closed-form
#'       cumulant generating function \eqn{K(t) = b[(1-c)^a - (1-ce^t)^a]}. Costs
#'       \eqn{O(1)} per observation. Inaccurate for small means with \eqn{c \to 1}.}
#'     \item{\code{"legacy"}}{the original raw, unscaled recursion. Reproduces
#'       pre-existing results bit-for-bit, including the underflow.}
#'   }
#' @param yswitch count above which \code{"hybrid"} switches to the saddlepoint.
#'   The exact branch then costs at most \code{yswitch^2}. Default 200.
#' @param order saddlepoint order, 1 or 2. Order 2 adds the Daniels correction
#'   \eqn{1 + \kappa_4/8 - 5\kappa_3^2/24}, improving the relative error from
#'   \eqn{O(1/y)} to \eqn{O(1/y^2)}. It falls back to order 1 if the correction
#'   factor is non-positive. Default 2.
#' @param eps relative severity-truncation tolerance used by \code{method = "trunc"}.
#'   Default 1e-12.
#' @param sptol largest \eqn{|corr|} at which \code{"hybrid"} will accept a
#'   saddlepoint evaluation, where \eqn{corr = \kappa_4/8 - 5\kappa_3^2/24} is the
#'   leading neglected term of the expansion. Above \code{sptol} the exact recursion
#'   is used instead. Default 1e-3.
#' @param ymax largest count at which \code{"hybrid"} will fall back to the exact
#'   \eqn{O(y^2)} recursion when \code{sptol} is exceeded. Above it the saddlepoint
#'   is kept regardless, because the recursion could not finish: at \eqn{y = 10^5} it
#'   is \eqn{\sim 10^{10}} operations for a single observation. An optimiser routinely
#'   transits regions of small \eqn{b} and large \eqn{y} where \eqn{|corr|} is large,
#'   so this bound is load-bearing rather than a formality; such regions are transient,
#'   since at a converged fit \eqn{b} grows with the fitted mean and drives
#'   \eqn{|corr|} back below \code{sptol}. Default 1000.
#'
#' @param nquad number of quadrature nodes used to compute the saddlepoint's
#'   normalising constant \eqn{S = \sum_k \hat p(k)}. Default 256.
#' @param normalize divide the saddlepoint pmf by \eqn{S = \sum_k \hat p(k)}, restoring a
#'   proper pmf. Opt-in, default \code{FALSE}. Where the expansion has converged this
#'   sharpens the log-pmf by 50-100x, but at the \eqn{c \to 1} boundary it buys nothing:
#'   the bias there is \eqn{y}-dependent, whereas \eqn{S} is a constant (at the fitted
#'   MEASLES parameters \eqn{S = 1.034}, i.e. 0.033 log units, against an error of 0.15).
#'   The quadrature is also expensive, since the integrand is a peak of width
#'   \eqn{1/\sqrt{K''(0)}} that narrows as \eqn{b} grows. It refines until convergence and
#'   declines to normalise rather than return a wrong constant.
#'
#' @return Invisibly, the previous settings, as a list suitable for passing back to
#'   \code{gpig_control()}.
#'
#' @details The accuracy of the saddlepoint is governed by the large parameter
#'   \eqn{b} -- the compound-Poisson rate is \eqn{b[1-(1-c)^a]} -- and not by
#'   \eqn{y}. At \eqn{b = 300} the second-order log-pmf error is around \eqn{10^{-6}}
#'   at every \eqn{y}, whereas at \eqn{b = 5} it grows with \eqn{y}. This is why
#'   \code{"hybrid"} guards the switch with \code{sptol} rather than trusting
#'   \code{yswitch} alone. In a GAM \eqn{b} is proportional to the fitted mean, so a
#'   large count is accompanied by a large \eqn{b} and the saddlepoint is at its most
#'   accurate exactly where the recursion is most expensive.
#'
#'   The saddlepoint mass function does not sum to exactly one; that residual is the
#'   usual price of an unnormalised saddlepoint. Use \code{method = "recursion"} when
#'   an exact likelihood is required.
#'
#' @examples
#' \dontrun{
#' old <- gpig_control("saddlepoint", order = 2)
#' fit <- egpd(list(lmu = y ~ s(t), logita = ~1, logitc = ~1), data = df, family = "gpig")
#' do.call(gpig_control, old)
#' }
#' @seealso \code{\link{egpd}}
#' @export
gpig_control <- function(method = c("hybrid", "recursion", "trunc", "saddlepoint", "legacy"),
                         yswitch = 200L, order = 2L, eps = 1e-12, sptol = 1e-3,
                         ymax = 1000L, nquad = 256L, normalize = FALSE) {
  method <- match.arg(method)
  if (!is.numeric(yswitch) || length(yswitch) != 1L || yswitch < 0)
    stop("'yswitch' must be a single non-negative integer", call. = FALSE)
  if (!(order %in% c(1L, 2L))) stop("'order' must be 1 or 2", call. = FALSE)
  if (!is.numeric(eps) || length(eps) != 1L || eps <= 0 || eps >= 1)
    stop("'eps' must be a single value in (0, 1)", call. = FALSE)
  if (!is.numeric(sptol) || length(sptol) != 1L || sptol <= 0)
    stop("'sptol' must be a single positive value", call. = FALSE)
  if (!is.numeric(ymax) || length(ymax) != 1L || ymax < 0)
    stop("'ymax' must be a single non-negative integer", call. = FALSE)
  if (!is.numeric(nquad) || length(nquad) != 1L || nquad < 16)
    stop("'nquad' must be a single integer >= 16", call. = FALSE)
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize))
    stop("'normalize' must be TRUE or FALSE", call. = FALSE)
  old <- gpig_get_opts()
  old$method <- names(GPIG_METHODS)[match(old$method, GPIG_METHODS)]
  gpig_set_opts(as.integer(GPIG_METHODS[[method]]), as.integer(yswitch),
                as.integer(order), as.numeric(eps), as.numeric(sptol),
                as.integer(ymax), as.integer(nquad), as.integer(normalize))
  old$normalize <- as.logical(old$normalize)
  invisible(old)
}

## Apply gpig.args around a fit; returns the previous settings for on.exit restore.
## A no-op (returning NULL) when the family is not a GPIG variant or no args given.
.gpig_apply_control <- function(family, gpig.args) {
  if (!length(gpig.args)) return(NULL)
  if (!family %in% c("gpig", "gpignat", "zigpig", "zigpignat")) {
    warning("'gpig.args' ignored for family '", family, "'", call. = FALSE)
    return(NULL)
  }
  do.call(gpig_control, gpig.args)
}
