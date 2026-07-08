## Discrete-GPD peaks-over-threshold (POT) fitting.
## Unlike the GAM backend egpd(family = "degpd"), whose shape uses a (bounded)
## log link xi = exp(eta) and is therefore constrained to xi >= 0, this fits the
## discrete EGPD tail by direct MLE through fitegpd(), which uses an IDENTITY
## link on xi -- so a bounded/short tail (xi < 0, finite upper endpoint) can be
## estimated.

#' Discrete-GPD peaks-over-threshold fit (unconstrained shape)
#'
#' Fit a discrete (E)GPD tail to the exceedances of a non-negative integer
#' series over a threshold, by maximum likelihood through \code{\link{fitegpd}}.
#' The shape parameter \eqn{\xi} uses an identity link, so it is unconstrained
#' and may be negative (a bounded/short tail with finite upper endpoint), unlike
#' the GAM backend \code{egpd(family = "degpd")} whose log-link shape is forced
#' \eqn{\xi \ge 0}.
#'
#' The integer exceedances \eqn{Z = X - u} for \eqn{X \ge u} are modelled with
#' the \code{"degpd"} family. With \code{type = 1} and \code{kappa} fixed to 1
#' this is the plain discrete GPD (dGPD); leaving \code{kappa} free gives the
#' Naveau discrete EGPD tail.
#'
#' @param x numeric vector of counts (non-negative integers; NA dropped).
#' @param threshold integer threshold \eqn{u}; exceedances are \eqn{x - u} for
#'   \eqn{x \ge u}.
#' @param type integer 1-6, the Naveau G-transformation (default 1).
#' @param fix.kappa logical; if TRUE (default) fix \eqn{\kappa = 1} so the tail
#'   is the plain discrete GPD. Set FALSE to estimate the full discrete EGPD.
#' @param hessian logical; compute standard errors (default TRUE).
#' @param ... passed to \code{\link{fitegpd}} / \code{\link{optim}}.
#'
#' @return A \code{fitegpd} object (class \code{"degpd_pot"}) with extra fields:
#'   \code{threshold}, \code{n.exceed}, \code{prop.exceed}, \code{upper.endpoint}
#'   (\eqn{u - \sigma/\xi} when \eqn{\xi < 0}, else \code{Inf}), and \code{xi.ci}
#'   (95\% Wald interval for \eqn{\xi}).
#' @seealso \code{\link{degpd_pot_stability}}, \code{\link{fitegpd}}
#' @export
fit_degpd_pot <- function(x, threshold, type = 1, fix.kappa = TRUE,
                          hessian = TRUE, ...) {
  x <- x[!is.na(x)]
  if (any(x != round(x))) warning("non-integer values in 'x'; discrete GPD expects counts")
  exc <- x[x >= threshold] - threshold
  if (length(exc) < 10)
    warning(sprintf("only %d exceedances over threshold %s", length(exc), threshold))
  fix.arg <- if (fix.kappa) list(kappa = 1) else NULL
  fit <- distfit(exc, type = type, family = "degpd",
                 fix.arg = fix.arg, hessian = hessian, ...)
  xi  <- unname(fit$estimate["xi"]); sig <- unname(fit$estimate["sigma"])
  se  <- if (!is.null(fit$sd)) unname(fit$sd["xi"]) else NA_real_
  fit$threshold      <- threshold
  fit$n.exceed       <- length(exc)
  fit$prop.exceed    <- length(exc) / length(x)
  fit$upper.endpoint <- if (!is.na(xi) && xi < 0) threshold - sig / xi else Inf
  fit$xi.ci          <- if (!is.na(se)) c(xi - 1.96 * se, xi + 1.96 * se) else c(NA, NA)
  class(fit) <- c("degpd_pot", class(fit))
  fit
}

#' Threshold stability for the discrete-GPD POT shape
#'
#' Fit \code{\link{fit_degpd_pot}} across a range of thresholds and tabulate the
#' shape estimate, its 95\% interval, the number of exceedances and the implied
#' upper endpoint -- the standard threshold-stability diagnostic, here with an
#' unconstrained (possibly negative) \eqn{\xi}.
#'
#' @param x numeric vector of counts.
#' @param thresholds integer thresholds to scan. Default: the 75th-97th
#'   empirical percentiles of \code{x}.
#' @param ... passed to \code{\link{fit_degpd_pot}}.
#' @return A data.frame with columns \code{threshold, n_exceed, prop_exceed,
#'   sigma, xi, xi_lo, xi_hi, upper_endpoint}.
#' @seealso \code{\link{fit_degpd_pot}}
#' @export
degpd_pot_stability <- function(x, thresholds = NULL, ...) {
  x <- x[!is.na(x)]
  if (is.null(thresholds))
    thresholds <- unique(as.integer(round(quantile(x, seq(0.75, 0.97, by = 0.02)))))
  rows <- lapply(thresholds, function(u) {
    f <- tryCatch(fit_degpd_pot(x, u, hessian = TRUE, ...), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    data.frame(threshold = u, n_exceed = f$n.exceed, prop_exceed = f$prop.exceed,
               sigma = unname(f$estimate["sigma"]), xi = unname(f$estimate["xi"]),
               xi_lo = f$xi.ci[1], xi_hi = f$xi.ci[2],
               upper_endpoint = f$upper.endpoint)
  })
  do.call(rbind, rows)
}
