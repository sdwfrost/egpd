## Discrete-GPD peaks-over-threshold (POT) fitting.
## Unlike the GAM backend egpd(family = "degpd"), whose shape uses a (bounded)
## log link xi = exp(eta) and is therefore constrained to xi >= 0, this fits the
## discrete EGPD tail by direct MLE through fitegpd(), which uses an IDENTITY
## link on xi -- so a bounded/short tail (xi < 0, finite upper endpoint) can be
## estimated.

#' Discrete-GPD Peaks-Over-Threshold Fit and Threshold Stability
#'
#' Fits a discrete (extended) Generalized Pareto tail to integer exceedances
#' \eqn{Z = X - u} (\eqn{X \ge u}) by maximum likelihood through
#' \code{\link{distfit}}, parameterising the tail shape \eqn{\xi} on the
#' unconstrained identity scale so negative shapes (bounded tails) can be
#' estimated. \code{degpd_pot_stability} evaluates \code{fit_degpd_pot} across a
#' grid of integer thresholds as a threshold-stability diagnostic.
#'
#' The integer exceedances are modelled with the \code{"degpd"} family. With
#' \code{type = 1} and \code{kappa} fixed to 1 this is the plain discrete GPD
#' (dGPD) of Hitz, Davis and Samorodnitsky (2024); leaving \code{kappa} free
#' gives the Naveau discrete EGPD tail. Unlike the GAM backend
#' \code{egpd(family = "degpd")}, whose log-link shape is forced
#' \eqn{\xi \ge 0}, the shape here may be negative.
#'
#' For \eqn{\xi < 0} the exceedance \eqn{Z} takes the values
#' \eqn{0, 1, \ldots, \lceil -\sigma/\xi \rceil - 1}, so
#' \code{upper.endpoint} is the largest count with positive probability,
#' \eqn{u + \lceil -\sigma/\xi \rceil - 1} (not the continuous endpoint
#' \eqn{u - \sigma/\xi} of the underlying GPD).
#'
#' \code{degpd_pot_stability} skips thresholds without exceedances and
#' thresholds at which the fit fails (including a nonzero optimizer convergence
#' code), with a warning naming each one; it stops
#' with an error if no threshold can be fitted.
#' When the fitted endpoint leaves at most two support cells, free scale
#' and shape are not separately identifiable. The fit warns, retains one
#' representative and its likelihood, and withholds Wald inference.
#' The stability table keeps the threshold but reports \code{NA} for the
#' non-unique scale and shape.
#' If the observed-information calculation crosses a non-smooth change of
#' support at the largest exceedance, estimates are retained but no Wald
#' standard errors or intervals are reported.
#'
#' @param x numeric vector of counts (non-negative integers; \code{NA} values are
#'   dropped).
#' @param threshold single finite integer threshold \eqn{u}; exceedances are
#'   \eqn{x - u} for \eqn{x \ge u}.
#' @param type integer 1--6 selecting the G-transformation carrier (default
#'   \code{1}).
#' @param fix.kappa logical; if \code{TRUE} (default, valid only for
#'   \code{type = 1}), fix \eqn{\kappa = 1} so the exceedance model reduces to
#'   the plain discrete Generalized Pareto distribution (dGPD). Set \code{FALSE}
#'   to estimate the full discrete EGPD carrier (required when \code{type} is in
#'   \code{2:6}).
#' @param fix.arg optional named fixed parameters passed to \code{\link{distfit}}.
#'   With \code{fix.kappa = TRUE}, \code{kappa = 1} is included and must
#'   not conflict with a supplied value. Fixed scale and shape values
#'   are used when reporting the endpoint and stability table.
#' @param hessian logical; compute observed-information standard errors
#'   (default \code{TRUE}).
#' @param ... additional arguments passed to \code{\link{distfit}} /
#'   \code{\link[stats]{optim}} (for \code{fit_degpd_pot}) or to
#'   \code{fit_degpd_pot} (for \code{degpd_pot_stability}, e.g.
#'   \code{hessian = FALSE}).
#'
#' @return \code{fit_degpd_pot} returns a \code{\link{distfit}} object (with
#'   subclass \code{"degpd_pot"}) augmented with:
#'   \describe{
#'     \item{threshold}{the integer threshold \eqn{u}}
#'     \item{n.exceed}{number of exceedances \eqn{x \ge u}}
#'     \item{prop.exceed}{proportion of non-missing observations exceeding \eqn{u}}
#'     \item{upper.endpoint}{largest attainable count
#'       \eqn{u + \lceil -\sigma/\xi \rceil - 1} when \eqn{\xi < 0}, otherwise
#'       \code{Inf}}
#'     \item{xi.ci}{95\% Wald confidence interval for \eqn{\xi} (\code{NA}
#'       when ordinary Wald inference is unavailable)}
#'   }
#'
#'   \code{degpd_pot_stability} returns a \code{data.frame} with columns
#'   \code{threshold}, \code{n_exceed}, \code{prop_exceed}, \code{sigma},
#'   \code{xi}, \code{xi_lo}, \code{xi_hi}, and \code{upper_endpoint}, one row per
#'   fitted threshold.
#' @references Hitz, A. S., Davis, R. A. and Samorodnitsky, G. (2024). Discrete
#'   extremes. \emph{Journal of Data Science}, 22(4), 524--536.
#' @seealso \code{\link{distfit}}, \code{\link{ddiscegpd}}
#' @examples
#' set.seed(1)
#' x <- rdiscegpd(300, sigma = 4, xi = 0.2, kappa = 1, type = 1)
#' fit <- fit_degpd_pot(x, threshold = 3)
#' fit$estimate
#' degpd_pot_stability(x, thresholds = c(2L, 4L, 6L))
#' @export
fit_degpd_pot <- function(x, threshold, type = 1, fix.kappa = TRUE,
                          hessian = TRUE, ..., fix.arg = NULL) {
  if (!is.numeric(threshold) || is.logical(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold != floor(threshold)) {
    stop("'threshold' must be a single finite integer", call. = FALSE)
  }
  if (abs(threshold) <= .Machine$integer.max) {
    threshold <- as.integer(threshold)
  }
  if (!is.numeric(type) || is.logical(type) || length(type) != 1L ||
      is.na(type) || !(type %in% 1:6)) {
    stop("'type' must be an integer in 1:6", call. = FALSE)
  }
  type <- as.integer(type)
  if (!is.logical(fix.kappa) || length(fix.kappa) != 1L || is.na(fix.kappa)) {
    stop("'fix.kappa' must be TRUE or FALSE", call. = FALSE)
  }
  if (fix.kappa && type != 1L) {
    stop(
      "fix.kappa = TRUE reduces to the plain discrete GPD only for type = 1; ",
      "use fix.kappa = FALSE for type = 2..6",
      call. = FALSE
    )
  }

  x <- .pot_counts(x)
  exc <- as.numeric(x[x >= threshold]) - threshold
  if (length(exc) == 0L)
    stop("no exceedances over threshold ", threshold, call. = FALSE)
  if (length(exc) < 10)
    warning(sprintf("only %d exceedances over threshold %s", length(exc), threshold))
  if (fix.kappa) {
    if ("kappa" %in% names(fix.arg) && !isTRUE(fix.arg[["kappa"]] == 1))
      stop("'fix.kappa = TRUE' conflicts with 'fix.arg$kappa'; use fix.kappa = FALSE.",
           call. = FALSE)
    if (!"kappa" %in% names(fix.arg)) fix.arg <- c(as.list(fix.arg), list(kappa = 1))
  }
  fit <- distfit(exc, type = type, family = "degpd",
                 fix.arg = fix.arg, hessian = hessian, ...)
  xi <- .pot_parameter(fit, "xi"); sig <- .pot_parameter(fit, "sigma")
  se  <- if (!is.null(fit$sd)) unname(fit$sd["xi"]) else NA_real_
  fit$threshold      <- threshold
  fit$n.exceed       <- length(exc)
  fit$prop.exceed    <- length(exc) / length(x)
  ## largest count with positive probability: Z = X - u <= ceiling(-sigma/xi) - 1
  fit$upper.endpoint <- if (!is.na(xi) && xi < 0) threshold + .gpd_last_count(sig, xi) else Inf
  fit$xi.ci          <- if (!is.na(se)) c(xi - 1.96 * se, xi + 1.96 * se) else c(NA_real_, NA_real_)
  class(fit) <- c("degpd_pot", class(fit))
  fit
}

.pot_counts <- function(x) {
  if (!is.numeric(x) || !is.null(dim(x)))
    stop("'x' must be a numeric vector of non-negative integer counts.", call. = FALSE)
  x <- x[!is.na(x)]
  if (any(!is.finite(x) | x < 0 | x != floor(x)))
    stop("'x' must contain finite non-negative integer counts (or NA).", call. = FALSE)
  x
}

.pot_parameter <- function(fit, name) {
  if (name %in% names(fit$fix.arg)) return(unname(fit$fix.arg[[name]]))
  unname(fit$estimate[name])
}

#' @rdname fit_degpd_pot
#' @param thresholds integer vector of thresholds to scan in
#'   \code{degpd_pot_stability}. Defaults to the unique rounded empirical
#'   percentiles 75, 77, \ldots, 97 of \code{x}. Integer-valued thresholds above
#'   the 32-bit integer range are retained as doubles.
#' @export
degpd_pot_stability <- function(x, thresholds = NULL, ...) {
  x <- .pot_counts(x)
  if (is.null(thresholds))
    thresholds <- round(quantile(x, seq(0.75, 0.97, by = 0.02), names = FALSE))
  if (!is.numeric(thresholds) || is.logical(thresholds) || length(thresholds) == 0L ||
      any(!is.finite(thresholds)) || any(thresholds != floor(thresholds)))
    stop("'thresholds' must be a non-empty vector of finite integers", call. = FALSE)
  thresholds <- unique(if (all(abs(thresholds) <= .Machine$integer.max))
    as.integer(thresholds) else as.numeric(thresholds))
  failures <- character(0)
  rows <- lapply(thresholds, function(u) {
    if (!any(x >= u)) {
      failures <<- c(failures, sprintf("threshold %s: no exceedances", format(u)))
      return(NULL)
    }
    f <- tryCatch(withCallingHandlers(fit_degpd_pot(x, u, ...),
      distfit_nonconvergence = function(w) invokeRestart("muffleWarning")), error = identity)
    if (inherits(f, "error")) {
      failures <<- c(failures, sprintf("threshold %s: %s", format(u), conditionMessage(f)))
      return(NULL)
    }
    if (!isTRUE(f$convergence == 0L)) {
      failures <<- c(failures, sprintf("threshold %s: optimizer did not converge (code %s)",
                                       format(u), format(f$convergence)))
      return(NULL)
    }
    data.frame(threshold = u, n_exceed = f$n.exceed, prop_exceed = f$prop.exceed,
               sigma = if (isTRUE(f$support.ridge) && !"sigma" %in% names(f$fix.arg)) NA_real_ else .pot_parameter(f,"sigma"),
               xi = if (isTRUE(f$support.ridge) && !"xi" %in% names(f$fix.arg)) NA_real_ else .pot_parameter(f,"xi"),
               xi_lo = f$xi.ci[1], xi_hi = f$xi.ci[2],
               upper_endpoint = f$upper.endpoint)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    stop("no threshold could be fitted:\n  ", paste(failures, collapse = "\n  "),
         call. = FALSE)
  }
  if (length(failures)) {
    warning("skipped ", length(failures), " threshold(s):\n  ",
            paste(failures, collapse = "\n  "), call. = FALSE)
  }
  do.call(rbind, rows)
}
