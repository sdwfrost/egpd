## Blended generalized Pareto distribution utilities

#' Blended generalized Pareto distribution functions
#'
#' Distribution utilities for the blended generalized Pareto construction of
#' Majumder and Richards (2026): \code{\link{pbgpd}}, \code{\link{dbgpd}},
#' \code{\link{qbgpd}}, and \code{\link{rbgpd}}. A finite-support bulk model on
#' \eqn{[0,1]} is blended with a GP tail using the paper's multiplicative CDF
#' blend and Beta-CDF weighting function.
#'
#' The GP threshold and scale are derived from \code{pa} and \code{pb} so that
#' the blended CDF matches the bulk CDF at the start and end of the blending
#' interval. For \eqn{y \le a}, the distribution equals the bulk model; for
#' \eqn{y \ge b}, it equals the exact GP model. The Beta weighting-function
#' shapes must satisfy \eqn{c_1, c_2 > 3} (Majumder and Richards, p. 12), which
#' makes the second derivative of the log-density continuous; all four
#' functions reject smaller values.
#'
#' Because the paper notes that some bulk/weight combinations can produce a
#' negative density in the blending interval, these helpers numerically screen
#' the blended density over that interval and stop if the supplied specification
#' appears invalid.
#' Validation evaluates a fine spatial grid plus Beta-weight quantile knots,
#' followed by local density minimization across concentrated transitions,
#' finite-difference sign changes, and low-density subintervals.
#' This is not a monotonicity proof for arbitrary bulk APIs. The source
#' construction and negative-density caveat are in the original-author
#' preprint arXiv:2504.19994v2, Eqs. (4)--(6).
#'
#' @references
#' Majumder, R. and Richards, J. (2026). Semi-parametric bulk and tail
#' regression using spline-based neural networks. \emph{Extremes}, 29(2),
#' 245--277. \doi{10.1007/s10687-026-00533-y}
#'
#' @seealso \code{\link{pbgpd}}, \code{\link{dbgpd}}, \code{\link{qbgpd}},
#'   \code{\link{rbgpd}}
#'
#' @examples
#' ## Default Beta-bulk blended GP
#' pbgpd(c(0.2, 0.8, 1.5), xi = 0.2)
#' dbgpd(c(0.2, 0.8, 1.5), xi = 0.2)
#' qbgpd(c(0.5, 0.9, 0.99), xi = 0.2)
#'
#' ## Custom bulk distribution on [0, 1]
#' qbgpd(
#'   c(0.5, 0.9, 0.99),
#'   xi = -0.15,
#'   pbulk = pbeta,
#'   dbulk = dbeta,
#'   qbulk = qbeta,
#'   bulk.args = list(shape1 = 3, shape2 = 7),
#'   pa = 0.7,
#'   pb = 0.95
#' )
#' @name bgpd
NULL

.bgpd_bulk_call <- function(fun, x, bulk.args) {
  do.call(match.fun(fun), c(list(x), bulk.args))
}

.bgpd_exprel <- function(x) {
  v <- expm1(x) / x
  v[x == 0] <- 1
  v
}

## All calls are at or above the lower matching knot. Use the conditional
## GP there, avoiding an extrapolated location and two nearly equal powers.
.bgpd_tail_hazard <- function(x, spec) {
  xi <- spec$xi
  if (!is.null(spec$steep_L)) {
    ## A = log(1 + t*expm1(xi*L))/xi, t = (x-a)/(b-a).
    lt <- log(x - spec$a) - log(spec$b - spec$a)
    z <- xi * spec$steep_L
    if (is.finite(z)) return(.logspace_add(0, z + lt) / xi)
    return(ifelse(x == spec$a, 0, spec$steep_L + lt / xi))
  }
  t <- (x - spec$a) / spec$scale_a
  if (xi == 0) return(t)
  z <- xi * t
  out <- rep(Inf, length(t))
  live <- which(is.finite(x) & x < spec$endpoint)
  if (length(live)) {
    far <- live[z[live] < -0.5 & is.finite(spec$endpoint)]
    near <- setdiff(live, far)
    if (length(far))
      out[far] <- log((spec$endpoint - x[far]) / (spec$endpoint - spec$a)) / xi
    small <- near[abs(z[near]) < sqrt(.Machine$double.eps)]
    out[small] <- t[small] * (1 - z[small] / 2)
    rest <- setdiff(near, small)
    out[rest] <- log1p(z[rest]) / xi
    ov <- rest[z[rest] == Inf]
    if (length(ov)) out[ov] <- (log(xi) + log(x[ov] - spec$a) - log(spec$scale_a)) / xi
  }
  out
}

.bgpd_gp_cdf <- function(x, spec) {
  -expm1(log1p(-spec$pa) - .bgpd_tail_hazard(x, spec))
}

.bgpd_gp_pdf <- function(x, spec, log = FALSE) {
  A <- .bgpd_tail_hazard(x, spec)
  ld <- if (is.null(spec$steep_L)) {
    log1p(-spec$pa) - log(spec$scale_a) - (1 + spec$xi) * A
  } else {
    ## Density = survival / local scale. The scale at a can round to zero.
    ls <- log(spec$xi) + .logspace_add(log(x - spec$a),
                                        log(spec$b - spec$a) - spec$xi * spec$steep_L)
    log1p(-spec$pa) - A - ls
  }
  ld[x >= spec$endpoint] <- -Inf
  if (log) ld else exp(ld)
}

.bgpd_bulk_logdens <- function(x, spec) {
  if ("log" %in% names(formals(spec$dbulk)))
    return(.bgpd_bulk_call(spec$dbulk, x, c(spec$bulk.args, list(log = TRUE))))
  log(.bgpd_bulk_call(spec$dbulk, x, spec$bulk.args))
}

.bgpd_gp_quantile <- function(p, spec) {
  if (!is.null(spec$steep_L)) {
    h <- log1p((p - spec$pb) / (1 - p))
    return(spec$b + .qgpd_from_logsurv(-h, (spec$b - spec$a) * spec$xi, spec$xi))
  }
  h <- log1p((p - spec$pa) / (1 - p))
  spec$a + .qgpd_from_logsurv(-h, spec$scale_a, spec$xi)
}

.bgpd_weight <- function(y, a, b, c1, c2) {
  pbeta((y - a) / (b - a), shape1 = c1, shape2 = c2)
}

.bgpd_weight_prime <- function(y, a, b, c1, c2) {
  dbeta((y - a) / (b - a), shape1 = c1, shape2 = c2) / (b - a)
}

.bgpd_bridge_params <- function(xi, qbulk, bulk.args, pa, pb) {
  support <- .bgpd_bulk_call(qbulk, c(0, 1), bulk.args)
  if (length(support) != 2 || any(!is.finite(support))) {
    stop("The bulk quantile function must return finite values at probabilities 0 and 1.")
  }
  if (abs(support[1]) > 1e-8 || abs(support[2] - 1) > 1e-8) {
    stop("The bulk distribution must have support on [0, 1].")
  }

  ab <- .bgpd_bulk_call(qbulk, c(pa, pb), bulk.args)
  a <- ab[1]
  b <- ab[2]
  if (!is.finite(a) || !is.finite(b) || a < 0 || b > 1 || b <= a) {
    stop("The bulk quantile function must return 0 <= Q(pa) < Q(pb) <= 1.")
  }

  L <- log1p((pb - pa) / (1 - pb))
  ## For a steep positive bridge, expm1(xi*L) may overflow and even
  ## scale_a may underflow. Anchor the inverse at b instead: scale_b =
  ## (b-a)*xi/(1-exp(-xi*L)), whose denominator rounds to exactly one.
  if (xi * L > 700) return(list(a = a, b = b, steep_L = L, endpoint = Inf))
  scale_a <- ((b - a) / L) / .bgpd_exprel(xi * L)
  if (!is.finite(scale_a) || scale_a <= 0) {
    stop("The derived GP scale is not positive.")
  }
  endpoint <- if (xi < 0) a - scale_a / xi else Inf
  if (xi < 0 && endpoint <= b) {
    stop("The derived GP upper endpoint must be larger than the upper blending bound.")
  }

  list(a = a, b = b, scale_a = scale_a, endpoint = endpoint)
}

.bgpd_prepare <- function(
    xi,
    pbulk,
    dbulk,
    qbulk,
    bulk.args,
    pa,
    pb,
    c1,
    c2,
    validate = TRUE) {
  if (!is.numeric(xi) || length(xi) != 1 || !is.finite(xi)) {
    stop("`xi` must be a finite numeric scalar.")
  }
  if (is.null(bulk.args)) {
    bulk.args <- list()
  }
  if (!is.list(bulk.args)) {
    stop("`bulk.args` must be a list.")
  }
  if (!is.numeric(pa) || length(pa) != 1 || !is.finite(pa) || pa <= 0 || pa >= 1) {
    stop("`pa` must be a scalar in (0, 1).")
  }
  if (!is.numeric(pb) || length(pb) != 1 || !is.finite(pb) || pb <= 0 || pb >= 1) {
    stop("`pb` must be a scalar in (0, 1).")
  }
  if (pa >= pb) {
    stop("`pa` must be strictly smaller than `pb`.")
  }
  ## Majumder and Richards (p. 12) assume c1, c2 > 3 so that the second
  ## derivative of the log-density is continuous.
  if (!is.numeric(c1) || length(c1) != 1 || !is.finite(c1) || c1 <= 3) {
    stop("`c1` must be a finite scalar greater than 3.")
  }
  if (!is.numeric(c2) || length(c2) != 1 || !is.finite(c2) || c2 <= 3) {
    stop("`c2` must be a finite scalar greater than 3.")
  }

  params <- .bgpd_bridge_params(
    xi = xi,
    qbulk = qbulk,
    bulk.args = bulk.args,
    pa = pa,
    pb = pb
  )

  spec <- c(
    list(
      xi = xi,
      pbulk = match.fun(pbulk),
      dbulk = match.fun(dbulk),
      qbulk = match.fun(qbulk),
      bulk.args = bulk.args,
      pa = pa,
      pb = pb,
      c1 = c1,
      c2 = c2
    ),
    params
  )

  if (validate) {
    ## Resolve concentrated Beta transitions AND narrow multimodal bulk troughs:
    ## evaluate on a fine 2049-point spatial grid plus Beta-weight quantile knots,
    ## then run local minimization on concentrated weight intervals, local-minimum
    ## brackets where finite differences change sign, boundary cells, and low-density
    ## subintervals.
    beta_knots <- sort(unique(c(
      0, 1,
      qbeta(plogis(seq(-28, 28, length.out = 65L)), c1, c2)
    )))
    beta_knots <- beta_knots[is.finite(beta_knots) & beta_knots >= 0 & beta_knots <= 1]
    t_grid <- sort(unique(c(seq(0, 1, length.out = 2049L), beta_knots)))
    grid <- spec$a + (spec$b - spec$a) * t_grid
    dens <- .dbgpd_from_spec(grid, spec, log = FALSE)
    if (any(!is.finite(dens))) {
      stop("The blended GP density is not finite across the blending interval.")
    }
    if (min(dens) < -1e-8) {
      stop(
        "The blended GP density becomes negative in the blending interval. ",
        "Choose a different bulk distribution or blending hyperparameters."
      )
    }

    obj_fun <- function(t) {
      .dbgpd_from_spec(spec$a + (spec$b - spec$a) * t, spec, log = FALSE)
    }

    ## Coarse Beta-weight intervals (resolves concentrated transitions)
    for (i in seq_len(length(beta_knots) - 1L)) {
      if (beta_knots[i + 1L] - beta_knots[i] <= .Machine$double.eps) next
      opt <- optimize(obj_fun, interval = beta_knots[c(i, i + 1L)],
                      tol = .Machine$double.eps^0.75)
      if (!is.finite(opt$objective)) {
        stop("The blended GP density is not finite across the blending interval.")
      }
      if (opt$objective < -1e-8) {
        stop(
          "The blended GP density becomes negative in the blending interval. ",
          "Choose a different bulk distribution or blending hyperparameters."
        )
      }
    }

    ## Fine-grid candidate subintervals: sign changes in finite differences,
    ## endpoints, and low-density regions.
    n_g <- length(t_grid)
    dd <- diff(dens)
    low_thresh <- max(1e-4, 0.05 * max(dens))
    cand_left <- logical(n_g - 1L)
    cand_left[1L] <- TRUE
    cand_left[n_g - 1L] <- TRUE
    cand_left[pmin(dens[-n_g], dens[-1L]) <= low_thresh] <- TRUE
    if (n_g >= 3L) {
      trough_idx <- which(dd[-length(dd)] <= 0 & dd[-1L] >= 0)
      if (length(trough_idx) > 0L) {
        cand_left[trough_idx] <- TRUE
        cand_left[trough_idx + 1L] <- TRUE
      }
    }
    for (i in which(cand_left)) {
      if (t_grid[i + 1L] - t_grid[i] <= .Machine$double.eps) next
      opt <- optimize(obj_fun, interval = t_grid[c(i, i + 1L)],
                      tol = .Machine$double.eps^0.75)
      if (!is.finite(opt$objective)) {
        stop("The blended GP density is not finite across the blending interval.")
      }
      if (opt$objective < -1e-8) {
        stop(
          "The blended GP density becomes negative in the blending interval. ",
          "Choose a different bulk distribution or blending hyperparameters."
        )
      }
    }
  }

  spec
}

.pbgpd_from_spec <- function(q, spec) {
  if (!is.numeric(q)) {
    stop("`q` must be numeric.")
  }
  q <- as.numeric(q)
  out <- rep(NA_real_, length(q))

  idx_bulk <- !is.na(q) & q >= 0 & q <= spec$a
  idx_blend <- !is.na(q) & q > spec$a & q < spec$b
  idx_tail <- !is.na(q) & q >= spec$b

  if (any(idx_bulk)) {
    out[idx_bulk] <- .bgpd_bulk_call(spec$pbulk, q[idx_bulk], spec$bulk.args)
  }

  if (any(idx_blend)) {
    Fbulk <- .bgpd_bulk_call(spec$pbulk, q[idx_blend], spec$bulk.args)
    Fgp <- .bgpd_gp_cdf(q[idx_blend], spec)
    w <- .bgpd_weight(q[idx_blend], a = spec$a, b = spec$b, c1 = spec$c1, c2 = spec$c2)
    out[idx_blend] <- exp((1 - w) * log(Fbulk) + w * log(Fgp))
  }

  if (any(idx_tail)) {
    out[idx_tail] <- .bgpd_gp_cdf(q[idx_tail], spec)
  }

  if (spec$xi < 0) {
    out[q >= spec$endpoint] <- 1
  }
  out[!is.na(q) & q < 0] <- 0
  out
}

.dbgpd_from_spec <- function(x, spec, log = FALSE) {
  if (!is.numeric(x)) {
    stop("`x` must be numeric.")
  }
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))

  idx_bulk <- !is.na(x) & x >= 0 & x <= spec$a
  idx_blend <- !is.na(x) & x > spec$a & x < spec$b
  idx_tail <- !is.na(x) & x >= spec$b

  if (any(idx_bulk)) {
    out[idx_bulk] <- if (log) .bgpd_bulk_logdens(x[idx_bulk], spec) else
      .bgpd_bulk_call(spec$dbulk, x[idx_bulk], spec$bulk.args)
  }

  if (any(idx_blend)) {
    Fbulk <- .bgpd_bulk_call(spec$pbulk, x[idx_blend], spec$bulk.args)
    Fgp <- .bgpd_gp_cdf(x[idx_blend], spec)
    w <- .bgpd_weight(x[idx_blend], a = spec$a, b = spec$b, c1 = spec$c1, c2 = spec$c2)
    lH <- (1 - w) * log(Fbulk) + w * log(Fgp)
    ratio <- log(Fgp) - log(Fbulk)
    if (log) {
      ## H' = H * [(1-w)f/F + w*g/G + w' log(G/F)]. Only the last
      ## term is signed. Keep the two positive densities in log space.
      t <- (x[idx_blend] - spec$a) / (spec$b - spec$a)
      a <- pbeta(t, spec$c1, spec$c2, lower.tail = FALSE, log.p = TRUE) +
        .bgpd_bulk_logdens(x[idx_blend], spec) - log(Fbulk)
      b <- pbeta(t, spec$c1, spec$c2, log.p = TRUE) +
        .bgpd_gp_pdf(x[idx_blend], spec, log = TRUE) - log(Fgp)
      pos <- .logspace_add(a, b)
      adj <- dbeta(t, spec$c1, spec$c2, log = TRUE) - log(spec$b - spec$a) + log(abs(ratio))
      ld <- .logspace_add(pos, adj)
      sub <- which(ratio < 0)
      ld[sub] <- pos[sub] + .log1mexp(adj[sub] - pos[sub])
      neg <- sub[adj[sub] > pos[sub]]
      invalid <- any(lH[neg] + adj[neg] + .log1mexp(pos[neg]-adj[neg]) > log(1e-8))
      out[idx_blend] <- lH + ld
    } else {
      fbulk <- .bgpd_bulk_call(spec$dbulk, x[idx_blend], spec$bulk.args)
      fgp <- .bgpd_gp_pdf(x[idx_blend], spec)
      wprime <- .bgpd_weight_prime(x[idx_blend], spec$a, spec$b, spec$c1, spec$c2)
      dens <- exp(lH) * (wprime * ratio + w * (fgp / Fgp) + (1 - w) * (fbulk / Fbulk))
      invalid <- any(dens < -1e-8)
      out[idx_blend] <- pmax(dens, 0)
    }
    if (invalid) stop("The blended GP density becomes negative in the blending interval. ",
                      "Choose a different bulk distribution or blending hyperparameters.")
  }

  if (any(idx_tail)) {
    out[idx_tail] <- .bgpd_gp_pdf(x[idx_tail], spec, log = log)
  }

  zero <- if (log) -Inf else 0
  out[!is.na(x) & x < 0] <- zero
  if (spec$xi < 0) {
    out[!is.na(x) & x >= spec$endpoint] <- zero
  }
  out
}

.qbgpd_from_spec <- function(p, spec) {
  if (!is.numeric(p)) {
    stop("`p` must be numeric.")
  }
  p <- as.numeric(p)
  if (any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("Probabilities must lie in [0, 1].")
  }

  out <- rep(NA_real_, length(p))
  idx_bulk <- !is.na(p) & p <= spec$pa
  idx_blend <- !is.na(p) & p > spec$pa & p < spec$pb
  idx_tail <- !is.na(p) & p >= spec$pb & p < 1

  if (any(idx_bulk)) {
    out[idx_bulk] <- .bgpd_bulk_call(spec$qbulk, p[idx_bulk], spec$bulk.args)
  }

  if (any(idx_blend)) {
    out[idx_blend] <- vapply(p[idx_blend], function(target) {
      fa <- .pbgpd_from_spec(spec$a, spec) - target
      if (fa >= 0) return(spec$a)
      fb <- .pbgpd_from_spec(spec$b, spec) - target
      if (fb <= 0) return(spec$b)
      width <- spec$b - spec$a
      root <- uniroot(
        function(t) .pbgpd_from_spec(spec$a + width * t, spec) - target,
        interval = c(0, 1),
        f.lower = fa,
        f.upper = fb,
        tol = .Machine$double.eps
      )$root
      spec$a + width * root
    }, numeric(1))
  }

  if (any(idx_tail)) {
    out[idx_tail] <- .bgpd_gp_quantile(p[idx_tail], spec)
  }

  idx_one <- !is.na(p) & p == 1
  idx_zero <- !is.na(p) & p == 0
  if (spec$xi >= 0) {
    out[idx_one] <- Inf
  } else {
    out[idx_one] <- spec$endpoint
  }
  out[idx_zero] <- 0
  out
}

#' CDF of the blended generalized Pareto distribution
#'
#' @details Numerical validation uses spatial and Beta-weight quantile
#'   knots plus local density minimization between knots, so narrow blending
#'   transitions are not skipped. This is not a monotonicity proof for
#'   arbitrary user-supplied bulk functions. The source construction and
#'   negative-density caveat are in the original-author preprint
#'   arXiv:2504.19994v2, Eqs. (4)--(6).
#'
#' @param q quantiles.
#' @param xi GP tail shape parameter.
#' @param pbulk bulk CDF supported on \eqn{[0, 1]}.
#' @param dbulk bulk density supported on \eqn{[0, 1]}. This is used to
#'   validate the blended construction, even for \code{pbgpd()} and
#'   \code{qbgpd()}. An explicit \code{log} argument on this callback lets
#'   \code{dbgpd(log = TRUE)} retain bulk densities below the linear range.
#' @param qbulk bulk quantile function supported on \eqn{[0, 1]}.
#' @param bulk.args named list of arguments passed to \code{pbulk},
#'   \code{dbulk}, and \code{qbulk}. Defaults to a Beta bulk.
#' @param pa lower bulk quantile used to start blending.
#' @param pb upper bulk quantile used to end blending and recover the exact GP
#'   tail.
#' @param c1,c2 Beta weighting-function shape parameters, both greater than
#'   3: Majumder and Richards (p. 12) assume \eqn{c_1, c_2 > 3}, which makes
#'   the second derivative of the log-density continuous. Smaller values are
#'   rejected by all four functions.
#'
#' @return CDF values.
#' @seealso \code{\link{bgpd}} for the construction and its reference.
#' @export
pbgpd <- function(
    q,
    xi,
    pbulk = pbeta,
    dbulk = dbeta,
    qbulk = qbeta,
    bulk.args = list(shape1 = 2, shape2 = 5),
    pa = 0.8,
    pb = 0.99,
    c1 = 5,
    c2 = 5) {
  spec <- .bgpd_prepare(
    xi = xi,
    pbulk = pbulk,
    dbulk = dbulk,
    qbulk = qbulk,
    bulk.args = bulk.args,
    pa = pa,
    pb = pb,
    c1 = c1,
    c2 = c2
  )
  .pbgpd_from_spec(q, spec)
}

#' Density of the blended generalized Pareto distribution
#'
#' @inheritParams pbgpd
#' @param x values.
#' @param log logical; should log-densities be returned?
#'
#' @return Density values.
#' @seealso \code{\link{bgpd}} for the construction and its reference.
#' @export
dbgpd <- function(
    x,
    xi,
    pbulk = pbeta,
    dbulk = dbeta,
    qbulk = qbeta,
    bulk.args = list(shape1 = 2, shape2 = 5),
    pa = 0.8,
    pb = 0.99,
    c1 = 5,
    c2 = 5,
    log = FALSE) {
  spec <- .bgpd_prepare(
    xi = xi,
    pbulk = pbulk,
    dbulk = dbulk,
    qbulk = qbulk,
    bulk.args = bulk.args,
    pa = pa,
    pb = pb,
    c1 = c1,
    c2 = c2
  )
  .dbgpd_from_spec(x, spec, log = log)
}

#' Quantile function of the blended generalized Pareto distribution
#'
#' @inheritParams pbgpd
#' @param p probabilities.
#'
#' @return Quantile values.
#' @seealso \code{\link{bgpd}} for the construction and its reference.
#' @export
qbgpd <- function(
    p,
    xi,
    pbulk = pbeta,
    dbulk = dbeta,
    qbulk = qbeta,
    bulk.args = list(shape1 = 2, shape2 = 5),
    pa = 0.8,
    pb = 0.99,
    c1 = 5,
    c2 = 5) {
  spec <- .bgpd_prepare(
    xi = xi,
    pbulk = pbulk,
    dbulk = dbulk,
    qbulk = qbulk,
    bulk.args = bulk.args,
    pa = pa,
    pb = pb,
    c1 = c1,
    c2 = c2
  )
  .qbgpd_from_spec(p, spec)
}

#' Random generation from the blended generalized Pareto distribution
#'
#' @inheritParams pbgpd
#' @param n number of samples.
#' @param unifsamp optional uniforms to transform instead of drawing new ones.
#'
#' @return Random samples.
#' @seealso \code{\link{bgpd}} for the construction and its reference.
#' @export
rbgpd <- function(
    n,
    xi,
    pbulk = pbeta,
    dbulk = dbeta,
    qbulk = qbeta,
    bulk.args = list(shape1 = 2, shape2 = 5),
    pa = 0.8,
    pb = 0.99,
    c1 = 5,
    c2 = 5,
    unifsamp = NULL) {
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n < 0 || n != as.integer(n)) {
    stop("`n` must be a single non-negative integer.")
  }
  if (is.null(unifsamp)) {
    unifsamp <- runif(n)
  } else if (length(unifsamp) != n) {
    unifsamp <- rep(unifsamp, length.out = n)
  }

  qbgpd(
    p = unifsamp,
    xi = xi,
    pbulk = pbulk,
    dbulk = dbulk,
    qbulk = qbulk,
    bulk.args = bulk.args,
    pa = pa,
    pb = pb,
    c1 = c1,
    c2 = c2
  )
}
