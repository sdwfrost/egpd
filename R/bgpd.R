## Blended generalized Pareto distribution utilities

.bgpd_bulk_call <- function(fun, x, bulk.args) {
  do.call(match.fun(fun), c(list(x), bulk.args))
}

.bgpd_gp_cdf <- function(q, u, sigma, xi) {
  out <- rep(NA_real_, length(q))
  idx <- !is.na(q) & q >= u
  if (!any(idx)) {
    out[!is.na(q) & q < u] <- 0
    return(out)
  }

  if (abs(xi) < 1e-8) {
    out[idx] <- 1 - exp(-(q[idx] - u) / sigma)
  } else {
    tmp <- 1 + xi * (q[idx] - u) / sigma
    tmp <- pmax(tmp, 0)
    out[idx] <- 1 - tmp^(-1 / xi)
    if (xi < 0) {
      endpoint <- u - sigma / xi
      out[q >= endpoint] <- 1
    }
  }

  out[!is.na(q) & q < u] <- 0
  out
}

.bgpd_gp_pdf <- function(x, u, sigma, xi) {
  out <- rep(NA_real_, length(x))
  idx <- !is.na(x) & x >= u
  if (!any(idx)) {
    out[!is.na(x) & x < u] <- 0
    return(out)
  }

  if (abs(xi) < 1e-8) {
    out[idx] <- exp(-(x[idx] - u) / sigma) / sigma
  } else {
    tmp <- 1 + xi * (x[idx] - u) / sigma
    valid <- tmp > 0
    out_idx <- numeric(sum(idx))
    out_idx[valid] <- (1 / sigma) * tmp[valid]^(-1 / xi - 1)
    out[idx] <- out_idx
  }

  out[!is.na(x) & x < u] <- 0
  out
}

.bgpd_gp_quantile <- function(p, u, sigma, xi) {
  if (abs(xi) < 1e-8) {
    return(u - sigma * log1p(-p))
  }
  u + sigma * ((1 - p)^(-xi) - 1) / xi
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

  if (abs(xi) < 1e-8) {
    loga <- log1p(-pa)
    logb <- log1p(-pb)
    sigma <- (a - b) / (logb - loga)
    u <- a + sigma * loga
  } else {
    A <- (1 - pa)^(-xi)
    B <- (1 - pb)^(-xi)
    sigma <- xi * (a - b) / (A - B)
    u <- a - sigma * (A - 1) / xi
  }

  endpoint <- if (xi < 0) u - sigma / xi else Inf
  if (!is.finite(sigma) || sigma <= 0) {
    stop("The derived GP scale is not positive.")
  }
  if (xi < 0 && endpoint <= b) {
    stop("The derived GP upper endpoint must be larger than the upper blending bound.")
  }

  list(a = a, b = b, u = u, sigma = sigma, endpoint = endpoint)
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
  if (!is.numeric(pa) || length(pa) != 1 || pa <= 0 || pa >= 1) {
    stop("`pa` must be a scalar in (0, 1).")
  }
  if (!is.numeric(pb) || length(pb) != 1 || pb <= 0 || pb >= 1) {
    stop("`pb` must be a scalar in (0, 1).")
  }
  if (pa >= pb) {
    stop("`pa` must be strictly smaller than `pb`.")
  }
  if (!is.numeric(c1) || length(c1) != 1 || c1 <= 0) {
    stop("`c1` must be a positive scalar.")
  }
  if (!is.numeric(c2) || length(c2) != 1 || c2 <= 0) {
    stop("`c2` must be a positive scalar.")
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
    grid <- seq(spec$a, spec$b, length.out = 201L)
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
    Fgp <- .bgpd_gp_cdf(q[idx_blend], u = spec$u, sigma = spec$sigma, xi = spec$xi)
    w <- .bgpd_weight(q[idx_blend], a = spec$a, b = spec$b, c1 = spec$c1, c2 = spec$c2)
    out[idx_blend] <- exp((1 - w) * log(Fbulk) + w * log(Fgp))
  }

  if (any(idx_tail)) {
    out[idx_tail] <- .bgpd_gp_cdf(q[idx_tail], u = spec$u, sigma = spec$sigma, xi = spec$xi)
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
    out[idx_bulk] <- .bgpd_bulk_call(spec$dbulk, x[idx_bulk], spec$bulk.args)
  }

  if (any(idx_blend)) {
    Fbulk <- .bgpd_bulk_call(spec$pbulk, x[idx_blend], spec$bulk.args)
    fbulk <- .bgpd_bulk_call(spec$dbulk, x[idx_blend], spec$bulk.args)
    Fgp <- .bgpd_gp_cdf(x[idx_blend], u = spec$u, sigma = spec$sigma, xi = spec$xi)
    fgp <- .bgpd_gp_pdf(x[idx_blend], u = spec$u, sigma = spec$sigma, xi = spec$xi)
    w <- .bgpd_weight(x[idx_blend], a = spec$a, b = spec$b, c1 = spec$c1, c2 = spec$c2)
    wprime <- .bgpd_weight_prime(x[idx_blend], a = spec$a, b = spec$b, c1 = spec$c1, c2 = spec$c2)
    H <- exp((1 - w) * log(Fbulk) + w * log(Fgp))
    dens <- H * (
      wprime * log(Fgp) +
      w * (fgp / Fgp) -
      wprime * log(Fbulk) +
      (1 - w) * (fbulk / Fbulk)
    )
    if (any(dens < -1e-8)) {
      stop(
        "The blended GP density becomes negative in the blending interval. ",
        "Choose a different bulk distribution or blending hyperparameters."
      )
    }
    dens[dens < 0] <- 0
    out[idx_blend] <- dens
  }

  if (any(idx_tail)) {
    out[idx_tail] <- .bgpd_gp_pdf(x[idx_tail], u = spec$u, sigma = spec$sigma, xi = spec$xi)
  }

  out[!is.na(x) & x < 0] <- 0
  if (spec$xi < 0) {
    out[x >= spec$endpoint] <- 0
  }

  if (log) {
    out <- ifelse(out > 0, log(out), -Inf)
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
  idx_tail <- !is.na(p) & p >= spec$pb

  if (any(idx_bulk)) {
    out[idx_bulk] <- .bgpd_bulk_call(spec$qbulk, p[idx_bulk], spec$bulk.args)
  }

  if (any(idx_blend)) {
    out[idx_blend] <- vapply(p[idx_blend], function(target) {
      uniroot(
        function(y) .pbgpd_from_spec(y, spec) - target,
        interval = c(spec$a, spec$b),
        tol = 1e-12
      )$root
    }, numeric(1))
  }

  if (any(idx_tail)) {
    out[idx_tail] <- .bgpd_gp_quantile(p[idx_tail], u = spec$u, sigma = spec$sigma, xi = spec$xi)
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
#' @param q quantiles.
#' @param xi GP tail shape parameter.
#' @param pbulk bulk CDF supported on \eqn{[0, 1]}.
#' @param dbulk bulk density supported on \eqn{[0, 1]}. This is used to
#'   validate the blended construction, even for \code{pbgpd()} and
#'   \code{qbgpd()}.
#' @param qbulk bulk quantile function supported on \eqn{[0, 1]}.
#' @param bulk.args named list of arguments passed to \code{pbulk},
#'   \code{dbulk}, and \code{qbulk}. Defaults to a Beta bulk.
#' @param pa lower bulk quantile used to start blending.
#' @param pb upper bulk quantile used to end blending and recover the exact GP
#'   tail.
#' @param c1,c2 Beta weighting-function shape parameters.
#'
#' @return CDF values.
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
