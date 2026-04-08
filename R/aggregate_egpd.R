## Multiscale aggregated EGPD framework based on compound Poisson EGPD models

.aggregate_duration_terms <- function(log_d, degree) {
  vapply(0:degree, function(k) log_d^k, numeric(length(log_d)))
}

.aggregate_eval_poly <- function(coef, durations) {
  log_d <- log(durations)
  as.vector(.aggregate_duration_terms(log_d, length(coef) - 1L) %*% coef)
}

.aggregate_windows <- function(x, d, overlap = TRUE) {
  x <- as.numeric(x)
  if (!isTRUE(all.equal(d, as.integer(d))) || d < 1L) {
    stop("'d' must be a positive integer", call. = FALSE)
  }
  d <- as.integer(d)

  if (overlap) {
    if (length(x) < d) {
      return(numeric())
    }
    cs <- c(0, cumsum(replace(x, is.na(x), 0)))
    na_cs <- c(0, cumsum(is.na(x)))
    sums <- cs[(d + 1L):(length(x) + 1L)] - cs[seq_len(length(x) - d + 1L)]
    miss <- na_cs[(d + 1L):(length(x) + 1L)] - na_cs[seq_len(length(x) - d + 1L)]
    return(sums[miss == 0L])
  }

  n_block <- floor(length(x) / d)
  if (n_block == 0L) {
    return(numeric())
  }
  mat <- matrix(x[seq_len(n_block * d)], ncol = d, byrow = TRUE)
  row_sums <- rowSums(replace(mat, is.na(mat), 0))
  row_sums[rowSums(is.na(mat)) == 0L]
}

.prepare_aggregate_input <- function(data, durations, overlap = TRUE) {
  if (missing(durations) || is.null(durations) || length(durations) == 0L) {
    stop("'durations' must contain at least one aggregation scale", call. = FALSE)
  }

  if (!is.numeric(durations) || any(!is.finite(durations))) {
    stop("'durations' must be a numeric vector of positive integers", call. = FALSE)
  }
  if (any(durations != floor(durations))) {
    stop("'durations' must be positive integers", call. = FALSE)
  }
  durations <- sort(unique(as.integer(durations)))
  if (any(durations < 1L)) {
    stop("'durations' must be positive integers", call. = FALSE)
  }

  if (is.numeric(data)) {
    if (any(data < 0, na.rm = TRUE)) {
      stop("'data' must be non-negative", call. = FALSE)
    }
    full_samples <- lapply(durations, function(d) .aggregate_windows(data, d, overlap = overlap))
    source <- "series"
  } else if (is.list(data)) {
    full_samples <- data
    if (!is.null(names(full_samples)) && all(as.character(durations) %in% names(full_samples))) {
      full_samples <- full_samples[as.character(durations)]
    } else if (length(full_samples) != length(durations)) {
      stop("When 'data' is a list, it must either be named by duration or have the same length as 'durations'",
           call. = FALSE)
    }
    full_samples <- lapply(full_samples, function(x) {
      x <- as.numeric(x)
      x[is.finite(x)]
    })
    if (any(vapply(full_samples, function(x) any(x < 0), logical(1)))) {
      stop("Aggregated samples must be non-negative", call. = FALSE)
    }
    source <- "aggregates"
  } else {
    stop("'data' must be either a numeric vector or a list of aggregated samples",
         call. = FALSE)
  }

  positive_samples <- lapply(full_samples, function(x) x[x > 0])
  if (any(lengths(positive_samples) == 0L)) {
    bad <- durations[lengths(positive_samples) == 0L]
    stop("No positive aggregated values are available for durations: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }

  compressed_positive <- lapply(positive_samples, function(x) {
    vals <- sort(unique(x))
    counts <- tabulate(match(x, vals), nbins = length(vals))
    list(values = vals, counts = counts)
  })

  list(
    source = source,
    durations = durations,
    full_samples = full_samples,
    positive_samples = positive_samples,
    compressed_positive = compressed_positive,
    n_total = lengths(full_samples),
    n_positive = lengths(positive_samples)
  )
}

.aggregate_positive_logdens <- function(x, family, sigma, kappa, xi, lambda,
                                        cpegpd.h = 0.2, K = NULL) {
  log_pos_prob <- log(-expm1(-lambda))
  if (!is.finite(log_pos_prob)) {
    return(rep(-Inf, length(x)))
  }

  if (family == "cpegpd") {
    dcpegpd(x, lambda = lambda, sigma = sigma, xi = xi, kappa = kappa,
            type = 1, h = cpegpd.h, K = K, log = TRUE) - log_pos_prob
  } else {
    dcpdegpd(x, lambda = lambda, sigma = sigma, xi = xi, kappa = kappa,
             type = 1, K = K, log = TRUE) - log_pos_prob
  }
}

.aggregate_positive_quantile <- function(p, family, sigma, kappa, xi, lambda,
                                         cpegpd.h = 0.2, K = NULL) {
  p <- pmin(pmax(p, 0), 1)
  p_full <- exp(-lambda) + p * (1 - exp(-lambda))
  if (family == "cpegpd") {
    qcpegpd(p_full, lambda = lambda, sigma = sigma, xi = xi, kappa = kappa,
            type = 1, h = cpegpd.h, K = K)
  } else {
    qcpdegpd(p_full, lambda = lambda, sigma = sigma, xi = xi, kappa = kappa,
             type = 1, K = K)
  }
}

.aggregate_feasible_regression <- function(durations, values, degree, slope = 1e-3) {
  x <- log(durations)
  X <- .aggregate_duration_terms(x, degree)
  fit <- lm.fit(X, log(pmax(values, 1e-8)))
  coef <- fit$coefficients
  coef[!is.finite(coef)] <- 0

  pred <- as.vector(X %*% coef)
  if (any(diff(pred) < 0)) {
    coef <- c(pred[1] - slope * x[1], rep(0, degree))
    if (degree >= 1L) {
      coef[2] <- max((pred[length(pred)] - pred[1]) / max(diff(range(x)), 1), slope)
    }
  }

  coef
}

.aggregate_auto_start <- function(prep, family, p, q, cpegpd.h = 0.2) {
  dur_fits <- lapply(prep$full_samples, function(x) {
    tryCatch(
      fitegpd(x, type = 1, family = family, hessian = FALSE,
              cpegpd.h = cpegpd.h, optim.method = "Nelder-Mead"),
      error = function(e) NULL
    )
  })

  sigma_hat <- vapply(seq_along(dur_fits), function(i) {
    fit <- dur_fits[[i]]
    if (!is.null(fit) && "sigma" %in% names(fit$estimate)) {
      return(unname(fit$estimate["sigma"]))
    }
    max(mean(prep$positive_samples[[i]]), 1e-2)
  }, numeric(1))

  lambda_hat <- vapply(seq_along(dur_fits), function(i) {
    fit <- dur_fits[[i]]
    if (!is.null(fit) && "lambda" %in% names(fit$estimate)) {
      return(unname(fit$estimate["lambda"]))
    }
    p0 <- mean(prep$full_samples[[i]] == 0)
    max(-log(pmax(p0, 1e-6)), 1e-2)
  }, numeric(1))

  kappa_hat <- vapply(dur_fits, function(fit) {
    if (!is.null(fit) && "kappa" %in% names(fit$estimate)) {
      unname(fit$estimate["kappa"])
    } else {
      0.3
    }
  }, numeric(1))

  xi_hat <- vapply(dur_fits, function(fit) {
    if (!is.null(fit) && "xi" %in% names(fit$estimate)) {
      unname(fit$estimate["xi"])
    } else {
      0.2
    }
  }, numeric(1))

  list(
    sigma_coef = .aggregate_feasible_regression(prep$durations, sigma_hat, p),
    lambda_coef = .aggregate_feasible_regression(prep$durations, lambda_hat, q),
    kappa = median(kappa_hat, na.rm = TRUE),
    xi = median(xi_hat, na.rm = TRUE)
  )
}

.bounded_to_theta <- function(x, bounds) {
  lo <- bounds[1]
  hi <- bounds[2]
  x <- pmin(pmax(x, lo + 1e-8), hi - 1e-8)
  qlogis((x - lo) / (hi - lo))
}

.theta_to_bounded <- function(theta, bounds) {
  lo <- bounds[1]
  hi <- bounds[2]
  lo + plogis(theta) * (hi - lo)
}

.normalize_aggregate_start <- function(start, p, q) {
  req <- c("sigma_coef", "lambda_coef", "kappa", "xi")
  if (!is.list(start) || !all(req %in% names(start))) {
    stop("'start' must be a list with entries 'sigma_coef', 'lambda_coef', 'kappa', and 'xi'",
         call. = FALSE)
  }
  if (length(start$sigma_coef) != p + 1L) {
    stop("'start$sigma_coef' must have length p + 1", call. = FALSE)
  }
  if (length(start$lambda_coef) != q + 1L) {
    stop("'start$lambda_coef' must have length q + 1", call. = FALSE)
  }
  if (!is.numeric(start$kappa) || length(start$kappa) != 1L || !is.finite(start$kappa) || start$kappa <= 0) {
    stop("'start$kappa' must be a positive number", call. = FALSE)
  }
  if (!is.numeric(start$xi) || length(start$xi) != 1L || !is.finite(start$xi) || start$xi <= 0) {
    stop("'start$xi' must be a positive number", call. = FALSE)
  }

  list(
    sigma_coef = as.numeric(start$sigma_coef),
    lambda_coef = as.numeric(start$lambda_coef),
    kappa = as.numeric(start$kappa),
    xi = as.numeric(start$xi)
  )
}

.aggregate_pack_theta <- function(start, kappa_bounds, xi_bounds) {
  c(start$sigma_coef,
    start$lambda_coef,
    .bounded_to_theta(start$kappa, kappa_bounds),
    .bounded_to_theta(start$xi, xi_bounds))
}

.aggregate_unpack_theta <- function(theta, p, q, kappa_bounds, xi_bounds) {
  s_idx <- seq_len(p + 1L)
  l_idx <- (max(s_idx) + 1L):(max(s_idx) + q + 1L)
  list(
    sigma_coef = theta[s_idx],
    lambda_coef = theta[l_idx],
    kappa = .theta_to_bounded(theta[length(theta) - 1L], kappa_bounds),
    xi = .theta_to_bounded(theta[length(theta)], xi_bounds)
  )
}

.aggregate_parameter_frame <- function(object, durations = object$durations) {
  sigma <- exp(.aggregate_eval_poly(object$sigma_coef, durations))
  lambda <- exp(.aggregate_eval_poly(object$lambda_coef, durations))
  data.frame(
    duration = durations,
    sigma = sigma,
    lambda = lambda,
    kappa = rep(object$kappa, length(durations)),
    xi = rep(object$xi, length(durations))
  )
}

.aggregate_nll <- function(theta, prep, family, p, q, cpegpd.h = 0.2,
                           kappa_bounds, xi_bounds) {
  pars <- .aggregate_unpack_theta(theta, p, q, kappa_bounds, xi_bounds)
  sigma_d <- exp(.aggregate_eval_poly(pars$sigma_coef, prep$durations))
  lambda_d <- exp(.aggregate_eval_poly(pars$lambda_coef, prep$durations))

  if (any(!is.finite(sigma_d)) || any(!is.finite(lambda_d)) ||
      !is.finite(pars$kappa) || !is.finite(pars$xi)) {
    return(1e20)
  }

  ll <- 0
  for (i in seq_along(prep$durations)) {
    comp <- prep$compressed_positive[[i]]
    K <- if (family == "cpegpd") {
      max(ceiling(max(prep$full_samples[[i]], 0) / cpegpd.h) + 50L, 100L)
    } else {
      max(max(prep$full_samples[[i]], 0), 100L) + 50L
    }
    ld <- .aggregate_positive_logdens(
      comp$values,
      family = family,
      sigma = sigma_d[i],
      kappa = pars$kappa,
      xi = pars$xi,
      lambda = lambda_d[i],
      cpegpd.h = cpegpd.h,
      K = K
    )
    if (any(!is.finite(ld))) {
      return(1e20)
    }
    ll <- ll + sum(comp$counts * ld)
  }

  -ll
}

.aggregate_penalized_nll <- function(theta, prep, family, p, q, cpegpd.h = 0.2,
                                     kappa_bounds, xi_bounds,
                                     penalty_scale = 1e8) {
  pars <- .aggregate_unpack_theta(theta, p, q, kappa_bounds, xi_bounds)
  sigma_lp <- .aggregate_eval_poly(pars$sigma_coef, prep$durations)
  lambda_lp <- .aggregate_eval_poly(pars$lambda_coef, prep$durations)
  sigma_diff <- diff(sigma_lp)
  lambda_diff <- diff(lambda_lp)

  if (any(!is.finite(sigma_diff)) || any(!is.finite(lambda_diff))) {
    return(1e20)
  }

  violations <- c(sigma_diff, lambda_diff)
  if (any(violations < 0)) {
    neg <- pmin(violations, 0)
    penalty <- 1e12 + penalty_scale * sum(neg^2)
  } else {
    penalty <- 0
  }

  .aggregate_nll(theta, prep = prep, family = family, p = p, q = q,
                 cpegpd.h = cpegpd.h,
                 kappa_bounds = kappa_bounds,
                 xi_bounds = xi_bounds) + penalty
}

#' Fit the multiscale aggregated EGPD framework
#'
#' Fits the positive-part compound Poisson EGPD framework of Ailliot,
#' Gaetan, and Naveau (2026) across several aggregation durations. The
#' model assumes that positive aggregates at duration \code{d} follow a
#' compound Poisson EGPD law with shared shape parameters and
#' duration-specific scale and Poisson mean:
#' \deqn{A_d \sim EGPD(\sigma_d, \kappa, \xi, \lambda_d)}
#' where \eqn{\log \sigma_d} and \eqn{\log \lambda_d} are polynomials in
#' \eqn{\log d}.
#'
#' @param data either a numeric vector of non-negative observations at the
#'   finest time scale, or a list of non-negative aggregated samples for each
#'   duration.
#' @param durations positive integer aggregation durations.
#' @param family character string, either \code{"cpegpd"} for continuous
#'   aggregates or \code{"cpdegpd"} for discrete aggregates.
#' @param p degree of the polynomial used for \eqn{\log \sigma_d}.
#' @param q degree of the polynomial used for \eqn{\log \lambda_d}.
#' @param overlap logical; if \code{TRUE}, rolling sums are used when
#'   \code{data} is a finest-scale series.
#' @param cpegpd.h discretization step used by the Panjer recursion for the
#'   continuous compound family.
#' @param start optional named list with entries \code{sigma_coef},
#'   \code{lambda_coef}, \code{kappa}, and \code{xi}.
#' @param kappa_bounds numeric vector of length two giving lower and upper
#'   bounds for \eqn{\kappa}.
#' @param xi_bounds numeric vector of length two giving lower and upper bounds
#'   for \eqn{\xi}.
#' @param optim.method optimization method passed to \code{\link{optim}}.
#' @param control optional optimization control list.
#' @param keep_data logical; should aggregated samples be stored in the fitted
#'   object?
#'
#' @return An object of class \code{"aggregate_egpd"}.
#' @export
fit_aggregate_egpd <- function(data, durations,
                               family = c("cpegpd", "cpdegpd"),
                               p = 2L, q = 2L,
                               overlap = TRUE,
                               cpegpd.h = 0.2,
                               start = NULL,
                               kappa_bounds = c(0.05, 5),
                               xi_bounds = c(0.01, 1.5),
                               optim.method = "Nelder-Mead",
                               control = list(),
                               keep_data = FALSE) {
  cl <- match.call()
  family <- match.arg(family)
  p_in <- p
  q_in <- q
  if (!is.numeric(p_in) || length(p_in) != 1L || !is.finite(p_in) ||
      p_in < 0L || p_in != floor(p_in)) {
    stop("'p' must be a non-negative integer", call. = FALSE)
  }
  if (!is.numeric(q_in) || length(q_in) != 1L || !is.finite(q_in) ||
      q_in < 0L || q_in != floor(q_in)) {
    stop("'q' must be a non-negative integer", call. = FALSE)
  }
  p <- as.integer(p_in)
  q <- as.integer(q_in)
  if (!is.numeric(kappa_bounds) || length(kappa_bounds) != 2L ||
      any(!is.finite(kappa_bounds)) || kappa_bounds[1] <= 0 ||
      kappa_bounds[1] >= kappa_bounds[2]) {
    stop("'kappa_bounds' must contain two increasing positive values", call. = FALSE)
  }
  if (!is.numeric(xi_bounds) || length(xi_bounds) != 2L ||
      any(!is.finite(xi_bounds)) || xi_bounds[1] <= 0 ||
      xi_bounds[1] >= xi_bounds[2]) {
    stop("'xi_bounds' must contain two increasing positive values", call. = FALSE)
  }

  prep <- .prepare_aggregate_input(data, durations, overlap = overlap)
  if (length(prep$durations) < 2L) {
    stop("At least two durations are required to fit the multiscale framework",
         call. = FALSE)
  }
  if (family == "cpdegpd") {
    discrete_ok <- all(vapply(prep$full_samples, function(x) all(x == floor(x)), logical(1)))
    if (!discrete_ok) {
      stop("cpdegpd requires integer-valued aggregated samples", call. = FALSE)
    }
  }

  start_vals <- if (is.null(start)) {
    .aggregate_auto_start(prep, family = family, p = p, q = q,
                          cpegpd.h = cpegpd.h)
  } else {
    .normalize_aggregate_start(start, p = p, q = q)
  }

  theta0 <- .aggregate_pack_theta(start_vals, kappa_bounds = kappa_bounds,
                                  xi_bounds = xi_bounds)
  opt_control <- modifyList(list(maxit = 500L), control)

  opt <- optim(
    par = theta0,
    fn = .aggregate_penalized_nll,
    prep = prep,
    family = family,
    p = p,
    q = q,
    cpegpd.h = cpegpd.h,
    kappa_bounds = kappa_bounds,
    xi_bounds = xi_bounds,
    method = optim.method,
    control = opt_control,
    hessian = FALSE
  )

  pars <- .aggregate_unpack_theta(opt$par, p = p, q = q,
                                  kappa_bounds = kappa_bounds,
                                  xi_bounds = xi_bounds)
  sigma_names <- paste0("s", 0:p)
  lambda_names <- paste0("l", 0:q)
  names(pars$sigma_coef) <- sigma_names
  names(pars$lambda_coef) <- lambda_names

  object <- structure(list(
    sigma_coef = pars$sigma_coef,
    lambda_coef = pars$lambda_coef,
    kappa = pars$kappa,
    xi = pars$xi,
    durations = prep$durations,
    n_total = prep$n_total,
    n_positive = prep$n_positive,
    family = family,
    p = p,
    q = q,
    overlap = overlap,
    cpegpd.h = if (family == "cpegpd") cpegpd.h else NULL,
    kappa_bounds = kappa_bounds,
    xi_bounds = xi_bounds,
    composite_loglik = -opt$value,
    convergence = opt$convergence,
    message = opt$message,
    optim = opt,
    call = cl,
    data_source = prep$source,
    data = if (isTRUE(keep_data)) prep$full_samples else NULL
  ), class = "aggregate_egpd")

  object$fitted <- .aggregate_parameter_frame(object)
  object
}

#' @export
coef.aggregate_egpd <- function(object, ...) {
  c(object$sigma_coef, object$lambda_coef,
    kappa = object$kappa, xi = object$xi)
}

#' @export
predict.aggregate_egpd <- function(object, durations = object$durations,
                                   type = c("parameter", "quantile"),
                                   prob = NULL, ...) {
  type <- match.arg(type)
  if (any(!is.finite(durations)) || any(durations != floor(durations))) {
    stop("'durations' must be positive integers", call. = FALSE)
  }
  durations <- sort(unique(as.integer(durations)))
  if (any(durations < 1L)) {
    stop("'durations' must be positive integers", call. = FALSE)
  }

  pars <- .aggregate_parameter_frame(object, durations = durations)
  if (type == "parameter") {
    return(pars)
  }

  if (is.null(prob)) {
    stop("'prob' must be supplied when type = 'quantile'", call. = FALSE)
  }
  prob <- as.numeric(prob)
  if (any(!is.finite(prob)) || any(prob < 0 | prob > 1)) {
    stop("'prob' must contain probabilities in [0, 1]", call. = FALSE)
  }

  out <- matrix(NA_real_, nrow = nrow(pars), ncol = length(prob),
                dimnames = list(as.character(durations), as.character(prob)))
  for (i in seq_len(nrow(pars))) {
    out[i, ] <- .aggregate_positive_quantile(
      p = prob,
      family = object$family,
      sigma = pars$sigma[i],
      kappa = object$kappa,
      xi = object$xi,
      lambda = pars$lambda[i],
      cpegpd.h = if (object$family == "cpegpd") object$cpegpd.h else 0.2
    )
  }
  out
}

#' @export
print.aggregate_egpd <- function(x, ...) {
  cat("Multiscale aggregated ", if (x$family == "cpegpd") "EGPD" else "discrete EGPD",
      " fit\n", sep = "")
  cat("Durations:", paste(x$durations, collapse = ", "), "\n")
  cat("Composite log-likelihood:", round(x$composite_loglik, 2), "\n")
  cat("Convergence:", if (x$convergence == 0) "successful" else paste("code", x$convergence), "\n")
  invisible(x)
}

#' @export
summary.aggregate_egpd <- function(object, ...) {
  structure(list(
    coefficients = coef(object),
    fitted = object$fitted,
    n_total = object$n_total,
    n_positive = object$n_positive,
    family = object$family,
    composite_loglik = object$composite_loglik,
    convergence = object$convergence,
    call = object$call
  ), class = "summary.aggregate_egpd")
}

#' @export
print.summary.aggregate_egpd <- function(x, digits = 4, ...) {
  cat("Multiscale aggregated ", if (x$family == "cpegpd") "EGPD" else "discrete EGPD",
      " fit\n", sep = "")
  cat("Composite log-likelihood:", round(x$composite_loglik, digits), "\n")
  cat("Convergence:", if (x$convergence == 0) "successful" else paste("code", x$convergence), "\n\n")
  cat("Coefficients:\n")
  print(round(x$coefficients, digits))
  cat("\nDuration-specific parameters:\n")
  print(round(x$fitted, digits), row.names = FALSE)
  invisible(x)
}

#' @export
plot.aggregate_egpd <- function(x, ...) {
  pars <- predict(x, type = "parameter")
  op <- par(mfrow = c(1, 2))
  on.exit(par(op), add = TRUE)

  plot(pars$duration, pars$sigma, type = "b", log = "x",
       xlab = "Duration", ylab = expression(sigma[d]),
       main = expression("Estimated " * sigma[d]), ...)
  plot(pars$duration, pars$lambda, type = "b", log = "x",
       xlab = "Duration", ylab = expression(lambda[d]),
       main = expression("Estimated " * lambda[d]), ...)
  invisible(x)
}
