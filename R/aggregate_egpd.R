## Multiscale aggregated EGPD framework based on compound Poisson EGPD models

.aggregate_duration_terms <- function(log_d, degree) {
  vapply(0:degree, function(k) log_d^k, numeric(length(log_d)))
}

.aggregate_eval_poly <- function(coef, durations) {
  log_d <- log(durations)
  as.vector(.aggregate_duration_terms(log_d, length(coef) - 1L) %*% coef)
}

## Integer grid units round(x / grid) when every observed value lies on the
## lattice {0, grid, 2 grid, ...} (up to floating-point error); NULL otherwise.
.aggregate_grid_units <- function(x, grid) {
  if (is.null(grid)) return(NULL)
  obs <- x[!is.na(x)]
  if (!length(obs) || any(!is.finite(obs))) return(NULL)
  scaled <- obs / grid
  units <- round(scaled)
  ## Zero is a separate atom, never a rounded positive gauge reading.
  if (any(!is.finite(scaled)) || any(units == 0 & obs != 0) ||
      any(abs(scaled - units) > .compound_grid_tol * pmax(1, abs(scaled))))
    return(NULL)
  if (sum(abs(units)) >= 2^52) return(NULL)
  out <- rep(NA_real_, length(x))
  out[!is.na(x)] <- units
  out
}

## Window sums without differencing a running total of non-integer values
## (its absolute error grows with the series length); integer-valued input is
## summed exactly. Windows containing NA are dropped.
.aggregate_window_sums <- function(x, d, overlap) {
  n <- length(x)
  if (overlap) {
    if (n < d) {
      return(numeric())
    }
    if (d == 1L) {
      return(x[!is.na(x)])
    }
    obs <- x[!is.na(x)]
    if (all(obs == round(obs)) && sum(abs(obs)) < 2^52) {
      cs <- c(0, cumsum(replace(x, is.na(x), 0)))
      na_cs <- c(0, cumsum(is.na(x)))
      idx <- seq_len(n - d + 1L)
      sums <- cs[idx + d] - cs[idx]
      return(sums[na_cs[idx + d] - na_cs[idx] == 0L])
    }
    sums <- as.numeric(stats::filter(x, rep(1, d), method = "convolution",
                                     sides = 1L))[d:n]
    return(sums[!is.na(sums)])
  }

  n_block <- floor(n / d)
  if (n_block == 0L) {
    return(numeric())
  }
  mat <- matrix(x[seq_len(n_block * d)], ncol = d, byrow = TRUE)
  row_sums <- rowSums(replace(mat, is.na(mat), 0))
  row_sums[rowSums(is.na(mat)) == 0L]
}

.aggregate_windows <- function(x, d, overlap = TRUE, grid = NULL) {
  x <- as.numeric(x)
  if (!isTRUE(all.equal(d, as.integer(d))) || d < 1L) {
    stop("'d' must be a positive integer", call. = FALSE)
  }
  d <- as.integer(d)

  ## On the gauge lattice, sum integer units exactly and return k * grid,
  ## the same representation as the lattice values of dcpegpd()/qcpegpd().
  units <- .aggregate_grid_units(x, grid)
  if (!is.null(units)) {
    return(.aggregate_window_sums(units, d, overlap) * grid)
  }
  .aggregate_window_sums(x, d, overlap)
}

.prepare_aggregate_input <- function(data, durations, overlap = TRUE,
                                     grid = NULL) {
  if (missing(durations) || is.null(durations) || length(durations) == 0L) {
    stop("'durations' must contain at least one aggregation scale", call. = FALSE)
  }

  if (!is.numeric(durations) || any(!is.finite(durations))) {
    stop("'durations' must be a numeric vector of positive integers", call. = FALSE)
  }
  if (any(durations != floor(durations) | durations > .Machine$integer.max)) {
    stop("'durations' must be positive integers", call. = FALSE)
  }
  input_durations <- durations
  durations <- sort(unique(as.integer(durations)))
  if (any(durations < 1L)) {
    stop("'durations' must be positive integers", call. = FALSE)
  }

  if (is.numeric(data)) {
    if (any(data < 0, na.rm = TRUE)) {
      stop("'data' must be non-negative", call. = FALSE)
    }
    full_samples <- lapply(durations, function(d)
      .aggregate_windows(data, d, overlap = overlap, grid = grid))
    source <- "series"
  } else if (is.list(data)) {
    full_samples <- data
    if (!is.null(names(full_samples))) {
      if (anyDuplicated(names(full_samples)) ||
          !all(as.character(durations) %in% names(full_samples)))
        stop("Named aggregate samples must be uniquely named by duration",
             call. = FALSE)
      full_samples <- full_samples[as.character(durations)]
    } else if (length(full_samples) != length(durations)) {
      stop("When 'data' is a list, it must either be named by duration or have the same length as 'durations'",
           call. = FALSE)
    } else {
      if (anyDuplicated(input_durations))
        stop("Unnamed aggregate samples require distinct durations", call. = FALSE)
      full_samples <- full_samples[order(input_durations)]
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

.aggregate_logpos <- function(family, sigma, kappa, xi, lambda) {
  lr <- log(lambda)
  if (family == "cpdegpd")
    lr <- lr + pegpd(1, sigma = sigma, xi = xi, kappa = kappa, type = 1,
                      lower.tail = FALSE, log.p = TRUE)
  .compound_log_nonzero(lr)
}

.aggregate_positive_logdens <- function(x, family, sigma, kappa, xi, lambda,
                                        cpegpd.h = 0.2, K = NULL) {
  out <- if (family == "cpegpd") {
    dcpegpd(x, lambda = lambda, sigma = sigma, xi = xi, kappa = kappa,
            type = 1, h = cpegpd.h, K = K, log = TRUE)
  } else {
    dcpdegpd(x, lambda = lambda, sigma = sigma, xi = xi, kappa = kappa,
             type = 1, K = K, log = TRUE)
  }
  log_pos_prob <- .aggregate_logpos(family, sigma, kappa, xi, lambda)
  if (!is.finite(log_pos_prob)) return(rep(-Inf, length(x)))
  out <- out - log_pos_prob
  out[!is.na(x) & x <= 0] <- -Inf
  out
}

## Invert on the conditional scale; transforming p to the full CDF can lose
## it entirely when the effective rate of positive jumps is small.
.aggregate_positive_quantile <- function(p, family, sigma, kappa, xi, lambda,
                                         cpegpd.h = 0.2, K = NULL) {
  if (any(!is.finite(p) | p < 0 | p > 1))
    stop("'p' must contain probabilities in [0, 1]")
  h <- if (family == "cpegpd") cpegpd.h else 1
  K <- .compound_grid_size(if (is.null(K)) 5000 else K, h = h)
  .validate_compound_params(lambda = lambda, sigma = sigma, xi = xi,
                            kappa = kappa, type = 1)
  .compound_quantiles(p, lambda, sigma, xi, kappa, NA, NA, 1L, K,
                      h = if (family == "cpegpd") h, positive = TRUE)
}

.aggregate_feasible_regression <- function(durations, values, degree, slope = 1e-3) {
  x <- log(durations)
  y <- log(pmax(values, 1e-8))
  X <- .aggregate_duration_terms(x, degree)
  fit <- lm.fit(X, y)
  coef <- fit$coefficients
  coef[!is.finite(coef)] <- 0

  if (degree >= 1L && any(coef[-1L] < 0)) {
    init_b <- c(mean(y), rep(slope, degree))
    nnls <- tryCatch(
      optim(
        par = init_b,
        fn = function(b) sum((y - as.vector(X %*% b))^2),
        gr = function(b) -2 * as.vector(crossprod(X, y - as.vector(X %*% b))),
        method = "L-BFGS-B",
        lower = c(-Inf, rep(0, degree))
      )$par,
      error = function(e) NULL
    )
    if (!is.null(nnls) && all(is.finite(nnls))) {
      coef <- c(nnls[1L], pmax(nnls[-1L], 0))
    } else {
      coef <- c(y[1L], rep(0, degree))
      coef[2L] <- max((y[length(y)] - y[1L]) / max(diff(range(x)), 1), 0)
    }
  }

  coef
}

## E[X] / sigma for the type-1 EGPD (G(u) = u^kappa), 0 <= xi < 1:
## X = H^{-1}(V) with V ~ Beta(kappa, 1) gives kappa B(kappa, 1 - xi).
.egpd1_mean_factor <- function(kappa, xi) {
  if (abs(xi) < 1e-8) return(digamma(kappa + 1) - digamma(1))
  (kappa * beta(kappa, 1 - xi) - 1) / xi
}

## Keep automatic starts off the flat ends of the logit-bounded scale.
.aggregate_interior <- function(x, bounds, margin = 0.02) {
  width <- bounds[2] - bounds[1]
  v <- min(max(x, bounds[1] + margin * width), bounds[2] - margin * width)
  if (v <= bounds[1] || v >= bounds[2]) v <- bounds[1] + width / 2
  if (v <= bounds[1] || v >= bounds[2])
    stop("Shape bounds have no representable interior.", call. = FALSE)
  v
}

## Moment start for one duration: lambda from the dry fraction (if informative),
## sigma from the mean positive aggregate E[A | A > 0] = lambda E[X] / (1 - e^-lambda).
.aggregate_moment_start <- function(full, pos, kappa, xi) {
  p0 <- mean(full == 0)
  lambda <- if (p0 > 0 && p0 < 1) -log(p0) else 1
  per_event <- mean(pos) * (-expm1(-lambda)) / lambda
  ## the mean is infinite for xi >= 1; a capped shape only sets the scale
  c(sigma = max(per_event / .egpd1_mean_factor(kappa, min(xi, 0.9)), 1e-8),
    lambda = lambda)
}

## Reject absurd single-duration fits (e.g. kappa -> 1e9 with sigma -> 0)
## before they seed the composite optimisation.
.aggregate_start_ok <- function(est, pos, kappa_bounds, xi_bounds) {
  nm <- c("sigma", "kappa", "xi", "lambda")
  if (is.null(est) || !all(nm %in% names(est))) return(FALSE)
  est <- est[nm]
  if (any(!is.finite(est))) return(FALSE)
  unname(
    est["kappa"] > kappa_bounds[1] && est["kappa"] < kappa_bounds[2] &&
      est["xi"] > xi_bounds[1] && est["xi"] < xi_bounds[2] &&
      est["sigma"] > 1e-3 * stats::median(pos) && est["sigma"] < 1e3 * max(pos) &&
      est["lambda"] > 1e-3 && est["lambda"] < 1e3
  )
}

.aggregate_moment_candidate <- function(prep, p, q, kappa, xi) {
  per <- vapply(seq_along(prep$durations), function(i)
    .aggregate_moment_start(prep$full_samples[[i]], prep$positive_samples[[i]],
                            kappa, xi), numeric(2))
  list(
    sigma_coef = .aggregate_feasible_regression(prep$durations, per[1L, ], p),
    lambda_coef = .aggregate_feasible_regression(prep$durations, per[2L, ], q),
    kappa = kappa,
    xi = xi
  )
}

.aggregate_auto_start <- function(prep, family, p, q, cpegpd.h = 0.2,
                                  kappa_bounds = c(0.05, 5),
                                  xi_bounds = c(0.01, 1.5)) {
  est <- lapply(prep$full_samples, function(x) {
    fit <- tryCatch(
      withCallingHandlers(distfit(x, type = 1, family = family, hessian = FALSE,
                                  cpegpd.h = cpegpd.h, optim.method = "Nelder-Mead"),
                          distfit_nonconvergence = function(w) invokeRestart("muffleWarning")),
      error = function(e) NULL
    )
    if (!isTRUE(fit$convergence == 0L)) NULL else fit$estimate
  })
  ok <- vapply(seq_along(est), function(i)
    .aggregate_start_ok(est[[i]], prep$positive_samples[[i]],
                        kappa_bounds, xi_bounds), logical(1))

  kappa0 <- if (any(ok)) stats::median(vapply(est[ok], function(e) e[["kappa"]], 0)) else 1
  xi0 <- if (any(ok)) stats::median(vapply(est[ok], function(e) e[["xi"]], 0)) else 0.2
  kappa0 <- .aggregate_interior(kappa0, kappa_bounds)
  xi0 <- .aggregate_interior(xi0, xi_bounds)

  per <- vapply(seq_along(est), function(i) {
    if (ok[i]) {
      unname(est[[i]][c("sigma", "lambda")])
    } else {
      unname(.aggregate_moment_start(prep$full_samples[[i]],
                                     prep$positive_samples[[i]], kappa0, xi0))
    }
  }, numeric(2))

  list(
    sigma_coef = .aggregate_feasible_regression(prep$durations, per[1L, ], p),
    lambda_coef = .aggregate_feasible_regression(prep$durations, per[2L, ], q),
    kappa = kappa0,
    xi = xi0
  )
}

## Candidate automatic starts: per-duration fits (sanitised) plus moment starts
## over a few shapes; the composite likelihood picks the best.
.aggregate_start_candidates <- function(prep, family, p, q, cpegpd.h,
                                        kappa_bounds, xi_bounds) {
  auto <- .aggregate_auto_start(prep, family = family, p = p, q = q,
                                cpegpd.h = cpegpd.h,
                                kappa_bounds = kappa_bounds,
                                xi_bounds = xi_bounds)
  xi_m <- .aggregate_interior(0.2, xi_bounds)
  moments <- lapply(c(0.5, 1, 2), function(k)
    .aggregate_moment_candidate(prep, p, q,
                                .aggregate_interior(k, kappa_bounds), xi_m))
  c(list(auto), moments)
}

.bounded_to_theta <- function(x, bounds) {
  log(x - bounds[1]) - log(bounds[2] - x)
}

.theta_to_bounded <- function(theta, bounds) {
  lo <- bounds[1]
  hi <- bounds[2]
  if (theta > 0) hi - plogis(theta, lower.tail = FALSE) * (hi - lo)
  else lo + exp(plogis(theta, log.p = TRUE) + log(hi - lo))
}

## Monotonicity constraint for a degree-`degree` polynomial P in t = log d:
## M %*% coef[-1] >= 0. With upto = Inf this is s_1, ..., s_p >= 0 (sufficient
## for P' >= 0 at every d >= 1). With finite upto = D it is the vector of
## Bernstein coefficients of P'(t) on [0, log D], so P' >= 0 on [1, D]
## (necessary and sufficient for degree <= 2, sufficient otherwise).
.aggregate_monotone_matrix <- function(degree, upto = Inf) {
  if (degree < 1L) return(matrix(0, 0L, 0L))
  if (!is.finite(upto)) return(diag(1, degree))
  tmax <- log(upto)
  m <- degree - 1L
  M <- matrix(0, degree, degree)
  for (j in 0:m) {
    for (k in 0:j) {
      M[j + 1L, k + 1L] <- choose(j, k) / choose(m, k) * (k + 1) * tmax^k
    }
  }
  M
}

.aggregate_slope_index <- function(p, q) {
  list(s = if (p >= 1L) 2:(p + 1L) else integer(),
       l = if (q >= 1L) (p + 3L):(p + q + 2L) else integer())
}

## Map theta <-> z, where z holds the constrained combinations M %*% slopes,
## so that L-BFGS-B can impose the constraint as the box z >= 0.
.aggregate_theta_to_z <- function(theta, p, q, Ms, Ml) {
  idx <- .aggregate_slope_index(p, q)
  z <- theta
  if (p >= 1L) z[idx$s] <- as.vector(Ms %*% theta[idx$s])
  if (q >= 1L) z[idx$l] <- as.vector(Ml %*% theta[idx$l])
  z
}

.aggregate_z_to_theta <- function(z, p, q, Ms, Ml) {
  idx <- .aggregate_slope_index(p, q)
  theta <- z
  if (p >= 1L) theta[idx$s] <- solve(Ms, z[idx$s])
  if (q >= 1L) theta[idx$l] <- solve(Ml, z[idx$l])
  theta
}

.aggregate_project_monotone <- function(slopes, M) {
  if (!length(slopes)) return(slopes)
  if (identical(M, diag(1, nrow(M)))) return(pmax(slopes, 0))
  solve(M, pmax(as.vector(M %*% slopes), 0))
}

## The "cpegpd" likelihood is the law of gauge-discretised sums. Warn when
## positive aggregates are off the h-lattice and h is coarse relative to them,
## where the upper-bin model is biased for continuous data (about h / 2 per event).
.aggregate_offgrid_warning <- function(prep, h) {
  vals <- unlist(prep$positive_samples, use.names = FALSE)
  if (!length(vals) || !is.null(.aggregate_grid_units(vals, h)))
    return(invisible(FALSE))
  med <- stats::median(prep$positive_samples[[1L]])
  if (h > 0.01 * med) {
    warning(sprintf(paste0(
      "Positive aggregates are not on the cpegpd.h = %g lattice. ",
      "family = \"cpegpd\" is the law of gauge-discretised sums (each event ",
      "amount rounded up to a multiple of h; Ailliot et al., Appendix F); for ",
      "continuous data with h > 1%% of the median positive value (%.3g) the fit ",
      "is biased by about h/2 per event. Use gauge-resolution data or a smaller ",
      "cpegpd.h."), h, med), call. = FALSE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

.normalize_aggregate_start <- function(start, p, q,
                                       kappa_bounds = c(0.05, 5),
                                       xi_bounds = c(0.01, 1.5),
                                       Ms = .aggregate_monotone_matrix(p),
                                       Ml = .aggregate_monotone_matrix(q)) {
  req <- c("sigma_coef", "lambda_coef", "kappa", "xi")
  if (!is.list(start) || !all(req %in% names(start))) {
    stop("'start' must be a list with entries 'sigma_coef', 'lambda_coef', 'kappa', and 'xi'",
         call. = FALSE)
  }
  if (!is.numeric(start$sigma_coef) || length(start$sigma_coef) != p + 1L ||
      any(!is.finite(start$sigma_coef))) {
    stop("'start$sigma_coef' must be a finite numeric vector of length p + 1",
         call. = FALSE)
  }
  if (!is.numeric(start$lambda_coef) || length(start$lambda_coef) != q + 1L ||
      any(!is.finite(start$lambda_coef))) {
    stop("'start$lambda_coef' must be a finite numeric vector of length q + 1",
         call. = FALSE)
  }
  if ((p >= 1L && any(Ms %*% start$sigma_coef[-1L] < -1e-10)) ||
      (q >= 1L && any(Ml %*% start$lambda_coef[-1L] < -1e-10))) {
    stop("'start$sigma_coef[-1]' and 'start$lambda_coef[-1]' must satisfy the non-negative monotonicity constraint (see 'monotone_upto') so fitted curves are non-decreasing",
         call. = FALSE)
  }
  if (!is.numeric(start$kappa) || length(start$kappa) != 1L || !is.finite(start$kappa) || start$kappa <= 0) {
    stop("'start$kappa' must be a positive number", call. = FALSE)
  }
  if (!is.numeric(start$xi) || length(start$xi) != 1L || !is.finite(start$xi) || start$xi <= 0) {
    stop("'start$xi' must be a positive number", call. = FALSE)
  }
  if (start$kappa <= kappa_bounds[1] || start$kappa >= kappa_bounds[2]) {
    stop("'start$kappa' must lie strictly inside 'kappa_bounds'", call. = FALSE)
  }
  if (start$xi <= xi_bounds[1] || start$xi >= xi_bounds[2]) {
    stop("'start$xi' must lie strictly inside 'xi_bounds'", call. = FALSE)
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

  if (any(!is.finite(sigma_d) | sigma_d <= 0) ||
      any(!is.finite(lambda_d) | lambda_d <= 0) ||
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
                                     penalty_scale = 1e5,
                                     Ms = .aggregate_monotone_matrix(p),
                                     Ml = .aggregate_monotone_matrix(q)) {
  pars <- .aggregate_unpack_theta(theta, p, q, kappa_bounds, xi_bounds)
  slope_coefs <- c(pars$sigma_coef[-1L], pars$lambda_coef[-1L])
  if (any(!is.finite(slope_coefs))) {
    return(1e20)
  }

  ## Smooth C^1 quadratic hinge penalty on negative constrained combinations
  ## M %*% slopes (the slopes themselves when monotone_upto = Inf), so that
  ## sigma_d and lambda_d are non-decreasing without a step cliff at zero.
  cons <- c(if (p >= 1L) as.vector(Ms %*% pars$sigma_coef[-1L]),
            if (q >= 1L) as.vector(Ml %*% pars$lambda_coef[-1L]))
  neg <- pmin(cons, 0)
  penalty <- penalty_scale * sum(neg^2)

  nll <- .aggregate_nll(theta, prep = prep, family = family, p = p, q = q,
                        cpegpd.h = cpegpd.h,
                        kappa_bounds = kappa_bounds,
                        xi_bounds = xi_bounds)
  if (!is.finite(nll) || nll >= 1e20) {
    return(1e20)
  }
  nll + penalty
}

## Project the slopes of theta onto the monotonicity constraint (penalised
## runs may end slightly outside it).
.aggregate_project_theta <- function(theta, p, q, Ms, Ml, kappa_bounds,
                                     xi_bounds) {
  pars <- .aggregate_unpack_theta(theta, p = p, q = q,
                                  kappa_bounds = kappa_bounds,
                                  xi_bounds = xi_bounds)
  if (p >= 1L)
    pars$sigma_coef[-1L] <- .aggregate_project_monotone(pars$sigma_coef[-1L], Ms)
  if (q >= 1L)
    pars$lambda_coef[-1L] <- .aggregate_project_monotone(pars$lambda_coef[-1L], Ml)
  ## Only the slopes were projected. Keep the original shape coordinates
  ## even if their response-scale values rounded onto a bound.
  theta[seq_len(p + q + 2L)] <- c(pars$sigma_coef, pars$lambda_coef)
  list(pars = pars, theta = theta)
}

## One optim() run from theta0. L-BFGS-B imposes the monotonicity constraint
## as box bounds on z = (s0, Ms %*% slopes, l0, Ml %*% slopes, ...) (with
## monotone_upto = Inf, z is theta itself); other methods minimise the
## penalised objective.
.aggregate_optim_run <- function(theta0, method, control, nll_args, Ms, Ml,
                                 monotone_upto) {
  p <- nll_args$p
  q <- nll_args$q
  if (method == "L-BFGS-B") {
    ## Bound the two shapes directly. Differencing logits near a bound
    ## can move the response parameter by less than numerical precision.
    ix <- tail(seq_along(theta0), 2L)
    kb <- nll_args$kappa_bounds; xb <- nll_args$xi_bounds
    lower_vec <- c(-Inf, rep(0, p), -Inf, rep(0, q), kb[1L], xb[1L])
    upper_vec <- c(rep(Inf, p + q + 2L), kb[2L], xb[2L])
    to_theta <- function(z) {
      if (is.finite(monotone_upto)) z <- .aggregate_z_to_theta(z, p, q, Ms, Ml)
      ## A finite-difference probe can round one ulp past a box face.
      ## Map it to the face, not to the log of a negative distance.
      z[ix] <- pmin(upper_vec[ix], pmax(lower_vec[ix], z[ix]))
      z[ix] <- c(.bounded_to_theta(z[ix[1L]], kb), .bounded_to_theta(z[ix[2L]], xb))
      z
    }
    z0 <- if (is.finite(monotone_upto)) .aggregate_theta_to_z(theta0, p, q, Ms, Ml) else theta0
    z0[ix] <- c(.theta_to_bounded(theta0[ix[1L]], kb), .theta_to_bounded(theta0[ix[2L]], xb))
    opt <- optim(
      par = pmin(pmax(z0, lower_vec), upper_vec),
      fn = function(z) do.call(.aggregate_nll, c(list(to_theta(z)), nll_args)),
      method = "L-BFGS-B",
      lower = lower_vec,
      upper = upper_vec,
      control = control,
      hessian = FALSE
    )
    opt$par <- to_theta(opt$par)
  } else {
    opt <- optim(
      par = theta0,
      fn = .aggregate_penalized_nll,
      prep = nll_args$prep,
      family = nll_args$family,
      p = p,
      q = q,
      cpegpd.h = nll_args$cpegpd.h,
      kappa_bounds = nll_args$kappa_bounds,
      xi_bounds = nll_args$xi_bounds,
      Ms = Ms,
      Ml = Ml,
      method = method,
      control = control,
      hessian = FALSE
    )
  }
  opt
}

.aggregate_near_bound <- function(kappa, xi, kb, xb) {
  min((kappa - kb[1])/diff(kb), (kb[2] - kappa)/diff(kb),
      (xi - xb[1])/diff(xb), (xb[2] - xi)/diff(xb)) < 1e-4
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
#' \eqn{\log d}. The fit maximises a composite likelihood across the
#' aggregation durations.
#'
#' The positive likelihood conditions on \eqn{A_d>0}, using the actual
#' zero atom \eqn{p_0=\exp[-\lambda_d(1-f_0)]} of the severity/grid law.
#' Zero aggregates do not enter the composite likelihood (Ailliot et al.,
#' Eq. (4.1)), so \eqn{\exp(-\lambda_d)} is a model quantity and need not
#' match the observed dry fraction.
#' For \code{family = "cpegpd"} (Ailliot et al. 2026, Appendix F,
#' Eqs. (F.2)--(F.3)) the aggregates are modelled as gauge-discretised
#' compound sums \eqn{\sum_i h\lceil X_i/h\rceil} with \eqn{h} =
#' \code{cpegpd.h} (see \code{\link{dcpegpd}}): positive severities are
#' discretized onto \eqn{\{h,2h,\ldots,Kh\}} with \eqn{f_0=0}, preserving the
#' exact zero atom \eqn{p_0=\exp(-\lambda_d)}, and positive observations
#' \eqn{x>0} map to bin \eqn{k=\lceil x/h\rceil}. This is exact for data
#' recorded on the \eqn{h} lattice (e.g. 0.2 mm tipping-bucket gauges, or
#' \code{rcpegpd(..., h = cpegpd.h)} draws). A warning is issued when the
#' positive aggregates are off the lattice and \code{cpegpd.h} exceeds 1\% of
#' their median at the finest duration, because the upper-bin law then
#' overstates continuous sums by about \eqn{h/2} per event.
#' Quantiles invert the distribution conditional on a positive sum, using its
#' survivor for upper probabilities. \code{K} limits the search, not the law:
#' probability beyond the grid is retained. The endpoints are the first
#' positive grid point and infinity. Insufficient grids fail explicitly.
#' See the original-author preprint arXiv:2601.08350v1, Definition 2.5
#' and Appendix F (F.2)--(F.3); the displayed (F.1) omits a zero-atom
#' subtraction and is not used as a normalization identity.
#'
#' A finest-scale series is aggregated with exact window sums: on the
#' \code{cpegpd.h} lattice the sums are formed in integer units, so every
#' window is scored at its own bin.
#'
#' \strong{Monotonicity.} Ailliot et al. (Sec. 3, Prop. 3.1) use non-decreasing
#' \eqn{\sigma_d} and \eqn{\lambda_d} to order return levels for positive jumps.
#' Gauge rounding preserves that ordering for \code{"cpegpd"}. For \code{"cpdegpd"}
#' the floored jumps themselves can be zero: the \emph{unconditional} aggregate
#' is stochastically ordered, but its quantiles \emph{conditional on a positive
#' total} need not be. By default (\code{monotone_upto = Inf})
#' all slope coefficients \eqn{s_1,\ldots,s_p} and \eqn{l_1,\ldots,l_q} are
#' constrained to be non-negative: a sufficient condition that keeps the
#' parameter curves non-decreasing for every \eqn{d \ge 1}.
#' A finite \code{monotone_upto = D} imposes the paper's
#' condition only on \eqn{[1, D]}, through non-negative Bernstein
#' coefficients of \eqn{d\log\sigma_d/d\log d} and
#' \eqn{d\log\lambda_d/d\log d} on \eqn{[0, \log D]} (necessary and
#' sufficient for \eqn{p, q \le 2}); this admits concave increasing curves,
#' and \code{predict()} warns when it extrapolates beyond \eqn{D}. L-BFGS-B
#' runs (including both final stages of the default \code{"auto"}) impose
#' the constraint as box bounds; other methods use a smooth quadratic penalty
#' and their result is projected onto the constraint.
#'
#' \strong{Starting values.} Without \code{start}, single-duration
#' \code{distfit()} fits are screened (shape parameters inside the bounds,
#' scale on the data scale) and replaced by moment starts where they fail;
#' this candidate and moment starts for \eqn{\kappa \in \{0.5, 1, 2\}} are
#' compared by composite likelihood, and starts are kept away from the
#' flat ends of the logit-bounded \eqn{\kappa} and \eqn{\xi} scales.
#'
#' \strong{Optimisation.} The positive-part composite likelihood can be very
#' flat along combinations of \eqn{\kappa}, \eqn{\sigma_d} and
#' \eqn{\lambda_d} and can have more than one local optimum, where a single
#' quasi-Newton run may stop early. The default \code{optim.method = "auto"}
#' therefore runs L-BFGS-B and, from the same start, BFGS followed by an
#' L-BFGS-B polish, and keeps the higher composite likelihood; both runs are
#' reported in \code{optim.runs}. The box optimizer differences the two shapes
#' on their natural scales, including inward steps at their bounds. Other
#' requested methods use the logit shape coordinates; a converged fit within
#' \eqn{10^{-4}} of either end of a shape interval (as a fraction of its width)
#' is checked with a box-constrained polish. Nelder-Mead is not recommended: on these
#' problems it often stops short of the optimum, even when it reports
#' convergence. A fit whose (selected) run did not converge gives a
#' warning with the \code{optim} convergence code and message, and the code
#' is stored in \code{convergence}.
#'
#' \strong{Uncertainty.} The fit returns point estimates only; use
#' \code{\link{bootstrap_aggregate_egpd}} for the moving-block bootstrap of
#' Ailliot et al. (Sec. 4.1).
#'
#' @param data either a numeric vector of non-negative observations at the
#'   finest time scale, or a list of non-negative aggregated samples for each
#'   duration.
#' @param durations positive integer aggregation durations. Unnamed sample
#'   lists are paired in input order before sorting; named lists match by duration.
#'   In \code{predict()}, the durations at which to predict (default: the
#'   fitted durations).
#' @param family character string, either \code{"cpegpd"} for continuous
#'   (gauge-resolution) aggregates or \code{"cpdegpd"} for discrete aggregates.
#' @param p degree of the polynomial used for \eqn{\log \sigma_d}. Both
#'   degrees must be less than the number of distinct fitted durations.
#' @param q degree of the polynomial used for \eqn{\log \lambda_d}.
#' @param overlap logical; if \code{TRUE}, rolling sums are used when
#'   \code{data} is a finest-scale series.
#' @param cpegpd.h gauge resolution \eqn{h} of the continuous compound family
#'   (discretization step of the Panjer recursion).
#' @param start optional named list with entries \code{sigma_coef},
#'   \code{lambda_coef}, \code{kappa}, and \code{xi}; \code{kappa} and
#'   \code{xi} must lie strictly inside their bounds and the slopes must
#'   satisfy the monotonicity constraint.
#' @param kappa_bounds numeric vector of length two giving lower and upper
#'   bounds for \eqn{\kappa}.
#' @param xi_bounds numeric vector of length two giving lower and upper bounds
#'   for \eqn{\xi}.
#' @param optim.method \code{"auto"} (default; see Details) or an
#'   \code{\link{optim}} method. \code{"L-BFGS-B"} imposes the
#'   monotonicity constraint as box bounds, other methods through a penalty.
#' @param control list of \code{\link{optim}} control settings, used for every
#'   run. Defaults: \code{maxit = 500} and, for \code{"auto"} and the
#'   gradient-based methods, finite-difference steps \code{ndeps = 1e-7}.
#'   \code{maxit = 0} evaluates the likelihood at the start without
#'   optimising.
#' @param keep_data logical; should the aggregated samples (and a finest-scale
#'   input series, for \code{\link{bootstrap_aggregate_egpd}}) be stored in
#'   the fitted object?
#' @param monotone_upto \code{Inf} (default) or a duration \eqn{D \ge}
#'   \code{max(durations)} up to which \eqn{\sigma_d} and \eqn{\lambda_d} are
#'   constrained to be non-decreasing (see Details).
#'
#' @return \code{fit_aggregate_egpd()} returns an object of class
#'   \code{"aggregate_egpd"}: a list whose components include the polynomial
#'   coefficients \code{sigma_coef} (\eqn{s_0, \ldots, s_p} of
#'   \eqn{\log \sigma_d}) and \code{lambda_coef} (\eqn{l_0, \ldots, l_q} of
#'   \eqn{\log \lambda_d}) in \eqn{\log d}, the shared shape parameters
#'   \code{kappa} and \code{xi}, \code{fitted} (the duration-specific
#'   parameters at the fitted \code{durations}), the maximised composite
#'   log-likelihood \code{composite_loglik}, \code{convergence},
#'   \code{message} and \code{optim} (the \code{\link{optim}} result of the
#'   selected run), \code{optim.runs} (for \code{"auto"}, a data frame with
#'   the composite log-likelihood and convergence code of each run and the
#'   one selected; otherwise \code{NULL}), the per-duration sample
#'   sizes \code{n_total} and \code{n_positive}, and, with
#'   \code{keep_data = TRUE}, the aggregated samples \code{data} and a
#'   finest-scale input \code{series}.
#'
#'   \code{coef()} returns the named vector of \code{sigma_coef},
#'   \code{lambda_coef}, \code{kappa} and \code{xi}; \code{summary()} returns
#'   an object of class \code{"summary.aggregate_egpd"}.
#'
#'   \code{predict(type = "parameter")} returns a data frame with columns
#'   \code{duration}, \code{sigma}, \code{lambda}, \code{kappa} and \code{xi}.
#'
#'   \code{predict(type = "quantile")} returns a matrix of positive-part
#'   quantiles, with rows corresponding to durations and columns to
#'   probabilities.
#'
#' @references
#' Ailliot, P., Gaetan, C. and Naveau, P. (2026). A parsimonious tail
#' compliant multiscale statistical model for aggregated rainfall.
#' \emph{Advances in Water Resources}, 208, 105216.
#' \doi{10.1016/j.advwatres.2026.105216}. Section, equation and appendix
#' numbers on this page refer to the original-author preprint
#' arXiv:2601.08350v1.
#'
#' @seealso \code{\link{bootstrap_aggregate_egpd}}, \code{\link{dcpegpd}},
#'   \code{\link{dcpdegpd}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' durs <- c(1, 2, 4, 8)
#' ## gauge totals on the 0.2 lattice, as modelled by family = "cpegpd"
#' agg <- setNames(lapply(durs, function(d) {
#'   rcpegpd(200, lambda = exp(log(0.8) + 0.25 * log(d)),
#'           sigma = exp(log(0.6) + 0.12 * log(d)),
#'           xi = 0.2, kappa = 0.35, type = 1, h = 0.2)
#' }), durs)
#'
#' fit <- fit_aggregate_egpd(
#'   agg,
#'   durations = durs,
#'   family = "cpegpd",
#'   p = 1,
#'   q = 1
#' )
#'
#' summary(fit)
#' fit$optim.runs
#' predict(fit, type = "parameter")
#' predict(fit, type = "quantile", prob = c(0.5, 0.9))
#' plot(fit)
#' }
#' @export
fit_aggregate_egpd <- function(data, durations,
                               family = c("cpegpd", "cpdegpd"),
                               p = 2L, q = 2L,
                               overlap = TRUE,
                               cpegpd.h = 0.2,
                               start = NULL,
                               kappa_bounds = c(0.05, 5),
                               xi_bounds = c(0.01, 1.5),
                               optim.method = "auto",
                               control = list(),
                               keep_data = FALSE,
                               monotone_upto = Inf) {
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
  if (!is.character(optim.method) || length(optim.method) != 1L ||
      is.na(optim.method))
    stop("'optim.method' must be \"auto\" or a single optim() method name",
         call. = FALSE)
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

  if (family == "cpegpd") {
    if (!is.numeric(cpegpd.h) || length(cpegpd.h) != 1L ||
        !is.finite(cpegpd.h) || cpegpd.h <= 0)
      stop("'cpegpd.h' must be positive and finite", call. = FALSE)
  }
  prep <- .prepare_aggregate_input(
    data, durations, overlap = overlap,
    grid = if (family == "cpegpd") cpegpd.h else NULL
  )
  if (length(prep$durations) < 2L) {
    stop("At least two durations are required to fit the multiscale framework",
         call. = FALSE)
  }
  if (max(p_in, q_in) >= length(prep$durations)) {
    stop(sprintf(paste0("Cannot identify the duration curves from %d distinct durations: ",
                        "both polynomial degrees ('p' and 'q') must be at most %d."),
                 length(prep$durations), length(prep$durations) - 1L), call. = FALSE)
  }
  p <- as.integer(p_in)
  q <- as.integer(q_in)
  if (family == "cpdegpd") {
    discrete_ok <- all(vapply(prep$full_samples, function(x) all(x == floor(x)), logical(1)))
    if (!discrete_ok) {
      stop("cpdegpd requires integer-valued aggregated samples", call. = FALSE)
    }
  } else {
    .aggregate_offgrid_warning(prep, cpegpd.h)
  }
  if (!is.numeric(monotone_upto) || length(monotone_upto) != 1L ||
      is.na(monotone_upto) || monotone_upto < max(prep$durations)) {
    stop("'monotone_upto' must be Inf or a duration >= max(durations)",
         call. = FALSE)
  }
  Ms <- .aggregate_monotone_matrix(p, monotone_upto)
  Ml <- .aggregate_monotone_matrix(q, monotone_upto)
  nll_args <- list(prep = prep, family = family, p = p, q = q,
                   cpegpd.h = cpegpd.h, kappa_bounds = kappa_bounds,
                   xi_bounds = xi_bounds)

  if (is.null(start)) {
    candidates <- .aggregate_start_candidates(prep, family = family, p = p,
                                              q = q, cpegpd.h = cpegpd.h,
                                              kappa_bounds = kappa_bounds,
                                              xi_bounds = xi_bounds)
    cand_nll <- vapply(candidates, function(cand) {
      do.call(.aggregate_nll, c(list(.aggregate_pack_theta(
        cand, kappa_bounds = kappa_bounds, xi_bounds = xi_bounds)), nll_args))
    }, numeric(1))
    start_vals <- candidates[[which.min(cand_nll)]]
  } else {
    start_vals <- .normalize_aggregate_start(start, p = p, q = q,
                                             kappa_bounds = kappa_bounds,
                                             xi_bounds = xi_bounds,
                                             Ms = Ms, Ml = Ml)
  }

  theta0 <- .aggregate_pack_theta(start_vals, kappa_bounds = kappa_bounds,
                                  xi_bounds = xi_bounds)
  default_ctrl <- list(maxit = 500L)
  if (optim.method %in% c("auto", "BFGS", "CG", "L-BFGS-B") && is.null(control$ndeps)) {
    default_ctrl$ndeps <- rep(1e-7, length(theta0))
  }
  opt_control <- modifyList(default_ctrl, control)

  zero_iter <- !is.null(opt_control$maxit) &&
    length(opt_control$maxit) == 1L &&
    !is.na(opt_control$maxit) &&
    as.integer(opt_control$maxit) == 0L

  run <- function(theta, method)
    .aggregate_optim_run(theta, method, opt_control, nll_args, Ms, Ml,
                         monotone_upto)
  project <- function(theta)
    .aggregate_project_theta(theta, p, q, Ms, Ml, kappa_bounds, xi_bounds)$theta
  optim_runs <- NULL
  run_label <- optim.method

  if (zero_iter) {
    val0 <- do.call(.aggregate_nll, c(list(theta0), nll_args))
    opt <- list(
      par = theta0,
      value = val0,
      counts = c(`function` = 1L, gradient = NA_integer_),
      convergence = 0L,
      message = NULL
    )
  } else if (optim.method == "auto") {
    ## The positive-part composite likelihood has flat ridges and can have
    ## several local optima in (kappa, sigma_d, lambda_d), where one
    ## quasi-Newton run may stop early: run L-BFGS-B and BFGS (then an
    ## L-BFGS-B polish) from the same start and keep the better optimum.
    r_box <- run(theta0, "L-BFGS-B")
    r_bfgs <- run(theta0, "BFGS")
    r_pol <- if (is.finite(r_bfgs$value) && r_bfgs$value < 1e20) {
      run(project(r_bfgs$par), "L-BFGS-B")
    } else {
      r_bfgs
    }
    runs <- list(r_box, r_pol)
    nll_runs <- vapply(runs, function(r) {
      if (!is.finite(r$value) || r$value >= 1e20) return(Inf)
      do.call(.aggregate_nll, c(list(project(r$par)), nll_args))
    }, numeric(1))
    best <- which.min(nll_runs)
    opt <- runs[[best]]
    labels <- c("L-BFGS-B", "BFGS + L-BFGS-B")
    run_label <- labels[best]
    optim_runs <- data.frame(
      run = labels,
      composite_loglik = -nll_runs,
      convergence = c(r_box$convergence, r_pol$convergence),
      selected = seq_along(runs) == best
    )
  } else {
    opt <- run(theta0, optim.method)
    if (!identical(optim.method, "L-BFGS-B") && isTRUE(opt$convergence == 0L)) {
      pars <- .aggregate_unpack_theta(opt$par, p, q, kappa_bounds, xi_bounds)
      if (.aggregate_near_bound(pars$kappa, pars$xi, kappa_bounds, xi_bounds)) {
        ## Verify a putative logit-boundary optimum using inward differences.
        refined <- run(project(opt$par), "L-BFGS-B")
        val <- do.call(.aggregate_nll, c(list(project(opt$par)), nll_args))
        if (refined$value <= val) {
          opt <- refined
          run_label <- paste(optim.method, "+ L-BFGS-B")
        }
      }
    }
  }

  if (!is.finite(opt$value) || opt$value >= 1e20)
    stop("No finite positive-part likelihood; check aggregate support and parameters.",
         call. = FALSE)
  if (!zero_iter && opt$convergence != 0L) {
    msg <- if (is.character(opt$message) && length(opt$message) == 1L &&
               nzchar(opt$message)) paste0(": ", opt$message) else ""
    warning(sprintf(paste0(
      "fit_aggregate_egpd() did not converge (%s convergence code %d%s); ",
      "increase control$maxit or try another optim.method."),
      run_label, as.integer(opt$convergence), msg), call. = FALSE)
  }

  proj <- .aggregate_project_theta(opt$par, p, q, Ms, Ml, kappa_bounds, xi_bounds)
  pars <- proj$pars
  theta_final <- proj$theta
  nll_final <- do.call(.aggregate_nll, c(list(theta_final), nll_args))
  if (!is.finite(nll_final) || nll_final >= 1e20)
    stop("No finite positive-part likelihood; check aggregate support and parameters.",
         call. = FALSE)

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
    monotone_upto = monotone_upto,
    optim.method = optim.method,
    control = control,
    composite_loglik = -nll_final,
    convergence = opt$convergence,
    message = opt$message,
    optim = opt,
    optim.runs = optim_runs,
    call = cl,
    data_source = prep$source,
    data = if (isTRUE(keep_data)) prep$full_samples else NULL,
    series = if (isTRUE(keep_data) && prep$source == "series")
      as.numeric(data) else NULL
  ), class = "aggregate_egpd")

  object$fitted <- .aggregate_parameter_frame(object)
  object
}

#' @rdname fit_aggregate_egpd
#' @param object an object of class \code{"aggregate_egpd"}.
#' @param ... further graphical arguments passed to \code{\link{plot}} by the
#'   \code{plot} method; ignored by the other methods.
#' @export
coef.aggregate_egpd <- function(object, ...) {
  c(object$sigma_coef, object$lambda_coef,
    kappa = object$kappa, xi = object$xi)
}

#' @rdname fit_aggregate_egpd
#' @param type either \code{"parameter"} or \code{"quantile"}.
#' @param prob probabilities used when \code{type = "quantile"}.
#' @param K number of Panjer grid bins used when \code{type = "quantile"}
#'   (default \code{5000L}).
#' @export
predict.aggregate_egpd <- function(object, durations = object$durations,
                                   type = c("parameter", "quantile"),
                                   prob = NULL, K = 5000L, ...) {
  type <- match.arg(type)
  if (!is.numeric(durations) || length(durations) == 0L ||
      any(!is.finite(durations)) ||
      any(durations != floor(durations) | durations < 1 |
          durations > .Machine$integer.max)) {
    stop("'durations' must be positive integers", call. = FALSE)
  }
  durations <- as.integer(durations)
  upto <- if (is.null(object$monotone_upto)) Inf else object$monotone_upto
  if (any(durations > upto)) {
    warning("Durations beyond monotone_upto = ", upto, ": sigma_d and lambda_d ",
            "are only constrained to be non-decreasing up to that duration, so ",
            "extrapolated return levels may cross.", call. = FALSE)
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
  cpegpd_h <- if (object$family == "cpegpd") object$cpegpd.h else 1
  K_grid <- .compound_grid_size(if (is.null(K)) 5000L else K, h = cpegpd_h)

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
      cpegpd.h = if (object$family == "cpegpd") object$cpegpd.h else 0.2,
      K = K_grid
    )
  }
  out
}

#' @rdname fit_aggregate_egpd
#' @param x an object of class \code{"aggregate_egpd"} or
#'   \code{"summary.aggregate_egpd"}.
#' @export
print.aggregate_egpd <- function(x, ...) {
  cat("Multiscale aggregated ", if (x$family == "cpegpd") "EGPD" else "discrete EGPD",
      " fit\n", sep = "")
  cat("Durations:", paste(x$durations, collapse = ", "), "\n")
  cat("Composite log-likelihood:", round(x$composite_loglik, 2), "\n")
  cat("Convergence:", if (x$convergence == 0) "successful" else paste("code", x$convergence), "\n")
  invisible(x)
}

#' @rdname fit_aggregate_egpd
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

#' @rdname fit_aggregate_egpd
#' @param digits number of digits used when printing summaries.
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

#' @rdname fit_aggregate_egpd
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

## Moving-block bootstrap indices (Kunsch 1989): ceiling(n / b) blocks of b
## consecutive time steps with uniform starts in 1..(n - b + 1), truncated to n.
.aggregate_block_indices <- function(n, block_length) {
  starts <- sample.int(n - block_length + 1L, ceiling(n / block_length),
                       replace = TRUE)
  idx <- as.vector(outer(seq_len(block_length) - 1L, starts, "+"))
  idx[seq_len(n)]
}

#' Moving-block bootstrap for the multiscale aggregated EGPD
#'
#' Uncertainty for \code{\link{fit_aggregate_egpd}} by the block bootstrap
#' used by Ailliot, Gaetan and Naveau (2026, Sec. 4.1): the finest-scale
#' series is resampled in moving blocks of \code{block_length} consecutive
#' time steps (Kunsch 1989), aggregated again at every duration, and the
#' composite likelihood is re-maximised with the settings of the original
#' fit, starting from its estimate. If the original or refitted shape is
#' near a bound, the resample is also fitted from its own automatic start.
#' A converged fit is preferred and, among equally converged fits, the
#' higher composite likelihood is used. Blocks preserve the serial dependence
#' that the working-independence composite likelihood ignores, so the
#' percentile intervals account for it.
#'
#' @param object an \code{"aggregate_egpd"} fit to a finest-scale series.
#' @param data the finest-scale numeric series used for \code{object}. May
#'   be omitted when the fit was made with \code{keep_data = TRUE}.
#' @param nboot number of bootstrap replicates (Ailliot et al. used 500).
#' @param block_length block length in finest-scale time steps, at least the
#'   longest duration and at most the series length. Choose it to cover the
#'   dependence range (Ailliot et al. used two weeks of 6-minute steps).
#' @param level confidence level of the percentile intervals.
#' @param seed optional integer seed for the block draws.
#' @param mc.cores number of parallel workers for the refits (forced to 1 on
#'   Windows and capped at 2 when \code{_R_CHECK_LIMIT_CORES_} is set, as by
#'   \code{R CMD check --as-cran}). Results do not depend on it.
#'
#' @return An object of class \code{"aggregate_egpd_boot"} with components
#'   \code{t0} (original coefficients), \code{t} (\code{nboot} by
#'   coefficient matrix; failed refits are \code{NA}), \code{se} (bootstrap
#'   standard errors), \code{ci} (percentile intervals), \code{parameters}
#'   (\eqn{\sigma_d} and \eqn{\lambda_d} with percentile intervals at the
#'   fitted durations), \code{level}, \code{block_length}, \code{nboot},
#'   \code{n_failed} and \code{n_nonconverged}. Failed refits and refits that
#'   did not converge are reported with warnings.
#'
#' @references
#' Kunsch, H. R. (1989). The jackknife and the bootstrap for general
#' stationary observations. \emph{The Annals of Statistics}, 17(3), 1217--1241.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- rcpdegpd(300, lambda = 0.6, sigma = 1.5, xi = 0.15, kappa = 0.8)
#' fit <- fit_aggregate_egpd(x, durations = c(1, 2), family = "cpdegpd",
#'                           p = 1, q = 0, keep_data = TRUE)
#' boot <- bootstrap_aggregate_egpd(fit, nboot = 5, block_length = 10, seed = 1)
#' boot
#' }
#' @export
bootstrap_aggregate_egpd <- function(object, data = NULL, nboot = 200L,
                                     block_length, level = 0.95,
                                     seed = NULL, mc.cores = 1L) {
  if (!inherits(object, "aggregate_egpd"))
    stop("'object' must be an \"aggregate_egpd\" fit", call. = FALSE)
  if (!identical(object$data_source, "series"))
    stop("The block bootstrap resamples the finest-scale series; refit with a ",
         "numeric series rather than a list of aggregates.", call. = FALSE)
  if (is.null(data)) data <- object$series
  if (is.null(data))
    stop("Supply the finest-scale series in 'data' (or fit with keep_data = TRUE).",
         call. = FALSE)
  if (!is.numeric(data) || is.matrix(data))
    stop("'data' must be the numeric finest-scale series", call. = FALSE)
  n <- length(data)
  if (!is.numeric(nboot) || length(nboot) != 1L || !is.finite(nboot) ||
      nboot < 2 || nboot != floor(nboot))
    stop("'nboot' must be an integer >= 2", call. = FALSE)
  if (missing(block_length) || !is.numeric(block_length) ||
      length(block_length) != 1L || !is.finite(block_length) ||
      block_length != floor(block_length) ||
      block_length < max(object$durations) || block_length > n)
    stop("'block_length' must be an integer between max(durations) and length(data)",
         call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1)
    stop("'level' must lie in (0, 1)", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  nboot <- as.integer(nboot)
  block_length <- as.integer(block_length)

  indices <- lapply(seq_len(nboot), function(i)
    .aggregate_block_indices(n, block_length))
  start <- list(
    sigma_coef = unname(object$sigma_coef),
    lambda_coef = unname(object$lambda_coef),
    kappa = .aggregate_interior(object$kappa, object$kappa_bounds, 1e-6),
    xi = .aggregate_interior(object$xi, object$xi_bounds, 1e-6)
  )
  near <- function(f) .aggregate_near_bound(f$kappa, f$xi, object$kappa_bounds, object$xi_bounds)
  refit <- function(idx) {
    one <- function(st) tryCatch(
      withCallingHandlers(
        fit_aggregate_egpd(
          data[idx], durations = object$durations, family = object$family,
          p = object$p, q = object$q, overlap = object$overlap,
          cpegpd.h = if (is.null(object$cpegpd.h)) 0.2 else object$cpegpd.h,
          start = st, kappa_bounds = object$kappa_bounds,
          xi_bounds = object$xi_bounds,
          optim.method = if (is.null(object$optim.method)) "auto" else object$optim.method,
          control = if (is.null(object$control)) list() else object$control,
          monotone_upto = if (is.null(object$monotone_upto)) Inf else object$monotone_upto
        ),
        warning = function(w) invokeRestart("muffleWarning")
      ),
      error = function(e) NULL
    )
    fit <- one(start)
    if (near(object) || (!is.null(fit) && near(fit))) {
      ## A resample is a different objective. A boundary solution (and
      ## its nuisance coefficients) need not be a good starting basin.
      other <- one(NULL)
      if (!is.null(other) && (is.null(fit) ||
          (other$convergence == 0L && fit$convergence != 0L) ||
          ((other$convergence == 0L) == (fit$convergence == 0L) &&
             other$composite_loglik > fit$composite_loglik))) fit <- other
    }
    if (is.null(fit)) return(NULL)
    list(coef = coef(fit), nonconverged = isTRUE(fit$convergence != 0L))
  }
  mc.cores <- .sanitize_mc_cores(mc.cores)
  reps <- if (mc.cores > 1L) {
    parallel::mclapply(indices, refit, mc.cores = mc.cores)
  } else {
    lapply(indices, refit)
  }

  t0 <- coef(object)
  tmat <- matrix(NA_real_, nboot, length(t0), dimnames = list(NULL, names(t0)))
  failed <- vapply(reps, function(r) is.null(r) || inherits(r, "try-error"), logical(1))
  for (i in which(!failed)) tmat[i, ] <- reps[[i]]$coef
  n_nonconv <- sum(vapply(reps[!failed], function(r) r$nonconverged, logical(1)))
  if (all(failed))
    stop("All bootstrap refits failed.", call. = FALSE)
  if (any(failed))
    warning(sum(failed), " of ", nboot, " bootstrap refits failed and are omitted.",
            call. = FALSE)
  if (n_nonconv > 0L)
    warning(n_nonconv, " of ", sum(!failed), " bootstrap refits did not converge ",
            "(kept; see 'n_nonconverged').", call. = FALSE)

  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  pct <- paste0(format(100 * probs, trim = TRUE, scientific = FALSE), "%")
  ci <- t(apply(tmat, 2L, stats::quantile, probs = probs, na.rm = TRUE,
                names = FALSE))
  dimnames(ci) <- list(names(t0), pct)

  dur <- object$durations
  sig <- apply(tmat, 1L, function(b) exp(.aggregate_eval_poly(b[seq_along(object$sigma_coef)], dur)))
  lam <- apply(tmat, 1L, function(b) exp(.aggregate_eval_poly(
    b[length(object$sigma_coef) + seq_along(object$lambda_coef)], dur)))
  sig <- matrix(sig, nrow = length(dur))
  lam <- matrix(lam, nrow = length(dur))
  q_rows <- function(m) t(apply(m, 1L, stats::quantile, probs = probs,
                                na.rm = TRUE, names = FALSE))
  sq <- q_rows(sig)
  lq <- q_rows(lam)
  parameters <- data.frame(duration = dur,
                           sigma = object$fitted$sigma,
                           sigma_lower = sq[, 1L], sigma_upper = sq[, 2L],
                           lambda = object$fitted$lambda,
                           lambda_lower = lq[, 1L], lambda_upper = lq[, 2L])

  structure(list(
    t0 = t0,
    t = tmat,
    se = apply(tmat, 2L, stats::sd, na.rm = TRUE),
    ci = ci,
    parameters = parameters,
    level = level,
    block_length = block_length,
    nboot = nboot,
    n_failed = sum(failed),
    n_nonconverged = n_nonconv,
    call = match.call()
  ), class = "aggregate_egpd_boot")
}

#' @export
print.aggregate_egpd_boot <- function(x, digits = 4, ...) {
  cat("Moving-block bootstrap for the multiscale aggregated EGPD\n")
  cat("Replicates:", x$nboot - x$n_failed, "of", x$nboot,
      "| block length:", x$block_length,
      "| not converged:", x$n_nonconverged, "\n\n")
  tab <- cbind(Estimate = x$t0, `Boot. SE` = x$se, x$ci)
  print(round(tab, digits))
  cat("\nDuration-specific parameters (", format(100 * x$level), "% percentile intervals):\n",
      sep = "")
  print(round(x$parameters, digits), row.names = FALSE)
  invisible(x)
}
