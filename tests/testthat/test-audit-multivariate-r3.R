## Round 3 regression tests: compound/aggregate EGPD (D1, D2, D2a, L2, L3),
## neural NPE importance sampling and metadata (D3, L5, L6), blended GP (L1)
## and exact multivariate references (test gaps).

## ---- independent references ------------------------------------------------

## Skellam(rate, rate) pmf of Delta = T1 - T2 for iid Poisson(rate) generators.
skellam_pmf <- function(D, rate) besselI(2 * rate, abs(D), expon.scaled = TRUE)

## Exact bivariate homogeneous MDGPD log-likelihood (AKN Definition 2.1 with
## rectangle inclusion-exclusion over the Skellam spectral law).
mdgpd2_prep <- function(Y) {
  key <- paste(Y[, 1], Y[, 2])
  u <- !duplicated(key)
  list(U = Y[u, , drop = FALSE],
       cnt = as.vector(table(factor(key, levels = key[u]))))
}
mdgpd2_loglik <- function(sigma, xi, rate, prep) {
  Dmax <- ceiling(2 * rate + 14 * sqrt(2 * rate) + 6)
  Dl <- -Dmax:Dmax
  pD <- skellam_pmf(Dl, rate)
  keep <- pD > 1e-18
  Dl <- Dl[keep]
  pD <- pD[keep] / sum(pD[keep])
  S1 <- pmin(0, Dl)
  S2 <- pmin(0, -Dl)
  hk <- function(k) {
    if (abs(xi) < 1e-12) return(k / sigma)
    a <- 1 + xi * k / sigma
    out <- log(pmax(a, 0)) / xi
    out[a <= 0] <- if (xi > 0) -Inf else Inf
    out
  }
  Fc <- function(a, b) {
    m <- pmax(outer(-hk(a), S1, "+"), outer(-hk(b), S2, "+"))
    1 - as.vector(pmin(exp(m), 1) %*% pD)
  }
  U <- prep$U
  p <- Fc(U[, 1], U[, 2]) - Fc(U[, 1] - 1, U[, 2]) - Fc(U[, 1], U[, 2] - 1) +
    Fc(U[, 1] - 1, U[, 2] - 1)
  if (any(!(p > 0))) return(-Inf)
  sum(prep$cnt * log(p))
}

## Gauge series of the Ailliot et al. Appendix-G type: iid finest-scale
## compound Poisson sums of EGPD amounts, each rounded up to the 0.2 grid.
gauge_series <- function(n, seed, sigma = 0.8, xi = 0.2, kappa = 0.35,
                         lambda = 0.7, h = 0.2) {
  set.seed(seed)
  N <- rpois(n, lambda)
  vapply(N, function(k) {
    if (k == 0) return(0)
    sum(h * ceiling(regpd(k, sigma = sigma, xi = xi, kappa = kappa, type = 1) / h))
  }, 0)
}

## Exact window sums in integer grid units (independent of the package).
exact_window_units <- function(x, d, h = 0.2) {
  u <- round(x / h)
  as.numeric(stats::filter(u, rep(1, d), sides = 1))[d:length(u)]
}

## ---- D2: series windows and lattice lookup ---------------------------------

test_that("D2: gauge-series windows land in the right Panjer bin", {
  ## closed form for lambda = 1, GPD severity (kappa = 1), h = 0.2:
  ## P(S_h = 2h) = e^-l (l f2 + l^2 f1^2 / 2),
  ## P(S_h = 3h) = e^-l (l f3 + l^2 f1 f2 + l^3 f1^3 / 6)
  Fg <- function(x) 1 - (1 + 0.1 * x)^(-10)
  f <- diff(Fg(c(0, 0.2, 0.4, 0.6)))
  ref <- exp(-1) * c(f[2] + f[1]^2 / 2, f[3] + f[1] * f[2] + f[1]^3 / 6)

  x <- c(rep(0.2, 2000), 0.4, 0.6)
  prep <- egpd:::.prepare_aggregate_input(x, c(1, 2))
  w <- tail(prep$full_samples[[1]], 2)
  expect_equal(dcpegpd(w, lambda = 1, sigma = 1, xi = 0.1, kappa = 1), ref,
               tolerance = 1e-12)
  ## running-total errors of the old windows (0.6 + 2.3e-14) stay in bin 3
  expect_equal(dcpegpd(c(0.39999999999997726, 0.60000000000002274), lambda = 1,
                       sigma = 1, xi = 0.1, kappa = 1), ref, tolerance = 1e-12)
  expect_equal(pcpegpd(0.60000000000002274, lambda = 1, sigma = 1, xi = 0.1,
                       kappa = 1),
               pcpegpd(0.6, lambda = 1, sigma = 1, xi = 0.1, kappa = 1))

  ## random 0.2-grid series: every positive window maps to its exact bin
  set.seed(42)
  xs <- round(rexp(3000) * rbinom(3000, 1, 0.5) / 0.2) * 0.2
  durs <- c(1L, 2L, 4L, 8L)
  for (grid in list(NULL, 0.2)) {
    prep <- if (is.null(grid)) {
      egpd:::.prepare_aggregate_input(xs, durs)
    } else {
      egpd:::.prepare_aggregate_input(xs, durs, grid = grid)
    }
    for (i in seq_along(durs)) {
      units <- exact_window_units(xs, durs[i])
      w <- prep$full_samples[[i]]
      expect_equal(length(w), length(units))
      bins <- egpd:::.compound_lookup(w[w > 0], as.numeric(0:20000), h = 0.2)
      expect_identical(bins, units[units > 0])
    }
  }
  ## on the declared grid the windows are exactly k * h
  wg <- egpd:::.prepare_aggregate_input(xs, durs, grid = 0.2)$full_samples[[4]]
  expect_identical(wg, exact_window_units(xs, 8L) * 0.2)
})

test_that("D2: composite log-likelihood is the same from series and list input", {
  xg <- gauge_series(1500, seed = 11)
  durs <- c(1L, 2L, 4L, 8L)
  lst <- setNames(lapply(durs, function(d) exact_window_units(xg, d) * 0.2), durs)
  st <- list(sigma_coef = c(log(0.8), 0), lambda_coef = c(log(0.7), 1),
             kappa = 0.35, xi = 0.2)
  f_series <- fit_aggregate_egpd(xg, durations = durs, family = "cpegpd",
                                 p = 1, q = 1, start = st,
                                 control = list(maxit = 0))
  f_list <- fit_aggregate_egpd(lst, durations = durs, family = "cpegpd",
                               p = 1, q = 1, start = st,
                               control = list(maxit = 0))
  expect_equal(f_series$composite_loglik, f_list$composite_loglik,
               tolerance = 1e-12)
  ## non-overlapping blocks use the same exact arithmetic
  prep <- egpd:::.prepare_aggregate_input(xg, durs, overlap = FALSE, grid = 0.2)
  blocks <- matrix(round(xg / 0.2)[seq_len(4 * floor(1500 / 4))], ncol = 4,
                   byrow = TRUE)
  expect_identical(prep$full_samples[[3]], rowSums(blocks) * 0.2)
})

## ---- D2a: automatic starts --------------------------------------------------

test_that("D2a: absurd per-duration fits do not seed the composite fit", {
  xg <- gauge_series(600, seed = 5)
  prep <- egpd:::.prepare_aggregate_input(xg, c(1L, 2L, 4L))
  local_mocked_bindings(
    distfit = function(x, ...) {
      list(estimate = c(sigma = 2.6e-7, xi = 0.6, kappa = 1e9, lambda = 0.7))
    },
    .package = "egpd"
  )
  st <- egpd:::.aggregate_auto_start(prep, "cpegpd", 1L, 1L)
  expect_gt(st$kappa, 0.05)
  expect_lt(st$kappa, 5)
  theta <- egpd:::.aggregate_pack_theta(st, c(0.05, 5), c(0.01, 1.5))
  ## off the flat logit ends (|qlogis| ~ 18.4 at the old clamp)
  expect_true(all(abs(theta[5:6]) < 5))
  ## moment fallback: sigma on the scale of the data, not 1e-7
  expect_gt(exp(st$sigma_coef[1]), 0.05)
  expect_lt(exp(st$sigma_coef[1]), 20)
})

test_that("D2a: default fit on series input matches list input and the truth-start fit", {
  skip_on_cran()
  xg <- gauge_series(3000, seed = 11)
  durs <- c(1L, 2L, 4L, 8L)
  lst <- setNames(lapply(durs, function(d) exact_window_units(xg, d) * 0.2), durs)
  f_s <- fit_aggregate_egpd(xg, durations = durs, family = "cpegpd", p = 1, q = 1)
  f_l <- fit_aggregate_egpd(lst, durations = durs, family = "cpegpd", p = 1, q = 1)
  expect_equal(coef(f_s), coef(f_l), tolerance = 1e-10)
  expect_equal(f_s$convergence, 0)
  f_t <- fit_aggregate_egpd(xg, durations = durs, family = "cpegpd", p = 1, q = 1,
                            start = list(sigma_coef = c(log(0.8), 0),
                                         lambda_coef = c(log(0.7), 1),
                                         kappa = 0.35, xi = 0.2),
                            optim.method = "L-BFGS-B")
  expect_lt(abs(f_t$composite_loglik - f_s$composite_loglik), 1)
  expect_lt(abs(f_s$kappa - 0.35), 0.15)
  expect_lt(abs(f_s$xi - 0.2), 0.15)
  expect_lt(abs(f_s$lambda_coef[["l1"]] - 1), 0.2)
})

test_that("D2a: user starts outside the shape bounds are rejected", {
  x <- list(`1` = c(0, 0.2, 0.4, 1), `2` = c(0.2, 0.6, 1.4, 0))
  st <- list(sigma_coef = c(0, 0), lambda_coef = c(0, 0.5), kappa = 10, xi = 0.2)
  expect_error(fit_aggregate_egpd(x, durations = 1:2, p = 1, q = 1, start = st),
               "kappa_bounds")
  st$kappa <- 1
  st$xi <- 2
  expect_error(fit_aggregate_egpd(x, durations = 1:2, p = 1, q = 1, start = st),
               "xi_bounds")
})

## ---- D1: rcpegpd simulates the d/p/q law -----------------------------------

test_that("D1: rcpegpd draws follow pcpegpd (gauge-discretised sums)", {
  set.seed(3)
  for (cs in list(list(lam = 10, n = 20000L), list(lam = 100, n = 4000L),
                  list(lam = 800, n = 1500L))) {
    S <- rcpegpd(cs$n, lambda = cs$lam, sigma = 1, xi = 0.1, kappa = 1)
    expect_true(all(abs(S / 0.2 - round(S / 0.2)) < 1e-9))
    grid <- sort(unique(S))
    Fhat <- pcpegpd(grid, lambda = cs$lam, sigma = 1, xi = 0.1, kappa = 1,
                    K = ceiling(max(S) / 0.2) + 10)
    expect_lt(max(abs(stats::ecdf(S)(grid) - Fhat)), 1.63 / sqrt(cs$n))
    med <- stats::median(S)
    expect_lt(abs(pcpegpd(med, lambda = cs$lam, sigma = 1, xi = 0.1, kappa = 1,
                          K = ceiling(max(S) / 0.2) + 10) - 0.5), 0.03)
  }
})

test_that("D1: means of S_h and S match closed forms; h = 0 is the continuous sum", {
  ## E[S_h] = lambda * h * sum_{j >= 0} P(X > jh) for GPD(1, 0.1) amounts
  j <- 0:200000
  ESh <- 10 * 0.2 * sum((1 + 0.1 * j * 0.2)^(-10))
  x <- (0:2000) * 0.2
  pmf <- dcpegpd(x, lambda = 10, sigma = 1, xi = 0.1, kappa = 1, K = 2000L)
  expect_equal(sum(x * pmf), ESh, tolerance = 1e-6)

  set.seed(31)
  Sh <- rcpegpd(20000, lambda = 10, sigma = 1, xi = 0.1, kappa = 1)
  sd_S <- sqrt(10 * 2 / (0.9 * 0.8))
  expect_lt(abs(mean(Sh) - ESh), 4 * sd_S / sqrt(20000))
  S0 <- rcpegpd(20000, lambda = 10, sigma = 1, xi = 0.1, kappa = 1, h = 0)
  expect_lt(abs(mean(S0) - 10 / 0.9), 4 * sd_S / sqrt(20000))
  expect_gt(mean(abs(S0 / 0.2 - round(S0 / 0.2)) > 1e-6), 0.5)
  expect_error(rcpegpd(5, lambda = 1, sigma = 1, xi = 0.1, kappa = 1, h = -0.2),
               "'h'")
  expect_error(rcpegpd(5, lambda = 1, sigma = 1, xi = 0.1, kappa = 1, h = NA),
               "'h'")
})

test_that("D1: rcpdegpd draws follow pcpdegpd", {
  set.seed(9)
  S <- rcpdegpd(20000, lambda = 3, sigma = 2, xi = 0.15, kappa = 0.8)
  grid <- 0:max(S)
  expect_lt(max(abs(stats::ecdf(S)(grid) -
                      pcpdegpd(grid, lambda = 3, sigma = 2, xi = 0.15, kappa = 0.8))),
            1.63 / sqrt(20000))
})

test_that("D1: off-grid continuous data with a coarse h trigger a warning", {
  set.seed(21)
  x <- rcpegpd(400, lambda = 0.7, sigma = 0.8, xi = 0.2, kappa = 0.35, h = 0)
  st <- list(sigma_coef = c(log(0.8), 0), lambda_coef = c(log(0.7), 1),
             kappa = 0.35, xi = 0.2)
  expect_warning(fit_aggregate_egpd(x, durations = c(1L, 2L), p = 1, q = 1,
                                    start = st, control = list(maxit = 0)),
                 "gauge-discretised")
  xg <- rcpegpd(400, lambda = 0.7, sigma = 0.8, xi = 0.2, kappa = 0.35)
  expect_no_warning(fit_aggregate_egpd(xg, durations = c(1L, 2L), p = 1, q = 1,
                                       start = st, control = list(maxit = 0)))
})

test_that("compound d/p/q/r: missing required parameters error; supplied NA gives NA", {
  dpq <- list(dcpegpd, pcpegpd, qcpegpd, dcpdegpd, pcpdegpd, qcpdegpd)
  for (f in dpq) {
    expect_error(f(c(0.5, 0.7), lambda = 1, sigma = 1, xi = 0.1), "`kappa' missing")
    expect_error(f(0.5, lambda = 1, kappa = 1), "`sigma', `xi' missing")
    expect_error(f(0.5, sigma = 1, xi = 0.1, kappa = 1), "`lambda' missing")
    expect_error(f(0.5, lambda = 1, sigma = 1, xi = 0.1, kappa = 1, type = 4),
                 "`delta' missing")
    expect_error(f(0.5, lambda = 1, sigma = 1, xi = 0.1, type = 6),
                 "`kappa', `delta', `prob' missing")
    expect_identical(f(c(0.5, 0.7), lambda = 1, sigma = 1, xi = 0.1, kappa = NA),
                     c(NA_real_, NA_real_))
    expect_identical(f(c(0.5, 0.7), lambda = NA, sigma = 1, xi = 0.1, kappa = 1),
                     c(NA_real_, NA_real_))
    expect_identical(f(0.5, lambda = 1, sigma = 1, xi = NaN, kappa = 1), NaN)
    ## carrier parameters the type does not use may stay NA
    expect_false(anyNA(f(0.5, lambda = 1, sigma = 1, xi = 0.1, kappa = 1,
                         delta = NA, prob = NA)))
    expect_error(f(0.5, lambda = 1, sigma = -1, xi = 0.1, kappa = 1), "'sigma'")
  }
  for (rf in list(rcpegpd, rcpdegpd)) {
    expect_error(rf(3, lambda = 1, sigma = 1, xi = 0.1), "`kappa' missing")
    expect_error(rf(3, sigma = 1, xi = 0.1, kappa = 1), "`lambda' missing")
    expect_warning(r <- rf(3, lambda = 1, sigma = 1, xi = 0.1, kappa = NA),
                   "NAs produced")
    expect_true(length(r) == 3L && all(is.na(r)))
  }
  expect_type(suppressWarnings(rcpdegpd(2, lambda = NA, sigma = 1, xi = 0.1,
                                        kappa = 1)), "integer")
})

## ---- L2: moving-block bootstrap ---------------------------------------------

test_that("L2: block indices are moving blocks of consecutive time steps", {
  set.seed(4)
  idx <- egpd:::.aggregate_block_indices(103L, 10L)
  expect_length(idx, 103L)
  starts <- idx[seq(1, 103, by = 10)]
  expect_true(all(starts >= 1 & starts <= 94))
  for (s in seq(1, 101, by = 10)) {
    run <- idx[s:min(s + 9, 103)]
    expect_identical(run, run[1] + seq_along(run) - 1L)
  }
})

test_that("L2: bootstrap SE of a block-resampled mean matches the exact MBB variance", {
  ## With n = k b, the moving-block bootstrap mean has variance
  ## Var*(B) / (k b^2), B the block sum with a uniform start over n - b + 1.
  set.seed(8)
  x <- as.numeric(arima.sim(list(ar = 0.6), 240) + 3)
  x[x < 0] <- 0
  b <- 12L
  k <- 240L / b
  starts <- seq_len(240L - b + 1L)
  B <- vapply(starts, function(s) sum(x[s:(s + b - 1L)]), 0)
  exact_sd <- sqrt(mean((B - mean(B))^2) / (k * b^2))

  fake <- structure(list(data_source = "series", series = x,
                         durations = c(1L, 2L), sigma_coef = c(s0 = 0),
                         lambda_coef = c(l0 = 0), kappa = 1, xi = 0.2,
                         kappa_bounds = c(0.05, 5), xi_bounds = c(0.01, 1.5),
                         fitted = data.frame(sigma = c(1, 1), lambda = c(1, 1))),
                    class = "aggregate_egpd")
  local_mocked_bindings(
    fit_aggregate_egpd = function(data, ...) {
      structure(list(sigma_coef = c(s0 = mean(data)), lambda_coef = c(l0 = 0),
                     kappa = 1, xi = 0.2), class = "aggregate_egpd")
    },
    .package = "egpd"
  )
  boot <- bootstrap_aggregate_egpd(fake, nboot = 4000L, block_length = b, seed = 1)
  expect_equal(unname(boot$se["s0"]), exact_sd, tolerance = 0.05)
  expect_lt(abs(mean(boot$t[, "s0"]) - mean(B) / b), 4 * exact_sd / sqrt(4000))
  expect_error(bootstrap_aggregate_egpd(fake, nboot = 10, block_length = 1),
               "block_length")
})

test_that("L2: bootstrap_aggregate_egpd refits real series reproducibly", {
  skip_on_cran()
  set.seed(1)
  x <- rcpdegpd(300, lambda = 0.8, sigma = 1.5, xi = 0.15, kappa = 0.8)
  fit <- suppressWarnings(  # convergence is not under test here
    fit_aggregate_egpd(x, durations = c(1L, 2L), family = "cpdegpd",
                       p = 1, q = 1, keep_data = TRUE,
                       control = list(maxit = 100)))
  b1 <- bootstrap_aggregate_egpd(fit, nboot = 3, block_length = 10, seed = 3)
  b2 <- bootstrap_aggregate_egpd(fit, data = x, nboot = 3, block_length = 10,
                                 seed = 3)
  expect_identical(b1$t, b2$t)
  expect_equal(dim(b1$ci), c(6L, 2L))
  expect_true(all(b1$parameters$sigma_lower <= b1$parameters$sigma_upper))
  expect_output(print(b1), "Moving-block bootstrap")
  lst_fit <- fit_aggregate_egpd(list(`1` = x, `2` = x[-1] + x[-300]),
                                durations = c(1L, 2L), family = "cpdegpd",
                                p = 1, q = 1, control = list(maxit = 0),
                                start = list(sigma_coef = c(0.4, 0),
                                             lambda_coef = c(-0.2, 0.5),
                                             kappa = 0.8, xi = 0.15))
  expect_error(bootstrap_aggregate_egpd(lst_fit, nboot = 3, block_length = 10),
               "series")
})

## ---- L3: monotonicity on a duration range ---------------------------------

test_that("L3: Bernstein constraint is exact for quadratics on [1, D]", {
  D <- 8
  M <- egpd:::.aggregate_monotone_matrix(2L, D)
  t <- seq(0, log(D), length.out = 2001)
  set.seed(6)
  for (i in 1:200) {
    s <- runif(2, -1, 1)
    deriv_ok <- all(s[1] + 2 * s[2] * t >= 0)
    expect_identical(all(M %*% s >= 0), deriv_ok)
  }
  ## cubic: non-negative Bernstein coefficients imply a non-negative derivative
  M3 <- egpd:::.aggregate_monotone_matrix(3L, D)
  for (i in 1:200) {
    s <- runif(3, -1, 1)
    if (all(M3 %*% s >= 0))
      expect_true(all(s[1] + 2 * s[2] * t + 3 * s[3] * t^2 >= -1e-12))
  }
  expect_identical(egpd:::.aggregate_monotone_matrix(2L, Inf), diag(1, 2))
})

test_that("L3: monotone_upto admits concave increasing sigma_d and warns beyond it", {
  skip_on_cran()
  durs <- c(1L, 2L, 4L, 8L)
  set.seed(12)
  ## log sigma_d = 0.5 log d - 0.12 (log d)^2: increasing on [1, 8] only
  samples <- setNames(lapply(durs, function(d) {
    rcpdegpd(250, sigma = exp(0.5 * log(d) - 0.12 * log(d)^2), xi = 0.15,
             kappa = 0.8, lambda = exp(log(0.9) + 0.3 * log(d)))
  }), durs)
  st <- list(sigma_coef = c(0, 0.2, 0), lambda_coef = c(log(0.9), 0.3),
             kappa = 0.8, xi = 0.15)
  ## Nelder-Mead (penalised path); convergence is not under test here
  f_all <- suppressWarnings(fit_aggregate_egpd(samples, durs, family = "cpdegpd",
                                               p = 2, q = 1, start = st,
                                               optim.method = "Nelder-Mead"))
  f_rng <- suppressWarnings(fit_aggregate_egpd(samples, durs, family = "cpdegpd",
                                               p = 2, q = 1, start = st,
                                               monotone_upto = 8,
                                               optim.method = "Nelder-Mead"))
  expect_gte(f_all$sigma_coef[["s2"]], 0)
  expect_lt(f_rng$sigma_coef[["s2"]], -0.05)
  expect_gt(f_rng$composite_loglik, f_all$composite_loglik)
  M <- egpd:::.aggregate_monotone_matrix(2L, 8)
  expect_true(all(M %*% f_rng$sigma_coef[-1] >= 0))
  pars <- predict(f_rng, durations = 1:8)
  expect_true(all(diff(pars$sigma) >= -1e-10))
  expect_warning(predict(f_rng, durations = 16L), "monotone_upto")
  ## L-BFGS-B imposes the same constraint as box bounds in z = M %*% slopes
  f_box <- suppressWarnings(fit_aggregate_egpd(samples, durs, family = "cpdegpd",
                                               p = 2, q = 1, start = st,
                                               monotone_upto = 8,
                                               optim.method = "L-BFGS-B",
                                               control = list(maxit = 3)))
  expect_true(all(M %*% f_box$sigma_coef[-1] >= -1e-12))
  expect_error(fit_aggregate_egpd(samples, durs, family = "cpdegpd", p = 2, q = 1,
                                  start = st, monotone_upto = 4),
               "monotone_upto")
})

## ---- Default optimiser of fit_aggregate_egpd -------------------------------

## Gain from re-optimising at the estimate with other methods: an independent
## check that the default fit is at a local optimum.
aggregate_polish_gain <- function(fit, data, durs) {
  st <- list(sigma_coef = unname(fit$sigma_coef),
             lambda_coef = unname(fit$lambda_coef),
             kappa = egpd:::.aggregate_interior(fit$kappa, fit$kappa_bounds, 1e-6),
             xi = egpd:::.aggregate_interior(fit$xi, fit$xi_bounds, 1e-6))
  ll <- vapply(c("Nelder-Mead", "BFGS"), function(m) suppressWarnings(
    fit_aggregate_egpd(data, durs, family = fit$family, p = fit$p, q = fit$q,
                       start = st, optim.method = m,
                       control = list(maxit = 3000, reltol = 1e-12))
  )$composite_loglik, 0)
  max(ll) - fit$composite_loglik
}

test_that("default optimiser converges on the documented examples", {
  skip_on_cran()
  durs <- c(1, 2, 4, 8)
  set.seed(1)  # previous help-page example (Nelder-Mead stopped at maxit)
  agg_round <- setNames(lapply(durs, function(d) {
    vapply(rpois(200, exp(log(0.8) + 0.25 * log(d))), function(n_events) {
      severity <- regpd(n_events, sigma = exp(log(0.6) + 0.12 * log(d)),
                        xi = 0.2, kappa = 0.35, type = 1)
      sum(0.2 * round(severity / 0.2))
    }, numeric(1))
  }), durs)
  set.seed(1)  # current help-page example
  agg_rcp <- setNames(lapply(durs, function(d) {
    rcpegpd(200, lambda = exp(log(0.8) + 0.25 * log(d)),
            sigma = exp(log(0.6) + 0.12 * log(d)),
            xi = 0.2, kappa = 0.35, type = 1, h = 0.2)
  }), durs)
  for (agg in list(agg_round, agg_rcp)) {
    expect_no_warning(
      fit <- fit_aggregate_egpd(agg, durations = durs, family = "cpegpd",
                                p = 1, q = 1))
    expect_equal(fit$convergence, 0)
    expect_identical(fit$optim.method, "auto")
    expect_equal(nrow(fit$optim.runs), 2L)
    expect_lt(aggregate_polish_gain(fit, agg, durs), 1e-3)
  }
})

test_that("default optimiser reaches the better of two local optima", {
  skip_on_cran()
  ## Positive-part likelihood with two local optima (kappa about 0.7 and 3.9):
  ## a single L-BFGS-B run from the automatic start stops at the lower one,
  ## Nelder-Mead (old default) short of either.
  durs <- c(1L, 2L, 4L, 8L)
  set.seed(23)
  agg <- setNames(lapply(durs, function(d)
    rcpegpd(200, lambda = 3 * d^0.5, sigma = 0.3 * d^0.1, xi = 0.15,
            kappa = 1, h = 0.2)), durs)
  fit <- fit_aggregate_egpd(agg, durations = durs, family = "cpegpd",
                            p = 1, q = 1)
  ## reference: L-BFGS-B from every automatic start candidate
  prep <- egpd:::.prepare_aggregate_input(agg, durs, grid = 0.2)
  cands <- egpd:::.aggregate_start_candidates(prep, "cpegpd", 1L, 1L, 0.2,
                                              c(0.05, 5), c(0.01, 1.5))
  ref <- max(vapply(cands, function(st) suppressWarnings(
    fit_aggregate_egpd(agg, durations = durs, family = "cpegpd", p = 1, q = 1,
                       start = st, optim.method = "L-BFGS-B")
  )$composite_loglik, 0))
  expect_gt(fit$composite_loglik, ref - 0.01)
  expect_equal(fit$convergence, 0)
})

test_that("non-convergence is reported with the method, code and message", {
  set.seed(1)
  x <- rcpdegpd(300, lambda = 0.6, sigma = 1.5, xi = 0.15, kappa = 0.8)
  expect_warning(
    f1 <- fit_aggregate_egpd(x, durations = c(1L, 2L), family = "cpdegpd",
                             p = 1, q = 1, keep_data = TRUE,
                             control = list(maxit = 1)),
    "did not converge \\((L-BFGS-B|BFGS \\+ L-BFGS-B) convergence code 1: ")
  expect_equal(f1$convergence, 1)
  expect_equal(f1$optim.runs$convergence, c(1L, 1L))
  expect_output(print(f1), "Convergence: code 1")
  expect_warning(
    fit_aggregate_egpd(x, durations = c(1L, 2L), family = "cpdegpd",
                       p = 1, q = 1, optim.method = "Nelder-Mead",
                       control = list(maxit = 5)),
    "did not converge \\(Nelder-Mead convergence code 1\\)")
  expect_error(fit_aggregate_egpd(x, durations = c(1L, 2L), family = "cpdegpd",
                                  optim.method = c("BFGS", "CG")),
               "optim.method")
  ## bootstrap refits inherit the settings and report non-convergence
  expect_warning(
    b <- bootstrap_aggregate_egpd(f1, nboot = 2, block_length = 10, seed = 1),
    "2 of 2 bootstrap refits did not converge")
  expect_identical(b$n_nonconverged, 2L)
})

## ---- L1: blended GP weight shapes -------------------------------------------

test_that("L1: bgpd requires c1, c2 > 3 in d/p/q/r", {
  fns <- list(function(...) pbgpd(0.5, ...), function(...) dbgpd(0.5, ...),
              function(...) qbgpd(0.5, ...), function(...) rbgpd(3, ...))
  for (fn in fns) {
    expect_error(fn(xi = 0.1, c1 = 3, c2 = 5), "`c1`")
    expect_error(fn(xi = 0.1, c1 = 5, c2 = 1), "`c2`")
    expect_no_error(fn(xi = 0.1, c1 = 3.5, c2 = 5))
  }
  expect_error(rbgpd(0, xi = 0.1, c1 = 2, c2 = 5), "`c1`")
})

## ---- D3: adaptive importance sampling for the NPE ---------------------------

test_that("D3: network-scale log prior matches the training prior sampler", {
  for (fam in c("mdgpd", "zimdgpd", "begpd", "bzidegpd")) {
    box <- egpd:::.neural_prior_box(fam)
    set.seed(2)
    draws <- egpd:::.sample_theta_prior(fam, 20000L)
    for (j in seq_along(box$lower)) {
      dens <- function(t) {
        th <- matrix((box$lower + box$upper) / 2, length(box$lower), length(t))
        th[j, ] <- t
        others <- egpd:::.neural_log_prior(th[, 1, drop = FALSE], box) -
          switch(box$transform[j],
                 log = th[j, 1] - log(diff(box$bounds[j, ])),
                 identity = -log(diff(box$bounds[j, ])),
                 logit = plogis(th[j, 1], log.p = TRUE) +
                   plogis(-th[j, 1], log.p = TRUE) - log(diff(box$bounds[j, ])))
        exp(egpd:::.neural_log_prior(th, box) - others)
      }
      expect_equal(integrate(dens, box$lower[j], box$upper[j])$value, 1,
                   tolerance = 1e-6)
      m1 <- integrate(function(t) t * dens(t), box$lower[j], box$upper[j])$value
      expect_lt(abs(mean(draws[j, ]) - m1), 4 * sd(draws[j, ]) / sqrt(20000))
    }
    outside <- matrix(box$upper + 0.1, ncol = 1)
    expect_identical(egpd:::.neural_log_prior(outside, box), -Inf)
  }
})

test_that("D3: importance sampling recovers an analytic Gaussian x prior posterior", {
  ## log r = N(c, S) kernel on the network scale; with the log-uniform prior
  ## (density e^theta in log-transformed coordinates, logistic for logit ones)
  ## the posterior is N(c + S g, S) for g the prior log-density gradient.
  cases <- list(
    list(fam = "mdgpd", c = c(log(2), 0.1, log(0.5)), sd = c(0.15, 0.05, 0.15)),
    list(fam = "bzidegpd", c = c(log(c(2, 1, 0.2, 5, 5, 0.2)), qlogis(0.3)),
         sd = rep(0.1, 7))
  )
  for (cs in cases) {
    box <- egpd:::.neural_prior_box(cs$fam)
    g <- ifelse(box$transform == "log", 1,
                ifelse(box$transform == "logit", 1 - 2 * plogis(cs$c), 0))
    lr <- function(theta) -0.5 * colSums(((theta - cs$c) / cs$sd)^2)
    set.seed(10)
    res <- egpd:::.npe_importance_sample(lr, cs$fam, 2L, 2000L)
    expect_gte(res$diagnostics$ess, 1000)
    expect_gt(res$diagnostics$unique, 1000)
    post_mean <- rowMeans(res$theta)
    expect_true(all(abs(post_mean - (cs$c + cs$sd^2 * g)) < 4 * cs$sd / sqrt(1000)))
    expect_equal(apply(res$theta, 1, sd), cs$sd, tolerance = 0.12)
  }
})

test_that("D3: exact MDGPD likelihood ratio at n = 1000 gives a non-degenerate posterior", {
  skip_on_cran()
  set.seed(7)
  Y <- rmdgpd(1000, sigma = 2, xi = 0.2, rate = 0.5)
  prep <- mdgpd2_prep(Y)
  path <- test_path("r3-exact-mdgpd.bson")
  file.create(path)
  on.exit(unlink(c(path, paste0(path, ".metadata.rds"))), add = TRUE)
  egpd:::.write_neural_metadata(path, "mdgpd", "npe", 2L)
  local_mocked_bindings(
    .check_julia_deps = function() NULL,
    .init_julia_begpd = function() NULL,
    .init_npe_begpd = function(...) "mock",
    .neural_loadstate = function(est, ...) est,
    .neural_logratio = function(estimator, Z, theta) {
      apply(theta, 2, function(t) mdgpd2_loglik(exp(t[1]), t[2], exp(t[3]), prep))
    },
    .package = "egpd"
  )
  set.seed(11)
  fit <- egpd:::.distfit_neuralbayes(Y, "mdgpd", model.path = path,
                                     nsamples = 1000L,
                                     npe.control = list(K0 = 5000L))
  diag <- fit$importance_sampling
  ## the 20000-point prior-grid SIR of round 2 had ESS ~ 1 and 1-3 distinct draws
  expect_gte(diag$ess, 1000)
  expect_gt(diag$unique, 600)
  post_sd <- apply(fit$posterior_samples, 1, sd)
  ## brute-force posterior on a 81 x 41 x 41 grid around the mode (independent
  ## of the sampler): mean (2.0139, 0.2021, 0.4993), sd (0.0112, 0.0047, 0.0250)
  expect_equal(unname(post_sd), c(0.0112, 0.0047, 0.0250), tolerance = 0.2)
  expect_equal(unname(rowMeans(fit$posterior_samples)), c(2.0139, 0.2021, 0.4993),
               tolerance = 0.002)
  ## xi and rate agree with the Laplace approximation (0.0045, 0.026); the
  ## likelihood has a kink in sigma, where the Laplace sd (0.004) is too small
  expect_equal(unname(post_sd[2:3]), c(0.0045, 0.026), tolerance = 0.25)
  expect_output(print(fit), "Importance sampling: ESS")
  expect_output(print(summary(fit)), "distinct of 1000 posterior draws")
})

test_that("D3: low effective sample size and few distinct draws are reported", {
  path <- test_path("r3-lowess.bson")
  file.create(path)
  on.exit(unlink(c(path, paste0(path, ".metadata.rds"))), add = TRUE)
  egpd:::.write_neural_metadata(path, "mdgpd", "npe", 2L)
  local_mocked_bindings(
    .check_julia_deps = function() NULL,
    .init_julia_begpd = function() NULL,
    .init_npe_begpd = function(...) "mock",
    .neural_loadstate = function(est, ...) est,
    .neural_logratio = function(estimator, Z, theta) {
      -0.5 * colSums(((theta - c(log(2), 0.1, log(0.5))) / 1e-6)^2)
    },
    .package = "egpd"
  )
  set.seed(1)
  w <- character()
  fit <- withCallingHandlers(
    egpd:::.distfit_neuralbayes(matrix(1:20, ncol = 2), "mdgpd",
                                model.path = path, nsamples = 500L,
                                npe.control = list(max_stages = 0L)),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("effective sample size of only 1 ", w, fixed = TRUE)))
  expect_true(any(grepl("Only 1 distinct parameter vectors among 500", w,
                        fixed = TRUE)))
  expect_lt(fit$importance_sampling$ess, 200)
})

test_that("D3: log-ratios beyond the Float32 sigmoid range are used on the log scale", {
  ## A well-trained network at large n gives log r far above ~17, where a
  ## Float32 classifier probability rounds to 1; egpd_logratio returns the
  ## logit, and the sampler only uses differences of log r.
  centre <- c(log(2), 0.1, log(0.5))
  sds <- c(0.15, 0.05, 0.15)
  path <- test_path("r3-saturating.bson")
  file.create(path)
  on.exit(unlink(c(path, paste0(path, ".metadata.rds"))), add = TRUE)
  egpd:::.write_neural_metadata(path, "mdgpd", "npe", 2L)
  lr_value <- function(theta) 500 - 0.5 * colSums(((theta - centre) / sds)^2)
  local_mocked_bindings(
    .check_julia_deps = function() NULL,
    .init_julia_begpd = function() NULL,
    .init_npe_begpd = function(...) "mock",
    .neural_loadstate = function(est, ...) est,
    .neural_logratio = function(estimator, Z, theta) lr_value(theta),
    .package = "egpd"
  )
  set.seed(6)
  fit <- egpd:::.distfit_neuralbayes(matrix(1:20, ncol = 2), "mdgpd",
                                     model.path = path, nsamples = 500L)
  expect_gte(fit$importance_sampling$ess, 500)
  ## posterior N(centre + sds^2 * (1, 0, 1), sds^2) on the network scale
  theta_draws <- rbind(log(fit$posterior_samples[1, ]), fit$posterior_samples[2, ],
                       log(fit$posterior_samples[3, ]))
  expect_true(all(abs(rowMeans(theta_draws) - (centre + sds^2 * c(1, 0, 1))) <
                    4 * sds / sqrt(500)))
  ## the pre-fix v0.1.3 path returned +Inf (saturated sigmoid): explicit error
  lr_value <- function(theta) rep(Inf, ncol(theta))
  expect_error(egpd:::.distfit_neuralbayes(matrix(1:20, ncol = 2), "mdgpd",
                                           model.path = path, nsamples = 50L),
               "finite or -Inf")
})

test_that("parallel workers are capped at 2 under the R CMD check core limit", {
  skip_on_os("windows")
  old <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else
    Sys.setenv("_R_CHECK_LIMIT_CORES_" = old), add = TRUE)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
  expect_lte(egpd:::.sanitize_mc_cores(NA), 2L)
  expect_lte(egpd:::.sanitize_mc_cores(NULL), 2L)
  expect_identical(egpd:::.sanitize_mc_cores(8L), 2L)
  expect_identical(egpd:::.sanitize_mc_cores(1L), 1L)
  ## parallel refuses more than 2 processes while the limit is set
  set.seed(1)
  expect_identical(egpd:::.neural_simulate(1:6, function(i) i, mc.cores = NA),
                   as.list(1:6))
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
  expect_identical(egpd:::.sanitize_mc_cores(8L), 8L)
  Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  expect_identical(egpd:::.sanitize_mc_cores(3L), 3L)
  for (fn in list(train_begpd, train_bdegpd, train_mdgpd))
    expect_null(formals(fn)$mc.cores)
})

## ---- L5 / L6: neural metadata, version and range checks, savepath ----------

test_that("L5: training range, data support and version mismatches warn", {
  path <- test_path("r3-meta.bson")
  file.create(path)
  on.exit(unlink(c(path, paste0(path, ".metadata.rds"))), add = TRUE)
  local_mocked_bindings(
    .neural_backend_versions = function() c(julia = "1.12.7",
      NeuralEstimators.jl = "0.1.0", Flux.jl = "0.14.25",
      NeuralEstimators.R = NA_character_, egpd = "1.0"),
    .package = "egpd"
  )
  egpd:::.write_neural_metadata(path, "mdgpd", "nbe", 2L,
                                training = list(m = c(500L, 1000L)))
  md <- readRDS(paste0(path, ".metadata.rds"))
  expect_identical(md$training$m, c(500L, 1000L))
  expect_identical(unname(md$versions["NeuralEstimators.jl"]), "0.1.0")
  expect_silent(egpd:::.validate_neural_model(path, "mdgpd", "nbe", 2L))

  local_mocked_bindings(
    .check_julia_deps = function() NULL,
    .init_julia_begpd = function() NULL,
    .init_nbe_begpd = function(...) "mock",
    .neural_loadstate = function(est, ...) est,
    .julia_let = function(...) matrix(c(log(2), 0.1, log(0.5)), ncol = 1),
    .neural_backend_versions = function() c(julia = "1.12.7",
      NeuralEstimators.jl = "0.3.0", Flux.jl = "0.14.25",
      NeuralEstimators.R = NA_character_, egpd = "1.0"),
    .package = "egpd"
  )
  x <- rbind(matrix(1:20, ncol = 2), c(-1, 0))
  w <- character()
  withCallingHandlers(
    egpd:::.distfit_neuralbayes(x, "mdgpd", model.path = path, estimator = "nbe"),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("outside the training sample-size range", w)))
  expect_true(any(grepl("1 of 11 rows", w)))
  expect_true(any(grepl("NeuralEstimators.jl 0.1.0 but 0.3.0", w)))
})

test_that("L5: support checks follow each simulator", {
  expect_warning(egpd:::.neural_check_support(cbind(c(0.5, -1), c(1, 2)), "begpd"),
                 "1 of 2 rows")
  expect_warning(egpd:::.neural_check_support(cbind(c(1.5, 1), c(1, 2)), "bdegpd"),
                 "non-negative integers")
  expect_silent(egpd:::.neural_check_support(cbind(c(0, 3), c(0, -2)), "zimdgpd"))
  expect_warning(egpd:::.neural_check_support(cbind(c(0, 3), c(0, -2)), "mdgpd"),
                 "strictly positive maximum")
  set.seed(3)
  expect_silent(egpd:::.neural_check_support(rmdgpd(200, 2, -0.2, rate = 1, d = 3),
                                             "mdgpd"))
  expect_silent(egpd:::.neural_check_support(
    rbzidegpd(200, 2, 1, 0.1, 5, 5, 0.2, pi0 = 0.3), "bzidegpd"))
})

test_that("L6: trainers default to a per-user directory, not inst/models", {
  for (fn in list(train_begpd, train_bdegpd, train_mdgpd))
    expect_null(formals(fn)$savepath)
  userdir <- file.path(tempdir(), "r3-user-data")
  old <- Sys.getenv("R_USER_DATA_DIR", unset = NA)
  Sys.setenv(R_USER_DATA_DIR = userdir)
  on.exit({
    if (is.na(old)) Sys.unsetenv("R_USER_DATA_DIR") else Sys.setenv(R_USER_DATA_DIR = old)
    unlink(userdir, recursive = TRUE)
  }, add = TRUE)
  local_mocked_bindings(
    .check_julia_deps = function() NULL,
    .init_julia_begpd = function() NULL,
    .julia_let = function(...) NULL,
    .init_npe_begpd = function(d = 6L, n = 2L) list(d = d, n = n),
    .init_nbe_begpd = function(d = 6L, n = 2L) list(d = d, n = n),
    .neural_train = function(estimator, ...) estimator,
    .neural_savestate = function(estimator, path) file.create(path),
    .package = "egpd"
  )
  expect_message(
    paths <- train_mdgpd(estimator = "npe", K = 4L, m = c(5L, 9L), epochs = 1L,
                         mc.cores = 1L, verbose = FALSE),
    "Saving trained neural models to")
  expect_true(startsWith(normalizePath(paths$npe),
                         normalizePath(tools::R_user_dir("egpd", "data"))))
  md <- readRDS(paste0(paths$npe, ".metadata.rds"))
  expect_identical(md$training$m, c(5L, 9L))
  expect_true(egpd:::.neural_contract_matches(md, "mdgpd", "npe", 2L))
  expect_error(train_mdgpd(savepath = c("a", "b"), K = 4L, m = 5L, epochs = 1L,
                           verbose = FALSE), "savepath")
})

## ---- Test gaps: exact MDGPD references (AKN Def. 2.1, Cor. 2.2, 2.3) ---------

## P(M <= k) = 1 - E[1 ^ exp(max(S - h(k)))], S = T - max(T), T iid Poisson(rate),
## evaluated exactly by summing over the (truncated) Poisson lattice.
mdgpd_exact_cdf <- function(k, sigma, xi, rate, d, Tmax = 30L) {
  sigma <- rep_len(sigma, d)
  xi <- rep_len(xi, d)
  h <- vapply(seq_len(d), function(j) {
    if (xi[j] == 0) return(k[j] / sigma[j])
    a <- 1 + xi[j] * k[j] / sigma[j]
    if (a <= 0) return(if (xi[j] < 0) Inf else -Inf)
    log(a) / xi[j]
  }, 0)
  Tg <- as.matrix(expand.grid(rep(list(0:Tmax), d)))
  w <- apply(matrix(dpois(Tg, rate), ncol = d), 1, prod)
  S <- Tg - apply(Tg, 1, max)
  m <- apply(sweep(S, 2, h, "-"), 1, max)
  1 - sum(w * pmin(1, exp(m))) / sum(w)
}

test_that("Definition 2.1 CDF matches the exact Poisson-generator expectation (d = 2)", {
  ## d = 3 is covered with the same exact expectation in test-audit-multivariate-r2.R
  set.seed(202)
  n <- 100000L
  y2 <- rmdgpd(n, sigma = c(3, 1.5), xi = c(-0.2, 0.25), rate = 0.8, d = 2L)
  for (k in list(c(1L, 1L), c(3L, 0L), c(5L, 4L), c(0L, 2L), c(9L, 12L))) {
    ref <- mdgpd_exact_cdf(k, c(3, 1.5), c(-0.2, 0.25), 0.8, 2L, Tmax = 25L)
    emp <- mean(y2[, 1] <= k[1] & y2[, 2] <= k[2])
    expect_lt(abs(emp - ref), 5 * sqrt(ref * (1 - ref) / n) + 1e-4)
  }
})

test_that("Corollary 2.2: marginal pmf of N1 for the standard bivariate MDGPD", {
  ## P(N1 = n) = (1 - e^-1) e^(1-n) [E(e^D 1(D < min(0, n))) + P(D >= 0) 1(n > 0)]
  set.seed(22)
  n <- 200000L
  for (rate in c(0.5, 2)) {
    N <- rmdgpd(n, sigma = 1, xi = 0, rate = rate)
    D <- -60:60
    pD <- skellam_pmf(D, rate)
    for (k in -4:6) {
      ref <- (1 - exp(-1)) * exp(1 - k) *
        (sum(exp(D) * pD * (D < min(0, k))) + sum(pD[D >= 0]) * (k > 0))
      emp <- mean(N[, 1] == k)
      expect_lt(abs(emp - ref), 5 * sqrt(ref * (1 - ref) / n) + 1e-4)
    }
  }
})

test_that("Corollary 2.3: threshold stability with exact spectral law (Eq. 10)", {
  set.seed(23)
  rate <- 0.7
  N <- rmdgpd(300000L, sigma = 1, xi = 0, rate = rate)
  D <- -50:50
  pD <- skellam_pmf(D, rate)
  S <- cbind(pmin(0, D), pmin(0, -D))
  for (m in list(c(1L, 2L), c(3L, 0L))) {
    exc <- N[, 1] > m[1] | N[, 2] > m[2]
    X <- sweep(N[exc, , drop = FALSE], 2, m)
    G <- apply(X, 1, max)
    for (j in 1:4)
      expect_lt(abs(mean(G <= j) - (1 - exp(-j))), 5 * sqrt(0.25 / length(G)))
    ## Eq. (10): P(S_m = l) = E[e^max(U) 1(U - max(U) = l)] / E[e^max(U)], U = S - m
    U <- sweep(S, 2, m)
    mU <- apply(U, 1, max)
    lab <- paste(U[, 1] - mU, U[, 2] - mU)
    ref <- tapply(pD * exp(mU), lab, sum) / sum(pD * exp(mU))
    emp_lab <- paste(X[, 1] - G, X[, 2] - G)
    top <- names(sort(ref, decreasing = TRUE))[1:5]
    for (l in top) {
      emp <- mean(emp_lab == l)
      expect_lt(abs(emp - ref[[l]]), 5 * sqrt(ref[[l]] * (1 - ref[[l]]) / length(G)))
      ## spectral vector independent of the geometric maximum
      expect_lt(abs(mean(emp_lab[G == 1] == l) - ref[[l]]),
                5 * sqrt(ref[[l]] * (1 - ref[[l]]) / sum(G == 1)))
    }
  }
})

## ---- Test gaps: BZIDEGPD joint zero from the BEGPD construction -------------

## P(floor(Y1) = 0, floor(Y2) = 0) for BEGPD by numerical integration of the
## radial/angular construction: Y = R (A, 1 - A), A = (1 - w) V + w U.
begpd_p00 <- function(kappa, sigma, xi, thL, thU, thw) {
  H <- function(x) 1 - (1 + xi * x / sigma)^(-1 / xi)
  r_of_u <- function(u) sigma / xi * ((1 - u^(1 / kappa))^(-xi) - 1)
  u1 <- H(1)^kappa
  u2 <- H(2)^kappa
  FA <- function(a, w) {
    if (w <= 0) return(pbeta(a, thL, thL))
    if (w >= 1) return(pbeta(a, thU, thU))
    integrate(function(v) dbeta(v, thL, thL) *
                pbeta(pmin(pmax((a - (1 - w) * v) / w, 0), 1), thU, thU),
              0, 1, rel.tol = 1e-9)$value
  }
  mid <- integrate(Vectorize(function(u) {
    r <- r_of_u(u)
    w <- pbeta((u - thw) / (1 - 2 * thw), 3, 3)
    FA(1 / r, w) - FA(1 - 1 / r, w)
  }), u1, u2, rel.tol = 1e-8)$value
  u1 + mid
}

test_that("rbdegpd/rbzidegpd joint-zero rate matches the exact BEGPD integral", {
  p00 <- begpd_p00(kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  set.seed(77)
  n <- 100000L
  Y <- rbdegpd(n, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  expect_lt(abs(mean(Y[, 1] == 0 & Y[, 2] == 0) - p00),
            5 * sqrt(p00 * (1 - p00) / n))
  pz <- 0.3 + 0.7 * p00
  Yz <- rbzidegpd(n, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5,
                  thw = 0.2, pi0 = 0.3)
  expect_lt(abs(mean(Yz[, 1] == 0 & Yz[, 2] == 0) - pz), 5 * sqrt(pz * (1 - pz) / n))
})

## ---- Test gaps: Julia helpers against the v0.1 / v0.2.1 / v0.3.0 APIs -------

test_that("Julia helpers follow the real NeuralEstimators.jl v0.1.0/0.1.3/0.1.4, v0.2.1 and v0.3.0 signatures", {
  skip_on_cran()
  julia_bin <- Sys.which("julia")
  if (!nzchar(julia_bin) && file.exists("~/.juliaup/bin/julia"))
    julia_bin <- path.expand("~/.juliaup/bin/julia")
  skip_if(!nzchar(julia_bin), "Julia executable not available")
  arch <- normalizePath(system.file("julia", "begpd_architecture.jl",
                                    package = "egpd", mustWork = TRUE),
                        winslash = "/", mustWork = TRUE)
  ## Base-only mock modules reproducing the constructor/keyword/return
  ## contracts of the upstream sources (v0.1.0 Estimators.jl/inference.jl;
  ## v0.1.3 and v0.1.4 RatioEstimator call methods copied with ASCII names;
  ## v0.2.1 and v0.3.0 Estimators/RatioEstimator.jl and DeepSet).
  harness <- paste0('arch = "', arch, '"\n', r"(
using Statistics: mean
module FluxMock
  export Chain, Dense, LayerNorm, SkipConnection, relu
  relu(x) = max(zero(x), x)
  struct Dense; in_d::Int; out_d::Int; act; end
  Dense(i::Integer, o::Integer, act = identity) = Dense(Int(i), Int(o), act)
  struct LayerNorm; d::Int; end
  struct Chain; layers::Tuple; end
  Chain(xs...) = Chain(xs)
  struct SkipConnection; layer; conn; end
end
mocklr(Z, th) = vec(sum(Float32.(th); dims = 1)) .+ Float32(mean(Z))

module Era1   # NeuralEstimators.jl 0.1.0
  using ..FluxMock
  module NeuralEstimators
    using Statistics: mean
    export samplesize
    samplesize(Z) = Float32(size(Z, 2))
    struct DeepSet; psi; phi; a; S; end
    DeepSet(psi, phi, a::Function = mean; S = nothing) = DeepSet(psi, phi, a, S)
    function (d::DeepSet)(tup::Tuple)          # (Z, theta) -> 1 x G logits
      Z, th = tup
      @assert Z isa AbstractMatrix
      reshape(vec(sum(th; dims = 1)) .+ Float32(mean(Z)), 1, :)
    end
    struct PointEstimator; net; end
    struct Ensemble; ests; end
    struct RatioEstimator; deepset::DeepSet; end
    function (est::RatioEstimator)(Zth::Tuple; classifier::Bool = false)
      c = 1 ./ (1 .+ exp.(-est.deepset(Zth)))
      classifier ? c : c ./ (1 .- c)
    end
    estimateinbatches(est, Z; batchsize::Integer = 32, use_gpu::Bool = true) = zeros(Float32, 3, 1)
    sampleposterior(est::RatioEstimator, Z, N::Integer = 1000; theta_grid = nothing, kwargs...) = theta_grid[:, 1:N]
  end
  using .NeuralEstimators
end
Base.include(Era1, arch)

module Era13   # NeuralEstimators.jl 0.1.3: field `network`, sigmoid call with `classifier`
  using ..FluxMock
  module NeuralEstimators
    using Statistics: mean
    export samplesize
    samplesize(Z) = Float32(size(Z, 2))
    struct DeepSet; psi; phi; a; S; end
    DeepSet(psi, phi, a::Function = mean; S = nothing) = DeepSet(psi, phi, a, S)
    function (d::DeepSet)(tup::Tuple)   # logits large enough to saturate a Float32 sigmoid
      Z, th = tup
      reshape(vec(sum(th; dims = 1)) .+ Float32(mean(Z)) .+ 40f0, 1, :)
    end
    abstract type NeuralEstimator end
    sigmoid(x) = 1 ./ (1 .+ exp.(-x))
    struct RatioEstimator{N} <: NeuralEstimator
      network::N
    end
    function (estimator::RatioEstimator)(Z, theta; kwargs...)
      estimator((Z, theta); kwargs...)
    end
    function (estimator::RatioEstimator)(Ztheta::Tuple; classifier::Bool = false)
      c = sigmoid(estimator.network(Ztheta))
      if typeof(c) <: AbstractVector
        c = reduce(vcat, c)
      end
      classifier ? c : c ./ (1 .- c)
    end
    estimate(estimator, z, x = nothing; batchsize::Integer = 32, use_gpu::Bool = true,
             kwargs...) = zeros(Float32, 3, 1)
  end
  using .NeuralEstimators
end
Base.include(Era13, arch)

module Era14   # NeuralEstimators.jl 0.1.4: field `network`, call returns the logit, no keywords
  using ..FluxMock
  module NeuralEstimators
    using Statistics: mean
    export samplesize
    samplesize(Z) = Float32(size(Z, 2))
    struct DeepSet; psi; phi; a; S; end
    DeepSet(psi, phi, a::Function = mean; S = nothing) = DeepSet(psi, phi, a, S)
    function (d::DeepSet)(tup::Tuple)
      Z, th = tup
      reshape(vec(sum(th; dims = 1)) .+ Float32(mean(Z)) .+ 40f0, 1, :)
    end
    abstract type NeuralEstimator end
    struct RatioEstimator{N} <: NeuralEstimator
      network::N
    end
    function (estimator::RatioEstimator)(Z, theta; kwargs...)
      estimator((Z, theta); kwargs...)
    end
    function (estimator::RatioEstimator)(Ztheta::Tuple)
      logr = estimator.network(Ztheta)
      if typeof(logr) <: AbstractVector
        logr = reduce(vcat, logr)
      end
      return logr
    end
    estimate(estimator, z, x = nothing; batchsize::Integer = 32, use_gpu::Bool = true,
             kwargs...) = zeros(Float32, 3, 1)
  end
  using .NeuralEstimators
end
Base.include(Era14, arch)

module Era2   # NeuralEstimators.jl 0.2.1
  using ..FluxMock
  module NeuralEstimators
    using Statistics: mean
    export samplesize
    samplesize(Z) = Float32(size(Z, 2))
    struct DeepSet; psi; phi; a; S; end
    DeepSet(psi, phi, a::Function = mean; S = nothing) = DeepSet(psi, phi, a, S)
    struct RatioEstimator; summary_network; summary_network_theta; inference_network; end
    RatioEstimator(summary_network, num_parameters::Integer, num_summaries::Integer;
                   num_summaries_theta::Integer = 2num_parameters, kwargs...) =
      RatioEstimator(summary_network, (num_parameters, num_summaries_theta),
                     (num_summaries + num_summaries_theta, 1))
    RatioEstimator(summary_network, num_parameters::Integer; num_summaries::Integer, kwargs...) =
      RatioEstimator(summary_network, num_parameters, num_summaries; kwargs...)
    function logratio(estimator::RatioEstimator, Z; grid, kwargs...)   # K x G
      Zs = Z isa AbstractVector ? Z : [Z]
      reduce(vcat, [permutedims(vec(sum(grid; dims = 1)) .+ Float32(mean(z))) for z in Zs])
    end
    sampleposterior(estimator::RatioEstimator, Z; grid, N::Integer = 1000,
                    logprior::Function = t -> 0f0, kwargs...) = grid[:, 1:N]
    estimate(estimator, z, args...; kwargs...) = ones(Float32, 3, 1)
  end
  using .NeuralEstimators
end
Base.include(Era2, arch)

module Era3   # NeuralEstimators.jl 0.3.0
  using ..FluxMock
  module NeuralEstimators
    using Statistics: mean
    export samplesize
    samplesize(Z) = Float32(size(Z, 2))
    struct DeepSet; psi; phi; a; S; end
    DeepSet(psi, phi, a::Function = mean; condition_on_sample_size::Bool = false) =
      DeepSet(psi, phi, a, condition_on_sample_size ? :logsamplesize : nothing)
    struct RatioEstimator; summary_network; summary_network_theta; inference_network; end
    RatioEstimator(num_parameters::Integer, summary_network = identity; num_summaries::Integer,
                   num_summaries_theta::Integer = 2num_parameters, kwargs...) =
      RatioEstimator(summary_network, (num_parameters, num_summaries_theta),
                     (num_summaries + num_summaries_theta, 1))
    RatioEstimator(summary_network, num_parameters::Integer; kwargs...) =
      RatioEstimator(num_parameters, summary_network; kwargs...)
    function logratio(estimator::RatioEstimator, Z; grid, kwargs...)
      Zs = Z isa AbstractVector ? Z : [Z]
      reduce(vcat, [permutedims(vec(sum(grid; dims = 1)) .+ Float32(mean(z))) for z in Zs])
    end
    estimate(estimator, z, args...; kwargs...) = fill(2f0, 3, 1)
  end
  using .NeuralEstimators
end
Base.include(Era3, arch)

Z = rand(2, 10); th = randn(3, 11)
r1 = Era1.initialize_ratio_network(2, 3)
@assert r1 isa Era1.NeuralEstimators.RatioEstimator
@assert r1.deepset.phi.layers[1].in_d == 128 + 1 + 3
r2 = Era2.initialize_ratio_network(2, 3)
@assert r2.summary_network isa Era2.NeuralEstimators.DeepSet
@assert r2.summary_network.S === Era2._begpd_expert_summary
@assert r2.summary_network_theta == (3, 6) && r2.inference_network == (134, 1)
r3 = Era3.initialize_ratio_network(2, 3)
@assert r3.summary_network.S === :logsamplesize
@assert r3.summary_network_theta == (3, 6) && r3.inference_network == (134, 1)
r13 = Era13.initialize_ratio_network(2, 3)
r14 = Era14.initialize_ratio_network(2, 3)
@assert r13 isa Era13.NeuralEstimators.RatioEstimator && hasproperty(r13, :network)
@assert r14 isa Era14.NeuralEstimators.RatioEstimator && hasproperty(r14, :network)
for (mod, r, off) in ((Era1, r1, 0f0), (Era13, r13, 40f0), (Era14, r14, 40f0),
                      (Era2, r2, 0f0), (Era3, r3, 0f0))
  lr = mod.egpd_logratio(r, Z, th; batchsize = 4)
  @assert lr isa Vector{Float64} && length(lr) == 11 && all(isfinite, lr)
  @assert maximum(abs.(lr .- (mocklr(Z, th) .+ off))) < 1e-4
  @assert length(mod.egpd_logratio(r, Z, th[:, 1:1])) == 1
end
## what the helper avoids: the v0.1.3 classifier output saturates to 1 (log r = Inf)
## and the v0.1.4 call method has no `classifier` keyword
@assert all(r13((Float32.(Z), Float32.(th)); classifier = true) .== 1f0)
@assert try r14((Float32.(Z), Float32.(th)); classifier = true); false catch e; e isa MethodError end
@assert size(Era1.egpd_estimate(Era1.initialize_point_network(2, 3), zeros(Float32, 2, 10))) == (3, 1)
@assert Era2.egpd_estimate(r2, zeros(Float32, 2, 10)) == ones(Float32, 3, 1)
@assert Era3.egpd_estimate(r3, zeros(Float32, 2, 10)) == fill(2f0, 3, 1)
println("ARCH_ALL_ERAS_OK")
)")
  out <- system2(julia_bin, c("--startup-file=no", "-e", shQuote(harness)),
                 stdout = TRUE, stderr = TRUE)
  expect_true(any(grepl("ARCH_ALL_ERAS_OK", out, fixed = TRUE)),
              info = paste(out, collapse = "\n"))
})
