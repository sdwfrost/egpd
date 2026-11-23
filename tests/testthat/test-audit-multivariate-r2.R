## Round 2 audit regression tests for multivariate, neural-Bayes/Julia,
## compound EGPD, blended GP, and aggregate EGPD

## Round 3: exact Definition 2.1 CDF, 1 - E[1 ^ exp(max(S - h(k)))] with
## S = T - max(T), T iid Poisson(rate), summed over the truncated Poisson
## lattice (previously estimated from the RNG stream that rmdgpd() uses).
mdgpd_def21_cdf <- function(k, sigma, xi, rate, d, Tmax = 16L) {
  sigma <- rep_len(sigma, d)
  xi <- rep_len(xi, d)
  k <- rep_len(k, d)
  h <- numeric(d)
  for (j in seq_len(d)) {
    if (xi[j] == 0) {
      h[j] <- k[j] / sigma[j]
    } else {
      arg <- 1 + xi[j] * k[j] / sigma[j]
      if (arg <= 0) {
        h[j] <- if (xi[j] < 0) Inf else -Inf
      } else {
        h[j] <- log1p(xi[j] * k[j] / sigma[j]) / xi[j]
      }
    }
  }
  Tg <- as.matrix(expand.grid(rep(list(0:Tmax), d)))
  w <- apply(matrix(dpois(Tg, rate), ncol = d), 1L, prod)
  shifted <- sweep(Tg - apply(Tg, 1L, max), 2L, h, "-")
  max_shift <- apply(shifted, 1L, max)
  1 - sum(w * pmin(1, exp(max_shift))) / sum(w)
}

test_that("MDGPD supports xi < 0, rate = 0, vector (sigma, xi), and matches Definition 2.1", {
  ## Exact upper endpoint when -sigma / xi is an integer vs non-integer
  set.seed(101)
  y_int_ep <- rmdgpd(30000L, sigma = 2, xi = -0.5, rate = 0.6, d = 2L)
  expect_true(all(is.finite(y_int_ep)))
  expect_type(y_int_ep, "integer")
  expect_true(all(y_int_ep <= 4L))
  expect_true(any(y_int_ep == 4L))
  expect_true(all(apply(y_int_ep, 1L, max) >= 1L))

  ## Conditional positive margin P(M_j > k | M_j > 0) = (1 + xi * k / sigma)_+^(-1/xi)
  pos1 <- y_int_ep[y_int_ep[, 1] > 0L, 1]
  expect_equal(mean(pos1 > 2L), (1 - 0.5 * 2 / 2)^(1 / 0.5), tolerance = 0.012)
  expect_equal(mean(pos1 > 4L), 0)

  ## Heterogeneous vector (sigma, xi) including Aka-Kratz-Naveau Table 3 negative shapes
  ## and mixed-sign shapes in d = 3, verified against Definition 2.1 CDF
  set.seed(202)
  n <- 60000L
  rate <- 0.7
  sigma_v <- c(3.2, 4.0, 2.5)
  xi_v <- c(-0.18, -0.25, 0.15)
  y_vec <- rmdgpd(n, sigma = sigma_v, xi = xi_v, rate = rate, d = 3L)
  expect_equal(dim(y_vec), c(n, 3L))
  expect_true(all(y_vec[, 1] <= ceiling(-sigma_v[1] / xi_v[1])))
  expect_true(all(y_vec[, 2] <= ceiling(-sigma_v[2] / xi_v[2])))

  for (k_test in list(c(3L, 4L, 2L), c(8L, 16L, 5L), c(0L, 1L, 1L))) {
    emp_cdf <- mean(y_vec[, 1] <= k_test[1] &
                    y_vec[, 2] <= k_test[2] &
                    y_vec[, 3] <= k_test[3])
    theo_cdf <- mdgpd_def21_cdf(k_test, sigma_v, xi_v, rate, 3L)
    expect_lt(abs(emp_cdf - theo_cdf), 5 * sqrt(theo_cdf * (1 - theo_cdf) / n) + 1e-4)
  }

  ## rate = 0 gives exact comonotonic limit S == 0
  set.seed(303)
  y_como <- rmdgpd(500L, sigma = 2.5, xi = -0.2, rate = 0, d = 4L)
  expect_identical(y_como[, 1], y_como[, 2])
  expect_identical(y_como[, 1], y_como[, 3])
  expect_identical(y_como[, 1], y_como[, 4])
  expect_true(all(y_como >= 1L))

  ## Zero-inflated heterogeneous negative-xi MDGPD
  set.seed(304)
  y_zi <- rzimdgpd(5000L, sigma = c(3, 5), xi = c(-0.2, 0.1),
                   rate = 0.4, pi0 = 0.25, d = 2L)
  expect_equal(mean(y_zi[, 1] == 0L & y_zi[, 2] == 0L), 0.25, tolerance = 0.02)

  ## Taylor continuity for |xi| near 0
  set.seed(305); y0 <- rmdgpd(200L, sigma = 3, xi = 0, rate = 0.5)
  set.seed(305); y_tiny_neg <- rmdgpd(200L, sigma = 3, xi = -1e-12, rate = 0.5)
  set.seed(305); y_tiny_pos <- rmdgpd(200L, sigma = 3, xi = 1e-12, rate = 0.5)
  expect_identical(y0, y_tiny_neg)
  expect_identical(y0, y_tiny_pos)

  ## Invalid vector lengths and non-finite entries rejected
  expect_error(rmdgpd(10, sigma = c(1, 2, 3), xi = 0.1, rate = 1, d = 2L), "'sigma'")
  expect_error(rmdgpd(10, sigma = c(1, -2), xi = 0.1, rate = 1, d = 2L), "'sigma'")
  expect_error(rmdgpd(10, sigma = 2, xi = c(0.1, 0.2, 0.3), rate = 1, d = 2L), "'xi'")
  expect_error(rmdgpd(10, sigma = 2, xi = c(0.1, NA), rate = 1, d = 2L), "'xi'")
})

test_that("log-space Panjer handles lambda = 800 and 5000 without underflow", {
  ## Exact analytical oracle: when severity is degenerate at 1 (f = c(0, 1, 0, ...)),
  ## S ~ Poisson(lambda) so Panjer must match dpois(0:K, lambda) to machine precision!
  K_800 <- 1200L
  f_unit <- c(0, 1, rep(0, K_800 - 1L))
  g_800 <- egpd:::.panjer_poisson(800, f_unit)
  expect_equal(sum(g_800), 1, tolerance = 1e-10)
  expect_equal(g_800, dpois(0:K_800, lambda = 800), tolerance = 1e-12)

  lg_800 <- egpd:::.panjer_poisson(800, f_unit, log_out = TRUE)
  expect_equal(lg_800[1L], -800, tolerance = 1e-12)
  expect_equal(lg_800, dpois(0:K_800, lambda = 800, log = TRUE), tolerance = 1e-10)

  ## Compound Poisson-EGPD at lambda = 800 sums to 1 and preserves log P(S = 0) = -800
  pmf_800 <- dcpegpd((0:2500) * 0.5, lambda = 800, sigma = 0.5, xi = 0.05,
                     kappa = 1.2, type = 1, h = 0.5, K = 2500L)
  expect_true(all(is.finite(pmf_800) & pmf_800 >= 0))
  expect_equal(sum(pmf_800), 1, tolerance = 1e-4)
  expect_equal(dcpegpd(0, lambda = 5000, sigma = 1, xi = 0.1, kappa = 1,
                       type = 1, h = 0.5, K = 10L, log = TRUE),
               -5000, tolerance = 1e-12)
})

test_that("cpegpd realizes Ailliot Eq. F.2 exact zero atom and validates all 8 d/p/q/r functions", {
  lam <- 1.7
  sig <- 1.4
  xi <- 0.15
  kap <- 0.8
  h <- 0.25

  ## Exact zero atom P(S = 0) = exp(-lambda) regardless of h
  for (h_val in c(0.05, 0.25, 1.0)) {
    expect_equal(dcpegpd(0, lambda = lam, sigma = sig, xi = xi,
                         kappa = kap, type = 1, h = h_val),
                 exp(-lam), tolerance = 1e-13)
    expect_equal(pcpegpd(0, lambda = lam, sigma = sig, xi = xi,
                         kappa = kap, type = 1, h = h_val),
                 exp(-lam), tolerance = 1e-13)
  }

  ## Positive values in ((k-1)h, kh] map to bin k = ceiling(x / h) >= 1 (never bin 0)
  f1 <- pegpd(h, sigma = sig, xi = xi, kappa = kap, type = 1)
  expect_equal(dcpegpd(c(1e-12, 0.5 * h, h), lambda = lam, sigma = sig,
                       xi = xi, kappa = kap, type = 1, h = h),
               rep(exp(-lam) * lam * f1, 3), tolerance = 1e-12)
  ## Floating-point grid multiple 3 * 0.1 does not jump to bin 4
  expect_equal(dcpegpd(3 * 0.1, lambda = lam, sigma = sig, xi = xi,
                       kappa = kap, type = 1, h = 0.1),
               dcpegpd(0.3, lambda = lam, sigma = sig, xi = xi,
                       kappa = kap, type = 1, h = 0.1))

  ## Upfront validation in all 8 compound functions (including tiny lambda where N == 0)
  fns_x <- list(dcpegpd, pcpegpd, dcpdegpd, pcpdegpd)
  for (fn in fns_x) {
    expect_error(fn(1, lambda = -1, sigma = 1, xi = 0.1, kappa = 1, type = 1), "'lambda'")
    expect_error(fn(1, lambda = 1, sigma = 0, xi = 0.1, kappa = 1, type = 1), "'sigma'")
    ## Round 3: supplied NA/NaN parameters now give NA/NaN (see the r3 contract
    ## test); other invalid values still fail upfront.
    expect_error(fn(1, lambda = 1, sigma = 1, xi = Inf, kappa = 1, type = 1), "'xi'")
    expect_error(fn(1, lambda = 1, sigma = 1, xi = 0.1, kappa = -1, type = 1), "'kappa'")
    expect_error(fn(1, lambda = 1, sigma = 1, xi = 0.1, delta = -1, type = 4), "'delta'")
    expect_error(fn(1, lambda = 1, sigma = 1, xi = 0.1, kappa = 1, delta = 2,
                    prob = 1.5, type = 6), "'prob'")
  }
  for (qfn in list(qcpegpd, qcpdegpd)) {
    expect_error(qfn(0.5, lambda = 0, sigma = 1, xi = 0.1, kappa = 1, type = 1), "'lambda'")
    expect_error(qfn(0.5, lambda = 1, sigma = -1, xi = 0.1, kappa = 1, type = 1), "'sigma'")
  }
  ## Even when lambda = 1e-12 so rpois draws 0, invalid severity params still fail
  for (rfn in list(rcpegpd, rcpdegpd)) {
    expect_error(rfn(5, lambda = 1e-12, sigma = -1, xi = 0.1, kappa = 1, type = 1), "'sigma'")
    expect_error(rfn(5, lambda = 1e-12, sigma = 1, xi = Inf, kappa = 1, type = 1), "'xi'")
    expect_error(rfn(5, lambda = 1e-12, sigma = 1, xi = 0.1, kappa = 0, type = 1), "'kappa'")
    expect_error(rfn(-1, lambda = 1, sigma = 1, xi = 0.1, kappa = 1, type = 1), "'n'")
  }
  expect_error(dcpegpd(1, lambda = 1, sigma = 1, xi = 0.1, kappa = 1, h = 0), "'h'")
  expect_error(dcpegpd(1, lambda = 1, sigma = 1, xi = 0.1, kappa = 1, K = 0), "'K'")

  ## rcpdegpd returns integer vector and guards 32-bit overflow
  set.seed(12)
  draws <- rcpdegpd(50, lambda = 2, sigma = 3, xi = 0.1, kappa = 1.2, type = 1)
  expect_type(draws, "integer")
  expect_identical(rcpdegpd(0, lambda = 2, sigma = 3, xi = 0.1, kappa = 1.2, type = 1),
                   integer(0))
  set.seed(13)
  expect_error(rcpdegpd(20, lambda = 5, sigma = 1e10, xi = 0.2, kappa = 1, type = 1),
               "integer range")
})

test_that("qbgpd handles endpoint roundoff and .bgpd_prepare catches narrow multimodal troughs", {
  ## Perturbed bulk CDF where pbulk(qbulk(pa)) > pa and Fgp(b) > pb-eps boundary roundoff
  pa <- 0.8
  pb <- 0.99
  pbulk_pert <- function(q) {
    val <- pbeta(q, 2, 5)
    ## Perturb slightly upward near qbeta(pa, 2, 5) within numerical tolerance
    a0 <- qbeta(pa, 2, 5)
    val + 1e-12 * exp(-((q - a0) / 1e-5)^2)
  }
  dbulk_pert <- function(x) {
    a0 <- qbeta(pa, 2, 5)
    dbeta(x, 2, 5) - 1e-12 * (2 * (x - a0) / 1e-10) * exp(-((x - a0) / 1e-5)^2)
  }
  qbulk_pert <- function(p) qbeta(p, 2, 5)

  probs_edge <- c(
    pa,
    pa + .Machine$double.eps,
    pa + 1e-13,
    0.5 * (pa + pb),
    pb - 1e-13,
    pb - .Machine$double.eps,
    pb
  )
  qs <- qbgpd(probs_edge, xi = 0.15, pbulk = pbulk_pert, dbulk = dbulk_pert,
              qbulk = qbulk_pert, bulk.args = list(), pa = pa, pb = pb)
  expect_true(all(is.finite(qs)))
  expect_true(all(diff(qs) >= -1e-15))

  ## Multimodal 4-component Beta mixture bulk on [0, 1] creating negative blended
  ## density troughs in the blending interval; must be rejected by .bgpd_prepare()
  w_mix <- c(0.70, 0.14, 0.14, 0.02)
  s1 <- c(2, 547 * 8, 553 * 8, 20)
  s2 <- c(20, 453 * 8, 447 * 8, 2)
  pbulk_mm <- function(q) {
    w_mix[1] * pbeta(q, s1[1], s2[1]) + w_mix[2] * pbeta(q, s1[2], s2[2]) +
      w_mix[3] * pbeta(q, s1[3], s2[3]) + w_mix[4] * pbeta(q, s1[4], s2[4])
  }
  dbulk_mm <- function(x) {
    w_mix[1] * dbeta(x, s1[1], s2[1]) + w_mix[2] * dbeta(x, s1[2], s2[2]) +
      w_mix[3] * dbeta(x, s1[3], s2[3]) + w_mix[4] * dbeta(x, s1[4], s2[4])
  }
  qbulk_mm <- function(p) {
    vapply(p, function(pi) {
      if (pi <= 0) return(0)
      if (pi >= 1) return(1)
      uniroot(function(x) pbulk_mm(x) - pi, c(0, 1), tol = 1e-14)$root
    }, numeric(1))
  }
  expect_true(all(dbulk_mm(seq(0.001, 0.999, length.out = 2000L)) > 0))
  expect_error(
    pbgpd(0.5, xi = -0.8, pbulk = pbulk_mm, dbulk = dbulk_mm,
          qbulk = qbulk_mm, bulk.args = list(), pa = 0.1, pb = 0.995,
          c1 = 25, c2 = 25),
    "negative"
  )

  ## Check the *bgpd Rd files have no duplicated Rd aliases: the source man/
  ## under devtools, the installed Rd database under R CMD check
  man_dir <- test_path("../../man")
  rd_list <- if (dir.exists(man_dir)) {
    lapply(Sys.glob(file.path(man_dir, "*bgpd.Rd")), tools::parse_Rd)
  } else {
    rd_db <- tools::Rd_db("egpd")
    rd_db[grepl("bgpd\\.Rd$", names(rd_db))]
  }
  expect_gt(length(rd_list), 0L)
  aliases <- unlist(lapply(rd_list, function(rd) {
    unlist(lapply(rd, function(el) {
      if (identical(attr(el, "Rd_tag"), "\\alias")) as.character(el) else character()
    }))
  }))
  expect_equal(anyDuplicated(aliases), 0L)
})

test_that("fit_aggregate_egpd handles BFGS zero slope, p=q=2 extrapolation, K forwarding, order, and maxit=0", {
  skip_on_cran()
  set.seed(77)
  durs <- c(1L, 2L, 4L, 8L)
  ## 1. BFGS starting at exact zero slope on positive-slope data
  samples_inc <- setNames(lapply(durs, function(d) {
    rcpdegpd(200L,
             sigma = exp(log(1.5) + 0.18 * log(d)),
             xi = 0.15, kappa = 0.5,
             lambda = exp(log(1.0) + 0.25 * log(d)),
             type = 1)
  }), durs)

  zero_slope_start <- list(
    sigma_coef = c(log(1.5), 0),
    lambda_coef = c(log(1.0), 0),
    kappa = 0.5,
    xi = 0.15
  )
  fit_bfgs <- fit_aggregate_egpd(
    samples_inc,
    durations = durs,
    family = "cpdegpd",
    p = 1L, q = 1L,
    start = zero_slope_start,
    optim.method = "BFGS",
    control = list(maxit = 150L)
  )
  expect_equal(fit_bfgs$convergence, 0L)
  expect_gt(fit_bfgs$sigma_coef["s1"], 0.05)
  expect_gte(fit_bfgs$lambda_coef["l1"], 0)
  ## composite_loglik equals exact unpenalized negative NLL
  theta_hat <- egpd:::.aggregate_pack_theta(
    list(sigma_coef = unname(fit_bfgs$sigma_coef),
         lambda_coef = unname(fit_bfgs$lambda_coef),
         kappa = fit_bfgs$kappa, xi = fit_bfgs$xi),
    fit_bfgs$kappa_bounds, fit_bfgs$xi_bounds
  )
  prep <- egpd:::.prepare_aggregate_input(samples_inc, durs)
  expect_equal(
    fit_bfgs$composite_loglik,
    -egpd:::.aggregate_nll(theta_hat, prep, "cpdegpd", 1L, 1L, 0.2,
                           fit_bfgs$kappa_bounds, fit_bfgs$xi_bounds),
    tolerance = 1e-10
  )

  ## 2. BFGS at zero slope on constant-parameter data stays non-negative and finite
  set.seed(78)
  samples_const <- setNames(lapply(durs, function(d) {
    rcpdegpd(200L, sigma = 1.8, xi = 0.15, kappa = 0.6, lambda = 1.2, type = 1)
  }), durs)
  fit_const <- fit_aggregate_egpd(
    samples_const,
    durations = durs,
    family = "cpdegpd",
    p = 1L, q = 1L,
    start = list(sigma_coef = c(log(1.8), 0), lambda_coef = c(log(1.2), 0),
                 kappa = 0.6, xi = 0.15),
    optim.method = "BFGS",
    control = list(maxit = 150L)
  )
  expect_gte(fit_const$sigma_coef["s1"], 0)
  expect_gte(fit_const$lambda_coef["l1"], 0)
  expect_true(is.finite(fit_const$composite_loglik))

  ## 3. Quadratic p = 2, q = 2 extrapolation non-crossing up to d = 128
  fit_quad <- fit_aggregate_egpd(
    samples_inc,
    durations = durs,
    family = "cpdegpd",
    p = 2L, q = 2L,
    start = list(sigma_coef = c(log(1.4), 0.1, 0.01),
                 lambda_coef = c(log(0.9), 0.2, 0.01),
                 kappa = 0.5, xi = 0.15),
    optim.method = "BFGS",
    control = list(maxit = 150L)
  )
  expect_true(all(fit_quad$sigma_coef[-1L] >= 0))
  expect_true(all(fit_quad$lambda_coef[-1L] >= 0))
  extrap_durs <- c(1L, 2L, 4L, 8L, 16L, 32L, 64L, 128L)
  extrap_pars <- predict(fit_quad, durations = extrap_durs, type = "parameter")
  expect_true(all(diff(extrap_pars$sigma) >= -1e-12))
  expect_true(all(diff(extrap_pars$lambda) >= -1e-12))
  extrap_q <- predict(fit_quad, durations = extrap_durs, type = "quantile",
                      prob = c(0.5, 0.9, 0.99), K = 8000L)
  for (col in seq_len(ncol(extrap_q))) {
    expect_true(all(diff(extrap_q[, col]) >= 0))
  }

  ## 4. Non-monotone start rejected; unconverged fit warns; maxit = 0 safe
  expect_error(
    fit_aggregate_egpd(samples_inc, durations = durs, family = "cpdegpd",
                       p = 1L, q = 1L,
                       start = list(sigma_coef = c(0.4, -0.1),
                                    lambda_coef = c(0, 0.1),
                                    kappa = 0.5, xi = 0.15)),
    "non-negative"
  )
  expect_warning(
    fit_aggregate_egpd(samples_inc, durations = durs, family = "cpdegpd",
                       p = 1L, q = 1L, start = zero_slope_start,
                       optim.method = "Nelder-Mead",
                       control = list(maxit = 1L)),
    "did not converge"
  )
  fit_zero_nm <- fit_aggregate_egpd(
    samples_inc, durations = durs, family = "cpdegpd",
    p = 1L, q = 1L, start = zero_slope_start,
    optim.method = "Nelder-Mead",
    control = list(maxit = 0L)
  )
  expect_equal(unname(fit_zero_nm$sigma_coef), zero_slope_start$sigma_coef)
  expect_equal(unname(fit_zero_nm$lambda_coef), zero_slope_start$lambda_coef)
  expect_true(is.finite(fit_zero_nm$composite_loglik))

  ## 5. predict preserves exact order and duplicates of durations and forwards K
  req_durs <- c(8L, 1L, 8L, 2L)
  p_ord <- predict(fit_bfgs, durations = req_durs, type = "parameter")
  expect_identical(p_ord$duration, req_durs)
  expect_identical(p_ord$sigma[1L], p_ord$sigma[3L])
  q_ord <- predict(fit_bfgs, durations = req_durs, type = "quantile",
                   prob = c(0.5, 0.9), K = 5000L)
  expect_identical(rownames(q_ord), as.character(req_durs))
  expect_identical(q_ord[1L, ], q_ord[3L, ])
  expect_error(
    predict(fit_bfgs, durations = c(8L), type = "quantile", prob = 0.99, K = 2L),
    "increase 'K'"
  )
})

test_that("rbegpd, rbdegpd, and rbzidegpd reject non-finite parameters and avoid large-kappa cancellation", {
  expect_error(rbegpd(10, kappa = Inf, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2), "'kappa'")
  expect_error(rbegpd(10, kappa = 2, sigma = NaN, xi = 0.1, thL = 5, thU = 5, thw = 0.2), "'sigma'")
  expect_error(rbegpd(10, kappa = 2, sigma = 1, xi = NA_real_, thL = 5, thU = 5, thw = 0.2), "'xi'")
  expect_error(rbegpd(10, kappa = 2, sigma = 1, xi = 0.1, thL = Inf, thU = 5, thw = 0.2), "'thL'")
  expect_error(rbegpd(10, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = NaN, thw = 0.2), "'thU'")
  expect_error(rbegpd(10, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = NA_real_), "'thw'")
  expect_error(rbdegpd(10, kappa = NA_real_, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2), "'kappa'")
  expect_error(rbzidegpd(10, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2, pi0 = NA_real_), "'pi0'")
  expect_error(rbzidegpd(10, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2, pi0 = NaN), "'pi0'")

  ## At kappa = 1e12, naive 1 - U^(1/kappa) cancels to 0 (giving R = 0 for xi = 0 and R = Inf for xi > 0)
  set.seed(99)
  Y_zero <- rbegpd(50, kappa = 1e12, sigma = 2, xi = 0, thL = 5, thU = 5, thw = 0.2)
  Y_pos  <- rbegpd(50, kappa = 1e12, sigma = 2, xi = 0.2, thL = 5, thU = 5, thw = 0.2)
  expect_true(all(is.finite(Y_zero) & Y_zero > 0))
  expect_true(all(is.finite(Y_pos) & Y_pos > 0))
})

test_that("stale bundled .bson files are removed, detectCores NA is guarded, and mocked neural S3/plots restore par()", {
  ## 1. No bundled .bson files remain in inst/models
  models_dir <- system.file("models", package = "egpd")
  bson <- if (nzchar(models_dir)) Sys.glob(file.path(models_dir, "*.bson")) else character()
  expect_length(bson, 0L)
  for (fam in c("begpd", "bdegpd", "bzidegpd", "mdgpd", "zimdgpd")) {
    expect_error(
      distfit(matrix(1:20, ncol = 2), family = fam, method = "neuralbayes"),
      "train_begpd\\(\\)|train_bdegpd\\(\\)|train_mdgpd\\(\\)"
    )
  }

  ## 2. Guard detectCores() returning NA
  expect_identical(egpd:::.sanitize_mc_cores(NA_integer_),
                   if (.Platform$OS.type == "windows") 1L else egpd:::.sanitize_mc_cores(NULL))
  local_mocked_bindings(
    .sanitize_mc_cores = egpd:::.sanitize_mc_cores,
    .package = "egpd"
  )
  res_sim <- egpd:::.neural_simulate(1:4, function(i) i * 2L, mc.cores = NA_integer_)
  expect_identical(res_sim, as.list(c(2L, 4L, 6L, 8L)))

  ## 3. Mocked neural fits, S3 methods, SIR grid scaling, negative-xi MDGPD, and diagnostic plots
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  make_mock_fit <- function(family, estimator = "npe", data_dim = 2L, nsamples = 25L) {
    path <- test_path(paste0("r2-mock-", family, "-", data_dim, "d-", estimator, ".bson"))
    file.create(path)
    on.exit(unlink(c(path, paste0(path, ".metadata.rds"))), add = TRUE)
    egpd:::.write_neural_metadata(path, family, estimator, data_dim)

    raw_vec <- switch(family,
      "begpd" = , "bdegpd" = log(c(2, 1.2, 0.1, 5, 5, 0.2)),
      "bzidegpd" = c(log(c(2, 1.2, 0.1, 5, 5, 0.2)), qlogis(0.3)),
      "mdgpd" = c(log(2.5), -0.22, log(0.6)),
      "zimdgpd" = c(log(2.5), -0.18, log(0.6), qlogis(0.25))
    )
    dat <- switch(family,
      "begpd" = rbegpd(80, 2, 1.2, 0.1, 5, 5, 0.2),
      "bdegpd" = rbdegpd(80, 2, 1.2, 0.1, 5, 5, 0.2),
      "bzidegpd" = rbzidegpd(80, 2, 1.2, 0.1, 5, 5, 0.2, 0.3),
      "mdgpd" = rmdgpd(80, 2.5, -0.22, rate = 0.6, d = data_dim),
      "zimdgpd" = rzimdgpd(80, 2.5, -0.18, rate = 0.6, pi0 = 0.25, d = data_dim)
    )

    ## Round 3 (D3): NPE draws come from adaptive importance sampling on the
    ## learned log-ratio (mocked as a narrow Gaussian bump), not a prior grid.
    local_mocked_bindings(
      .check_julia_deps = function() NULL,
      .init_julia_begpd = function() NULL,
      .init_npe_begpd = function(...) "mock_npe",
      .init_nbe_begpd = function(...) "mock_nbe",
      .neural_loadstate = function(est, ...) est,
      .neural_logratio = function(estimator, Z, theta) {
        -0.5 * colSums(((theta - raw_vec) / 0.05)^2)
      },
      .julia_let = function(...) matrix(raw_vec, ncol = 1L),
      .package = "egpd"
    )
    fit <- distfit(dat, family = family, method = "neuralbayes",
                   model.path = path, estimator = estimator, nsamples = nsamples)
    list(fit = fit)
  }

  ## Many posterior draws: importance sampling keeps them distinct
  big <- make_mock_fit("mdgpd", "npe", data_dim = 2L, nsamples = 2500L)
  expect_gte(big$fit$importance_sampling$ess, 1000)
  expect_gt(big$fit$importance_sampling$unique, 1000)
  expect_lt(big$fit$estimate["xi"], 0)

  ## Exercise S3 methods & plots across families, estimators, and d = 2, 3
  cases <- list(
    list("begpd", "npe", 2L),
    list("begpd", "nbe", 2L),
    list("bdegpd", "npe", 2L),
    list("bzidegpd", "nbe", 2L),
    list("mdgpd", "npe", 2L),
    list("mdgpd", "npe", 3L),
    list("zimdgpd", "nbe", 3L)
  )
  for (cs in cases) {
    res <- make_mock_fit(cs[[1]], cs[[2]], cs[[3]], nsamples = 30L)
    fit <- res$fit
    expect_output(print(fit))
    s <- summary(fit)
    expect_s3_class(s, "summary.distfit")
    expect_output(print(s))
    expect_equal(nobs(fit), 80L)
    expect_warning(ll <- logLik(fit), "not available")
    expect_true(is.na(as.numeric(ll)))
    if (cs[[2]] == "npe") {
      expect_equal(nrow(vcov(fit)), fit$npar)
      expect_equal(dim(confint(fit)), c(fit$npar, 2L))
    } else {
      expect_error(vcov(fit), "not available")
      expect_error(confint(fit), "not available")
    }

    ## Verify par("mfrow"), par("mar"), and par("oma") are restored cleanly!
    par(mfrow = c(1L, 3L), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(1, 1, 1, 1))
    before_par <- par(c("mfrow", "mar", "oma"))
    expect_no_error(plot(fit))
    after_par <- par(c("mfrow", "mar", "oma"))
    expect_identical(after_par$mfrow, before_par$mfrow)
    expect_equal(after_par$mar, before_par$mar)
    expect_equal(after_par$oma, before_par$oma)
  }
})
