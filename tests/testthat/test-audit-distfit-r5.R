## Round-5 regressions for the Bernstein fit, the Julia set-up, the compound
## Poisson-DEGPD and distfit(fix.arg) (verify-C clusters VC5-02 to VC5-08).
## Expected values are closed forms, finite sums or independent implementations,
## never the package compared with itself.

## ---- VC5-04: Bernstein density when kappa S is not small --------------------

## Independent Bernstein-EGPD log-density from base R. log v = kappa log(1 - S)
## is -kappa S below S = e^-40 (log(1 - S) = -S to double precision there), and
## the basis is v^(j - 1) (1 - v)^(m - j) with 1 - v = -expm1(log v).
bern_logdens_ind <- function(x, sigma, xi, kappa, w) {
  m <- length(w)
  t <- 1 + xi * x / sigma
  lS <- if (xi == 0) -x / sigma else -log(t) / xi
  lv <- ifelse(lS < -40, -exp(log(kappa) + lS), kappa * log1p(-exp(lS)))
  v <- exp(lv)
  omv <- -expm1(lv)
  b <- 0
  for (j in seq_len(m)) b <- b + w[j] * m * choose(m - 1, j - 1) * v^(j - 1) * omv^(m - j)
  log(b) + log(kappa) + lv * (kappa - 1) / kappa + lS - log(t) - log(sigma)
}
## the same as a negative log-likelihood of (log kappa, log sigma, xi, alpha_1..m)
bern_nll_ind <- function(th, x, m) {
  a <- th[4:(3 + m)]
  w <- exp(a - max(a))
  -sum(bern_logdens_ind(x, exp(th[2]), th[3], exp(th[1]), w / sum(w)))
}

test_that("VC5-04: Bernstein density is right where S < e^-40 but kappa S is not small", {
  ## m = 2, w = (1/2, 1/2): b(v) = 1 (partition of unity), so the density is
  ## kappa H^(kappa - 1) h; for xi = 0, sigma = 1: kappa exp((kappa - 1) log1p(-e^-x)) e^-x
  dens <- function(x, kappa) egpd:::.bernstein_full_density(x, 1, 0, kappa, c(.5, .5), 2)
  closed <- function(x, kappa) exp(log(kappa) + (kappa - 1) * log1p(-exp(-x)) - x)
  expect_equal(dens(45, 1e20), 0.1635205317, tolerance = 1e-9)     # was 0.477
  x <- c(38, 41, 45, 50, 55)
  expect_equal(dens(x, 1e20), closed(x, 1e20), tolerance = 1e-12)
  ## continuous across S = e^-40 (x = 40) for kappa = e^40, where kappa S = 1:
  ## both sides are e^(-1) (was 0.368 and 0.503)
  xs <- 40 + c(-1e-7, 1e-7)
  expect_equal(dens(xs, exp(40)), closed(xs, exp(40)), tolerance = 1e-12)
  expect_equal(dens(xs, exp(40)), rep(exp(-1), 2), tolerance = 1e-6)
  ## a density: total mass 1 (m = 8, w = e_1, kappa = e^41; it was 278)
  f <- function(z) egpd:::.bernstein_full_density(z, 1, 0, exp(41), c(1, rep(0, 7)), 8)
  expect_equal(integrate(f, 30, 60, rel.tol = 1e-10, subdivisions = 1000)$value, 1,
               tolerance = 1e-8)
  ## S subnormal (x = 709.5): b(v) = 1.4 - 0.8 v for w = (0.7, 0.3), v = exp(-kappa S)
  ## and kappa S = e^-0.5, so log f = log b - kappa S + log(kappa S)
  kS <- exp(-0.5)
  expect_equal(egpd:::.bernstein_full_density(709.5, 1, 0, exp(709), c(.7, .3), 2),
               (1.4 - 0.8 * exp(-kS)) * exp(-kS) * kS, tolerance = 1e-12)
})

test_that("VC5-04: the Bernstein gradient is that of the value", {
  for (cc in list(list(th = c(log(1e20), 0, 0, 0, 0), x = 45, m = 2),
                  list(th = c(log(1e20), 0, 0.01, 0.3, -0.2, 0), x = c(45, 46), m = 3),
                  list(th = c(log(3), log(1.3), 0.1, 0.2, -0.1, 0.3, 0),
                       x = c(0.4, 1.7, 3.2, 6.5), m = 4))) {
    expect_equal(egpd:::.bernstein_nll(cc$th, cc$x, cc$m, NULL),
                 bern_nll_ind(cc$th, cc$x, cc$m), tolerance = 1e-12)
    expect_equal(egpd:::.bernstein_nll(cc$th, cc$x, cc$m, NULL, deriv = TRUE),
                 numDeriv::grad(function(t) bern_nll_ind(t, cc$x, cc$m), cc$th),
                 tolerance = 1e-7)
  }
  ## S subnormal: m = 2, xi = 0 in closed form; with c = kappa S, v = exp(-c) and
  ## b = 2 w1 (1 - v) + 2 w2 v, d ll / d log kappa = 1 - c + 2 (w1 - w2) v c / b and
  ## d ll / d log sigma = (x / sigma) d ll / d log kappa - 1
  th <- c(709, 0, 0, log(.7 / .3), 0)
  c0 <- exp(-0.5)
  v <- exp(-c0)
  dk <- 1 - c0 + 2 * 0.4 * v * c0 / (1.4 - 0.8 * v)
  expect_equal(egpd:::.bernstein_nll(th, 709.5, 2, NULL, deriv = TRUE)[1:2],
               -c(dk, 709.5 * dk - 1), tolerance = 1e-10)
})

test_that("VC5-04: a large-kappa Bernstein fit reports the likelihood of its estimate", {
  set.seed(1)
  x <- 45 + rexp(100)
  fb <- suppressWarnings(distfit(x, method = "bernstein"))
  e <- fb$estimate
  ## the fit reaches kappa far above 1e15, where S at the data is below e^-37
  expect_gt(e[["kappa"]], 1e15)
  expect_equal(fb$loglik, sum(bern_logdens_ind(x, e[["sigma"]], e[["xi"]], e[["kappa"]],
                                               fb$bernstein.weights)),
               tolerance = 1e-10)                          # was 0.224 too high
})

test_that("VC5-04: the Bernstein CDF keeps kappa S where H rounds to 1", {
  ## m = 2, w = (1/2, 1/2): B(v) = v, so F(x) = H^kappa = exp(kappa log1p(-e^-x))
  expect_equal(egpd:::.bernstein_full_cdf(c(40, 45, 50), 1, 0, 1e20, c(.5, .5), 2),
               exp(1e20 * log1p(-exp(-c(40, 45, 50)))), tolerance = 1e-12)   # was 1
  p <- c(0.1, 0.5, 0.9)
  q <- egpd:::.bernstein_full_quantile(p, 1, 0, 1e20, c(.5, .5), 2)
  expect_equal(q, -log(-expm1(log(p) / 1e20)), tolerance = 1e-12)
  expect_equal(egpd:::.bernstein_full_cdf(q, 1, 0, 1e20, c(.5, .5), 2), p, tolerance = 1e-10)
})

## ---- VC5-02: Julia set-up -----------------------------------------------------

test_that("VC5-02: the Julia set-up loads Optimisers, which train() needs", {
  ## NeuralEstimators::train() runs `using Optimisers: Adam`, which fails unless
  ## Optimisers is a direct dependency of the Julia environment: the session set-up
  ## loads it (so that a missing package fails before any data are simulated) and
  ## names it in the install command
  code <- character()
  local_mocked_bindings(
    .check_julia_deps = function() invisible(NULL),
    .julia_eval = function(x, ...) {
      code <<- c(code, x)
      if (grepl("Optimisers", x))
        stop("ArgumentError: Package Optimisers not found in current path.")
      invisible(NULL)
    }
  )
  env <- egpd:::.begpd_env
  old <- env$julia_initialized
  env$julia_initialized <- FALSE
  on.exit(env$julia_initialized <- old, add = TRUE)
  expect_error(egpd:::.init_julia_begpd(),
               'Optimisers not found.*Pkg.add\\(\\[.*"Optimisers"\\]\\)')
  expect_match(code[1], "^using .*\\bOptimisers\\b")
  expect_false(env$julia_initialized)
})

## ---- VC5-05: qcpdegpd(1) when the grid cannot hold the severity --------------

test_that("VC5-05: qcpdegpd(1) is Inf whenever the severity can be positive", {
  ## xi = 0: the EGPD has unbounded support, so P(Y >= 1) = 1 - F(1) > 0 (here
  ## 1 - (1 - exp(-1 / sigma))^kappa, essentially 1) and S is unbounded for any
  ## lambda > 0, even when every grid mass P(Y = j), j <= K, underflows
  expect_identical(qcpdegpd(1, lambda = 1, sigma = 1e4, xi = 0, kappa = 100, K = 1), Inf)
  expect_identical(qcpdegpd(1, lambda = 1, sigma = 1e6, xi = 0, kappa = 1000), Inf)
  expect_identical(qcpdegpd(1, lambda = 1e-17, sigma = 1, xi = 0.1, kappa = 1), Inf)
  ## xi = -1, sigma = 0.5: X < 0.5, so Y = floor(X) = 0 and S = 0 almost surely
  expect_identical(qcpdegpd(1, lambda = 1, sigma = 0.5, xi = -1, kappa = 1), 0)
  ## the endpoint sigma / -xi = 1 is still degenerate (P(X >= 1) = 0)
  expect_identical(qcpdegpd(1, lambda = 2, sigma = 0.5, xi = -0.5, kappa = 1), 0)
})

## ---- VC5-06: positive part when P(Y >= 1) is tiny ----------------------------

test_that("VC5-06: CPDEGPD positive part with P(Y >= 1) below 1e-16", {
  ## kappa = 1, sigma = 0.01, xi = 0.02: P(Y >= y) = (1 + 2 y)^-50, so a = P(Y >= 1)
  ## = 3^-50 (1.4e-24), f1 = P(Y = 1) = 3^-50 - 5^-50, f2 = 5^-50 - 7^-50. With
  ## lambda = 1, P(S = 0) = exp(-a), P(S = 1) = exp(-a) f1 and
  ## P(S = 2) = exp(-a) (f2 + f1^2 / 2), exactly.
  a <- 3^-50
  f1 <- 3^-50 - 5^-50
  f2 <- 5^-50 - 7^-50
  lpos <- log(-expm1(-a))
  l1 <- -a + log(f1) - lpos
  l2 <- -a + log(f2 + f1^2 / 2) - lpos
  ## P(S = 1 | S > 0) = 1 - 8.1e-12, so the .5, .9 and .99 quantiles are all 1
  expect_gt(exp(l1), 0.99)
  expect_equal(egpd:::.aggregate_positive_quantile(c(.5, .9, .99), "cpdegpd", sigma = .01,
                                                  kappa = 1, xi = .02, lambda = 1, K = 50),
               c(1, 1, 1))                                 # was an error
  ld <- egpd:::.aggregate_positive_logdens(c(1, 2), "cpdegpd", sigma = .01, kappa = 1,
                                           xi = .02, lambda = 1, K = 50)
  expect_lt(abs(ld[1] - l1), 1e-13)                        # was -Inf
  expect_equal(ld[2], l2, tolerance = 1e-12)               # was -Inf
  expect_equal(dcpdegpd(0, lambda = 1, sigma = .01, xi = .02, kappa = 1, K = 50, log = TRUE),
               -a, tolerance = 1e-12)                      # was 0
})

test_that("VC5-06: CPDEGPD positive-part quantiles are right where 1 - P(Y = 0) is inexact", {
  ## xi = 1, kappa = 1: P(Y >= y) = sigma / (sigma + y). For sigma = 3e-16,
  ## 1 - P(Y = 0) rounds to 3.3e-16; S | S > 0 is a single event Y | Y >= 1 up to
  ## O(a) = 3e-16, with P(Y <= k | Y >= 1) = 1 - (sigma + 1) / (sigma + k + 1)
  p <- c(.45, .55, .62, .7, .78, .82, .86)
  sig <- 3e-16
  cdf1 <- function(k) 1 - (sig + 1) / (sig + k + 1)
  expected <- vapply(p, function(pp) which(cdf1(1:60) >= pp)[1], numeric(1))
  expect_equal(expected, c(1, 2, 2, 3, 4, 5, 7))
  expect_equal(egpd:::.aggregate_positive_quantile(p, "cpdegpd", sigma = sig, kappa = 1,
                                                  xi = 1, lambda = 1, K = 60),
               expected)                                   # was 1 2 3 4 7 11 22
})

## ---- VC5-08: fixed probabilities on the boundary -----------------------------

test_that("VC5-08: fix.arg rejects boundary values outside the family's parameter space", {
  ## GPIG: a in (0, 1], c in [0, 1]; ZIPIG, ZIGPIG and ZIBELL: pi in [0, 1).
  ## Outside, the range error (and no warnings) instead of a generic support error.
  set.seed(3)
  x <- rpois(80, 1.5)
  fit_err <- function(...) {
    w <- character()
    e <- withCallingHandlers(tryCatch(distfit(...), error = identity),
                             warning = function(cnd) {
                               w <<- c(w, conditionMessage(cnd))
                               invokeRestart("muffleWarning")
                             })
    list(error = if (inherits(e, "error")) conditionMessage(e) else "", warnings = w)
  }
  r <- fit_err(x, family = "gpig", fix.arg = list(a = 0))
  expect_match(r$error, "'fix.arg' value for a must be a single number in \\(0, 1\\]")
  expect_length(r$warnings, 0)
  for (fam in c("zipig", "zigpig", "zibell")) {
    r <- fit_err(rep(0, 10), family = fam, fix.arg = list(pi = 1))
    expect_match(r$error, "'fix.arg' value for pi must be a single number in \\[0, 1\\)",
                 info = fam)
    expect_length(r$warnings, 0)
  }
  ## the boundaries the families do allow are still accepted: a = 1 and c = 1 give
  ## the Poisson fit (GPIG a = 1 is Poisson with mean b c; with c = 1 the
  ## underdispersed sample takes a to 1), whose maximum is the closed form at mean(x)
  ll_pois <- sum(stats::dpois(x, mean(x), log = TRUE))
  f_a <- suppressWarnings(distfit(x, family = "gpig", fix.arg = list(a = 1)))
  expect_equal(f_a$loglik, ll_pois, tolerance = 1e-7)
  f_c <- suppressWarnings(distfit(x, family = "gpig", fix.arg = list(c = 1)))
  expect_equal(f_c$loglik, ll_pois, tolerance = 1e-7)
  ## a structural-zero probability of 1 is a model for ZIDEGPD (all zeros: P = 1)
  f_z <- suppressWarnings(distfit(rep(0, 10), family = "zidegpd", fix.arg = list(pi = 1)))
  expect_equal(f_z$loglik, 0)
})
