## Round 7 (native): likelihoods where exp() of a link-scale coefficient is subnormal
## (below DBL_MIN, i.e. a log-link coefficient below -708.4). exp() then keeps only a few
## significant bits while the coefficient itself is exact. Continuous EGPD carrier 5
## cancelled 0.5 log(kappa) against the log of its normaliser formed from exp(log kappa)
## (up to 0.2 per observation off, likelihoods far above the maximum, the 1e20 penalty
## below -744.1); the zero-truncated DEGPD formed log(1 - f(0)) and its cells from
## subnormal products (and its Hessian weight f0 / (1 - f0)^2 overflowed from log kappa
## ~ -355 down). References: own 400-bit Rmpfr implementations written from the
## model definitions (scripts in the round-7 notes); each limit agrees with the model at
## the subnormal coefficient to double precision.

r7_ld <- function(y, eta) {
  list(X = rep(list(matrix(1, length(y), 1)), length(eta)), y = cbind(y),
       idpars = seq_along(eta), dupid = 0L, duplicate = 0L,
       offsets = lapply(eta, rep, length(y)), censored = FALSE, xi.max = Inf)
}
r7_fns <- function(nm) get(nm, envir = asNamespace("egpd"))

test_that("continuous carrier 5 keeps its GPD limit for subnormal kappa", {
  f5 <- r7_fns(".egpd5fns")
  ## rows (y, log psi, xi); the same reference holds at every log kappa below
  ## (-708.39 and -708.40 straddle DBL_MIN, i.e. the switch)
  rows <- list(c(1.9, 0.1, 0.2), c(0, 0.1, 0.2), c(25, log(2), -0.03))
  ref <- c(1.87317917187294403, 0.1, 15.88993119283872879)
  for (i in seq_along(rows)) for (lk in c(-700, -708.39, -708.40, -742, -745, -800)) {
    q <- rows[[i]]
    expect_equal(f5$d0(c(0, 0, 0), r7_ld(q[1], c(q[2], q[3], lk))), ref[i], tolerance = 1e-13,
                 info = sprintf("row %d, log kappa %g", i, lk))
  }
  ## gradient and packed Hessian (log psi, xi, log kappa): the log-kappa entries are O(kappa)
  g1 <- c(-0.53517684210664529, 0.28763766106262639, 0, 1.1423821860165566,
          0.68465924536157197, 0, 0.096897338010231385, 0, 0)
  g2 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
  for (lk in c(-742, -800)) {
    expect_equal(as.numeric(f5$d120(c(0, 0, 0), r7_ld(1.9, c(0.1, 0.2, lk)))), g1, tolerance = 1e-12)
    expect_equal(as.numeric(f5$d120(c(0, 0, 0), r7_ld(0, c(0.1, 0.2, lk)))), g2, tolerance = 1e-12)
  }
})

test_that("a carrier-5 fit's likelihood cannot exceed its maximum at subnormal kappa", {
  set.seed(5); y <- regpd(150, sigma = 2, xi = 0.15, kappa = 1.5)
  fit <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "egpd", egpd.args = list(m = 5))
  nll <- function(th) fit$likfns$d0(th, fit$likdata)
  ## maximum of an independent double-precision likelihood (pchisq normaliser)
  nll_max <- 317.264046791304
  expect_equal(nll(fit$coefficients), nll_max, tolerance = 1e-10)
  ## at the GPD (kappa -> 0) optimum: 400-bit reference 319.03378593794719 for each log kappa
  th0 <- c(1.21241185156892, -0.0855199408081512)
  for (lk in c(-742, -742.75, -743, -745, -800)) {
    v <- nll(c(th0, lk))
    expect_equal(v, 319.03378593794719, tolerance = 1e-12, info = sprintf("log kappa %g", lk))
    expect_gt(v, nll_max)
  }
  o <- optim(th0, function(p) nll(c(p, -742.75)), method = "Nelder-Mead",
             control = list(reltol = 1e-12))
  expect_gt(o$value, 319.033785937 - 1e-9)
})

test_that("carrier 6 at y = 0 and DEGPD carrier 5 at a large count keep their limits", {
  ## kappa = 0 (log kappa -800) at y = 0 was the 1e20 penalty
  expect_equal(r7_fns(".egpd6fns")$d0(c(0, 0, 0), r7_ld(0, c(0.1, 0.2, -800))),
               -1.4060768630739156, tolerance = 1e-13)
  ## y = 1e6 (the big-count quadrature uses the continuous carrier-5 density)
  for (lk in c(-742, -745, -800)) {
    expect_equal(r7_fns(".degpd5fns")$d0(c(0, 0, 0), r7_ld(1e6, c(log(3e5), log(0.3), lk))),
                 15.615176619397726, tolerance = 1e-13)
    expect_equal(r7_fns(".degpd5idfns")$d0(c(0, 0, 0), r7_ld(1e6, c(log(3e5), 0.3, lk))),
                 15.615176619397726, tolerance = 1e-13)
  }
})

test_that("zero-truncated DEGPD keeps its limits for subnormal kappa and shape", {
  zt <- r7_fns(".ztdegpd1fns"); ztid <- r7_fns(".ztdegpd1idfns")
  ## y = 3, sigma = 3, xi = 0.2: the kappa -> 0 limit
  for (lk in c(-700, -708.39, -708.40, -742, -745, -800)) {
    expect_equal(zt$d0(c(0, 0, 0), r7_ld(3, c(log(3), log(0.2), lk))), 2.1657596022085537,
                 tolerance = 1e-12, info = sprintf("log kappa %g", lk))
    expect_equal(ztid$d0(c(0, 0, 0), r7_ld(3, c(log(3), 0.2, lk))), 2.1657596022085537,
                 tolerance = 1e-12, info = sprintf("identity, log kappa %g", lk))
  }
  ## kappa = 1.5, xi subnormal (log link): the exponential limit
  for (lx in c(-700, -742, -744, -800))
    expect_equal(zt$d0(c(0, 0, 0), r7_ld(3, c(log(3), lx, log(1.5)))), 1.8815318292433059,
                 tolerance = 1e-13, info = sprintf("log xi %g", lx))
  ## gradients and packed Hessians (log sigma, xi or log xi, log kappa)
  g_log <- c(-0.0014077829456734287, 0.051974177155983084, 0, 0.30125385084393758,
             0.010585836459244682, 0, 0.046648065486387655, 0, 0)
  g_id <- c(-0.0014077829456734291, 0.25987088577991541, 0, 0.30125385084393763,
            0.052929182296223412, 0, -0.13315279173988565, 0, 0)
  g_xi <- c(0.33037614525427139, 0, -0.094755384195667577, 1.0343571893590495, 0,
            0.51277801028776093, 0, 0, 0.15230402849713756)
  for (lk in c(-400, -600, -742, -800)) {
    expect_equal(as.numeric(zt$d120(c(0, 0, 0), r7_ld(3, c(log(3), log(0.2), lk)))), g_log, tolerance = 1e-11,
                 info = sprintf("log kappa %g", lk))
    expect_equal(as.numeric(ztid$d120(c(0, 0, 0), r7_ld(3, c(log(3), 0.2, lk)))), g_id, tolerance = 1e-11,
                 info = sprintf("identity, log kappa %g", lk))
  }
  expect_equal(as.numeric(zt$d120(c(0, 0, 0), r7_ld(3, c(log(3), -744, log(1.5))))), g_xi, tolerance = 1e-12)
})

test_that("a zero-truncated DEGPD fit's likelihood cannot exceed its maximum at subnormal kappa", {
  set.seed(5); cnt <- rdiscegpd(150, sigma = 3, xi = 0.15, kappa = 1.5, type = 1)
  fit <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = cnt[cnt > 0]), family = "ztdegpd",
              degpd.args = list(link = "identity"))
  nll <- function(th) fit$likfns$d0(th, fit$likdata)
  ## maximum of an independent double-precision likelihood
  nll_max <- 319.82056173184441
  expect_equal(nll(fit$coefficients), nll_max, tolerance = 1e-10)
  ## at the optimum of the kappa -> 0 limit: 400-bit reference for each log kappa
  lim <- c(2.3634349060504674, -0.44309756191373595)
  for (lk in c(-739.8, -742, -800)) {
    v <- nll(c(lim, lk))
    expect_equal(v, 320.75836465656448, tolerance = 1e-11, info = sprintf("log kappa %g", lk))
    expect_gt(v, nll_max)
  }
  o <- optim(lim, function(p) { v <- nll(c(p, -739.8)); if (is.finite(v)) v else 1e20 },
             method = "Nelder-Mead", control = list(reltol = 1e-12))
  expect_gt(o$value, 320.758364656 - 1e-8)
})
