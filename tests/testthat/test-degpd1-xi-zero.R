## The DEGPD model-1 density is built from (1 + xi*y/sigma)^(1/xi). Forming the base
## as `1 + z` and then raising it throws away z once z < eps: at xi = 1e-15 the base
## rounds to exactly 1, both survival terms collapse to the same value, and the
## log-likelihood becomes -Inf. Errors were already ~3e-2 at xi = 1e-13 and changed
## sign by 1e-14.
##
## This is not a corner case: the *log* link forces xi = exp(eta) > 0, so on light-
## or bounded-tailed data the optimiser drives eta -> -Inf straight into that region
## and then optimises a corrupted objective. `degpd_pow1p()` in src/degpd.cpp routes
## these powers through log1p.

## exact reference, independent of the package: F(x) = H(x)^kappa with H the GPD cdf,
## discretised as P(Y = y) = F(y+1) - F(y)
ref_ldegpd1 <- function(y, sigma, xi, kappa) {
  H <- function(x) if (abs(xi) < 1e-300) -expm1(-x / sigma)
                   else -expm1(-log1p(xi * x / sigma) / xi)
  log(H(y + 1)^kappa - H(y)^kappa)
}

## the internal (C++) negative log-likelihood at a single observation
d0_at <- function(fit, y, sigma, xi, kappa) {
  ld <- fit$likdata
  ld$X <- list(matrix(1, 1, 1), matrix(1, 1, 1), matrix(1, 1, 1))
  ld$y <- matrix(y, 1, 1)
  ld$nobs <- 1L
  ld$idpars <- c(1, 2, 3)
  ld$offsets <- list(numeric(0), numeric(0), numeric(0))
  -fit$likfns$d0(c(log(sigma), log(xi), log(kappa)), ld)
}

test_that("the log-link DEGPD1 density stays accurate as xi approaches zero", {
  set.seed(11)
  d <- data.frame(cases = rnbinom(60, size = 2, mu = 20))
  fit <- egpd(list(lsigma = cases ~ 1, lxi = ~1, lkappa = ~1), data = d,
              family = "degpd", degpd.args = list(m = 1), restarts = FALSE, trace = 0)

  for (xi in 10^-(4:16)) {
    got  <- d0_at(fit, 7, 21, xi, 2.11)
    want <- ref_ldegpd1(7, 21, xi, 2.11)
    expect_true(is.finite(got),
                info = sprintf("xi = %.0e gave a non-finite log-density", xi))
    expect_equal(got, want, tolerance = 1e-10,
                 info = sprintf("xi = %.0e", xi))
  }
})

test_that("the internal density agrees with dDEGPD1 across the whole xi range", {
  set.seed(11)
  d <- data.frame(cases = rnbinom(60, size = 2, mu = 20))
  fit <- egpd(list(lsigma = cases ~ 1, lxi = ~1, lkappa = ~1), data = d,
              family = "degpd", degpd.args = list(m = 1), restarts = FALSE, trace = 0)

  for (xi in c(0.5, 0.1, 1e-3, 1e-6, 1e-9, 1e-12, 1e-15)) {
    for (y in c(0, 1, 7, 40)) {
      expect_equal(d0_at(fit, y, 21, xi, 2.11),
                   dDEGPD1(y, mu = 21, sigma = xi, nu = 2.11, log = TRUE),
                   tolerance = 1e-6,
                   info = sprintf("xi = %.0e, y = %d", xi, y))
    }
  }
})

test_that("degpd1d12 derivatives still match the objective at tiny xi", {
  set.seed(11)
  d <- data.frame(cases = rnbinom(80, size = 3, mu = 15))
  fit <- egpd(list(lsigma = cases ~ 1, lxi = ~1, lkappa = ~1), data = d,
              family = "degpd", degpd.args = list(m = 1), restarts = FALSE, trace = 0)
  ld <- fit$likdata
  f <- function(b) fit$likfns$d0(b, ld)

  for (xi in c(0.3, 1e-3, 1e-9, 1e-15)) {
    b <- c(log(15), log(xi), log(1.5))
    analytic <- colSums(fit$likfns$d120(b, ld)[, 1:3, drop = FALSE])
    numeric <- vapply(1:3, function(k) {
      h <- 1e-6 * max(1, abs(b[k]))
      bp <- bm <- b; bp[k] <- b[k] + h; bm[k] <- b[k] - h
      (f(bp) - f(bm)) / (2 * h)
    }, 0)
    expect_equal(analytic, numeric, tolerance = 1e-5,
                 info = sprintf("xi = %.0e", xi))
  }
})

test_that("a log-link fit reports the log-likelihood its own density implies", {
  ## the end-to-end symptom: with smooths the fitted xi is often driven to ~0, and
  ## $logLik then disagreed with sum(dDEGPD1(...)) at the fitted parameters
  set.seed(7)
  n <- 300
  d <- data.frame(t = seq(0, 10, length.out = n))
  d$season <- (seq_len(n) %% 52) / 52
  d$cases <- rnbinom(n, size = 2, mu = 20 + 15 * sin(2 * pi * d$season))

  fit <- egpd(list(lsigma = cases ~ s(t, k = 10) + s(season, bs = "cc", k = 6),
                   lxi = ~1, lkappa = ~1),
              data = d, family = "degpd", degpd.args = list(m = 1),
              restarts = FALSE, trace = 0, knots = list(season = c(0, 1)))

  pa <- as.data.frame(predict(fit, type = "response"))
  manual <- sum(dDEGPD1(d$cases, mu = pa$scale, sigma = pa$shape, nu = pa$kappa,
                        log = TRUE))
  expect_equal(as.numeric(fit$logLik), manual, tolerance = 1e-4)
})
