## Regression tests for fitting-level fixes and features: logLik/AIC population,
## the predict() `xi` alias, offset-based parameter fixing, ZIDEGPD<-DEGPD nesting,
## and the bounded xi.max shape link.

kn <- list(season = c(0.5, 53.5))
sform <- cases ~ s(season, bs = "cc", k = 6)

fit_degpd <- function(d, m = 1, ...) suppressWarnings(suppressMessages(
  egpd(list(lsigma = sform, lxi = ~ 1, lkappa = ~ 1),
       data = d, family = "degpd", degpd.args = list(m = m, ...), knots = kn)))

mkdata <- function(seed = 1, n = 200, size = 2) {
  set.seed(seed); ss <- ((seq_len(n) - 1) %% 52) + 1
  data.frame(cases = rnbinom(n, mu = exp(1 + 0.6 * cos(2 * pi * ss / 52)), size = size),
             season = ss)
}

test_that("logLik / AIC / BIC are finite on a converged DEGPD fit", {
  f <- fit_degpd(mkdata())
  expect_true(is.finite(f$logLik))
  expect_true(is.finite(AIC(f)))
  expect_true(is.finite(BIC(f)))
})

test_that("predict(type='response') exposes an `xi` alias equal to `shape`", {
  d <- mkdata(); f <- fit_degpd(d)
  pr <- as.data.frame(predict(f, d[1, , drop = FALSE], type = "response"))
  expect_true("xi" %in% names(pr))
  expect_equal(pr$xi, pr$shape)
})

test_that("intercept-free offset pins the shape parameter", {
  d <- mkdata(); d$zoff <- log(0.15)
  f <- suppressWarnings(suppressMessages(egpd(
    list(lsigma = sform, lxi = ~ offset(zoff) - 1, lkappa = ~ 1),
    data = d, family = "degpd", degpd.args = list(m = 1), knots = kn)))
  xi <- as.data.frame(predict(f, d[1, , drop = FALSE], type = "response"))$shape[1]
  expect_equal(xi, 0.15, tolerance = 1e-6)
})

test_that("ZIDEGPD warm-start respects nesting (logLik >= nested DEGPD)", {
  ## Intercept-only models: no smoothing parameters, so ZIDEGPD strictly nests DEGPD
  ## (pi -> 0). The warm-start guarantees the richer fit cannot sit below the simpler.
  d <- mkdata(seed = 2); d$cases[sample(nrow(d), 40)] <- 0
  fd <- suppressWarnings(suppressMessages(egpd(
    list(cases ~ 1, ~ 1, ~ 1), data = d, family = "degpd", degpd.args = list(m = 1))))
  fz <- suppressWarnings(suppressMessages(egpd(
    list(cases ~ 1, ~ 1, ~ 1, ~ 1), data = d, family = "zidegpd", zidegpd.args = list(m = 1))))
  expect_true(is.finite(fz$logLik))
  expect_gte(fz$logLik, fd$logLik - 1e-4)
})

test_that("bounded xi.max caps the estimated tail index", {
  set.seed(3); n <- 80; ss <- ((seq_len(n) - 1) %% 52) + 1
  y <- rpois(n, 1); y[c(15, 40, 60)] <- c(120, 260, 480)   # heavy upper tail
  d <- data.frame(cases = y, season = ss)
  f <- fit_degpd(d, m = 1, xi.max = 0.5)
  xi <- as.data.frame(predict(f, d[1, , drop = FALSE], type = "response"))$shape[1]
  expect_lte(xi, 0.5 + 1e-6)
})
