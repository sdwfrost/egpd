## Round-7 audit regressions (comppareto): a CompPareto fit whose default start walked out
## onto the plateau beyond the data (every observation below the fitted splice point) is
## refitted from other starts. Expected log-likelihoods are the maxima of an independent
## implementation of the documented law (stats::dlnorm/plnorm only), found by Nelder-Mead
## and BFGS from 38 starts (theta at 14 quantiles from 0.3 to 0.97 with alpha at the Lomax
## tail estimate above it, and theta at 6 quantiles with alpha 0.5, 1, 2 or 5) and polished
## over the splice points at the neighbouring observations, where the likelihood has a kink.

## log-likelihood of the composite lognormal-Pareto law at (meanlog, log sdlog, log alpha,
## log theta): body g(y) / ((1 + phi) G(theta)) for y <= theta, tail
## phi / (1 + phi) alpha (2 theta)^alpha / (y + theta)^(alpha + 1) above it, with
## phi = 2 theta g(theta) / (alpha G(theta))
cp_lnorm_ll <- function(par, y) {
  ml <- par[1]; sl <- exp(par[2]); al <- exp(par[3]); th <- exp(par[4])
  lG <- stats::plnorm(th, ml, sl, log.p = TRUE)
  lphi <- log(2 * th) + stats::dlnorm(th, ml, sl, log = TRUE) - log(al) - lG
  l1p <- log1p(exp(lphi))
  body <- y <= th
  sum(stats::dlnorm(y[body], ml, sl, log = TRUE) - l1p - lG) +
    sum(lphi - l1p + log(al) + al * log(2 * th) - (al + 1) * log(y[!body] + th))
}

test_that("r7: a CompPareto fit stranded beyond the data is refitted to the maximum", {
  skip_on_cran()
  set.seed(11)
  yc <- regpd(400, sigma = 2, xi = 0.15, kappa = 1.5)
  d <- data.frame(y = yc)
  fit <- suppressWarnings(egpd(list(y ~ 1, ~1, ~1, ~1), d, family = "comppareto"))
  b <- unname(fit$coefficients)
  expect_equal(fit$convergence, 0)
  ## the independent maximum is at theta = 3.3889, an observation (a kink of the
  ## likelihood), which the coefficient search approaches to within 1.4e-4
  expect_equal(as.numeric(logLik(fit)), -804.280573579216, tolerance = 1e-6)
  expect_equal(cp_lnorm_ll(b, yc), as.numeric(logLik(fit)), tolerance = 1e-12)
  ## both sides of the splice point are observed (113 or 114 observations above it, by the
  ## side of the kink at 3.3889 where the search stops), and the Hessian has full rank
  nt <- sum(yc > exp(b[4]))
  expect_gte(nt, 100)
  expect_lte(nt, 130)
  expect_equal(unname(fit$hessian.rank), c(4, 4))
  ## the default start alone ends on the plateau, 8.73 below: theta = 433 > max(y) = 22.19
  f0 <- suppressWarnings(egpd(list(y ~ 1, ~1, ~1, ~1), d, family = "comppareto",
                              restarts = FALSE))
  expect_equal(sum(yc > exp(f0$coefficients[4])), 0)
  expect_equal(as.numeric(logLik(f0)), cp_lnorm_ll(unname(f0$coefficients), yc), tolerance = 1e-12)
  expect_lt(as.numeric(logLik(f0)), -813)
})

test_that("r7: CompPareto restarts improve a stranded fit and retain an ordinary good fit", {
  ## a second stranded default fit (log-likelihood -320.37, no observation above theta)
  set.seed(11504)
  y <- regpd(150, sigma = 2, xi = 0.15, kappa = 1.5)
  fit <- suppressWarnings(egpd(list(y ~ 1, ~1, ~1, ~1), data.frame(y = y), family = "comppareto"))
  expect_equal(fit$convergence, 0)
  expect_equal(as.numeric(logLik(fit)), -311.168396345415, tolerance = 3e-6)
  expect_equal(cp_lnorm_ll(unname(fit$coefficients), y), as.numeric(logLik(fit)), tolerance = 1e-12)
  expect_gt(sum(y > exp(fit$coefficients[4])), 0)
  ## a default fit that reaches the maximum is returned as it is, bit for bit
  set.seed(21)
  y1 <- rcomppareto(180, spec = "lnorm", meanlog = 0.2, sdlog = 0.35, alpha = 1.5, theta = 2)
  d1 <- data.frame(y = y1)
  f1 <- egpd(list(y ~ 1, ~1, ~1, ~1), d1, family = "comppareto")
  f0 <- egpd(list(y ~ 1, ~1, ~1, ~1), d1, family = "comppareto", restarts = FALSE)
  expect_equal(as.numeric(logLik(f1)), -377.913755112192, tolerance = 1e-9)
  expect_identical(f1$coefficients, f0$coefficients)
  expect_identical(f1$Vp, f0$Vp)
  expect_identical(f1$convergence.report, f0$convergence.report)
})

test_that("r7: CompPareto restarts also run when the response has a missing value", {
  ## the stranded seed-11504 data with one NA row: .setup.data() drops the row, and the
  ## restart must still run (it was skipped, returning the stranded -320.37 fit)
  set.seed(11504)
  y <- regpd(150, sigma = 2, xi = 0.15, kappa = 1.5)
  fit <- suppressWarnings(egpd(list(y ~ 1, ~1, ~1, ~1), data.frame(y = c(y, NA)),
                               family = "comppareto"))
  expect_equal(fit$nobs, 150)
  expect_equal(as.numeric(logLik(fit)), -311.168396345415, tolerance = 3e-6)
  expect_gt(sum(y > exp(fit$coefficients[4])), 0)
})
