## distfit() and an intercept-only egpd() GAM are two independent routes to the
## same MLE: distfit() runs Nelder-Mead on a transformed parameter vector, while
## egpd() runs penalised Newton steps inside a REML outer loop. With no covariates
## and no smooths there is nothing to penalise, so the two must agree -- on the
## parameters, and on the likelihood those parameters imply.
##
## The likelihood assertions here are deliberately made against a hand-computed
## sum over the family's own density rather than against the two objects' stored
## values. That is what makes them a regression test for the logLik constant bug:
## egpd() used to subtract the REML normalising constant .5 * Mp * log(2 * pi)
## from $logLik, which left the parameters correct but the reported likelihood
## short by ~0.92 per unpenalised coefficient (and by more as Mp grew with model
## complexity, so it did not cancel in an AIC comparison).

## predict(type = "response") appends an `xi` alias to the DEGPD parameters, so
## take the leading npar columns rather than everything returned.
gam_pars <- function(fit, npar) {
  p <- vapply(predict(fit, type = "response"), function(z) z[1], numeric(1))
  unname(p[seq_len(npar)])
}

fit_gam <- function(y, family, ...) suppressWarnings(suppressMessages(
  egpd(list(y ~ 1, ~ 1, ~ 1), data = data.frame(y = y), family = family,
       trace = -1, ...)))

# ---- DEGPD (model 1, identity link on xi) -----------------------------------

test_that("distfit and an intercept-only GAM agree on DEGPD, xi > 0", {
  skip_on_cran()
  set.seed(1)
  y <- rDEGPD1(2000, mu = 1.2, sigma = 0.25, nu = 1.5)

  fd <- distfit(y, type = 1, family = "degpd")
  fg <- fit_gam(y, "degpd", degpd.args = list(m = 1, link = "identity"),
                restarts = FALSE)

  pd <- unname(fd$estimate[c("sigma", "xi", "kappa")])
  pg <- gam_pars(fg, 3)
  expect_equal(pd, pg, tolerance = 5e-3)

  ## and the two parameter vectors imply the same likelihood
  ll <- function(p) sum(dDEGPD1(y, mu = p[1], sigma = p[2], nu = p[3], log = TRUE))
  expect_equal(ll(pd), ll(pg), tolerance = 1e-6)
})

test_that("distfit and an intercept-only GAM agree on DEGPD, xi < 0", {
  skip_on_cran()
  set.seed(1)
  y <- rDEGPD1(2000, mu = 1.2, sigma = -0.15, nu = 1.5)

  fd <- distfit(y, type = 1, family = "degpd")
  fg <- fit_gam(y, "degpd", degpd.args = list(m = 1, link = "identity"),
                restarts = FALSE)

  pd <- unname(fd$estimate[c("sigma", "xi", "kappa")])
  pg <- gam_pars(fg, 3)
  ## the bounded-tail case is the looser of the two: sigma and xi trade off more
  ## strongly when the support is finite
  expect_equal(pd, pg, tolerance = 5e-3)

  ll <- function(p) sum(dDEGPD1(y, mu = p[1], sigma = p[2], nu = p[3], log = TRUE))
  expect_equal(ll(pd), ll(pg), tolerance = 1e-6)
})

# ---- generalized Waring ------------------------------------------------------

test_that("distfit and an intercept-only GAM agree on the generalized Waring", {
  skip_on_cran()
  set.seed(1)
  y <- rgwaring(2000, mu = 4, k = 2, rho = 3)

  fd <- distfit(y, family = "gwaring")
  fg <- fit_gam(y, "gwaring", restarts = FALSE)

  pd <- unname(fd$estimate[c("mu", "k", "rho")])
  pg <- gam_pars(fg, 3)
  expect_equal(pd, pg, tolerance = 1e-2)

  ll <- function(p) sum(dgwaring(y, mu = p[1], k = p[2], rho = p[3], log = TRUE))
  expect_equal(ll(pd), ll(pg), tolerance = 1e-6)
})

test_that("neither optimiser is materially beaten on the generalized Waring", {
  skip_on_cran()
  ## k and rho trade off along a flat ridge, so the two optimisers do not always
  ## stop at the same point even though both are consistent. Across seeds the gap
  ## is one-sided: distfit's Nelder-Mead never lands materially below the GAM's
  ## Newton steps, while the reverse does happen -- egpd()'s multi-start covers
  ## DEGPD model 1 only, so the Waring fit gets a single start. Seed 101 is the
  ## worst case seen (~0.5 nats over 2000 observations); this test pins the
  ## direction and magnitude so a regression either way is visible.
  ll_of <- function(y, p) sum(dgwaring(y, mu = p[1], k = p[2], rho = p[3], log = TRUE))
  gaps <- vapply(c(1, 7, 42, 101, 2024), function(s) {
    set.seed(s)
    y <- rgwaring(2000, mu = 4, k = 2, rho = 3)
    fd <- distfit(y, family = "gwaring")
    fg <- fit_gam(y, "gwaring", restarts = FALSE)
    ll_of(y, unname(fd$estimate[c("mu", "k", "rho")])) - ll_of(y, gam_pars(fg, 3))
  }, numeric(1))

  ## the GAM never beats distfit by more than round-off
  expect_true(all(gaps > -1e-3))
  ## and distfit's advantage stays within the known ridge tolerance
  expect_true(all(gaps < 1))
})

# ---- the reported likelihood is the plain log-likelihood ---------------------

test_that("egpd() reports an unpenalised logLik with no REML constant", {
  skip_on_cran()
  set.seed(1)
  y <- rDEGPD1(2000, mu = 1.2, sigma = 0.25, nu = 1.5)
  fg <- fit_gam(y, "degpd", degpd.args = list(m = 1, link = "identity"),
                restarts = FALSE)
  p <- gam_pars(fg, 3)
  manual <- sum(dDEGPD1(y, mu = p[1], sigma = p[2], nu = p[3], log = TRUE))

  expect_equal(as.numeric(logLik(fg)), manual, tolerance = 1e-6)
  expect_equal(fg$AIC, -2 * manual + 2 * fg$df, tolerance = 1e-6)
  expect_equal(fg$BIC, -2 * manual + log(length(y)) * fg$df, tolerance = 1e-6)
})

test_that("distfit and egpd() report the same likelihood on the same fit", {
  skip_on_cran()
  set.seed(1)
  y <- rgwaring(2000, mu = 4, k = 2, rho = 3)
  fd <- distfit(y, family = "gwaring")
  fg <- fit_gam(y, "gwaring", restarts = FALSE)

  ## each object's stored value must equal a hand computation from the density,
  ## which is what makes the two directly comparable
  ld <- sum(dgwaring(y, mu = fd$estimate[["mu"]], k = fd$estimate[["k"]],
                     rho = fd$estimate[["rho"]], log = TRUE))
  pg <- gam_pars(fg, 3)
  lg <- sum(dgwaring(y, mu = pg[1], k = pg[2], rho = pg[3], log = TRUE))

  expect_equal(as.numeric(logLik(fd)), ld, tolerance = 1e-6)
  expect_equal(as.numeric(logLik(fg)), lg, tolerance = 1e-6)
  expect_equal(as.numeric(logLik(fd)), as.numeric(logLik(fg)), tolerance = 1e-3)
})

test_that("the REML constant does not creep back in as models gain smooths", {
  skip_on_cran()
  ## Mp counts one per linear predictor plus each smooth's null-space dimension,
  ## so a model with a smooth is exactly where the old offset was largest. Fitting
  ## one and checking against the density catches a reintroduction that an
  ## intercept-only test could miss.
  set.seed(3)
  n <- 1200; tt <- seq(0, 1, length.out = n)
  y <- rDEGPD1(n, mu = exp(0.2 + sin(2 * pi * tt)), sigma = 0.2, nu = 1.4)
  fg <- suppressWarnings(suppressMessages(
    egpd(list(y ~ s(t, k = 8), ~ 1, ~ 1), data = data.frame(y = y, t = tt),
         family = "degpd", degpd.args = list(m = 1, link = "identity"),
         trace = -1, restarts = FALSE)))
  p <- predict(fg, type = "response")
  manual <- sum(dDEGPD1(y, mu = p[[1]], sigma = p[[2]], nu = p[[3]], log = TRUE))

  expect_equal(as.numeric(logLik(fg)), manual, tolerance = 1e-6)
  ## df exceeds the number of linear predictors, confirming the smooth is active
  ## and so Mp > 3 -- the case the old code got most wrong
  expect_gt(fg$df, 3)
})
