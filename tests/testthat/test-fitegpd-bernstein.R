## Tests for fitegpd() Bernstein polynomial fitting

test_that("fitegpd Bernstein method converges on type 1 data", {
  set.seed(42)
  x <- regpd(500, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
  fit <- fitegpd(x, type = 1, method = "bernstein", bernstein.m = 8)

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$convergence, 0)
  expect_equal(fit$method, "bernstein")
  expect_equal(fit$family, "egpd")
  expect_equal(fit$bernstein.m, 8)
  expect_true(!is.null(fit$bernstein.weights))
  expect_equal(length(fit$bernstein.weights), 8)
  expect_true(abs(sum(fit$bernstein.weights) - 1) < 1e-10)
  expect_true(all(fit$bernstein.weights >= 0))
})

test_that("fitegpd Bernstein weights sum to 1", {
  set.seed(123)
  x <- regpd(300, sigma = 1.5, xi = 0.2, kappa = 2, type = 1)
  fit <- fitegpd(x, type = 1, method = "bernstein", bernstein.m = 6)

  expect_equal(sum(fit$bernstein.weights), 1, tolerance = 1e-10)
})

test_that("fitegpd Bernstein S3 methods work", {
  set.seed(42)
  x <- regpd(300, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
  fit <- fitegpd(x, type = 1, method = "bernstein", bernstein.m = 8)

  ## coef
  expect_named(coef(fit), c("sigma", "xi", "kappa"))

  ## logLik
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")

  ## npar should account for Bernstein degrees of freedom
  ## 3 params (sigma, xi, kappa) + (m-1) Bernstein weights = 3 + 7 = 10
  expect_equal(fit$npar, 10)

  ## AIC
  expect_true(is.numeric(AIC(fit)))

  ## summary
  s <- summary(fit)
  expect_s3_class(s, "summary.fitegpd")
  expect_output(print(s))
})

test_that("fitegpd Bernstein plot does not error", {
  set.seed(42)
  x <- regpd(200, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
  fit <- fitegpd(x, type = 1, method = "bernstein", bernstein.m = 6)

  expect_no_error(plot(fit))
})

test_that("fitegpd Bernstein with different degrees converges", {
  set.seed(42)
  x <- regpd(300, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)

  fit4 <- fitegpd(x, type = 1, method = "bernstein", bernstein.m = 4)
  fit12 <- fitegpd(x, type = 1, method = "bernstein", bernstein.m = 12)

  expect_equal(fit4$convergence, 0)
  expect_equal(fit12$convergence, 0)
  expect_equal(fit4$bernstein.m, 4)
  expect_equal(fit12$bernstein.m, 12)
})

test_that("fitegpd Bernstein rejects non-egpd family", {
  x <- c(1, 2, 3, 4, 5)
  expect_error(
    fitegpd(x, type = 1, method = "bernstein", family = "degpd"),
    "Bernstein method is only available"
  )
})

test_that("AIC comparison between MLE and Bernstein works", {
  set.seed(42)
  x <- regpd(300, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)

  fit_mle <- fitegpd(x, type = 1, method = "mle")
  fit_bern <- fitegpd(x, type = 1, method = "bernstein", bernstein.m = 8)

  ## Both should produce finite AIC
  expect_true(is.finite(AIC(fit_mle)))
  expect_true(is.finite(AIC(fit_bern)))
})
