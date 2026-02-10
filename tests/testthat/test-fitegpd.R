## Tests for fitegpd() MLE fitting

test_that("fitegpd type 1 continuous EGPD recovers parameters", {
  set.seed(42)
  true_sigma <- 2; true_xi <- 0.1; true_kappa <- 1.5
  x <- regpd(500, sigma = true_sigma, xi = true_xi, kappa = true_kappa, type = 1)
  fit <- fitegpd(x, type = 1, family = "egpd")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$convergence, 0)
  expect_equal(fit$family, "egpd")
  expect_equal(fit$type, 1)
  expect_equal(fit$method, "mle")
  expect_equal(fit$n, 500)
  expect_equal(fit$npar, 3)
  expect_true(all(c("sigma", "xi", "kappa") %in% names(fit$estimate)))

  ## Parameter recovery: within 3 SEs of true values
  expect_true(abs(fit$estimate["sigma"] - true_sigma) < 3 * fit$sd["sigma"])
  expect_true(abs(fit$estimate["kappa"] - true_kappa) < 3 * fit$sd["kappa"])
})

test_that("fitegpd type 4 continuous EGPD converges", {
  set.seed(123)
  x <- regpd(300, sigma = 1.5, xi = 0.2, delta = 1.2, type = 4)
  fit <- fitegpd(x, type = 4, family = "egpd")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$convergence, 0)
  expect_true(all(c("sigma", "xi", "delta") %in% names(fit$estimate)))
  expect_equal(fit$npar, 3)
})

test_that("fitegpd type 5 continuous EGPD converges", {
  set.seed(456)
  x <- regpd(300, sigma = 2, xi = 0.1, delta = 1, kappa = 2, type = 5)
  fit <- fitegpd(x, type = 5, family = "egpd")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$convergence, 0)
  expect_true(all(c("sigma", "xi", "delta", "kappa") %in% names(fit$estimate)))
  expect_equal(fit$npar, 4)
})

test_that("fitegpd discrete EGPD type 1 converges", {
  set.seed(789)
  x <- rdiscegpd(500, sigma = 3, xi = 0.1, kappa = 1.5, type = 1)
  fit <- fitegpd(x, type = 1, family = "degpd")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$convergence, 0)
  expect_equal(fit$family, "degpd")
})

test_that("fitegpd zero-inflated discrete EGPD type 1 converges", {
  set.seed(101)
  x <- rzidiscegpd(500, pi = 0.3, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
  fit <- fitegpd(x, type = 1, family = "zidegpd")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$convergence, 0)
  expect_equal(fit$family, "zidegpd")
  expect_true("pi" %in% names(fit$estimate))
})

test_that("fitegpd zero-inflated continuous EGPD type 1 converges", {
  set.seed(202)
  x <- rziegpd(500, pi = 0.2, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
  fit <- fitegpd(x, type = 1, family = "ziegpd")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$convergence, 0)
  expect_equal(fit$family, "ziegpd")
  expect_true("pi" %in% names(fit$estimate))
})

test_that("fitegpd fix.arg works correctly", {
  set.seed(42)
  x <- regpd(300, sigma = 2, xi = 0.05, kappa = 1.5, type = 1)
  fit <- fitegpd(x, type = 1, fix.arg = list(xi = 0.05))

  expect_equal(fit$convergence, 0)
  expect_false("xi" %in% names(fit$estimate))
  expect_equal(fit$npar, 2)
  expect_equal(fit$fix.arg, list(xi = 0.05))
})

test_that("fitegpd S3 methods work", {
  set.seed(42)
  x <- regpd(200, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
  fit <- fitegpd(x, type = 1)

  ## coef
  expect_named(coef(fit), c("sigma", "xi", "kappa"))

  ## vcov
  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 3)

  ## logLik
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "df"), 3)

  ## nobs
  expect_equal(nobs(fit), 200)

  ## AIC (uses logLik method)
  expect_true(is.numeric(AIC(fit)))

  ## confint
  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 3)
  expect_equal(ncol(ci), 2)

  ## summary
  s <- summary(fit)
  expect_s3_class(s, "summary.fitegpd")

  ## print (should not error)
  expect_output(print(fit))
  expect_output(print(s))
})

test_that("fitegpd plot does not error", {
  set.seed(42)
  x <- regpd(200, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
  fit <- fitegpd(x, type = 1)

  expect_no_error(plot(fit))
})

test_that("fitegpd input validation works", {
  expect_error(fitegpd(1, type = 1), "'x' must be a numeric")
  expect_error(fitegpd(c(1, 2, 3), type = 7), "'type' must be an integer")
  expect_error(fitegpd(c(1, 2, 3), type = 1, method = "bernstein", family = "degpd"),
               "Bernstein method is only available")
  expect_error(fitegpd(c(1, 2, 3), type = 1, fix.arg = list(foo = 1)),
               "Unknown fixed parameters")
})

test_that("fitegpd with custom start values works", {
  set.seed(42)
  x <- regpd(300, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
  fit <- fitegpd(x, type = 1, start = list(sigma = 1.5, kappa = 1.2))

  expect_equal(fit$convergence, 0)
})

test_that("fitegpd type 6 continuous EGPD converges", {
  set.seed(999)
  x <- regpd(500, sigma = 2, xi = 0.1, kappa = 1.5, delta = 2, prob = 0.4, type = 6)
  fit <- fitegpd(x, type = 6, family = "egpd")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$convergence, 0)
  expect_true(all(c("sigma", "xi", "kappa", "delta", "prob") %in% names(fit$estimate)))
  expect_equal(fit$npar, 5)
})
