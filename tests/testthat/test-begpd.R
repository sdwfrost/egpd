## Tests for multivariate EGPD (BEGPD)

## ---- Pure R tests (no Julia needed) ----

test_that("rbegpd returns correct dimensions", {
  set.seed(42)
  Y <- rbegpd(100, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  expect_true(is.matrix(Y))
  expect_equal(nrow(Y), 100)
  expect_equal(ncol(Y), 2)
  expect_equal(colnames(Y), c("Y1", "Y2"))
})

test_that("rbegpd produces finite values", {
  set.seed(42)
  Y <- rbegpd(500, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  expect_true(all(is.finite(Y)))
})

test_that("rbegpd validates inputs", {
  expect_error(rbegpd(-1, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2),
               "'n' must be a positive integer")
  expect_error(rbegpd(100, kappa = -1, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2),
               "'kappa' must be a positive number")
  expect_error(rbegpd(100, kappa = 2, sigma = 0, xi = 0.1, thL = 5, thU = 5, thw = 0.2),
               "'sigma' must be a positive number")
  expect_error(rbegpd(100, kappa = 2, sigma = 1, xi = -0.1, thL = 5, thU = 5, thw = 0.2),
               "'xi' must be a positive number")
  expect_error(rbegpd(100, kappa = 2, sigma = 1, xi = 0.1, thL = -1, thU = 5, thw = 0.2),
               "'thL' must be a positive number")
  expect_error(rbegpd(100, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 0, thw = 0.2),
               "'thU' must be a positive number")
  expect_error(rbegpd(100, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.5),
               "'thw' must be a number in")
  expect_error(rbegpd(100, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0),
               "'thw' must be a number in")
})

test_that("rbegpd is reproducible with set.seed", {
  set.seed(1)
  Y1 <- rbegpd(50, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  set.seed(1)
  Y2 <- rbegpd(50, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  expect_identical(Y1, Y2)
})

test_that("fitegpd rejects begpd with method='mle'", {
  Y <- matrix(rnorm(200), ncol = 2)
  expect_error(fitegpd(Y, family = "begpd", method = "mle"),
               "family='begpd' requires method='neuralbayes'")
})

test_that("fitegpd rejects neuralbayes with non-multivariate family", {
  x <- rnorm(100)
  expect_error(fitegpd(x, family = "egpd", method = "neuralbayes"),
               "method='neuralbayes' requires a multivariate family")
})


## ---- Julia-dependent tests ----

test_that("NPE inference produces valid fitegpd object", {
  skip_if_not_installed("JuliaConnectoR")
  skip_if_not_installed("NeuralEstimators")
  skip_if(system.file("models", "NPE.bson", package = "egpd") == "",
          "Bundled NPE.bson model not found")

  set.seed(42)
  Y <- rbegpd(500, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  fit <- fitegpd(Y, family = "begpd", method = "neuralbayes", estimator = "npe")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$family, "begpd")
  expect_equal(fit$method, "neuralbayes")
  expect_equal(fit$estimator_type, "npe")
  expect_equal(fit$npar, 6)
  expect_equal(length(fit$estimate), 6)
  expect_true(all(c("kappa", "sigma", "xi", "thL", "thU", "thw") %in% names(fit$estimate)))
  expect_true(all(fit$estimate > 0))
  expect_true(!is.null(fit$posterior_samples))
  expect_equal(nrow(fit$posterior_samples), 6)
  expect_true(is.na(fit$loglik))
  expect_true(is.na(fit$aic))
})

test_that("NBE inference produces valid fitegpd object", {
  skip_if_not_installed("JuliaConnectoR")
  skip_if_not_installed("NeuralEstimators")
  skip_if(system.file("models", "NBE.bson", package = "egpd") == "",
          "Bundled NBE.bson model not found")

  set.seed(42)
  Y <- rbegpd(500, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  fit <- fitegpd(Y, family = "begpd", method = "neuralbayes", estimator = "nbe")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$estimator_type, "nbe")
  expect_null(fit$posterior_samples)
  expect_equal(length(fit$estimate), 6)
  expect_true(all(fit$estimate > 0))
})

test_that("S3 methods work for begpd NPE fit", {
  skip_if_not_installed("JuliaConnectoR")
  skip_if_not_installed("NeuralEstimators")
  skip_if(system.file("models", "NPE.bson", package = "egpd") == "",
          "Bundled NPE.bson model not found")

  set.seed(42)
  Y <- rbegpd(500, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  fit <- fitegpd(Y, family = "begpd", method = "neuralbayes", estimator = "npe")

  ## print
  expect_output(print(fit), "bivariate BEGPD")

  ## summary
  s <- summary(fit)
  expect_s3_class(s, "summary.fitegpd")
  expect_output(print(s), "Posterior summary")

  ## coef
  expect_named(coef(fit), c("kappa", "sigma", "xi", "thL", "thU", "thw"))

  ## vcov (computed from posterior samples)
  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 6)
  expect_equal(ncol(V), 6)

  ## confint (credible intervals from posterior)
  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 6)
  expect_equal(ncol(ci), 2)

  ## logLik returns NA with warning
  expect_warning(logLik(fit), "not available")

  ## nobs
  expect_equal(nobs(fit), 500)

  ## plot
  expect_no_error(plot(fit))
})
