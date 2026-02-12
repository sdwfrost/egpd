## Tests for bivariate discrete EGPD (BDEGPD / BZIDEGPD)

## ---- Pure R tests (no Julia needed) ----

test_that("rbdegpd returns correct dimensions", {
  set.seed(42)
  Y <- rbdegpd(100, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  expect_true(is.matrix(Y))
  expect_equal(nrow(Y), 100)
  expect_equal(ncol(Y), 2)
  expect_equal(colnames(Y), c("Y1", "Y2"))
})

test_that("rbdegpd produces non-negative integers", {
  set.seed(42)
  Y <- rbdegpd(500, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  expect_true(all(is.finite(Y)))
  expect_true(all(Y == floor(Y)))
  expect_true(is.integer(Y))
})

test_that("rbdegpd is reproducible with set.seed", {
  set.seed(1)
  Y1 <- rbdegpd(50, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  set.seed(1)
  Y2 <- rbdegpd(50, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  expect_identical(Y1, Y2)
})

test_that("rbdegpd validates inputs via rbegpd", {
  expect_error(rbdegpd(-1, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2),
               "'n' must be a positive integer")
  expect_error(rbdegpd(100, kappa = -1, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2),
               "'kappa' must be a positive number")
})

test_that("rbzidegpd returns correct dimensions", {
  set.seed(42)
  Y <- rbzidegpd(100, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2, pi0 = 0.3)
  expect_true(is.matrix(Y))
  expect_equal(nrow(Y), 100)
  expect_equal(ncol(Y), 2)
  expect_equal(colnames(Y), c("Y1", "Y2"))
})

test_that("rbzidegpd produces non-negative integers", {
  set.seed(42)
  Y <- rbzidegpd(500, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2, pi0 = 0.3)
  expect_true(all(is.finite(Y)))
  expect_true(all(Y == floor(Y)))
  expect_true(is.integer(Y))
})

test_that("rbzidegpd has correct zero-inflation proportion", {
  set.seed(42)
  n <- 10000
  pi0 <- 0.4
  Y <- rbzidegpd(n, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2, pi0 = pi0)
  ## Both columns zero simultaneously
  zero_rows <- Y[, 1] == 0 & Y[, 2] == 0
  ## Zero proportion should be at least pi0 (some natural zeros too)
  expect_true(mean(zero_rows) >= pi0 - 0.05)
  ## But close to pi0 + natural zero rate (allow generous tolerance)
  expect_true(mean(zero_rows) < 1)
})

test_that("rbzidegpd validates pi0", {
  expect_error(rbzidegpd(100, kappa = 2, sigma = 1, xi = 0.1,
                          thL = 5, thU = 5, thw = 0.2, pi0 = 0),
               "'pi0' must be a number in \\(0, 1\\)")
  expect_error(rbzidegpd(100, kappa = 2, sigma = 1, xi = 0.1,
                          thL = 5, thU = 5, thw = 0.2, pi0 = 1),
               "'pi0' must be a number in \\(0, 1\\)")
  expect_error(rbzidegpd(100, kappa = 2, sigma = 1, xi = 0.1,
                          thL = 5, thU = 5, thw = 0.2, pi0 = -0.1),
               "'pi0' must be a number in \\(0, 1\\)")
  expect_error(rbzidegpd(100, kappa = 2, sigma = 1, xi = 0.1,
                          thL = 5, thU = 5, thw = 0.2, pi0 = 1.5),
               "'pi0' must be a number in \\(0, 1\\)")
})

test_that("rbzidegpd is reproducible with set.seed", {
  set.seed(1)
  Y1 <- rbzidegpd(50, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2, pi0 = 0.3)
  set.seed(1)
  Y2 <- rbzidegpd(50, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2, pi0 = 0.3)
  expect_identical(Y1, Y2)
})

test_that("fitegpd rejects bdegpd with method='mle'", {
  Y <- matrix(as.integer(rpois(200, 5)), ncol = 2)
  expect_error(fitegpd(Y, family = "bdegpd", method = "mle"),
               "family='bdegpd' requires method='neuralbayes'")
})

test_that("fitegpd rejects bzidegpd with method='mle'", {
  Y <- matrix(as.integer(rpois(200, 5)), ncol = 2)
  expect_error(fitegpd(Y, family = "bzidegpd", method = "mle"),
               "family='bzidegpd' requires method='neuralbayes'")
})


## ---- Julia-dependent tests ----

test_that("NPE inference for bdegpd produces valid fitegpd object", {
  skip_if_not_installed("JuliaConnectoR")
  skip_if_not_installed("NeuralEstimators")
  skip_if(system.file("models", "BDEGPD_NPE.bson", package = "egpd") == "",
          "Bundled BDEGPD_NPE.bson model not found")

  set.seed(42)
  Y <- rbdegpd(500, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  fit <- fitegpd(Y, family = "bdegpd", method = "neuralbayes", estimator = "npe")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$family, "bdegpd")
  expect_equal(fit$method, "neuralbayes")
  expect_equal(fit$estimator_type, "npe")
  expect_equal(fit$npar, 6)
  expect_equal(length(fit$estimate), 6)
  expect_true(all(c("kappa", "sigma", "xi", "thL", "thU", "thw") %in% names(fit$estimate)))
  expect_true(all(fit$estimate > 0))
  expect_true(!is.null(fit$posterior_samples))
  expect_equal(nrow(fit$posterior_samples), 6)
})

test_that("NBE inference for bdegpd works", {
  skip_if_not_installed("JuliaConnectoR")
  skip_if_not_installed("NeuralEstimators")
  skip_if(system.file("models", "BDEGPD_NBE.bson", package = "egpd") == "",
          "Bundled BDEGPD_NBE.bson model not found")

  set.seed(42)
  Y <- rbdegpd(500, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  fit <- fitegpd(Y, family = "bdegpd", method = "neuralbayes", estimator = "nbe")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$estimator_type, "nbe")
  expect_null(fit$posterior_samples)
  expect_equal(length(fit$estimate), 6)
  expect_true(all(fit$estimate > 0))
})

test_that("NPE inference for bzidegpd produces valid fitegpd object with 7 params", {
  skip_if_not_installed("JuliaConnectoR")
  skip_if_not_installed("NeuralEstimators")
  skip_if(system.file("models", "BZIDEGPD_NPE.bson", package = "egpd") == "",
          "Bundled BZIDEGPD_NPE.bson model not found")

  set.seed(42)
  Y <- rbzidegpd(500, kappa = 2, sigma = 1, xi = 0.1,
                  thL = 5, thU = 5, thw = 0.2, pi0 = 0.3)
  fit <- fitegpd(Y, family = "bzidegpd", method = "neuralbayes", estimator = "npe")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$family, "bzidegpd")
  expect_equal(fit$npar, 7)
  expect_equal(length(fit$estimate), 7)
  expect_true("pi0" %in% names(fit$estimate))
  expect_true(fit$estimate["pi0"] > 0 && fit$estimate["pi0"] < 1)
  expect_true(!is.null(fit$posterior_samples))
  expect_equal(nrow(fit$posterior_samples), 7)
})

test_that("S3 methods work for bdegpd NPE fit", {
  skip_if_not_installed("JuliaConnectoR")
  skip_if_not_installed("NeuralEstimators")
  skip_if(system.file("models", "BDEGPD_NPE.bson", package = "egpd") == "",
          "Bundled BDEGPD_NPE.bson model not found")

  set.seed(42)
  Y <- rbdegpd(500, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  fit <- fitegpd(Y, family = "bdegpd", method = "neuralbayes", estimator = "npe")

  ## print
  expect_output(print(fit), "Experimental")

  ## summary
  s <- summary(fit)
  expect_s3_class(s, "summary.fitegpd")
  expect_output(print(s), "Posterior summary")

  ## coef
  expect_named(coef(fit), c("kappa", "sigma", "xi", "thL", "thU", "thw"))

  ## vcov
  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 6)

  ## confint
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
