## Tests for bivariate MDGPD (Aka-Kratz-Naveau) and ZIMDGPD

## ---- Pure R tests (no Julia needed) ----

test_that("rmdgpd returns correct dimensions", {
  set.seed(42)
  Y <- rmdgpd(100, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5)
  expect_true(is.matrix(Y))
  expect_equal(nrow(Y), 100)
  expect_equal(ncol(Y), 2)
  expect_equal(colnames(Y), c("Y1", "Y2"))
})

test_that("rmdgpd produces non-negative integers", {
  set.seed(42)
  Y <- rmdgpd(1000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5)
  expect_true(all(is.finite(Y)))
  expect_true(all(Y >= 0))
  expect_true(all(Y == floor(Y)))
  expect_true(is.integer(Y))
})

test_that("rmdgpd is reproducible with set.seed", {
  set.seed(1)
  Y1 <- rmdgpd(50, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5)
  set.seed(1)
  Y2 <- rmdgpd(50, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5)
  expect_identical(Y1, Y2)
})

test_that("rmdgpd validates inputs", {
  expect_error(rmdgpd(-1, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5),
               "'n' must be a positive integer")
  expect_error(rmdgpd(100, sigma = -1, xi = 0.2, lambda = 1, rho = 0.5),
               "'sigma' must be a positive number")
  expect_error(rmdgpd(100, sigma = 2, xi = -0.1, lambda = 1, rho = 0.5),
               "'xi' must be a non-negative number")
  expect_error(rmdgpd(100, sigma = 2, xi = 0.2, lambda = 0, rho = 0.5),
               "'lambda' must be a positive number")
  expect_error(rmdgpd(100, sigma = 2, xi = 0.2, lambda = 1, rho = -0.1),
               "'rho' must be a number in")
  expect_error(rmdgpd(100, sigma = 2, xi = 0.2, lambda = 1, rho = 1),
               "'rho' must be a number in")
})

test_that("rmdgpd with high rho gives strong positive dependence", {
  set.seed(42)
  Y <- rmdgpd(5000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.95)
  ## With high rho, Y1 and Y2 should be very similar
  cor_val <- cor(Y[, 1], Y[, 2])
  expect_true(cor_val > 0.5)
})

test_that("rmdgpd with rho = 0 gives weaker dependence", {
  set.seed(42)
  Y_high <- rmdgpd(5000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.95)
  Y_low  <- rmdgpd(5000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.0)
  ## Correlation should be higher with high rho
  expect_true(cor(Y_high[, 1], Y_high[, 2]) > cor(Y_low[, 1], Y_low[, 2]))
})

test_that("rmdgpd with xi = 0 works (geometric-like marginals)", {
  set.seed(42)
  Y <- rmdgpd(500, sigma = 2, xi = 0, lambda = 1, rho = 0.5)
  expect_true(all(is.finite(Y)))
  expect_true(all(Y >= 0))
  expect_true(is.integer(Y))
})

test_that("rmdgpd max component equals geometric", {
  ## The maximum of (N1, N2) in the standard MDGPD is Geom(1-e^{-1}).
  ## For large sigma and small xi, the non-standard transform is nearly

  ## linear, so we test that max(Y1, Y2) is approximately geometric.
  set.seed(42)
  n <- 10000
  Y <- rmdgpd(n, sigma = 1, xi = 0, lambda = 1, rho = 0.5)
  maxY <- pmax(Y[, 1], Y[, 2])
  ## Geometric mean = (1-p)/p where p = 1-e^{-1} ≈ 0.632
  ## So mean ≈ e^{-1}/(1-e^{-1}) ≈ 0.582
  ## But with sigma = 1 and xi = 0, M = floor(sigma*N) = floor(N) = N
  ## So max(M1, M2) = max(N1, N2) = G ~ Geom(1-e^{-1})
  expected_mean <- exp(-1) / (1 - exp(-1))
  expect_true(abs(mean(maxY) - expected_mean) < 0.1)
})

test_that("rzimdgpd returns correct dimensions", {
  set.seed(42)
  Y <- rzimdgpd(100, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = 0.3)
  expect_true(is.matrix(Y))
  expect_equal(nrow(Y), 100)
  expect_equal(ncol(Y), 2)
  expect_equal(colnames(Y), c("Y1", "Y2"))
})

test_that("rzimdgpd produces non-negative integers", {
  set.seed(42)
  Y <- rzimdgpd(500, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = 0.3)
  expect_true(all(is.finite(Y)))
  expect_true(all(Y >= 0))
  expect_true(all(Y == floor(Y)))
  expect_true(is.integer(Y))
})

test_that("rzimdgpd has correct zero-inflation proportion", {
  set.seed(42)
  n <- 10000
  pi0 <- 0.4
  Y <- rzimdgpd(n, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = pi0)
  zero_rows <- Y[, 1] == 0 & Y[, 2] == 0
  ## Zero proportion should be at least pi0 (some natural zeros too)
  expect_true(mean(zero_rows) >= pi0 - 0.05)
  expect_true(mean(zero_rows) < 1)
})

test_that("rzimdgpd validates pi0", {
  expect_error(rzimdgpd(100, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = 0),
               "'pi0' must be a number in \\(0, 1\\)")
  expect_error(rzimdgpd(100, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = 1),
               "'pi0' must be a number in \\(0, 1\\)")
  expect_error(rzimdgpd(100, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = -0.1),
               "'pi0' must be a number in \\(0, 1\\)")
})

test_that("rzimdgpd is reproducible with set.seed", {
  set.seed(1)
  Y1 <- rzimdgpd(50, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = 0.3)
  set.seed(1)
  Y2 <- rzimdgpd(50, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = 0.3)
  expect_identical(Y1, Y2)
})

## ---- d=3 (trivariate) tests ----

test_that("rmdgpd d=3 returns correct dimensions", {
  set.seed(42)
  Y <- rmdgpd(100, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, d = 3)
  expect_true(is.matrix(Y))
  expect_equal(nrow(Y), 100)
  expect_equal(ncol(Y), 3)
  expect_equal(colnames(Y), c("Y1", "Y2", "Y3"))
})

test_that("rmdgpd d=3 produces non-negative integers", {
  set.seed(42)
  Y <- rmdgpd(1000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, d = 3)
  expect_true(all(is.finite(Y)))
  expect_true(all(Y >= 0))
  expect_true(all(Y == floor(Y)))
  expect_true(is.integer(Y))
})

test_that("rmdgpd d=3 is reproducible with set.seed", {
  set.seed(1)
  Y1 <- rmdgpd(50, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, d = 3)
  set.seed(1)
  Y2 <- rmdgpd(50, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, d = 3)
  expect_identical(Y1, Y2)
})

test_that("rmdgpd d=3 with high rho gives strong positive dependence", {
  set.seed(42)
  Y <- rmdgpd(5000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.95, d = 3)
  ## All pairwise correlations should be positive
  expect_true(cor(Y[, 1], Y[, 2]) > 0.3)
  expect_true(cor(Y[, 1], Y[, 3]) > 0.3)
  expect_true(cor(Y[, 2], Y[, 3]) > 0.3)
})

test_that("rmdgpd d=3 dependence increases with rho", {
  set.seed(42)
  Y_high <- rmdgpd(5000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.95, d = 3)
  Y_low  <- rmdgpd(5000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.0, d = 3)
  ## Average pairwise correlation should be higher with high rho
  avg_cor_high <- mean(c(cor(Y_high[,1], Y_high[,2]),
                          cor(Y_high[,1], Y_high[,3]),
                          cor(Y_high[,2], Y_high[,3])))
  avg_cor_low  <- mean(c(cor(Y_low[,1], Y_low[,2]),
                          cor(Y_low[,1], Y_low[,3]),
                          cor(Y_low[,2], Y_low[,3])))
  expect_true(avg_cor_high > avg_cor_low)
})

test_that("rmdgpd d=3 with xi = 0 works", {
  set.seed(42)
  Y <- rmdgpd(500, sigma = 2, xi = 0, lambda = 1, rho = 0.5, d = 3)
  expect_true(all(is.finite(Y)))
  expect_true(all(Y >= 0))
  expect_true(is.integer(Y))
})

test_that("rmdgpd validates d parameter", {
  expect_error(rmdgpd(100, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, d = 1),
               "'d' must be an integer >= 2")
})

test_that("rmdgpd d=4 works (higher dimensions)", {
  set.seed(42)
  Y <- rmdgpd(100, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, d = 4)
  expect_true(is.matrix(Y))
  expect_equal(ncol(Y), 4)
  expect_equal(colnames(Y), c("Y1", "Y2", "Y3", "Y4"))
  expect_true(all(Y >= 0))
  expect_true(is.integer(Y))
})

test_that("rzimdgpd d=3 returns correct dimensions", {
  set.seed(42)
  Y <- rzimdgpd(100, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = 0.3, d = 3)
  expect_true(is.matrix(Y))
  expect_equal(nrow(Y), 100)
  expect_equal(ncol(Y), 3)
  expect_equal(colnames(Y), c("Y1", "Y2", "Y3"))
})

test_that("rzimdgpd d=3 has correct zero-inflation proportion", {
  set.seed(42)
  n <- 10000
  pi0 <- 0.4
  Y <- rzimdgpd(n, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = pi0, d = 3)
  zero_rows <- rowSums(Y) == 0
  ## Zero proportion should be at least pi0 (some natural zeros too)
  expect_true(mean(zero_rows) >= pi0 - 0.05)
  expect_true(mean(zero_rows) < 1)
})

test_that("rzimdgpd d=3 is reproducible with set.seed", {
  set.seed(1)
  Y1 <- rzimdgpd(50, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = 0.3, d = 3)
  set.seed(1)
  Y2 <- rzimdgpd(50, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = 0.3, d = 3)
  expect_identical(Y1, Y2)
})

test_that("fitegpd rejects mdgpd with method='mle'", {
  Y <- matrix(as.integer(rpois(200, 5)), ncol = 2)
  expect_error(fitegpd(Y, family = "mdgpd", method = "mle"),
               "family='mdgpd' requires method='neuralbayes'")
})

test_that("fitegpd rejects zimdgpd with method='mle'", {
  Y <- matrix(as.integer(rpois(200, 5)), ncol = 2)
  expect_error(fitegpd(Y, family = "zimdgpd", method = "mle"),
               "family='zimdgpd' requires method='neuralbayes'")
})


## ---- Julia-dependent tests ----

test_that("NPE inference for mdgpd produces valid fitegpd object", {
  skip_if_not_installed("JuliaConnectoR")
  skip_if_not_installed("NeuralEstimators")
  skip_if(system.file("models", "MDGPD_NPE.bson", package = "egpd") == "",
          "Bundled MDGPD_NPE.bson model not found")

  set.seed(42)
  Y <- rmdgpd(500, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5)
  fit <- fitegpd(Y, family = "mdgpd", method = "neuralbayes", estimator = "npe")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$family, "mdgpd")
  expect_equal(fit$method, "neuralbayes")
  expect_equal(fit$estimator_type, "npe")
  expect_equal(fit$npar, 4)
  expect_equal(length(fit$estimate), 4)
  expect_true(all(c("sigma", "xi", "lambda", "rho") %in% names(fit$estimate)))
  expect_true(fit$estimate["sigma"] > 0)
  expect_true(fit$estimate["xi"] > 0)
  expect_true(fit$estimate["lambda"] > 0)
  expect_true(fit$estimate["rho"] > 0 && fit$estimate["rho"] < 1)
  expect_true(!is.null(fit$posterior_samples))
  expect_equal(nrow(fit$posterior_samples), 4)
})

test_that("NBE inference for mdgpd works", {
  skip_if_not_installed("JuliaConnectoR")
  skip_if_not_installed("NeuralEstimators")
  skip_if(system.file("models", "MDGPD_NBE.bson", package = "egpd") == "",
          "Bundled MDGPD_NBE.bson model not found")

  set.seed(42)
  Y <- rmdgpd(500, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5)
  fit <- fitegpd(Y, family = "mdgpd", method = "neuralbayes", estimator = "nbe")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$estimator_type, "nbe")
  expect_null(fit$posterior_samples)
  expect_equal(length(fit$estimate), 4)
})

test_that("NPE inference for zimdgpd produces valid fitegpd object with 5 params", {
  skip_if_not_installed("JuliaConnectoR")
  skip_if_not_installed("NeuralEstimators")
  skip_if(system.file("models", "ZIMDGPD_NPE.bson", package = "egpd") == "",
          "Bundled ZIMDGPD_NPE.bson model not found")

  set.seed(42)
  Y <- rzimdgpd(500, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = 0.3)
  fit <- fitegpd(Y, family = "zimdgpd", method = "neuralbayes", estimator = "npe")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$family, "zimdgpd")
  expect_equal(fit$npar, 5)
  expect_equal(length(fit$estimate), 5)
  expect_true("pi0" %in% names(fit$estimate))
  expect_true(fit$estimate["pi0"] > 0 && fit$estimate["pi0"] < 1)
  expect_true(!is.null(fit$posterior_samples))
  expect_equal(nrow(fit$posterior_samples), 5)
})

test_that("S3 methods work for mdgpd NPE fit", {
  skip_if_not_installed("JuliaConnectoR")
  skip_if_not_installed("NeuralEstimators")
  skip_if(system.file("models", "MDGPD_NPE.bson", package = "egpd") == "",
          "Bundled MDGPD_NPE.bson model not found")

  set.seed(42)
  Y <- rmdgpd(500, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5)
  fit <- fitegpd(Y, family = "mdgpd", method = "neuralbayes", estimator = "npe")

  ## print
  expect_output(print(fit), "Experimental")

  ## summary
  s <- summary(fit)
  expect_s3_class(s, "summary.fitegpd")
  expect_output(print(s), "Posterior summary")

  ## coef
  expect_named(coef(fit), c("sigma", "xi", "lambda", "rho"))

  ## vcov
  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 4)

  ## confint
  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 4)
  expect_equal(ncol(ci), 2)

  ## logLik returns NA with warning
  expect_warning(logLik(fit), "not available")

  ## nobs
  expect_equal(nobs(fit), 500)

  ## plot
  expect_no_error(plot(fit))
})
