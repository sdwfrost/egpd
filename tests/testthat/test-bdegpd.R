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
  ## Round 3: the previous bounds (>= pi0 - 0.05, < 1) were tautological.
  ## With sigma = 20 and kappa = 10, P(Y1 + Y2 < 2) = H(2)^kappa = 5.8e-11,
  ## so BDEGPD rows are never (0, 0) and the joint-zero rate must equal pi0.
  ## (The case with natural zeros is checked against an exact integral in
  ## test-audit-multivariate-r3.R.)
  set.seed(42)
  n <- 10000
  pi0 <- 0.4
  Y <- rbzidegpd(n, kappa = 10, sigma = 20, xi = 0.1, thL = 5, thU = 5,
                 thw = 0.2, pi0 = pi0)
  zero_rows <- Y[, 1] == 0 & Y[, 2] == 0
  expect_lt(abs(mean(zero_rows) - pi0), 5 * sqrt(pi0 * (1 - pi0) / n))
  expect_true(all(Y[!zero_rows, 1] + Y[!zero_rows, 2] >= 1))
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

test_that("distfit rejects bdegpd with method='mle'", {
  Y <- matrix(as.integer(rpois(200, 5)), ncol = 2)
  expect_error(distfit(Y, family = "bdegpd", method = "mle"),
               "family='bdegpd' requires method='neuralbayes'")
})

test_that("distfit rejects bzidegpd with method='mle'", {
  Y <- matrix(as.integer(rpois(200, 5)), ncol = 2)
  expect_error(distfit(Y, family = "bzidegpd", method = "mle"),
               "family='bzidegpd' requires method='neuralbayes'")
})


## ---- Julia-dependent tests ----
## No neural weights are bundled (inst/models/COMPATIBILITY.dcf: "Bundled-Weights:
## none"), so these tests train tiny networks (2 epochs on 200 parameter draws)
## once per file and family and run inference with them. They check loading,
## estimation and the S3 methods end to end with real weights, not estimation
## accuracy, and skip unless Julia can be started through JuliaConnectoR and loads
## the required Julia packages (the R-side neural code is tested with mocks
## elsewhere).

skip_if_no_julia <- function() {
  skip_on_cran()
  skip_if_not_installed("JuliaConnectoR")
  skip_if_not_installed("NeuralEstimators")
  skip_if_not(isTRUE(tryCatch(JuliaConnectoR::juliaSetupOk(), error = function(e) FALSE)),
              "Julia cannot be started by JuliaConnectoR")
  ## egpd and NeuralEstimators::train() load these with `using`, which needs each
  ## to be a direct dependency of the Julia environment
  pkgs <- c("NeuralEstimators", "Flux", "BSON", "Optimisers")
  found <- tryCatch(JuliaConnectoR::juliaEval(paste0(
    "[Base.find_package(p) !== nothing for p in [",
    paste0('"', pkgs, '"', collapse = ", "), "]]")), error = function(e) NULL)
  if (!is.logical(found) || length(found) != length(pkgs) || anyNA(found))
    found <- rep(FALSE, length(pkgs))
  skip_if_not(all(found), paste("Julia packages not loadable:",
                                paste(pkgs[!found], collapse = ", ")))
}

tiny_bdegpd_models <- local({
  paths <- list()
  function(family) {
    if (is.null(paths[[family]]))
      paths[[family]] <<- suppressMessages(train_bdegpd(
        savepath = file.path(tempdir(), paste0("egpd-test-", family)), family = family,
        K = 200L, m = 300L:600L, epochs = 2L, stopping_epochs = 1L, mc.cores = 1L,
        seed = 1L, verbose = FALSE))
    paths[[family]]
  }
})

test_that("NPE inference for bdegpd produces valid distfit object", {
  skip_if_no_julia()

  set.seed(42)
  Y <- rbdegpd(500, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  fit <- distfit(Y, family = "bdegpd", method = "neuralbayes", estimator = "npe",
                 model.path = tiny_bdegpd_models("bdegpd")$npe)

  expect_s3_class(fit, "distfit")
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
  skip_if_no_julia()

  set.seed(42)
  Y <- rbdegpd(500, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  fit <- distfit(Y, family = "bdegpd", method = "neuralbayes", estimator = "nbe",
                 model.path = tiny_bdegpd_models("bdegpd")$nbe)

  expect_s3_class(fit, "distfit")
  expect_equal(fit$estimator_type, "nbe")
  expect_null(fit$posterior_samples)
  expect_equal(length(fit$estimate), 6)
  expect_true(all(fit$estimate > 0))
})

test_that("NPE inference for bzidegpd produces valid distfit object with 7 params", {
  skip_if_no_julia()

  set.seed(42)
  Y <- rbzidegpd(500, kappa = 2, sigma = 1, xi = 0.1,
                  thL = 5, thU = 5, thw = 0.2, pi0 = 0.3)
  fit <- distfit(Y, family = "bzidegpd", method = "neuralbayes", estimator = "npe",
                 model.path = tiny_bdegpd_models("bzidegpd")$npe)

  expect_s3_class(fit, "distfit")
  expect_equal(fit$family, "bzidegpd")
  expect_equal(fit$npar, 7)
  expect_equal(length(fit$estimate), 7)
  expect_true("pi0" %in% names(fit$estimate))
  expect_true(fit$estimate["pi0"] > 0 && fit$estimate["pi0"] < 1)
  expect_true(!is.null(fit$posterior_samples))
  expect_equal(nrow(fit$posterior_samples), 7)
})

test_that("S3 methods work for bdegpd NPE fit", {
  skip_if_no_julia()

  set.seed(42)
  Y <- rbdegpd(500, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
  fit <- distfit(Y, family = "bdegpd", method = "neuralbayes", estimator = "npe",
                 model.path = tiny_bdegpd_models("bdegpd")$npe)

  ## print
  expect_output(print(fit), "Experimental")

  ## summary
  s <- summary(fit)
  expect_s3_class(s, "summary.distfit")
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
