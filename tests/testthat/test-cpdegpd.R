## Tests for Compound Poisson-Discrete EGPD (cpdegpd)

test_that("rcpdegpd produces non-negative integer data with zeros", {
  set.seed(42)
  x <- rcpdegpd(500, sigma = 3, xi = 0.1, kappa = 1.5, lambda = 2, type = 1)

  expect_true(is.numeric(x))
  expect_equal(length(x), 500)
  expect_true(all(x >= 0))
  expect_true(all(x == floor(x)))
  expect_true(any(x == 0))
  expect_true(any(x > 0))
})

test_that("dcpdegpd sums to approximately 1 over integer grid", {
  K <- 300
  xgrid <- 0:K
  pmf <- dcpdegpd(xgrid, lambda = 2, sigma = 3, xi = 0.1, kappa = 1.5,
                   type = 1, K = K)

  expect_true(all(pmf >= 0))
  expect_true(abs(sum(pmf) - 1) < 0.01)
})

test_that("pcpdegpd is monotonically increasing and reaches ~1", {
  qseq <- 0:100
  cdf <- pcpdegpd(qseq, lambda = 2, sigma = 3, xi = 0.1, kappa = 1.5, type = 1)

  expect_true(all(diff(cdf) >= -1e-10))
  expect_true(cdf[1] > 0)
  expect_true(tail(cdf, 1) > 0.95)
})

test_that("pcpdegpd equals cumsum of dcpdegpd", {
  K <- 200
  xgrid <- 0:K
  pmf <- dcpdegpd(xgrid, lambda = 2, sigma = 3, xi = 0.1, kappa = 1.5,
                   type = 1, K = K)
  cdf_from_pmf <- cumsum(pmf)
  cdf_direct <- pcpdegpd(xgrid, lambda = 2, sigma = 3, xi = 0.1, kappa = 1.5,
                          type = 1, K = K)

  expect_equal(cdf_from_pmf, cdf_direct, tolerance = 1e-10)
})

test_that("qcpdegpd inverts pcpdegpd", {
  pvals <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  q <- qcpdegpd(pvals, lambda = 2, sigma = 3, xi = 0.1, kappa = 1.5, type = 1)

  expect_true(all(q >= 0))
  expect_true(all(q == floor(q)))
  expect_true(all(diff(q) >= 0))

  ## Check that P(S <= q) >= p for all p
  cdf_at_q <- pcpdegpd(q, lambda = 2, sigma = 3, xi = 0.1, kappa = 1.5, type = 1)
  expect_true(all(cdf_at_q >= pvals - 1e-10))
})

test_that("fitegpd with family='cpdegpd' type 1 converges on simulated data", {
  set.seed(42)
  true_sigma <- 3; true_xi <- 0.1; true_kappa <- 1.5; true_lambda <- 2
  x <- rcpdegpd(500, sigma = true_sigma, xi = true_xi, kappa = true_kappa,
                 lambda = true_lambda, type = 1)
  fit <- fitegpd(x, type = 1, family = "cpdegpd")

  expect_s3_class(fit, "fitegpd")
  expect_equal(fit$convergence, 0)
  expect_equal(fit$family, "cpdegpd")
  expect_equal(fit$type, 1)
  expect_equal(fit$npar, 4)
  expect_true(all(c("sigma", "xi", "kappa", "lambda") %in% names(fit$estimate)))
})

test_that("fitegpd cpdegpd parameter recovery", {
  set.seed(123)
  true_sigma <- 3; true_xi <- 0.1; true_kappa <- 1.5; true_lambda <- 2
  x <- rcpdegpd(800, sigma = true_sigma, xi = true_xi, kappa = true_kappa,
                 lambda = true_lambda, type = 1)
  fit <- fitegpd(x, type = 1, family = "cpdegpd")

  expect_true(fit$estimate["sigma"] > 0.5 && fit$estimate["sigma"] < 8)
  expect_true(fit$estimate["kappa"] > 0.3 && fit$estimate["kappa"] < 5)
  expect_true(fit$estimate["lambda"] > 0.5 && fit$estimate["lambda"] < 5)
})

test_that("fitegpd cpdegpd with fix.arg works", {
  set.seed(42)
  x <- rcpdegpd(300, sigma = 3, xi = 0.1, kappa = 1.5, lambda = 2, type = 1)
  fit <- fitegpd(x, type = 1, family = "cpdegpd", fix.arg = list(xi = 0.1))

  expect_equal(fit$convergence, 0)
  expect_false("xi" %in% names(fit$estimate))
  expect_equal(fit$npar, 3)
  expect_equal(fit$fix.arg, list(xi = 0.1))
})

test_that("fitegpd cpdegpd S3 methods work", {
  set.seed(42)
  x <- rcpdegpd(300, sigma = 3, xi = 0.1, kappa = 1.5, lambda = 2, type = 1)
  fit <- fitegpd(x, type = 1, family = "cpdegpd")

  ## summary
  s <- summary(fit)
  expect_s3_class(s, "summary.fitegpd")

  ## print (should not error)
  expect_output(print(fit))
  expect_output(print(s))

  ## AIC
  expect_true(is.numeric(AIC(fit)))

  ## confint
  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 4)
})

test_that("fitegpd cpdegpd plot does not error", {
  set.seed(42)
  x <- rcpdegpd(300, sigma = 3, xi = 0.1, kappa = 1.5, lambda = 2, type = 1)
  fit <- fitegpd(x, type = 1, family = "cpdegpd")

  expect_no_error(plot(fit))
})
