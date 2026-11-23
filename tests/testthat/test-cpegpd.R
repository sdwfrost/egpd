## Tests for Compound Poisson-EGPD (cpegpd)

test_that("rcpegpd produces non-negative data with zeros and positive values", {
  set.seed(42)
  x <- rcpegpd(500, sigma = 2, xi = 0.1, kappa = 1.5, lambda = 2, type = 1)

  expect_true(is.numeric(x))
  expect_equal(length(x), 500)
  expect_true(all(x >= 0))
  expect_true(any(x == 0))
  expect_true(any(x > 0))
})

test_that("dcpegpd sums to approximately 1 over full grid", {
  h <- 0.2
  K <- 500
  xgrid <- (0:K) * h
  pmf <- dcpegpd(xgrid, lambda = 2, sigma = 2, xi = 0.1, kappa = 1.5,
                  type = 1, h = h, K = K)

  expect_true(all(pmf >= 0))
  expect_true(abs(sum(pmf) - 1) < 0.05)
})

test_that("pcpegpd is monotonically increasing and reaches ~1", {
  qseq <- seq(0, 50, by = 0.5)
  cdf <- pcpegpd(qseq, lambda = 2, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)

  expect_true(all(diff(cdf) >= -1e-10))
  expect_true(cdf[1] > 0)
  expect_true(tail(cdf, 1) > 0.95)
})

test_that("qcpegpd inverts pcpegpd", {
  pvals <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  q <- qcpegpd(pvals, lambda = 2, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)

  ## generalised inverse on the lattice {0, h, 2h, ...}, h = 0.2: the smallest
  ## lattice point whose CDF reaches p
  expect_equal(q / 0.2, round(q / 0.2))
  cdf <- function(v) pcpegpd(v, lambda = 2, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
  expect_true(all(cdf(q) >= pvals))
  expect_true(all(cdf(q - 0.2) < pvals))
})

test_that("distfit with family='cpegpd' type 1 converges on simulated data", {
  set.seed(42)
  true_sigma <- 2; true_xi <- 0.1; true_kappa <- 1.5; true_lambda <- 2
  x <- rcpegpd(500, sigma = true_sigma, xi = true_xi, kappa = true_kappa,
                lambda = true_lambda, type = 1)
  fit <- distfit(x, type = 1, family = "cpegpd")

  expect_s3_class(fit, "distfit")
  expect_equal(fit$convergence, 0)
  expect_equal(fit$family, "cpegpd")
  expect_equal(fit$type, 1)
  expect_equal(fit$npar, 4)
  expect_true(all(c("sigma", "xi", "kappa", "lambda") %in% names(fit$estimate)))
  expect_true(!is.null(fit$cpegpd.h))
})

test_that("distfit cpegpd parameter recovery", {
  set.seed(123)
  true_sigma <- 2; true_xi <- 0.1; true_kappa <- 1.5; true_lambda <- 2
  x <- rcpegpd(800, sigma = true_sigma, xi = true_xi, kappa = true_kappa,
                lambda = true_lambda, type = 1)
  fit <- distfit(x, type = 1, family = "cpegpd")

  ## every parameter, xi included, within 4 standard errors of the truth
  truth <- c(sigma = true_sigma, xi = true_xi, kappa = true_kappa, lambda = true_lambda)
  expect_true(all(is.finite(fit$sd)))
  expect_true(all(abs(fit$estimate[names(truth)] - truth) < 4 * fit$sd[names(truth)]))
})

test_that("distfit cpegpd with fix.arg works", {
  set.seed(42)
  x <- rcpegpd(300, sigma = 2, xi = 0.1, kappa = 1.5, lambda = 2, type = 1)
  fit <- distfit(x, type = 1, family = "cpegpd", fix.arg = list(xi = 0.1))

  expect_equal(fit$convergence, 0)
  expect_false("xi" %in% names(fit$estimate))
  expect_equal(fit$npar, 3)
  expect_equal(fit$fix.arg, list(xi = 0.1))
})

test_that("distfit cpegpd S3 methods work", {
  set.seed(42)
  x <- rcpegpd(300, sigma = 2, xi = 0.1, kappa = 1.5, lambda = 2, type = 1)
  fit <- distfit(x, type = 1, family = "cpegpd")

  ## summary
  s <- summary(fit)
  expect_s3_class(s, "summary.distfit")
  expect_true(!is.null(s$cpegpd.h))

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

test_that("distfit cpegpd plot does not error", {
  set.seed(42)
  x <- rcpegpd(300, sigma = 2, xi = 0.1, kappa = 1.5, lambda = 2, type = 1)
  fit <- distfit(x, type = 1, family = "cpegpd")

  expect_no_error(plot(fit))
})
test_that("compound distribution queries preserve NA and NaN separately", {
  for (f in list(dcpdegpd, pcpdegpd, qcpdegpd, dcpegpd, pcpegpd, qcpegpd)) {
    run <- function(z, rate) f(z, lambda=rate, sigma=1.5, xi=0, kappa=1, K=50)
    expect_identical(run(c(0,NA_real_,NaN),.8)[2:3], c(NA_real_,NaN))
    expect_identical(run(c(0,NaN),NA_real_), c(NA_real_,NaN))
    expect_true(is.nan(run(NA_real_,NaN)))
  }
})

test_that("compound distribution diagnostics distinguish parameter length from domain", {
  for (f in list(dcpegpd,dcpdegpd)) {
    expect_error(f(0,lambda=c(5,5),sigma=1,xi=0,kappa=1), "scalar: lambda")
    expect_error(f(0,lambda=5,sigma=c(1,1),xi=0,kappa=1), "scalar: sigma")
    expect_error(f(0,lambda=-5,sigma=1,xi=0,kappa=1), "positive and finite")
    expect_length(f(c(0,1),lambda=5,sigma=1,xi=0,kappa=1),2L)
  }
})
