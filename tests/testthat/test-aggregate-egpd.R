multiscale_fixture <- function(family = c("cpegpd", "cpdegpd"), n = 300L) {
  family <- match.arg(family)
  durations <- c(1L, 2L, 4L, 8L)
  sigma_coef <- c(log(if (family == "cpegpd") 0.6 else 1.2), if (family == "cpegpd") 0.12 else 0.08)
  lambda_coef <- c(log(if (family == "cpegpd") 0.8 else 0.9), if (family == "cpegpd") 0.25 else 0.22)
  kappa <- if (family == "cpegpd") 0.35 else 0.4
  xi <- if (family == "cpegpd") 0.2 else 0.18

  generator <- if (family == "cpegpd") rcpegpd else rcpdegpd
  samples <- lapply(durations, function(d) {
    sigma <- exp(sigma_coef[1] + sigma_coef[2] * log(d))
    lambda <- exp(lambda_coef[1] + lambda_coef[2] * log(d))
    generator(n, sigma = sigma, xi = xi, kappa = kappa, lambda = lambda, type = 1)
  })
  names(samples) <- as.character(durations)

  list(
    durations = durations,
    samples = samples,
    start = list(
      sigma_coef = sigma_coef,
      lambda_coef = lambda_coef,
      kappa = kappa,
      xi = xi
    ),
    truth = list(kappa = kappa, xi = xi)
  )
}

test_that("fit_aggregate_egpd fits the continuous multiscale framework", {
  set.seed(2)
  fixture <- multiscale_fixture("cpegpd", n = 300L)
  fit <- fit_aggregate_egpd(
    fixture$samples,
    durations = fixture$durations,
    family = "cpegpd",
    p = 1,
    q = 1,
    start = fixture$start,
    optim.method = "BFGS",
    control = list(maxit = 200)
  )

  expect_s3_class(fit, "aggregate_egpd")
  expect_equal(fit$convergence, 0)

  pars <- predict(fit, type = "parameter")
  expect_equal(pars$duration, fixture$durations)
  expect_true(all(diff(pars$sigma) >= -1e-8))
  expect_true(all(diff(pars$lambda) >= -1e-8))
  expect_equal(fit$xi, fixture$truth$xi, tolerance = 0.05)

  qs <- predict(fit, type = "quantile", prob = c(0.5, 0.9))
  expect_equal(dim(qs), c(length(fixture$durations), 2L))
  expect_true(all(qs[, 2] >= qs[, 1]))
})

test_that("fit_aggregate_egpd fits the discrete multiscale framework", {
  set.seed(3)
  fixture <- multiscale_fixture("cpdegpd", n = 250L)
  fit <- fit_aggregate_egpd(
    fixture$samples,
    durations = fixture$durations,
    family = "cpdegpd",
    p = 1,
    q = 1,
    start = fixture$start,
    optim.method = "BFGS",
    control = list(maxit = 200)
  )

  expect_s3_class(fit, "aggregate_egpd")
  expect_equal(fit$family, "cpdegpd")
  expect_equal(fit$convergence, 0)

  pars <- predict(fit, type = "parameter")
  expect_equal(pars$duration, fixture$durations)
  expect_true(all(diff(pars$sigma) >= -1e-8))
  expect_true(all(diff(pars$lambda) >= -1e-8))

  qs <- predict(fit, type = "quantile", prob = c(0.5, 0.9))
  expect_true(all(qs == floor(qs)))
})

test_that("fit_aggregate_egpd accepts finest-scale series input", {
  set.seed(4)
  x <- rcpegpd(200, sigma = 0.8, xi = 0.2, kappa = 0.35, lambda = 0.7, type = 1)
  fit <- NULL

  expect_no_error({
    fit <- fit_aggregate_egpd(
      x,
      durations = c(1L, 2L, 4L),
      family = "cpegpd",
      p = 1,
      q = 1,
      start = list(
        sigma_coef = c(log(0.8), 0.05),
        lambda_coef = c(log(0.7), 0.05),
        kappa = 0.35,
        xi = 0.2
      ),
      optim.method = "BFGS",
      control = list(maxit = 50)
    )
  })

  expect_equal(predict(fit, type = "parameter")$duration, c(1L, 2L, 4L))
})
