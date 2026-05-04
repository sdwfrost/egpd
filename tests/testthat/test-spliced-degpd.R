test_that("fit_spliced_degpd fits and predicts a spliced count model", {
  skip_if_not_installed("qgam")

  set.seed(101)
  n <- 300
  x <- runif(n, -1, 1)
  y <- rpois(n, lambda = exp(1 + 0.8 * x)) +
    rbinom(n, 1, 0.12) * rdiscegpd(n, sigma = exp(0.8 + 0.5 * x), xi = 0.15, kappa = 1.3)
  df <- data.frame(y = y, x = x)

  fit <- suppressWarnings(
    fit_spliced_degpd(
      y ~ s(x, k = 6),
      data = df,
      bulk_probs = c(0.1, 0.5, 0.75, 0.9),
      threshold_prob = 0.75,
      tail_formula = list(lsigma = ~x, lxi = ~1, lkappa = ~1),
      keep_data = TRUE
    )
  )

  expect_s3_class(fit, "spliced_degpd")
  expect_s3_class(summary(fit), "summary.spliced_degpd")
  expect_equal(fit$threshold_prob, 0.75)
  expect_true(fit$n_exceed > 0)

  qhat <- predict(fit, type = "quantile", prob = c(0.25, 0.75, 0.95))
  expect_equal(dim(qhat), c(nrow(df), 3L))
  expect_true(all(qhat[[2]] >= qhat[[1]]))
  expect_true(all(qhat[[3]] >= qhat[[2]]))

  pars <- predict(fit, type = "parameter")
  expect_true(all(c("threshold", "threshold_quantile", "scale", "shape", "kappa") %in% names(pars)))
  expect_true(all(pars$threshold >= 0))

  cdf <- predict(fit, newdata = df[1:3, ], type = "cdf", at = 0:40)
  expect_equal(dim(cdf), c(3L, 41L))
  expect_true(all(as.matrix(cdf) >= 0))
  expect_true(all(as.matrix(cdf) <= 1))
  expect_true(all(apply(as.matrix(cdf), 1, function(z) all(diff(z) >= -1e-10))))

  pmf <- predict(fit, newdata = df[1:3, ], type = "pmf", at = 0:60)
  expect_equal(dim(pmf), c(3L, 61L))
  expect_true(all(as.matrix(pmf) >= -1e-12))
  expect_true(all(rowSums(as.matrix(pmf)) > 0.85))
  expect_true(all(rowSums(as.matrix(pmf)) <= 1 + 1e-8))
})

test_that("fit_spliced_degpd works without a tail model", {
  skip_if_not_installed("qgam")

  set.seed(102)
  x <- runif(150, -1, 1)
  y <- rpois(150, lambda = exp(0.6 + 0.4 * x))
  df <- data.frame(y = y, x = x)

  fit <- fit_spliced_degpd(
    y ~ x,
    data = df,
    bulk_probs = c(0.2, 0.5, 0.8),
    threshold_prob = 0.8,
    fit_tail = FALSE,
    keep_data = TRUE
  )

  expect_null(fit$tail_fit)

  bulk <- predict(fit, type = "bulk")
  expect_equal(dim(bulk), c(nrow(df), 3L))

  cdf <- predict(fit, newdata = df[1:2, ], type = "cdf", at = 0:20)
  expect_equal(dim(cdf), c(2L, 21L))
})
