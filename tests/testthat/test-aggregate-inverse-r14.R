## Freeze the duration-one laws rather than making an inversion test depend on
## the optimizer. predict() still takes its usual fitted-coefficient path.
r14_aggregate <- function(lambda, sigma, xi, family = "cpegpd") {
  structure(list(durations = 1L, family = family, cpegpd.h = 1,
                 lambda_coef = log(lambda), sigma_coef = log(sigma),
                 kappa = 1, xi = xi), class = "aggregate_egpd")
}

r14_q <- function(fit, p, K)
  as.numeric(predict(fit, type = "quantile", prob = p, K = K))

test_that("a rare large gauge jump is retained beyond the quantile grid", {
  lambda <- 1001
  sigma <- 1e-14
  fit <- r14_aggregate(lambda, sigma, 1)
  sj <- function(j) sigma / (sigma + j)
  den <- -expm1(-lambda)
  ## Condition on the total Poisson count. For N <= k count configurations
  ## with exactly one non-unit jump; for N > k every configuration crosses.
  lower <- function(k) {
    n <- seq_len(k)
    one <- sum(dpois(n, lambda) * n * exp((n - 1) * log1p(-sj(1))) * sj(k - n + 1))
    (ppois(k, lambda, lower.tail = FALSE) + one) / den
  }
  omitted <- ppois(1, lambda * sj(1), lower.tail = FALSE) / den
  p <- c(.95, 1 - 2^-45, 1 - 2^-47)
  k <- c(1053, 1357, 2411)
  for (i in seq_along(p)) {
    expect_gt(lower(k[i] - 1), 1 - p[i])
    expect_lte(lower(k[i]) + omitted, 1 - p[i])
  }
  for (K in c(2411, 3500, 4500))
    expect_equal(r14_q(fit, p, K), k)
  expect_error(r14_q(fit, p[3], 2410), "exceeds the compound search grid")
  expect_error(r14_q(fit, .95, 900), "exceeds the compound search grid")

  ## The unconditional law has the same brackets (exp(-1001) zero atom).
  ## Changing the gauge and scale together only changes the units.
  expect_equal(qcpegpd(p, lambda, kappa = 1, sigma = sigma / 4, xi = 1,
                       h = .25, K = 3500), k / 4)
})

test_that("the extreme light-tailed gauge quantile is 191, not 189", {
  lambda <- 100
  fit <- r14_aggregate(lambda, .005, .02)
  den <- -expm1(-lambda)
  ## S = N unless at least one jump exceeds one bin.
  err <- -expm1(-lambda * (1 + .02 / .005)^(-1 / .02)) / den
  expect_gt(ppois(190, lambda, lower.tail = FALSE) / den, 2^-52)
  expect_lte(ppois(191, lambda, lower.tail = FALSE) / den + err, 2^-52)
  for (K in c(191, 400, 800))
    expect_equal(r14_q(fit, c(.95, 1 - 2^-52), K), c(117, 191))
  expect_error(r14_q(fit, 1 - 2^-52, 190), "exceeds the compound search grid")
})

test_that("native discrete jumps use their own zero atom and positive law", {
  ## floor(Uniform(0, 1.5)) is Bernoulli(1/3): the compound law is exactly
  ## Poisson(.1), *not* Poisson(.3), and prediction conditions it on S > 0.
  fit <- r14_aggregate(.3, 1.5, -1, "cpdegpd")
  p <- c(.5, .95, 1 - 2^-52)
  tail <- ppois(0:15, .1, lower.tail = FALSE)
  inv <- function(den) vapply(p, function(x) which(tail / den <= 1 - x)[1L] - 1, 0)
  expect_equal(r14_q(fit, p, 10), inv(-expm1(-.1)))
  expect_error(r14_q(fit, p, 9), "exceeds the compound search grid")
  expect_equal(qcpdegpd(p, .3, kappa = 1, sigma = 1.5, xi = -1, K = 9), inv(1))

  ## Now the positive rate is below 1e-399. Given a positive total there is
  ## one jump except with probability < 1e-399; its survivor is
  ## ((b+1)/(b+k+1))^16. Both the normalizer and the off-grid mass underflow.
  b <- 1e-25
  fit <- r14_aggregate(1, b / 16, 1/16, "cpdegpd")
  expect_gt(((b + 1) / (b + 8))^16, 2^-50)
  expect_lt(((b + 1) / (b + 9))^16, 2^-50)
  expect_equal(r14_q(fit, 1 - 2^-50, 8), 8)
  expect_error(r14_q(fit, 1 - 2^-50, 7), "exceeds the compound search grid")
})
