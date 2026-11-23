test_that("a bounded power-count fit reaches the highest cell edge from either side", {
  n <- c(4, 12, 20, 28, 36, 80)
  x <- rep(17 + 0:5, n)
  ll6 <- function(k) sum(n * log(diff(((0:6)/6)^k)))
  ref <- optimize(ll6, c(2, 3), maximum = TRUE)
  for (start in list(list(sigma = 5.55, kappa = 1.5), list(sigma = 6.6, kappa = 3))) {
    expect_warning(f <- fit_degpd_pot(x, 17, fix.kappa = FALSE, fix.arg = list(xi = -1),
                         start = start, lower = c(0, 0), upper = rep(log(12), 2)),
                   "change in count support")
    expect_equal(f$convergence, 0)
    expect_identical(coef(f)[["sigma"]], 6)
    expect_equal(f$upper.endpoint, 22)
    expect_lt(abs(f$loglik - ref$objective), 1e-6)
    expect_true(all(is.na(f$sd)))
  }
  ## A scale bound that excludes 6 must still be honored.
  expect_warning(b <- distfit(x - 17, family = "degpd", fix.arg = list(xi = -1, kappa = 2),
                              lower = log(4), upper = log(5.97)), "optimizer bound")
  expect_equal(coef(b), c(sigma = 5.97))
})

test_that("a rounded conditional scale does not lose an observed count", {
  x <- c(0, 0, 0, 1)
  k <- 2^52
  s <- 1 + .Machine$double.eps
  lp0 <- -k*log1p(.Machine$double.eps)
  for (family in c("degpd", "zidegpd")) {
    pins <- list(xi = -1, kappa = k)
    if (family == "zidegpd") pins$pi <- 0
    f <- distfit(x, family = family, fix.arg = pins, hessian = FALSE)
    expect_identical(coef(f)[["sigma"]], s)
    expect_equal(f$loglik, 3*lp0 + log(-expm1(lp0)), tolerance = 1e-12)
  }
})
