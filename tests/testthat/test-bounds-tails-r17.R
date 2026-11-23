test_that("aggregate fits evaluate the supplied legal shape start", {
  d <- list(`1` = rep(0:4,6), `3` = rep(c(0,1,1,2,3,7),6))
  fit <- function(k, bounds) fit_aggregate_egpd(d, c(1L,3L), family = "cpdegpd", p = 0, q = 0,
    start = list(sigma_coef = log(1.7), lambda_coef = log(.85), kappa = k, xi = .2),
    kappa_bounds = bounds, control = list(maxit = 0L))
  k <- .9000000005
  narrow <- fit(k,c(.9,.900000001))
  expect_equal(narrow$kappa, k, tolerance = 1e-15)
  expect_equal(narrow$composite_loglik, fit(k,c(.1,2))$composite_loglik, tolerance = 1e-12)
  expect_equal(fit(2e-10,c(1e-12,5))$kappa/2e-10, 1, tolerance = 1e-12)
})

test_that("extreme lower probabilities invert the represented upper tail", {
  ## Independent gamma arrival-time integration under IG mixing brackets
  ## 1-1e-14 at 11190; the tolerance on p used to give 10773.
  p <- 1-1e-14
  expect_equal(qpig(c(.9,p),5,50), c(5,11190))
  expect_equal(qpig(1-p,5,50,lower.tail = FALSE),11190)
  expect_warning(q <- qpig(p,5,50,max.value = 11189), "not attained")
  expect_true(is.na(q))
  t <- 60
  S <- function(k) (t/(t+1))^k * t*(k+t^2)/((t-1)*(t+1)^2)
  s <- 2^-48
  k <- qprl(1-s,t)
  expect_true(S(k) <= s && s < S(k-1))
})
