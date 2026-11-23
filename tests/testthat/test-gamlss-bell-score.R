test_that("GAMLSS Bell scores agree with the log-PMF at small and ordinary means", {
  f <- BELL()
  mu <- c(1.1e-8, 1e-16, 1e-10, .03, 1, 4)
  y <- c(1, 1, 0, 0, 1, 3)
  h <- 1e-5
  fd <- (dBELL(y, mu * exp(h), log = TRUE) - dBELL(y, mu * exp(-h), log = TRUE)) / (2*h)
  sm <- f$dldm(y, mu)
  expect_equal(mu * sm, fd, tolerance = 1e-8)
  expect_equal(f$d2ldm2(y, mu), pmin(-sm^2, -1e-15))
  expect_identical(f$dldm(NA_real_, 1), NA_real_)
  expect_warning(bad <- f$dldm(c(-1, .5, 0), c(1, 1, 0)), "undefined")
  expect_true(all(is.nan(bad)))
})

test_that("ZIBELL weights Bell zeros and keeps the one-sided score at sigma = 0", {
  f <- ZIBELL()
  theta <- c(1.1e-8, 1e-16, .7, .7, 8)
  mu <- theta * exp(theta)
  y <- c(1, 0, 0, 3, 0)
  pi <- c(.25, .75, 0, 0, 0)
  sm <- f$dldm(y, mu, pi)
  sp <- f$dldd(y, mu, pi)
  ## At pi = 0 even a vanishingly rare zero is known to be a Bell zero.
  expect_equal(sm * mu, (y - mu) / (1 + theta) * c(1, .25, 1, 1, 1),
               tolerance = 1e-13)
  expect_equal(sp[1:4], c(-1 / .75, 1e-16, expm1(expm1(.7)), -1), tolerance = 1e-12)
  expect_equal(sp[2] / mu[2], 1, tolerance = 1e-13)
  expect_equal(f$d2ldm2(y, mu, pi), pmin(-sm^2, -1e-15))
  expect_equal(f$d2ldd2(y, mu, pi), pmin(-sp^2, -1e-15))
  expect_equal(f$d2ldmdd(y, mu, pi), -sm * sp)
  ## Non-negligible probabilities for both components; independent log-PMF differences.
  h <- 1e-5
  ll <- function(m, p) dZIBELL(0, mu = m, sigma = p, log = TRUE)
  expect_equal(f$dldm(0, 3, .2), (ll(3+h, .2) - ll(3-h, .2)) / (2*h), tolerance = 1e-8)
  expect_equal(f$dldd(0, 3, .2), (ll(3, .2+h) - ll(3, .2-h)) / (2*h), tolerance = 1e-8)
})

test_that("a Bell offset fit with a low-exposure positive count maximizes the likelihood", {
  skip_if_not_installed("gamlss")
  d <- data.frame(y = c(1, 0, 1, 1, 2, 1, 4, 0, 2, 5, 0, 2, 1),
                  z = c(log(1.1e-8) - .55, rep(0, 12)))
  ll <- function(b) sum(dBELL(d$y, mu = exp(b + d$z), log = TRUE))
  reference <- optimize(ll, c(0, 1), maximum = TRUE, tol = 1e-10)
  f <- gamlss::gamlss(y ~ offset(z), data = d, family = egpd::BELL(),
                      control = gamlss::gamlss.control(n.cyc = 200, mu.step = .3,
                                                       c.crit = 1e-9, trace = FALSE))
  expect_true(f$converged)
  expect_lt(abs(unname(coef(f)) - reference$maximum), 1e-5)
})
