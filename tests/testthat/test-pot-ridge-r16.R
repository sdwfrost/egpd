test_that("one- and two-cell fitted POT supports withhold separate Wald inference", {
  u <- 8
  z <- rep(0:1, c(45, 15))
  x <- c(1:5, u + z)
  for (xi in c(NA_real_, -.75, -2.4)) {
    start <- if (!is.na(xi)) list(xi = xi, sigma = -xi / (-expm1(-xi * log(.25))))
    expect_warning(f <- fit_degpd_pot(x, u, start = start), "not separately identifiable")
    expect_equal(as.numeric(logLik(f)), sum(dbinom(z, 1, .25, log = TRUE)), tolerance = 1e-8)
    expect_true(all(is.na(c(f$sd, f$xi.ci, confint(f)))))
    expect_error(vcov(f), "not separately identifiable")
    expect_equal(attr(logLik(f), "df"), 2)
  }
  expect_warning(st <- degpd_pot_stability(x, u, hessian = FALSE), "not separately identifiable")
  expect_true(all(is.na(st[c("sigma", "xi", "xi_lo", "xi_hi")])))
  expect_warning(zero <- fit_degpd_pot(x, u + 1, start = list(sigma = .45, xi = -1)),
                 "not separately identifiable")
  expect_true(all(is.na(zero$sd)))
  ## No ridge with a fixed shape and one estimable parameter, or with
  ## three fitted cells of masses (.4, .4, .2) and both parameters free.
  pinned <- distfit(z, family = "degpd", fix.arg = list(kappa = 1, xi = -1),
                    start = list(sigma = 4/3))
  expect_equal(unname(coef(pinned)), 4/3, tolerance = 1e-5)
  expect_gt(unname(pinned$sd), 0)
  three <- fit_degpd_pot(u + rep(0:2, c(24, 24, 12)), u, start = list(sigma = 2.35, xi = -.96))
  expect_equal(unname(coef(three)), c(2.5, -1), tolerance = .001)
  expect_true(all(is.finite(three$xi.ci)))
})
