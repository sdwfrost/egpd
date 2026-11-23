test_that("a fitted two-cell jump law identifies only its thinned mean", {
  y <- rep(0:3, c(20,50,24,6))
  mu <- mean(y)
  ## Floor of Uniform(0,s), 1 < s < 2: effective rate is lambda*(1-1/s).
  for (s in c(1.12,1.75)) {
    expect_warning(f <- distfit(y, family = "cpdegpd", fix.arg = list(xi = -1, kappa = 1),
                     start = list(sigma = s, lambda = mu/(1-1/s))), "Fitted jump support")
    expect_equal(as.numeric(logLik(f)), sum(dpois(y,mu,log = TRUE)), tolerance = 1e-8)
    expect_true(isTRUE(f$support.ridge))
    expect_true(all(is.na(c(f$sd,confint(f)))))
    expect_error(vcov(f), "not separately identifiable")
  }
  g <- distfit(y, family = "cpdegpd", fix.arg = list(xi = -1,kappa = 1,sigma = 1.5))
  expect_equal(unname(coef(g)), 3*mu, tolerance = 1e-5)
  expect_equal(unname(g$sd), 3*sqrt(mu/length(y)), tolerance = 1e-5)
})
