test_that("a power-carrier sup at the excluded uniform endpoint is not a regular MLE", {
  x <- (1:60 - .5)/12
  for (k in c(1, 2.5)) {
    fixed <- list(kappa = k, xi = -1)
    for (m in c("mle", "bernstein"))
      expect_error(distfit(x, method = m, bernstein.m = 1, fix.arg = fixed),
                   "Nonregular continuous endpoint")
    expect_error(distfit(c(0,x), family = "ziegpd", fix.arg = c(fixed, list(pi = .2))),
                 "Nonregular continuous endpoint")
    box <- distfit(x, fix.arg = fixed, optim.method = "Brent", lower = log(6), upper = log(8),
                   hessian = FALSE)
    expect_equal(unname(coef(box)), 6, tolerance = 1e-4)
  }
  ## With sigma fixed, kappa really is estimable on this same law.
  power <- distfit(x, fix.arg = list(sigma = 6, xi = -1), optim.method = "Brent",
                   lower = log(.2), upper = log(3), control = list(reltol = 1e-12))
  expect_equal(unname(coef(power)), -1/mean(log(x/6)), tolerance = 1e-5)
})

test_that("a fitted two-cell law does not identify a scale and a power", {
  x <- rep(0:1, c(45,15))
  for (s in c(1.1, 1.8)) {
    expect_warning(f <- distfit(x, family = "degpd", fix.arg = list(xi = -1),
                                start = list(sigma = s, kappa = -log(.75)/log(s))),
                   "not separately identifiable")
    expect_equal(as.numeric(logLik(f)), sum(dbinom(x,1,.25,log = TRUE)), tolerance = 1e-8)
    expect_true(all(is.na(c(f$sd, confint(f)))))
    expect_error(vcov(f), "not separately identifiable")
  }
  one <- distfit(x, family = "degpd", fix.arg = list(xi = -1, kappa = 1),
                 start = list(sigma = 4/3))
  expect_gt(unname(one$sd), 0)
})

test_that("bounds apply on the actual reduced parameter", {
  fixed <- list(sigma = 1.5, xi = -1)
  for (method in c("Brent","L-BFGS-B")) {
    x <- rep(0:1,c(20,60))
    expect_warning(f <- distfit(x, family = "zidegpd", fix.arg = fixed, optim.method = method,
                                lower = qlogis(.2), upper = qlogis(.4)), "optimizer bound")
    expect_equal(coef(f), c(p1 = .4))
    expect_equal(as.numeric(logLik(f)), sum(dbinom(x,1,.4,log = TRUE)))
    expect_true(all(is.na(confint(f))))
    x <- rep(0:4,15)
    p <- distfit(x, family = "cpdegpd", fix.arg = fixed, hessian = FALSE,
                 optim.method = method, lower = log(.5), upper = log(.9))
    expect_equal(coef(p), c(mu = .9))
    expect_equal(as.numeric(logLik(p)), sum(dpois(x,.9,log = TRUE)))
  }
  z <- distfit(rep(0,20), family = "cpdegpd", fix.arg = fixed, hessian = FALSE,
               lower = log(.3))
  expect_equal(coef(z), c(mu = .3))
  expect_error(distfit(0:1, family = "zidegpd", fix.arg = fixed, lower = 1, upper = 0),
               "increasing scalar bounds")
})
