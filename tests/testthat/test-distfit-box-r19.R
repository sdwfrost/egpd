test_that("Bernstein refinement stays within the requested working-coordinate box", {
  y <- 5*qexp(((1:120)-.5)/120)
  fit <- function(method) distfit(y, method = "bernstein", bernstein.m = 1L,
    fix.arg = list(xi = 0), start = list(kappa = 1, sigma = 1), optim.method = method,
    lower = rep(log(.2),2), upper = rep(log(1.5),2))
  ll <- function(k,s) sum(log(k)-log(s)-y/s+(k-1)*log(-expm1(-y/s)))
  for (m in c("Nelder-Mead","L-BFGS-B")) {
    expect_warning(f <- fit(m), "optimizer bound")
    expect_equal(unname(coef(f)), c(1.5,1.5), tolerance = 1e-6)
    expect_equal(as.numeric(logLik(f)), ll(1.5,1.5), tolerance = 1e-9)
    expect_true(all(is.na(c(f$sd, confint(f)))))
    expect_error(vcov(f), "optimizer bound")
  }
})

test_that("ordinary Poisson-limit fits flag active boxes but retain interior SEs", {
  y <- rep(0:4,c(6,12,24,12,6)); n <- length(y)
  for (m in c("Brent","L-BFGS-B")) {
    fit <- function(lo,hi) distfit(y, family = "gpig", fix.arg = list(a = 1, c = .5),
        start = list(b = 2), optim.method = m, lower = log(2*lo), upper = log(2*hi))
    for (box in list(c(.4,.9), c(2.6,4))) {
      expect_warning(f <- fit(box[1],box[2]), "optimizer bound")
      expected <- max(box[1],min(mean(y),box[2]))
      expect_equal(unname(coef(f))/2, expected, tolerance = 1e-7)
      expect_equal(as.numeric(logLik(f)), sum(dpois(y,expected,log = TRUE)), tolerance = 1e-9)
      expect_true(all(is.na(c(f$sd, confint(f)))))
      expect_error(vcov(f), "optimizer bound")
    }
    f <- fit(.4,3.8)
    expect_equal(unname(f$sd)/2, sqrt(mean(y)/n), tolerance = 1e-5)
  }
})
