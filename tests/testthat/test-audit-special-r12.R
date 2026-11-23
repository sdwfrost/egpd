test_that("a scalar extracted from coef() keeps the outer fixed-parameter name", {
  x <- rep(0:4, c(69, 26, 12, 8, 7))
  fit <- distfit(x, family = "gpig", fix.arg = list(a = c(a = 1), c = .4),
                 start = list(b = c(pilot = 2)), hessian = FALSE)
  expect_equal(unname(coef(fit)), mean(x) / .4, tolerance = 1e-5)
  expect_equal(as.numeric(logLik(fit)), sum(dpois(x, mean(x), log = TRUE)))

  p <- .32
  f <- distfit(x, family = "zibell", fix.arg = list(pi = c(old = p)))
  exact <- uniroot(function(t) {
    p0 <- exp(-expm1(t))
    sum(x)/t - exp(t) * (sum(x > 0) + sum(x == 0) * (1-p)*p0/(p+(1-p)*p0))
  }, c(.01, 2), tol = 1e-12)$root
  expect_equal(unname(coef(f)), exact, tolerance = 1e-5)
})

test_that("generated hurdle and splice responses do not replace user predictors", {
  ya <- rep(0:4, c(45,28,12,5,4)); yb <- rep(0:4, c(9,13,16,30,25))
  hd <- data.frame(y = c(ya,yb), .positive = rep(0:1, c(length(ya),length(yb))))
  hd$.positive. <- 7
  h <- fit_hurdle(list(y ~ 1, ~1, ~1), hd, zero.formula = ~log1p(.positive),
                  family = "degpd", degpd.args = list(m = 1, link = "identity"),
                  fixed = list(kappa = 1, xi = 0), trace = -1)
  expect_equal(predict(h, data.frame(.positive = 0:1))$pzero,
               c(mean(ya == 0), mean(yb == 0)), tolerance = 1e-7)
  skip_if_not_installed("qgam")
  set.seed(475)
  sd <- data.frame(exceedance = rep(c(0,4),220), v = runif(440))
  sd$exceedance. <- 10
  sd$y <- rgeom(nrow(sd), 1/(3+.7*sd$exceedance+sd$v))
  f <- fit_spliced_degpd(y ~ v, sd, bulk_probs = c(.25,.55), threshold_prob = .55,
                        tail_formula = ~sqrt(1+exceedance),
                        egpd_control = list(fix.arg = list(xi = 0, kappa = 1),
                                            control = list(reltol = 1e-12)))
  z <- sd$y - predict(f, sd, type = "parameter")$threshold
  ref <- vapply(c(0,4), function(g) 1/log1p(1/mean(z[z >= 0 & sd$exceedance == g])), 0)
  got <- predict(f, data.frame(exceedance=c(0,4),v=.3), type="parameter")$scale
  expect_equal(got, ref, tolerance=1e-4)
})

test_that("MDGPD ceilings respect the open lower endpoint of each component", {
  set.seed(61)
  n <- 3000L
  y <- rmdgpd(n, sigma=c(.5,1,1.5), xi=.5, d=3, rate=620)
  set.seed(61)
  T <- matrix(rpois(3*n,620), n)
  W <- T-apply(T,1,max)+rexp(n)
  ## Translate by the integer sigma/xi only *after* taking the ceiling.
  k <- 1:3
  ref <- sweep(ceiling(sweep(exp(.5*W), 2, k, `*`)), 2, k, `-`)
  expect_equal(unname(y), ref)
  expect_true(any(y < 0))
  ## 1/(the double nearest 1/3) rounds to 3, but the exact ratio is
  ## slightly above 3: -3 can occur for a sufficiently negative excess.
  set.seed(7)
  expect_equal(min(rmdgpd(400, 1, 1/3, rate=5000)), -3)
})
