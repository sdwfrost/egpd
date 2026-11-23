test_that("a tilted uniform carrier can attain zero with a regular nuisance parameter", {
  pos <- seq(.12, .28, length.out = 50)
  x <- c(rep(0, 10), pos)
  pins <- list(sigma = 1, xi = -1)
  ## On (0,1), log f = -k*(1-x)^2/2 - log Z(k). The score decreases
  ## strictly and is already negative at k=0 for these positive values.
  expect_warning(f <- distfit(x, type = 2, family = "ziegpd", fix.arg = pins),
                 "Carrier MLE.*natural boundary")
  expect_equal(coef(f), c(kappa = 0, pi = 1/6))
  expect_equal(as.numeric(logLik(f)), 10*log(1/6) + 50*log(5/6))
  expect_true(all(is.na(c(f$sd, confint(f)))))
  expect_error(vcov(f), "natural boundary")
  expect_warning(b <- distfit(pos, type = 2, fix.arg = pins, optim.method = "Brent",
                              lower = log(.2), upper = log(3)), "optimizer bound")
  expect_equal(coef(b), c(kappa = .2))
})

test_that("a three-cell compound mixture refits the intensity on its natural face", {
  set.seed(3)
  y <- qpois(runif(500), .9) + 2*qpois(runif(500), .5)
  pins <- list(sigma = 1.2, xi = -.5, kappa = .6, delta = 3)
  ## At weight zero the jump CDF uses only the second power. The rates
  ## of single and double jumps are lambda times their two cell masses.
  F <- (1 - (1 - (1:2)/2.4)^2)^3
  tab <- tabulate(y + 1)
  nll <- function(rate) {
    mass <- vapply(0:(length(tab)-1L), function(k) {
      j <- 0:floor(k/2)
      sum(dpois(j, rate*(1-F[2])) * dpois(k-2*j, rate*diff(F)))
    }, numeric(1))
    -sum(tab*log(mass))
  }
  ref <- optimize(nll, c(.5, 6), tol = 1e-10)
  expect_warning(f <- distfit(y, type = 6, family = "cpdegpd", fix.arg = pins,
                              start = list(prob = .5, lambda = 2)), "natural boundary")
  expect_equal(unname(coef(f)["prob"]), 0)
  expect_equal(unname(coef(f)["lambda"]), ref$minimum, tolerance = 2e-4)
  expect_equal(f$loglik, -ref$objective, tolerance = 1e-8)
  expect_true(all(is.na(f$sd)))
})

test_that("two-cell mixture weights distinguish included faces from nearby interior fits", {
  base <- list(sigma = 2, xi = -1, kappa = 1, delta = 2)
  y <- rep(1:4, 10)
  fit <- function(rate, ...) distfit(y, type = 6, family = "cpdegpd",
                                     fix.arg = c(base, list(lambda = rate)), ...)
  ## The aggregate is Poisson(rate*(3-p)/4). Only p is unknown.
  for (p in c(.0001, .9999)) {
    rate <- 4*mean(y)/(3-p)
    f <- fit(rate)
    expect_equal(unname(coef(f)), p, tolerance = 1e-12)
    expect_equal(unname(f$sd), 4/rate*sqrt(mean(y)/length(y)), tolerance = 1e-12)
  }
  expect_warning(f0 <- fit(3), "natural boundary")
  expect_warning(f1 <- fit(6), "natural boundary")
  expect_equal(c(unname(coef(f0)), unname(coef(f1))), c(0, 1))
  expect_true(all(is.na(c(f0$sd, f1$sd))))
  expect_warning(b <- fit(3, lower = qlogis(.1)), "optimizer bound")
  expect_equal(coef(b), c(prob = .1))

  ## Fixed zero inflation scales the Bernoulli success probability.
  z <- rep(0:1, c(45, 55))
  g <- distfit(z, type = 6, family = "zidegpd", fix.arg = c(base, list(pi = .2)))
  expect_equal(coef(g), c(prob = .25))
  expect_equal(unname(g$sd), sqrt(.55*.45/length(z))/.2)
})
