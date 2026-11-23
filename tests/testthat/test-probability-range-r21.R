test_that("nonzero subnormal shapes retain the GPD exponential limit", {
  tiny <- .Machine$double.xmin*.Machine$double.eps
  xi <- c(-tiny,0,tiny)
  p <- -expm1(-.25)
  expect_equal(pegpd(.25, kappa=1, sigma=1, xi=xi), rep(p,3))
  expect_equal(degpd_density(.25, kappa=1, sigma=1, xi=xi), rep(exp(-.25),3))
  expect_equal(qegpd(p, kappa=1, sigma=1, xi=xi), rep(.25,3))
  expect_equal(ddiscegpd(0, kappa=1, sigma=2, xi=xi), rep(-expm1(-.5),3))
})

test_that("GPD inverses avoid intermediate overflow when the answer is representable", {
  q <- qegpd(-c(log(2),16), sigma=c(1e308,1e-60), xi=c(2,50), kappa=1,
              lower.tail=FALSE, log.p=TRUE)
  want <- c(1.5e308,exp(800 + log(1e-60/50)))
  expect_equal(q/want, c(1,1), tolerance=2e-13)
  expect_equal(pegpd(q, sigma=c(1e308,1e-60), xi=c(2,50), kappa=1,
                     lower.tail=FALSE, log.p=TRUE), -c(log(2),16))
  expect_equal(qegpd(.5, sigma=.Machine$double.xmax, xi=2, kappa=1), Inf)
  expect_equal(qegpd(1, sigma=8, xi=-.5, kappa=1), 16)
})

test_that("standalone small-delta densities agree with the analytical derivative", {
  u <- c(.25,.5,.75); A <- -log1p(-u); d <- 1e-310
  # For small delta, G4 = 1-exp(-A)*(1+A); g4 = A to working precision.
  G <- -expm1(-A+log1p(A))
  expect_equal(d.G(u, type=4, delta=d), A, tolerance=1e-13)
  expect_equal(d.G(u, type=5, delta=d, kappa=3), 1.5*sqrt(G)*A, tolerance=1e-13)
  expect_equal(d.G(c(0,1), type=4, delta=d, log=TRUE), c(-Inf,-log(d)))
})

test_that("exact count boundaries distinguish adjacent upper-tail probabilities", {
  m <- c(4,14,27); s <- 2^(2-2*m); k <- 2^m-3
  Q <- function(p,lower) qdiscegpd(p, kappa=1, sigma=1, xi=.5, lower.tail=lower)
  expect_equal(Q(1-s,TRUE), k)
  expect_equal(Q(s,FALSE), k)
  expect_equal(Q(s*(1-2^-53),FALSE), k+1)
  expect_equal(Q(s*(1+2^-52),FALSE), k)
  expect_equal(qzidiscegpd(.5*s,pi=.5,kappa=1,sigma=1,xi=.5,lower.tail=FALSE), k)
  # With zero shift and zero bulk probability a splice has precisely this tail.
  qsplice <- egpd:::.spliced_degpd_row_quantile(1-s, -1, .75, .75, TRUE, 0, 0, 1, .5, 1)
  expect_equal(qsplice, k)
})
