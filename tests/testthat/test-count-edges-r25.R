test_that("a rounded endpoint does not erase the last floor or ceiling cell", {
  s <- 64; xi <- -64/3
  # For the binary shape supplied, 1 + 3*xi/s is exactly 2^-54, not zero.
  S <- (2^-54)^(-1/xi)
  expect_equal(ddiscegpd(3, sigma=s, xi=xi, kappa=1), S)
  expect_equal(pdiscegpd(2, sigma=s, xi=xi, kappa=1, lower.tail=FALSE), S)
  expect_equal(qdiscegpd(c(.9, 1), sigma=s, xi=xi, kappa=1), c(3,3))
  set.seed(25081); E <- rexp(200)
  set.seed(25081); y <- rmdgpd(200, sigma=s, xi=xi, rate=0)
  expect_identical(y[,1] > 3, E > -log(S))
  expect_equal(ddiscegpd(3, sigma=63, xi=-21, kappa=1), 0)
})

test_that("floor draws do not generate an impossible count at an integral endpoint", {
  set.seed(25081); u <- runif(200)
  want <- (64*log1p(-u) <= log(2/3)) + (64*log1p(-u) <= log(1/3))
  # X = 3*(1-(1-U)^64), always strictly below 3 at interior uniforms.
  set.seed(25081)
  expect_equal(rdiscegpd(200, sigma=192, xi=-64, kappa=1), want)
  expect_equal(rdiscegpd(3, sigma=3, xi=-1, kappa=1, unifsamp=c(.5,.75,1)), c(1,2,2))
})

test_that("small standardized power-carrier inputs survive scaling", {
  s <- 2^900; k <- 1/512; x <- 2^-636; p <- 1/8
  expect_equal(pegpd(x, sigma=s, xi=-1, kappa=k), p)
  expect_equal(degpd_density(x, sigma=s, xi=-1, kappa=k)/(k*p/x), 1)
  expect_equal(qegpd(log1p(-p), sigma=s, xi=-1, kappa=k, lower.tail=FALSE, log.p=TRUE)/x, 1)
  expect_equal(regpd(1, sigma=s, xi=-1, kappa=k, unifsamp=p)/x, 1)
  pm <- .25*p + .75*p^2
  expect_equal(qegpd(pm, sigma=s, xi=-1, prob=.25, kappa=k, delta=2*k, type=6)/x, 1)
})

test_that("the near-critical generalized Poisson survivor has an attained inverse", {
  lam <- 1-2^-40
  # Finite-sum reference evaluated at 192 bits for the exact binary inputs.
  S <- c(2.479869581281919e-13, 2.346723962929187e-13)
  expect_equal(pgpois(8:9, 1, lam, lower.tail=FALSE)/S, c(1,1), tolerance=1e-12)
  expect_equal(pgpois(8:9, 1, lam, lower.tail=FALSE, log.p=TRUE), log(S))
  expect_equal(qgpois(2.4e-13, 1, lam, lower.tail=FALSE, max.value=24), 9)
  expect_warning(q <- qgpois(2.4e-13, 1, lam, lower.tail=FALSE, max.value=8), "not attained")
  expect_true(is.na(q))
})
