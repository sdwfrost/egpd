test_that("bounded continuous and floored means resolve a narrow endpoint region", {
  s <- 2^26
  k <- 30000
  p <- list(sigma = s, xi = -1, kappa = k)
  ## X/s ~ Beta(k, 1). For the floor, sum the last terms of the CDF
  ## series: the omitted terms together are at most s*(1-n/s)^k.
  mc <- s * (k/(k + 1))
  n <- ceiling(40*s/k)
  md <- s - 1 - sum(exp(k * log1p(-seq_len(n)/s)))
  expect_lt(s * exp(k * log1p(-n/s)), 1e-8)
  expect_equal(egpd_bamlss(xi = "identity")$mu(p)/mc, 1, tolerance = 1e-10)
  floor_mean <- degpd_bamlss(xi = "identity")$mu(p)
  expect_lt(abs(floor_mean/md - 1), 1e-7)
  expect_gt(floor_mean, mc - 1)
  expect_lte(floor_mean, mc)

  ## A mixture of this power and a uniform, including its endpoint weights.
  p$prob <- c(0, .35, 1)
  p$delta <- 1
  expect_equal(egpd_bamlss(2)$mu(p)/(p$prob*mc + (1-p$prob)*s/2),
               rep(1, 3), tolerance = 1e-10)
  expect_equal(degpd_bamlss(2)$mu(p)/(p$prob*md + (1-p$prob)*(s-1)/2),
               rep(1, 3), tolerance = 1e-7)
})

test_that("bounded integration also resolves non-power carriers and short supports", {
  ## For the truncated-normal carrier and xi=-1, 1-X/s is half-normal
  ## truncated at 1. Its mean is the deficit from the scaled endpoint.
  s <- 2^26
  k <- 1e10
  p <- list(sigma = s, xi = -1, kappa = k)
  deficit <- sqrt(2/(pi*k)) * (-expm1(-k/2)) / (2*pnorm(sqrt(k)) - 1)
  mc <- s * (1 - deficit)
  expect_equal(egpd_bamlss(5)$mu(p)/mc, 1, tolerance = 1e-10)
  expect_lt(abs(degpd_bamlss(5)$mu(p) - mc), 1)

  ## With xi=-2 the endpoint is sigma/2; check integral and lattice
  ## convention when that endpoint is exactly an integer or is below one.
  p <- list(sigma = c(64, 1), xi = -2, kappa = 2)
  expect_equal(egpd_bamlss()$mu(p), p$sigma/2 * (1 - 1/6))
  j <- 1:31
  expect_equal(degpd_bamlss()$mu(p), c(sum(1 - (1 - sqrt(1 - j/32))^2), 0))
})

test_that("zero inflation keeps representable means even when the component overflows", {
  ## For kappa=1 the unit-scale mean is 1/(1-xi). At these scales the
  ## floor correction (at most 1-pi) is smaller than one ulp of the mean.
  p <- list(sigma = c(1e308, 1e308, 1e300, .Machine$double.xmax),
            xi = c(.5, .5, 1-2^-48, 1-2^-30), kappa = 1,
            pi = c(.5, .875, 1-2^-48, 1-2^-30))
  ref <- p$sigma * ((1-p$pi)/(1-p$xi))
  for (fam in list(ziegpd_bamlss(), zidegpd_bamlss())) {
    expect_equal(fam$mu(p)/ref, rep(1, 4), tolerance = 1e-10)
  }

  ## Carrier m=3, delta=1 has G(u)=u^2. Its survivor in GPD hazard t
  ## is 2*exp(-t)-exp(-2*t), so this is not just the GPD submodel.
  p <- list(sigma = 1.6e308, xi = c(-1e-12, 0, .5), delta = 1, pi = .75)
  ref <- (p$sigma*(1-p$pi)) * (2/(1-p$xi) - 1/(2-p$xi))
  for (fam in list(ziegpd_bamlss(3), zidegpd_bamlss(3))) {
    expect_equal(fam$mu(p)/ref, rep(1, 3), tolerance = 1e-10)
  }

  ## Even the unit-scale component can overflow. For the power carrier
  ## 1 + xi*E[X] = Gamma(1-xi)*Gamma(k+1)/Gamma(k+1-xi) at sigma=1;
  ## the gamma ratio is k^xi to relative O(1/k).
  rate <- 2^-28
  p <- list(sigma = 1, xi = 1-rate, kappa = 1e308, pi = 1-rate)
  ref <- (exp(lgamma(1+rate) + p$xi*log(p$kappa)) - rate)/p$xi
  for (fam in list(ziegpd_bamlss(), zidegpd_bamlss())) {
    expect_equal(fam$mu(p)/ref, 1, tolerance = 1e-10)
  }

  ## Ordinary-scale control: S(j)=(1+j)^-2, so the floor cannot be
  ## replaced by the continuous value before applying the mixture weight.
  p <- list(sigma = .5, xi = .5, kappa = 1, pi = .5)
  expect_equal(ziegpd_bamlss()$mu(p), .5)
  expect_equal(zidegpd_bamlss()$mu(p), .5*(pi^2/6 - 1), tolerance = 1e-7)
})
