## Test analytical derivatives against numerical derivatives
## for all gamlss family types

test_that("EGPD1 analytical derivatives match numerical", {
  y  <- c(0.5, 1.0, 2.0, 5.0, 10.0)
  mu <- rep(2.0, 5)
  sigma <- rep(0.3, 5)
  nu <- rep(1.5, 5)

  for (idx in 1:3) {
    ad <- egpd::.egpd_ad3(y, mu, sigma, nu, "egpd1", idx)
    nd <- egpd::.egpd_nd3(y, mu, sigma, nu, dEGPD1, idx)
    expect_equal(ad, nd, tolerance = 1e-4,
                 label = paste0("EGPD1 idx=", idx))
  }
})

test_that("EGPD3 analytical derivatives match numerical", {
  y  <- c(0.5, 1.0, 2.0, 5.0, 10.0)
  mu <- rep(2.0, 5)
  sigma <- rep(0.3, 5)
  nu <- rep(1.5, 5)

  for (idx in 1:3) {
    ad <- egpd::.egpd_ad3(y, mu, sigma, nu, "egpd3", idx)
    nd <- egpd::.egpd_nd3(y, mu, sigma, nu, dEGPD3, idx)
    expect_equal(ad, nd, tolerance = 1e-4,
                 label = paste0("EGPD3 idx=", idx))
  }
})

test_that("EGPD5 analytical derivatives match numerical", {
  y  <- c(0.5, 1.0, 2.0, 5.0, 10.0)
  mu <- rep(2.0, 5)
  sigma <- rep(0.3, 5)
  nu <- rep(2.0, 5)

  for (idx in 1:3) {
    ad <- egpd::.egpd_ad3(y, mu, sigma, nu, "egpd5", idx)
    nd <- egpd::.egpd_nd3(y, mu, sigma, nu, dEGPD5, idx)
    expect_equal(ad, nd, tolerance = 1e-4,
                 label = paste0("EGPD5 idx=", idx))
  }
})

test_that("EGPD6 analytical derivatives match numerical", {
  ## Use small y so the type-3 truncated beta density stays well above 0;
  ## large y causes density underflow in dEGPD6 and numerical derivs → 0.
  y  <- c(0.1, 0.3, 0.5, 0.8, 1.0)
  mu <- rep(2.0, 5)
  sigma <- rep(0.3, 5)
  nu <- rep(2.0, 5)

  for (idx in 1:3) {
    ad <- egpd::.egpd_ad3(y, mu, sigma, nu, "egpd6", idx)
    nd <- egpd::.egpd_nd3(y, mu, sigma, nu, dEGPD6, idx)
    expect_equal(ad, nd, tolerance = 1e-3,
                 label = paste0("EGPD6 idx=", idx))
  }
})

test_that("EGPD4 analytical derivatives match numerical", {
  y   <- c(0.5, 1.0, 2.0, 5.0, 10.0)
  mu  <- rep(2.0, 5)
  sigma <- rep(0.3, 5)
  nu  <- rep(1.5, 5)
  tau <- rep(2.0, 5)

  for (idx in 1:4) {
    ad <- egpd::.egpd_ad4(y, mu, sigma, nu, tau, "egpd4", idx)
    nd <- egpd::.egpd_nd4(y, mu, sigma, nu, tau, dEGPD4, idx)
    expect_equal(ad, nd, tolerance = 1e-3,
                 label = paste0("EGPD4 idx=", idx))
  }
})

test_that("DEGPD1 analytical derivatives match numerical", {
  y  <- c(0, 1, 2, 5, 10)
  mu <- rep(2.0, 5)
  sigma <- rep(0.3, 5)
  nu <- rep(1.5, 5)

  for (idx in 1:3) {
    ad <- egpd::.egpd_ad3(y, mu, sigma, nu, "degpd1", idx)
    nd <- egpd::.egpd_nd3(y, mu, sigma, nu, dDEGPD1, idx)
    expect_equal(ad, nd, tolerance = 1e-4,
                 label = paste0("DEGPD1 idx=", idx))
  }
})

test_that("DEGPD3 analytical derivatives match numerical for mu/sigma", {
  y  <- c(0, 1, 2, 5, 10)
  mu <- rep(2.0, 5)
  sigma <- rep(0.3, 5)
  nu <- rep(1.5, 5)

  for (idx in 1:2) {
    ad <- egpd::.egpd_ad3(y, mu, sigma, nu, "degpd3", idx)
    nd <- egpd::.egpd_nd3(y, mu, sigma, nu, dDEGPD3, idx)
    expect_equal(ad, nd, tolerance = 1e-4,
                 label = paste0("DEGPD3 idx=", idx))
  }
  # idx=3 falls back to numerical, so should match exactly
  ad3 <- egpd::.egpd_ad3(y, mu, sigma, nu, "degpd3", 3L)
  nd3 <- egpd::.egpd_nd3(y, mu, sigma, nu, dDEGPD3, 3L)
  expect_equal(ad3, nd3, tolerance = 1e-6,
               label = "DEGPD3 idx=3 (numerical fallback)")
})

test_that("DEGPD5 analytical derivatives match numerical", {
  y  <- c(0, 1, 2, 5, 10)
  mu <- rep(2.0, 5)
  sigma <- rep(0.3, 5)
  nu <- rep(2.0, 5)

  for (idx in 1:3) {
    ad <- egpd::.egpd_ad3(y, mu, sigma, nu, "degpd5", idx)
    nd <- egpd::.egpd_nd3(y, mu, sigma, nu, dDEGPD5, idx)
    expect_equal(ad, nd, tolerance = 1e-4,
                 label = paste0("DEGPD5 idx=", idx))
  }
})

test_that("DEGPD6 analytical derivatives match numerical for mu/sigma", {
  y  <- c(0, 1, 2, 5, 10)
  mu <- rep(2.0, 5)
  sigma <- rep(0.3, 5)
  nu <- rep(2.0, 5)

  for (idx in 1:2) {
    ad <- egpd::.egpd_ad3(y, mu, sigma, nu, "degpd6", idx)
    nd <- egpd::.egpd_nd3(y, mu, sigma, nu, dDEGPD6, idx)
    expect_equal(ad, nd, tolerance = 1e-4,
                 label = paste0("DEGPD6 idx=", idx))
  }
})

test_that("DEGPD4 analytical derivatives match numerical for mu/sigma", {
  y   <- c(0, 1, 2, 5, 10)
  mu  <- rep(2.0, 5)
  sigma <- rep(0.3, 5)
  nu  <- rep(1.5, 5)
  tau <- rep(2.0, 5)

  for (idx in 1:2) {
    ad <- egpd::.egpd_ad4(y, mu, sigma, nu, tau, "degpd4", idx)
    nd <- egpd::.egpd_nd4(y, mu, sigma, nu, tau, dDEGPD4, idx)
    expect_equal(ad, nd, tolerance = 1e-4,
                 label = paste0("DEGPD4 idx=", idx))
  }
})

test_that("ZIEGPD1 analytical derivatives match numerical", {
  y   <- c(0, 0, 0.5, 1.0, 5.0)
  mu  <- rep(2.0, 5)
  sigma <- rep(0.3, 5)
  nu  <- rep(1.5, 5)
  tau <- rep(0.2, 5)

  for (idx in 1:4) {
    ad <- egpd::.egpd_ad4(y, mu, sigma, nu, tau, "ziegpd1", idx)
    nd <- egpd::.egpd_nd4(y, mu, sigma, nu, tau, dZIEGPD1, idx)
    expect_equal(ad, nd, tolerance = 1e-4,
                 label = paste0("ZIEGPD1 idx=", idx))
  }
})

test_that("ZIDEGPD1 analytical derivatives match numerical", {
  y   <- c(0, 0, 1, 2, 5)
  mu  <- rep(2.0, 5)
  sigma <- rep(0.3, 5)
  nu  <- rep(1.5, 5)
  tau <- rep(0.2, 5)

  for (idx in 1:4) {
    ad <- egpd::.egpd_ad4(y, mu, sigma, nu, tau, "zidegpd1", idx)
    nd <- egpd::.egpd_nd4(y, mu, sigma, nu, tau, dZIDEGPD1, idx)
    expect_equal(ad, nd, tolerance = 1e-4,
                 label = paste0("ZIDEGPD1 idx=", idx))
  }
})

test_that("family constructors produce correct derivative structure", {
  fam3 <- EGPD1()
  expect_true(is.function(fam3$dldm))
  expect_true(is.function(fam3$dldd))
  expect_true(is.function(fam3$dldv))
  expect_true(is.function(fam3$d2ldm2))

  fam4 <- EGPD4()
  expect_true(is.function(fam4$dldm))
  expect_true(is.function(fam4$dldt))
  expect_true(is.function(fam4$d2ldt2))
})

test_that("second derivatives use -(d1)^2 formula", {
  y  <- c(0.5, 1.0, 2.0, 5.0)
  mu <- rep(2.0, 4)
  sigma <- rep(0.3, 4)
  nu <- rep(1.5, 4)

  fam <- EGPD1()
  d1 <- fam$dldm(y, mu, sigma, nu)
  d2 <- fam$d2ldm2(y, mu, sigma, nu)
  expected <- ifelse(-d1 * d1 < -1e-15, -d1 * d1, -1e-15)
  expect_equal(d2, expected)
})
