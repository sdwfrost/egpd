## Round-6 regression (gwwarn): a targeted warning where the generalized Waring tail is approximate.
## pgwaring(), qgwaring() and rgwaring() take tail probabilities beyond the pmf ladder from the
## mixture integral gw_lcdf(). Mapped against 512-bit Rmpfr sums and an independent quadrature
## (files/egpd-fix-r6/gwwarn), it is exact to about 1e-11 except for k > 1e15 (beyond the fitting
## bound) evaluated at counts above 1000. There P(Y > q) was off by up to about 1e-3, so the
## functions warn (class egpd_gwaring_accuracy), once per call, and return the values unchanged.
## P(Y <= q) is off by at most ~1e-3 P(Y > q) there, so calls that use only lower-tail
## probabilities on the probability scale (and rgwaring) warn only where P(Y > q) > 1e-8.

test_that("calls that use the approximate tail warn once, with class egpd_gwaring_accuracy (gwwarn)", {
  ## (mu, k, rho) = (1e-3, 1e20, 30): log S(1e6) by exact 512-bit summation is -46.2501996617448;
  ## the integral gives -46.24962, 5.8e-4 off in S -- hence the warning
  expect_warning(v <- pgwaring(1e6, 1e-3, 1e20, 30, lower.tail = FALSE, log.p = TRUE),
                 "k = 1e\\+20 > 1e15.*count of 1000000 > 1000.*up to about 1e-3.*approximate",
                 class = "egpd_gwaring_accuracy")
  expect_equal(v, -46.2501996617448, tolerance = 1e-4)
  ## the value is returned unchanged
  expect_identical(v, egpd:::.count_p(egpd:::.gwaring_fam, 1e6, list(mu = 1e-3, k = 1e20, rho = 30),
                                      FALSE, TRUE))
  ## the lower tail on the log scale, log P(Y <= q) ~ -P(Y > q), carries the same relative error
  expect_warning(pgwaring(1e6, 1e-3, 1e20, 30, log.p = TRUE), class = "egpd_gwaring_accuracy")
  expect_warning(pgwaring(2e3, 1, 1e16, 3, lower.tail = FALSE), class = "egpd_gwaring_accuracy")
  ## one warning for a call whose parameter rows use the integral more than once
  nw <- 0
  withCallingHandlers(pgwaring(c(2e3, 1e6), mu = c(1, 2), k = 1e20, rho = 3, lower.tail = FALSE),
                      egpd_gwaring_accuracy = function(w) {
                        nw <<- nw + 1
                        invokeRestart("muffleWarning")
                      })
  expect_identical(nw, 1)
  ## quantiles: the search extends past 1000 counts for a far target (k = 1e17)
  expect_warning(q <- qgwaring(5e-16, 1, 1e17, 3, lower.tail = FALSE),
                 class = "egpd_gwaring_accuracy")
  expect_true(is.finite(q) && q > 1000)
  ## lower-tail probabilities on the probability scale, their quantiles and random draws are off
  ## by at most ~1e-3 P(Y > q) at the evaluated count: they warn where that exceeds 1e-8. At
  ## (mu, k, rho) = (300, 1e16, 1e5), P(Y > 1200) = 5.3e-8, and the searches start at 4 mu = 1200
  expect_warning(pgwaring(1200, 300, 1e16, 1e5), class = "egpd_gwaring_accuracy")
  expect_warning(qgwaring(0.5, 300, 1e16, 1e5), class = "egpd_gwaring_accuracy")
  set.seed(1)
  expect_warning(rgwaring(3, mu = 300, k = 1e16, rho = 1e5), class = "egpd_gwaring_accuracy")
})

test_that("validated evaluations do not warn (gwwarn)", {
  ## inside the fitting bounds, incl. k = 1e15 and rho = 1e15, counts up to 1e5
  for (k in c(0.5, 2, 1e8, 1e15)) for (rho in c(1.5, 3, 1e3, 1e15)) for (mu in c(1e-3, 1, 30)) {
    expect_no_warning(pgwaring(c(0, 1, 10, 1000, 1e5), mu, k, rho, lower.tail = FALSE, log.p = TRUE))
    expect_no_warning(pgwaring(c(0, 1, 10, 1000, 1e5), mu, k, rho))
    expect_no_warning(qgwaring(c(0.1, 0.5, 0.9), mu, k, rho))
  }
  ## at k = 1e15, the bound, even for a count of 1e6
  expect_no_warning(pgwaring(1e6, 1e-3, 1e15, 30, lower.tail = FALSE, log.p = TRUE))
  ## the round-5 validation grid: k up to 1e20 with counts up to 1000
  for (k in c(1e13, 1e15, 1e16, 1e18, 1e20)) for (rho in c(1.5, 3, 30)) for (mu in c(1e-3, 1, 1e3))
    expect_no_warning(pgwaring(c(0, 1, 10, 1000), mu, k, rho, lower.tail = FALSE, log.p = TRUE))
  ## draws and zero-truncated quantiles at ordinary parameters
  set.seed(2)
  expect_no_warning(rgwaring(200, mu = 3, k = 2, rho = 3))
  expect_no_warning(egpd:::.zt_quantile("gwaring", list(mu = rep(3, 3), k = rep(2, 3), rho = rep(3, 3)),
                                        c(0.5, 0.9, 0.999)))
  ## a huge mean never reaches the integral (its tail stays above 1e-3 within any feasible
  ## ladder), and the integral's lower tail, approximate at k = 1e12 on its direct path, is used
  ## by no exported function: neither warns
  expect_no_warning(pgwaring(c(0, 10, 1000), 1e40, 1e20, 2, log.p = TRUE))
  expect_no_warning(pgwaring(c(0, 10, 1000), 1e13, 1e12, 30, log.p = TRUE))
  ## lower-tail probabilities on the probability scale where the integral's P(Y > q) is at most
  ## 1e-8 (here 5.5e-15, 7.3e-16 and 8.2e-21): accurate to ~1e-11 whatever its relative error
  expect_no_warning(pgwaring(2e3, 1, 1e16, 3))
  expect_no_warning(pgwaring(1e6, 1e-3, 1e20, 30))
  expect_no_warning(qgwaring(c(0.5, 0.999), 1e3, 1e20, 3))
  set.seed(3)
  expect_no_warning(rgwaring(3, mu = 1e3, k = 1e20, rho = 3))
})
