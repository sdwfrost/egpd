test_that("both parts of a splice use the evaluated count response", {
  skip_if_not_installed("qgam")
  set.seed(2350)
  d <- data.frame(x = runif(240, -1, 1))
  d$y <- rnbinom(nrow(d), size = 2.5, mu = 6 * exp(.3 * d$x))
  d$raw <- d$y^2
  d[["sqrt(raw)"]] <- d$y + 8
  d[["case count"]] <- d$y
  fit <- function(f, data) {
    set.seed(2351)
    fit_spliced_degpd(f, data, bulk_probs = c(.25, .5, .75), threshold_prob = .75,
                      egpd_control = list(fix.arg = list(xi = 0, kappa = 1)))
  }
  ref <- fit(y ~ x, d)
  th <- predict(ref, d, type = "parameter")$threshold
  e <- d$y[d$y >= th] - th[d$y >= th]
  ll <- sum(dgeom(e, prob = 1/(1 + mean(e)), log = TRUE))
  for (form in list(sqrt(raw) ~ x, `case count` ~ x)) {
    got <- fit(form, d)
    expect_equal(predict(got, d, type = "bulk"), predict(ref, d, type = "bulk"))
    expect_equal(got$n_exceed, length(e))
    expect_equal(got$tail_fit$loglik, ll, tolerance = 1e-8)
  }
  d$raw[3] <- NA
  expect_warning(m <- fit(sqrt(raw) ~ x, d), "1 incomplete bulk")
  omit <- fit(y ~ x, d[-3, ])
  expect_equal(m$n, nrow(d) - 1L)
  expect_equal(m$tail_fit$loglik, omit$tail_fit$loglik)
})

test_that("orthogonal tail polynomials use just the complete exceedances", {
  d <- data.frame(t = seq(-1, 1, length.out = 65), e = rep(0:4, 13))
  d$t[3] <- NA_real_
  fit <- function(f, data) egpd:::.fit_spliced_degpd_tail(
    list(lsigma = f, lxi = ~1, lkappa = ~1), data,
    fit_control = list(fix.arg = list(xi = 0, kappa = 1)))
  expect_warning(a <- fit(e ~ poly(t, 2), d), "1 incomplete tail.*64 exceedances")
  b <- fit(e ~ t + I(t^2), d[-3, ])
  expect_identical(a$nobs, 64L)
  expect_equal(a$loglik, b$loglik, tolerance = 1e-7)
  b <- fit(e ~ poly(t, 2), d[-3, ])
  expect_equal(a$coefficients, b$coefficients)
})

test_that("a small positive gauge value remains a positive window", {
  tiny <- 2^-28
  x <- rep(c(0, tiny, 0, 1, 1, 2), 20)
  windows <- list(`1` = x, `2` = x[-1] + head(x, -1))
  fit <- function(v) fit_aggregate_egpd(v, c(1, 2), family = "cpegpd", cpegpd.h = 1,
    p = 0, q = 0, start = list(sigma_coef = log(2), lambda_coef = log(1.2),
                               xi = .25, kappa = 1), control = list(maxit = 0))
  expect_warning(series <- fit(x), "not on the.*lattice")
  expect_warning(pre <- fit(windows), "not on the.*lattice")
  expect_equal(unname(series$n_positive), c(80, 119))
  expect_equal(series$composite_loglik, pre$composite_loglik, tolerance = 1e-12)
  # A positive reading stays in bin one even if x/h itself rounds to zero.
  smallest <- .Machine$double.xmin * .Machine$double.eps
  mass <- dcpegpd(c(0, smallest, 2), lambda = 1, sigma = 2, xi = 0, kappa = 1, h = 2)
  expect_equal(mass, exp(-1) * c(1, -expm1(-1), -expm1(-1)))
})
