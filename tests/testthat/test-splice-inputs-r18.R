test_that("without a tail only the requested bulk grid determines the distribution", {
  skip_if_not_installed("qgam")
  d <- expand.grid(x = seq(0,1,length.out = 15), j = 1:12)
  d$y <- 2 + floor((3+d$x)*qexp((d$j-.5)/12))
  fit <- function(tp) {
    set.seed(123)
    fit_spliced_degpd(y ~ x, d, bulk_probs = c(.2,.6), threshold_prob = tp, fit_tail = FALSE)
  }
  a <- fit(.9); b <- fit(.6); nd <- d[c(3,10),]
  expect_identical(a$bulk_probs, c(.2,.6))
  expect_equal(predict(a, nd, prob = c(.95,1)), predict(b, nd, prob = c(.95,1)))
  expect_true(all(is.na(predict(a, nd, type = "parameter")$threshold)))
})

test_that("incomplete spliced bulk rows do not reach the tail or disturb known forecasts", {
  skip_if_not_installed("qgam")
  d <- expand.grid(x = seq(0,1,length.out = 16), j = 1:20)
  d$y <- 2 + floor((3+d$x)*qexp((d$j-.5)/20))
  d$off <- .1*(d$j%%2)
  d$x[1:4] <- NA
  fit <- function(ctl) {
    set.seed(1806)
    fit_spliced_degpd(y ~ x + offset(off), d, bulk_probs = c(.2,.65), threshold_prob = .65,
                      keep_data = TRUE, egpd_control = ctl)
  }
  expect_warning(f <- fit(list(fix.arg = list(kappa = 1, xi = 0))), "4 incomplete bulk")
  expect_equal(f$n, nrow(d)-4)
  expect_identical(which(is.na(predict(f, type = "parameter")$threshold)), 1:4)
  nd <- d[8:11,]
  nd$off[3] <- NA
  for (type in c("parameter","pmf","cdf","quantile")) {
    pred <- function(data) predict(f, data, type = type, prob = .9, at = c(2,30))
    a <- pred(nd)
    expect_true(all(is.na(a[3,])))
    expect_equal(a[-3,,drop = FALSE], pred(nd[-3,,drop = FALSE]))
    expect_true(all(is.na(pred(nd[3,,drop = FALSE]))))
  }
})

test_that("unknown tail controls are errors, not ignored pins or starts", {
  d <- data.frame(e = rep(0:4,10))
  f <- list(lsigma = e ~ 1, lxi = ~1, lkappa = ~1)
  fit <- function(ctl) egpd:::.fit_spliced_degpd_tail(f,d,ctl)
  expect_error(fit(list(fixed = list(xi = 0))), "egpd_control.*fix.arg")
  expect_error(fit(list(start = list(sigma = 11))), "start.*beta")
  expect_error(fit(list(start = list(beta = NA_real_))), "finite coefficients")
  expect_error(fit(list(fix.arg = list(kappa = -1))), "kappa")
  good <- fit(list(start = list(beta = log(11)), fix.arg = list(xi = 0,kappa = 1),
                   control = list(maxit = 0)))
  expect_equal(exp(unname(good$coefficients)), 11)
})
