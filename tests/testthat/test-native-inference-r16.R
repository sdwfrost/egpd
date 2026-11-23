test_that("GPIG profiles restore and reuse the evaluator that produced the fit", {
  freq <- c(96, 44, 21, 12, 8, 4, 3, 1)
  j <- seq_along(freq) - 1
  d <- data.frame(y = rep(j, freq))
  cc <- .7
  ## a = .5: coefficients of the pgf in terms of K_{j-1/2}(b).
  ll <- function(b) sum(freq * (b*(sqrt(1-cc)-1) + log(b) - .5*log(pi) +
    (j-.5)*log(b/2) + log(besselK(b, abs(j-.5), expon.scaled = TRUE)) +
    j*log(cc) - lgamma(j+1)))
  old <- gpig_control("saddlepoint", order = 1L, normalize = TRUE)
  on.exit(do.call(gpig_control, old))
  ambient <- egpd:::gpig_get_opts()
  f <- egpd(list(y ~ 1, ~1, ~1), d, family = "gpignat", fixed = list(a = .5, c = cc),
            gpig.args = list(method = "recursion"))
  ci <- confint(f, "b", method = "profile")
  bh <- predict(f, type = "response")$b[1]
  expect_equal(as.numeric(logLik(f)), ll(bh), tolerance = 1e-10)
  expect_lt(max(abs(2*(ll(bh) - vapply(ci, ll, 0)) - qchisq(.95, 1))), .002)
  expect_identical(egpd:::gpig_get_opts(), ambient)
  inherited <- egpd(list(y ~ 1, ~1, ~1), d, family = "gpignat", fixed = list(a = .5, c = cc))
  reference <- confint(inherited, "b", method = "profile")
  gpig_control("recursion")
  expect_equal(confint(inherited, "b", method = "profile"), reference)
  current <- egpd:::gpig_get_opts()
  broken <- inherited
  broken$likfns$d0 <- function(...) stop("profile test error")
  expect_error(confint(broken, "b", method = "profile"), "profile test error")
  expect_identical(egpd:::gpig_get_opts(), current)
  f$gpig.control <- NULL
  expect_error(confint(f, "b", method = "profile"), "refit")
})

test_that("PRL response intervals and SEs describe the effective mean", {
  d <- data.frame(x = rep(0:1, each = 300))
  u <- rep(((1:300) - .5) / 300, 2)
  d$y <- floor(-ifelse(d$x == 0, 6, 18) * log1p(-u))
  f <- egpd(y ~ x, d, family = "prl")
  nd <- data.frame(x = c(-1, .35, 40))
  X <- cbind(1, nd$x)
  eta <- drop(X %*% coef(f))
  se <- sqrt(rowSums((X %*% vcov(f)) * X))
  clamp <- function(v) pmax(4.0001, pmin(1e12, exp(v)))
  pr <- predict(f, nd, type = "response", se.fit = TRUE)
  expect_equal(pr$fitted$mu, clamp(eta))
  expect_equal(pr$se.fit$mu, c(0, exp(eta[2])*se[2], 0))
  expect_equal(unname(confint(f, "mu", newdata = nd)),
               unname(cbind(clamp(eta - qnorm(.975)*se), clamp(eta + qnorm(.975)*se))))
  expect_equal(unname(confint(f, "mu", level = 0, newdata = nd)[, 1]), clamp(eta))
  expect_equal(unname(confint(f, "mu", level = 0)[1, ]), rep(unname(coef(f)[1]), 2))
  ## Pin down the convention at the two kinks without finding x by
  ## numerical inversion (exp(log(1e12)) is not exactly 1e12).
  for (bound in c(4.0001, 1e12)) {
    f[[1]]$coefficients[1] <- log(bound)
    expect_true(is.na(predict(f, data.frame(x = 0), type = "response", se.fit = TRUE)$se.fit$mu))
  }
})
