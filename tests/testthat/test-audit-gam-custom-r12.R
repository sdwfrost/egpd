## An evgam-style Gaussian location/log-variance likelihood. The callback sees
## the compact design, but y, offsets and the returned score rows are not compact.
r12_eta <- function(b, ld, i) {
  eta <- drop(ld$X[[i]] %*% b[ld$idpars == i])
  if (ld$duplicate == 1L) eta <- eta[ld$dupid + 1L]
  off <- ld$offsets[[i]]
  if (length(off)) eta + off else eta
}
r12_normal <- list(
  d0 = function(b, ld)
    -sum(dnorm(drop(ld$y), r12_eta(b, ld, 1), exp(r12_eta(b, ld, 2)/2), log = TRUE)),
  d120 = function(b, ld) {
    r <- drop(ld$y) - r12_eta(b, ld, 1)
    iv <- exp(-r12_eta(b, ld, 2))
    cbind(-r*iv, (1-r^2*iv)/2, iv, r*iv, r^2*iv/2)
  })

test_that("custom variance pins and formula offsets reach likelihoods in y order", {
  set.seed(82)
  d <- data.frame(x = sample(0:4, 240, TRUE), off = sample(c(-.75, 1.25), 240, TRUE))
  d$y <- 2.5 + .4*d$x + d$off + 1.6*rnorm(nrow(d))
  form <- list(location = y ~ x + offset(off), variance = ~ 1)
  fns <- c(r12_normal, list(unlink = list(NULL, exp)))
  ref <- lm(form[[1]], d)
  ll <- sum(dnorm(d$y, fitted(ref), 2, log = TRUE))
  for (compact in c(TRUE, FALSE)) {
    f <- egpd(form, d, family = "custom", custom.fns = fns, fixed = list(variance = 4),
              compact = compact)
    expect_identical(f$convergence, 0L)
    expect_identical(nrow(f$likdata$X[[1L]]) < nrow(d), compact)
    expect_equal(as.numeric(logLik(f)), ll)
    expect_equal(unname(predict(f)$location), unname(fitted(ref)))
    nd <- d[c(3L, 17L), , drop = FALSE]
    expect_equal(predict(f, nd)$variance, rep(log(4), 2))
    expect_equal(predict(f, nd, type = "response")$variance, c(4, 4))
    expect_equal(unname(confint(f, "variance", newdata = nd)), matrix(4, 2, 2))
    expect_equal(unname(simulate(f, 3, seed = 9, newdata = nd, type = "response")$variance),
                 matrix(4, 2, 3))
  }
})

test_that("only a link prefix, not any leading l, supplies a default custom link", {
  set.seed(82)
  d <- data.frame(y = 2.5 + 1.6*rnorm(240))
  form <- list(location = y ~ 1, logvariance = ~ 1)
  fit <- function(...) egpd(form, d, family = "custom", custom.fns = r12_normal, ...)
  v <- fit(fixed = list(variance = 4))
  expect_equal(as.numeric(logLik(v)), sum(dnorm(d$y, mean(d$y), 2, log = TRUE)))
  f <- fit(fixed = list(location = -1))
  expected <- mean((d$y + 1)^2)
  expect_equal(predict(f, type = "response")$location, rep(-1, nrow(d)))
  expect_equal(predict(f, type = "response")$variance[1], expected, tolerance = 1e-6)

  ## Report a unit-interval transformation of the log variance instead.
  for (link in c("logit", "probit", "cloglog")) {
    names(form)[2L] <- paste0(link, "p")
    f <- fit(fixed = list(p = .8))
    expect_equal(predict(f)[[2L]][1L], make.link(link)$linkfun(.8))
    expect_equal(predict(f, type = "response")$p[1L], .8)
  }
})

test_that("custom pins honor a decreasing, restricted-domain inverse link or error", {
  set.seed(82)
  d <- data.frame(y = 2.5 + 1.6*rnorm(240))
  ## Same likelihood, but report precision and restrict it to (0, 1).
  precision <- function(e) {
    if (any(e <= 0)) stop("eta must be positive")
    exp(-e)
  }
  fit <- function(u, p) egpd(list(location = y ~ 1, logprecision = ~ 1), d,
                             family = "custom", fixed = list(precision = p),
                             custom.fns = c(r12_normal, list(unlink = list(NULL, u))))
  f <- fit(precision, .25)
  expect_equal(predict(f)$logprecision[1], log(4))
  expect_equal(predict(f, type = "response")$precision[1], .25)
  expect_equal(as.numeric(logLik(f)), sum(dnorm(d$y, mean(d$y), 2, log = TRUE)))
  expect_error(fit(precision, 2), "fixed.*precision.*unlink")
  ## A sign change alone does not establish that a pin is attainable.
  expect_error(fit(function(e) ifelse(e < 0, e, 2 + e), 1), "fixed.*precision.*unlink")
  narrow <- function(e) {
    if (any(abs(e - log(4)) > 1e-7)) stop("eta outside domain")
    exp(-e)
  }
  attr(narrow, "link") <- function(p) -log(p)
  pinned <- fit(narrow, .25)
  expect_equal(as.numeric(logLik(pinned)), as.numeric(logLik(f)))
  expect_equal(predict(pinned, type = "response", se.fit = TRUE)$se.fit$precision[1], 0)
  attr(narrow, "link") <- function(p) -log(p) + 1e-8
  expect_error(fit(narrow, .25), "fixed.*precision.*unlink")
})
