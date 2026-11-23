test_that("a hurdle subsample is shared by the Bernoulli and positive likelihoods", {
  d <- data.frame(y = c(rep(0, 40), rep(1:7, 20)))
  d$y[2] <- NA
  complete <- d[!is.na(d$y), , drop = FALSE]
  set.seed(140)
  sampled <- complete[sort(sample.int(nrow(complete), 40)), , drop = FALSE]
  set.seed(140)
  expect_warning(h <- fit_hurdle(list(y ~ 1, ~1, ~1), d, family = "degpd", maxdata = 40,
    degpd.args = list(link = "identity"), fixed = list(xi = 0, kappa = 1), trace = -1),
    "omitted 1 observation")
  expect_true(h$converged)
  expect_identical(h$data, sampled)
  expect_equal(nobs(h), 40)
  expect_equal(nobs(h$positive), sum(sampled$y > 0))
  p <- predict(h)
  pos <- sampled$y > 0
  ## With xi=0, kappa=1, Y-1 | Y>0 is geometric.
  ll <- sum(log(p$pzero[!pos])) + sum(log1p(-p$pzero[pos]) +
    dgeom(sampled$y[pos] - 1, -expm1(-1/p$scale[pos]), log = TRUE))
  expect_equal(as.numeric(logLik(h)), ll, tolerance = 1e-10)
  expect_equal(AIC(h), 2*h$df - 2*ll, tolerance = 1e-10)
})

test_that("hurdle splits the materialized response, even on a later named formula", {
  d <- data.frame(raw = rep(1:5, 15))
  times <- 0L
  shifted <- function(x) { times <<- times + 1L; x - 1L }
  h <- fit_hurdle(list(lsigma = ~1, xi = shifted(raw) ~ 1, lkappa = ~1), d,
    family = "degpd", degpd.args = list(link = "identity"), fixed = list(xi = 0, kappa = 1))
  expect_equal(times, 1)
  expect_equal(h$n_zero, sum(d$raw == 1))
  expect_equal(as.numeric(h$positive$likdata$y), d$raw[d$raw > 1] - 1L)
  expect_equal(predict(h)$pzero, rep(mean(d$raw == 1), nrow(d)), tolerance = 1e-8)
})

test_that("bounded count quantiles reach the last count, even when the lower CDF rounds to one", {
  sig <- c(1, 7.53, 10)
  last <- ceiling(20*sig) - 1
  for (fn in list(qdiscegpd, function(p, ...) qzidiscegpd(p, pi = .3, ...))) {
    q <- function(p, ...) fn(p, sigma = sig, xi = -.05, kappa = 1, ...)
    expect_equal(q(1), last)
    expect_equal(q(0, lower.tail = FALSE), last)
    expect_equal(q(0, log.p = TRUE), last)
    expect_equal(q(-Inf, lower.tail = FALSE, log.p = TRUE), last)
  }
  expect_equal(qzidiscegpd(1, pi = 1, sigma = 10, xi = -.05, kappa = 1), 0)
  ## An interior representable probability near one in a broad geometric
  ## tail: find the first count with -(y+1)/sigma <= log(1-p).
  p <- 1 - 2^-53
  expect_equal(qdiscegpd(p, sigma = 100, xi = 0, kappa = 1), ceiling(-100*log1p(-p))-1)
  expect_equal(qzidiscegpd(p, pi = .3, sigma = 100, xi = 0, kappa = 1),
               ceiling(-100*(log1p(-p)-log(.7)))-1)
})

test_that("a count above the rounded CDF yields its finite randomized normal residual", {
  d <- data.frame(y = c(rep(0:3, 50), 150))
  f <- egpd(list(y ~ 1, ~1, ~1), d, family = "degpd", degpd.args = list(link = "identity"),
            fixed = list(xi = 0, kappa = 1), restarts = FALSE)
  s <- predict(f, type = "response")$scale[1]
  set.seed(66)
  w <- runif(nrow(d))
  ## The geometric survivor interpolated using the same uniform draw.
  ls <- -d$y/s + log1p(w * expm1(-1/s))
  expect_equal(rqresid(f, seed = 66), qnorm(ls, lower.tail = FALSE, log.p = TRUE),
               tolerance = 1e-10)
})

test_that("CI levels are probabilities and known values remain points even at level one", {
  d <- data.frame(y = 1:20)
  b <- distfit(d$y, family = "bell")
  for (bad in list(-.2, 1.01, NA_real_, c(.5,.8)))
    expect_error(confint(b, level = bad), "level")
  f <- egpd(list(y ~ 1, ~1, ~1), d, fixed = list(xi = 0, kappa = 1), restarts = FALSE)
  ci <- confint(f, level = 1, method = "profile")
  expect_equal(unname(ci), rbind(c(0, Inf), c(0, 0), c(1, 1)))
  expect_equal(unname(confint(f, level = 1, newdata = d[1, , drop = FALSE])), unname(ci))
  expect_true(all(is.finite(confint(f, level = 1 - 2^-53, method = "wald"))))
})

test_that("POT stability skips nonconverged optimizations and names their thresholds", {
  set.seed(501)
  x <- rgeom(300, .25)
  expect_error(degpd_pot_stability(x, c(3,5), hessian = FALSE, control = list(maxit = 1L)),
               "no threshold could be fitted:.*threshold 3.*threshold 5")
  fit <- fit_degpd_pot
  local_mocked_bindings(fit_degpd_pot = function(x, threshold, ...)
    fit(x, threshold, ..., control = list(maxit = if (threshold == 3) 1L else 1000L)),
    .package = "egpd")
  expect_warning(st <- degpd_pot_stability(x, c(3,5), hessian = FALSE),
                 "threshold 3: optimizer did not converge")
  expect_identical(st$threshold, 5L)
})
