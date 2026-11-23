test_that("compaction respects exact numeric and matrix-valued predictors", {
  d <- expand.grid(t = (-5:5)/8, g = c(-1, 1), i = 1:15)
  d$time <- 2^30 + d$t/65536
  d$M <- I(cbind(t = d$t, g = d$g))
  u <- qexp(((1:15)-.5)/15)
  d$y <- exp(.4 + .9*d$t - .3*d$g) * u[d$i]/mean(u)
  for (fo in list(y ~ I((time-1073741824)*65536) + g, y ~ M)) {
    f <- egpd(list(fo, ~1, ~1), d, family = "egpd", compact = TRUE,
              egpd.args = list(link = "identity"), fixed = list(xi = 0, kappa = 1),
              restarts = FALSE)
    ## Each cell's mean is exactly exp(.4 + .9*t - .3*g).
    expect_equal(unname(coef(f)[1:3]), c(.4, .9, -.3), tolerance = 1e-5)
    expect_equal(length(unique(f$compactid)), 22)
    scale <- exp(.4 + .9*d$t - .3*d$g)
    expect_equal(as.numeric(logLik(f)), sum(-log(scale) - d$y/scale), tolerance = 1e-9)
  }
})

test_that("spliced bulk quantiles include every offset", {
  skip_if_not_installed("qgam")
  d <- expand.grid(a = c(0,5), b = c(0,12), i = 1:35)
  d$x <- (d$i-17)/17
  d$y <- 20 + round(4*d$x) + d$a + d$b + rep(-2:2, length.out = nrow(d))
  fit <- function(fo) {
    set.seed(1701)
    fit_spliced_degpd(fo, d, bulk_probs = c(.25,.5,.75), threshold_prob = .75, fit_tail = FALSE)
  }
  nd <- data.frame(a = 0, b = c(0,12), x = 0)
  ref <- predict(fit(y ~ x + offset(a+b)), nd, type = "bulk")
  for (fo in list(y ~ x + offset(a) + offset(b), y ~ x + offset(b) + offset(a))) {
    got <- predict(fit(fo), nd, type = "bulk")
    expect_equal(got, ref, tolerance = 1e-7)
    expect_equal(as.numeric(got[2,] - got[1,]), rep(12,3))
  }
})

test_that("spliced forecasts preserve the locations of estimable rows", {
  skip_if_not_installed("qgam")
  d <- data.frame(g = factor(rep(c("A","B"), each = 90)), y = rep(3:7,36))
  d$y <- d$y + 60*(d$g == "B")
  set.seed(1701)
  f <- fit_spliced_degpd(y ~ g, d, bulk_probs = c(.25,.5,.75), threshold_prob = .75,
                         tail_formula = ~g, egpd_control = list(fix.arg = list(xi = 0, kappa = 1)))
  nd <- data.frame(g = factor(c("B","new","A","new")), row.names = LETTERS[1:4])
  ok <- nd$g != "new"
  for (type in c("bulk","quantile","pmf","cdf","parameter")) {
    pred <- function(x) predict(f, x, type = type, prob = c(.5,.9), at = c(6,65))
    expect_warning(mixed <- pred(nd), "new.*NA")
    expect_equal(mixed[ok, , drop = FALSE], pred(nd[ok, , drop = FALSE]))
    expect_true(all(is.na(mixed[!ok, , drop = FALSE])))
    expect_warning(none <- pred(nd[!ok, , drop = FALSE]), "new.*NA")
    expect_equal(none, mixed[!ok, , drop = FALSE])
  }
})
