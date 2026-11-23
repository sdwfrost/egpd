test_that("binding preserves each observation and its matrix and time columns", {
  d <- data.frame(id = letters[1:5], y = c(2,1,7,3,4))
  d$t <- as.POSIXct("2020-06-01", tz = "UTC") + 3600*(1:5)
  d$M <- I(cbind(u = 1:5, v = 11:15))
  z <- dfbind(list(a = d[1:2,], b = d[3:5,]))
  expect_identical(z$id, d$id)
  expect_identical(z$y, d$y)
  expect_identical(z$t, d$t)
  expect_equal(unname(z$M), unname(d$M))
})

test_that("hurdle selection uses the same evaluated rows for both components", {
  d <- data.frame(g = rep(c("a","b","other"), each = 36))
  d$y <- rep(c(0,1,2,4,0,2),18) + ifelse(d$g == "b", rep(0:2,36), 0)
  d$group <- factor(d$g, levels = c("a","b"))
  fit <- function(fo, gate)
    fit_hurdle(list(fo,~1,~1), d, zero.formula = gate, family = "degpd",
               degpd.args = list(link = "identity"), fixed = list(xi = 0, kappa = 1))
  for (gate in c(TRUE,FALSE)) {
    expect_warning(a <- if (gate) fit(y ~ 1, ~factor(g, levels = c("a","b"))) else
                   fit(y ~ factor(g, levels = c("a","b")), ~1), "omitted 36")
    expect_warning(b <- if (gate) fit(y ~ 1, ~group) else fit(y ~ group, ~1), "omitted 36")
    expect_true(a$converged)
    expect_equal(nobs(a), 72)
    expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-9)
  }
})
