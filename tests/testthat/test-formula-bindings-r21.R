test_that("a transformed count response never shadows a real predictor", {
  d <- expand.grid(g = 0:1, i = 1:32)
  d$y <- floor((2 + 6*d$g)*qexp((d$i-.5)/32))
  d$raw <- d$y^2
  d$.egpd_response <- d[["sqrt(raw)"]] <- d$g
  scale <- 1/log1p(1/as.numeric(tapply(d$y, d$g, mean)))
  exact <- sum(dgeom(d$y, -expm1(-1/scale[1+d$g]), log = TRUE))
  for (fo in list(sqrt(raw) ~ I(.egpd_response), sqrt(raw) ~ I(`sqrt(raw)`))) {
    f <- egpd(list(fo, ~1, ~1), d, family = "degpd", degpd.args = list(link = "identity"),
              fixed = list(xi = 0, kappa = 1), restarts = FALSE, trace = -1)
    expect_equal(as.numeric(logLik(f)), exact, tolerance = 1e-8)
    expect_equal(fitted(f, type = "response")$scale, scale[1+d$g], tolerance = 1e-6)
    expect_equal(predict(f, d[c(2,1),], type = "response")$scale, rev(scale), tolerance = 1e-6)
  }
})

test_that("pinning kappa preserves a same-named exposure at training and prediction", {
  d <- expand.grid(g = 0:1, i = 1:32)
  d$.fixedpar_lkappa <- log(1 + 3*d$g)
  d$y <- exp(d$.fixedpar_lkappa)*qexp((d$i-.5)/32)
  f <- egpd(list(y ~ offset(.fixedpar_lkappa), ~1, ~1), d,
            egpd.args = list(link = "identity"), fixed = list(kappa = 1, xi = 0),
            restarts = FALSE, trace = -1)
  s <- mean(d$y/exp(d$.fixedpar_lkappa))
  expect_equal(as.numeric(logLik(f)), sum(dexp(d$y, exp(-d$.fixedpar_lkappa)/s, log = TRUE)))
  p <- predict(f, data.frame(.fixedpar_lkappa = log(c(2,5))), type = "response")
  expect_equal(p$scale, s*c(2,5), tolerance = 1e-6)
  expect_equal(p$kappa, c(1,1))
})

test_that("list component names are not required or screened for missingness", {
  d <- expand.grid(x = c(0,.5), i = 1:32)
  d$y <- (1 + 6*d$x)*qexp((d$i-.5)/32)
  form <- local({cfg <- list(gain = 2); y ~ I(x * cfg$gain)})
  fit <- function(dd) egpd(list(form, ~1, ~1), dd, egpd.args = list(link = "identity"),
                           fixed = list(xi = 0, kappa = 1), restarts = FALSE, trace = -1)
  a <- fit(d)
  d$gain <- ifelse(d$i < 16, NA_real_, 11)
  b <- fit(d)
  expect_equal(nobs(b), nrow(d))
  expect_equal(coef(a), coef(b))
  p <- function(nd) predict(b, nd, type = "response")$scale
  ref <- as.numeric(tapply(d$y, d$x, mean))
  expect_equal(p(data.frame(x = c(0,.5))), ref, tolerance = 1e-6)
  expect_equal(p(data.frame(x = c(0,.5), gain = NA_real_)), ref, tolerance = 1e-6)
})
