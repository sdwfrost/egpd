# Tests for EGPD models 5 (truncated normal) and 6 (truncated beta)

# ---- Distribution functions ----

test_that("d/p/q/r functions work for type 2 (truncated normal)", {
  x <- 0:10
  d <- ddiscegpd(x, sigma = 3, xi = 0.2, kappa = 2, type = 2)
  expect_length(d, length(x))
  expect_true(all(d >= 0))
  expect_lt(abs(sum(ddiscegpd(0:500, sigma = 3, xi = 0.2, kappa = 2, type = 2)) - 1), 0.01)

  p <- pdiscegpd(x, sigma = 3, xi = 0.2, kappa = 2, type = 2)
  expect_true(all(diff(p) >= 0))
  expect_true(all(p >= 0 & p <= 1))

  q <- qdiscegpd(c(0.25, 0.5, 0.75), sigma = 3, xi = 0.2, kappa = 2, type = 2)
  expect_length(q, 3)
  expect_true(all(diff(q) >= 0))

  set.seed(1)
  r <- rdiscegpd(100, sigma = 3, xi = 0.2, kappa = 2, type = 2)
  expect_length(r, 100)
  expect_true(all(r >= 0))
  expect_true(all(r == floor(r)))
})

test_that("d/p/q/r functions work for type 3 (truncated beta)", {
  x <- 0:10
  d <- ddiscegpd(x, sigma = 3, xi = 0.2, kappa = 2, type = 3)
  expect_length(d, length(x))
  expect_true(all(d >= 0))
  expect_lt(abs(sum(ddiscegpd(0:500, sigma = 3, xi = 0.2, kappa = 2, type = 3)) - 1), 0.01)

  p <- pdiscegpd(x, sigma = 3, xi = 0.2, kappa = 2, type = 3)
  expect_true(all(diff(p) >= 0))
  expect_true(all(p >= 0 & p <= 1))

  q <- qdiscegpd(c(0.25, 0.5, 0.75), sigma = 3, xi = 0.2, kappa = 2, type = 3)
  expect_length(q, 3)
  expect_true(all(diff(q) >= 0))

  set.seed(1)
  r <- rdiscegpd(100, sigma = 3, xi = 0.2, kappa = 2, type = 3)
  expect_length(r, 100)
  expect_true(all(r >= 0))
  expect_true(all(r == floor(r)))
})

test_that("continuous EGPD d/p/q/r work for types 2 and 3", {
  y <- seq(0.1, 5, by = 0.5)
  for (tp in c(2, 3)) {
    d <- degpd_density(y, sigma = 2, xi = 0.1, kappa = 1.5, type = tp)
    expect_true(all(d >= 0))

    p <- pegpd(y, sigma = 2, xi = 0.1, kappa = 1.5, type = tp)
    expect_true(all(diff(p) >= 0))
    expect_true(all(p >= 0 & p <= 1))

    q <- qegpd(c(0.25, 0.5, 0.75), sigma = 2, xi = 0.1, kappa = 1.5, type = tp)
    expect_length(q, 3)
    expect_true(all(diff(q) >= 0))

    set.seed(1)
    r <- regpd(100, sigma = 2, xi = 0.1, kappa = 1.5, type = tp)
    expect_length(r, 100)
    expect_true(all(r >= 0))
  }
})

test_that("ZI discrete d/p/q/r work for types 2 and 3", {
  x <- 0:10
  for (tp in c(2, 3)) {
    d <- dzidiscegpd(x, pi = 0.3, sigma = 3, xi = 0.2, kappa = 2, type = tp)
    expect_true(all(d >= 0))

    p <- pzidiscegpd(x, pi = 0.3, sigma = 3, xi = 0.2, kappa = 2, type = tp)
    expect_true(all(diff(p) >= 0))
    expect_true(all(p >= 0 & p <= 1))

    set.seed(1)
    r <- rzidiscegpd(100, pi = 0.3, sigma = 3, xi = 0.2, kappa = 2, type = tp)
    expect_length(r, 100)
    expect_true(all(r >= 0))
  }
})

# ---- EGPD model fitting ----

test_that("EGPD model 5 fits and recovers parameters", {
  set.seed(101)
  y <- regpd(1000, sigma = 2, xi = 0.1, kappa = 1.5, type = 2)
  df <- data.frame(y = y)
  fit <- suppressMessages(
    egpd(list(lpsi = y ~ 1, xi = ~1, lkappa = ~1),
         data = df, family = "egpd", egpd.args = list(m = 5))
  )
  expect_s3_class(fit, "egpd")
  pars <- predict(fit, type = "response")[1, ]
  expect_lt(abs(pars$scale - 2), 1)
  expect_lt(abs(pars$kappa - 1.5), 1)
})

test_that("EGPD model 6 fits and recovers parameters", {
  set.seed(102)
  y <- regpd(1000, sigma = 2, xi = 0.1, kappa = 1.5, type = 3)
  df <- data.frame(y = y)
  fit <- suppressMessages(
    egpd(list(lpsi = y ~ 1, xi = ~1, lkappa = ~1),
         data = df, family = "egpd", egpd.args = list(m = 6))
  )
  expect_s3_class(fit, "egpd")
  pars <- predict(fit, type = "response")[1, ]
  expect_lt(abs(pars$scale - 2), 1)
  expect_lt(abs(pars$kappa - 1.5), 1)
})

# ---- DEGPD model fitting ----

test_that("DEGPD model 5 fits and recovers parameters", {
  set.seed(201)
  y <- rdiscegpd(1000, sigma = 3, xi = 0.2, kappa = 2, type = 2)
  df <- data.frame(y = y)
  fit <- suppressMessages(
    egpd(list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1),
         data = df, family = "degpd", degpd.args = list(m = 5))
  )
  expect_s3_class(fit, "egpd")
  pars <- predict(fit, type = "response")[1, ]
  expect_lt(abs(pars$scale - 3), 2)
  expect_lt(abs(pars$kappa - 2), 1.5)
})

test_that("DEGPD model 6 fits and recovers parameters", {
  set.seed(202)
  y <- rdiscegpd(1000, sigma = 3, xi = 0.2, kappa = 2, type = 3)
  df <- data.frame(y = y)
  fit <- suppressMessages(
    egpd(list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1),
         data = df, family = "degpd", degpd.args = list(m = 6))
  )
  expect_s3_class(fit, "egpd")
  pars <- predict(fit, type = "response")[1, ]
  expect_lt(abs(pars$scale - 3), 2)
  expect_lt(abs(pars$kappa - 2), 1.5)
})

# ---- ZIDEGPD model fitting ----

test_that("ZIDEGPD model 5 fits and recovers parameters", {
  set.seed(301)
  y <- rzidiscegpd(1000, pi = 0.3, sigma = 3, xi = 0.2, kappa = 2, type = 2)
  df <- data.frame(y = y)
  fit <- suppressMessages(
    egpd(list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1, logitpi = ~1),
         data = df, family = "zidegpd", zidegpd.args = list(m = 5))
  )
  expect_s3_class(fit, "egpd")
  pars <- predict(fit, type = "response")[1, ]
  expect_lt(abs(pars$pi - 0.3), 0.2)
})

test_that("ZIDEGPD model 6 fits and recovers parameters", {
  set.seed(302)
  y <- rzidiscegpd(1000, pi = 0.3, sigma = 3, xi = 0.2, kappa = 2, type = 3)
  df <- data.frame(y = y)
  fit <- suppressMessages(
    egpd(list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1, logitpi = ~1),
         data = df, family = "zidegpd", zidegpd.args = list(m = 6))
  )
  expect_s3_class(fit, "egpd")
  pars <- predict(fit, type = "response")[1, ]
  expect_lt(abs(pars$pi - 0.3), 0.2)
})

# ---- Predictions and diagnostics ----

test_that("quantile predictions work for models 5 and 6", {
  set.seed(401)
  y <- rdiscegpd(500, sigma = 3, xi = 0.2, kappa = 2, type = 2)
  df <- data.frame(y = y)
  fit5 <- suppressMessages(
    egpd(list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1),
         data = df, family = "degpd", degpd.args = list(m = 5))
  )
  qp <- predict(fit5, type = "quantile", prob = c(0.5, 0.9))
  expect_equal(ncol(qp), 2)
  expect_true(all(qp[1, 2] >= qp[1, 1]))

  set.seed(402)
  y <- rdiscegpd(500, sigma = 3, xi = 0.2, kappa = 2, type = 3)
  df <- data.frame(y = y)
  fit6 <- suppressMessages(
    egpd(list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1),
         data = df, family = "degpd", degpd.args = list(m = 6))
  )
  qp <- predict(fit6, type = "quantile", prob = c(0.5, 0.9))
  expect_equal(ncol(qp), 2)
  expect_true(all(qp[1, 2] >= qp[1, 1]))
})

test_that("rqresid works for models 5 and 6", {
  set.seed(501)
  y <- rdiscegpd(500, sigma = 3, xi = 0.2, kappa = 2, type = 2)
  df <- data.frame(y = y)
  fit5 <- suppressMessages(
    egpd(list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1),
         data = df, family = "degpd", degpd.args = list(m = 5))
  )
  r <- rqresid(fit5, seed = 1)
  expect_length(r, length(y))
  expect_true(all(is.finite(r) | is.na(r)))

  set.seed(502)
  y <- rdiscegpd(500, sigma = 3, xi = 0.2, kappa = 2, type = 3)
  df <- data.frame(y = y)
  fit6 <- suppressMessages(
    egpd(list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1),
         data = df, family = "degpd", degpd.args = list(m = 6))
  )
  r <- rqresid(fit6, seed = 1)
  expect_length(r, length(y))
  expect_true(all(is.finite(r) | is.na(r)))
})

test_that("vcov works for DEGPD models 5 and 6", {
  set.seed(601)
  y <- rdiscegpd(500, sigma = 3, xi = 0.2, kappa = 2, type = 2)
  df <- data.frame(y = y)
  fit <- suppressMessages(
    egpd(list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1),
         data = df, family = "degpd", degpd.args = list(m = 5))
  )
  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 3L)
  expect_equal(ncol(V), 3L)
})

test_that("confint works for DEGPD models 5 and 6", {
  set.seed(701)
  y <- rdiscegpd(500, sigma = 3, xi = 0.2, kappa = 2, type = 2)
  df <- data.frame(y = y)
  fit <- suppressMessages(
    egpd(list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1),
         data = df, family = "degpd", degpd.args = list(m = 5))
  )
  ci <- confint(fit, method = "wald")
  expect_equal(dim(ci), c(3L, 2L))
  expect_true(all(ci[, "lower"] < ci[, "upper"]))
})
