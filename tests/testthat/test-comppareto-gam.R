numeric_derivative_comp <- function(cdf_fun, x) {
  h <- 1e-6 * pmax(1, x)
  (cdf_fun(x + h) - cdf_fun(pmax(x - h, 0))) / (2 * h)
}

test_that("CompPareto p/q/d functions roundtrip for supported body distributions", {
  specs <- list(
    lnorm = list(meanlog = 0.4, sdlog = 0.5, alpha = 1.4, theta = 2),
    gamma = list(shape = 1.8, scale = 1.1, alpha = 1.4, theta = 2),
    weibull = list(shape = 1.5, scale = 2.2, alpha = 1.4, theta = 2),
    exp = list(rate = 0.8, alpha = 1.4, theta = 1.6)
  )
  probs <- c(0.15, 0.5, 0.85)

  for (spec in names(specs)) {
    pars <- specs[[spec]]
    q <- do.call(qcomppareto, c(list(p = probs, spec = spec), pars))
    p_back <- do.call(pcomppareto, c(list(q = q, spec = spec), pars))
    d_num <- numeric_derivative_comp(
      function(x) do.call(pcomppareto, c(list(q = x, spec = spec), pars)),
      q
    )
    d_val <- do.call(dcomppareto, c(list(x = q, spec = spec), pars))

    expect_equal(p_back, probs, tolerance = 1e-6, info = spec)
    expect_true(
      max(abs(d_num - d_val) / pmax(abs(d_val), 1e-8)) < 1e-4,
      info = spec
    )
  }
})

test_that("CompPareto distribution functions match fixed upstream reference values", {
  expect_equal(
    dcomppareto(5, spec = "lnorm", alpha = 1.7, theta = 3.2, meanlog = 0.4, sdlog = 0.7),
    0.040855346910105825,
    tolerance = 1e-12
  )
  expect_equal(
    pcomppareto(5, spec = "lnorm", alpha = 1.7, theta = 3.2, meanlog = 0.4, sdlog = 0.7),
    0.80293303255125403,
    tolerance = 1e-12
  )
  expect_equal(
    qcomppareto(0.9, spec = "lnorm", alpha = 1.7, theta = 3.2, meanlog = 0.4, sdlog = 0.7),
    9.0212634851351687,
    tolerance = 1e-10
  )

  expect_equal(
    dcomppareto(7, spec = "gamma", alpha = 2.3, theta = 4.5, shape = 1.8, scale = 1.1),
    0.019833069959575306,
    tolerance = 1e-12
  )
  expect_equal(
    pcomppareto(7, spec = "gamma", alpha = 2.3, theta = 4.5, shape = 1.8, scale = 1.1),
    0.9008346502021235,
    tolerance = 1e-12
  )
  expect_equal(
    qcomppareto(0.9, spec = "gamma", alpha = 2.3, theta = 4.5, shape = 1.8, scale = 1.1),
    6.9581686205196984,
    tolerance = 1e-10
  )

  expect_equal(
    dcomppareto(4, spec = "weibull", alpha = 1.4, theta = 2.7, shape = 1.3, scale = 2.2),
    0.073252460688395707,
    tolerance = 1e-12
  )
  expect_equal(
    pcomppareto(4, spec = "weibull", alpha = 1.4, theta = 2.7, shape = 1.3, scale = 2.2),
    0.64943465241982057,
    tolerance = 1e-12
  )
  expect_equal(
    qcomppareto(0.9, spec = "weibull", alpha = 1.4, theta = 2.7, shape = 1.3, scale = 2.2),
    13.713284836512656,
    tolerance = 1e-10
  )

  expect_equal(
    dcomppareto(3, spec = "exp", alpha = 1.1, theta = 1.9, rate = 0.8),
    0.074031763139675549,
    tolerance = 1e-12
  )
  expect_equal(
    pcomppareto(3, spec = "exp", alpha = 1.1, theta = 1.9, rate = 0.8),
    0.67022214601417263,
    tolerance = 1e-12
  )
  expect_equal(
    qcomppareto(0.9, spec = "exp", alpha = 1.1, theta = 1.9, rate = 0.8),
    12.597947516197779,
    tolerance = 1e-10
  )
})

test_that("CompPareto GAM fits an intercept-only lognormal-body model", {
  set.seed(21)
  y <- rcomppareto(
    180,
    spec = "lnorm",
    meanlog = 0.2,
    sdlog = 0.35,
    alpha = 1.5,
    theta = 2
  )
  df <- data.frame(y = y)

  fit <- suppressMessages(
    egpd(
      list(meanlog = y ~ 1, logsdlog = ~1, logalpha = ~1, logtheta = ~1),
      data = df,
      family = "comppareto",
      comppareto.args = list(spec = "lnorm")
    )
  )

  expect_s3_class(fit, "egpd")
  expect_equal(fit$family, "comppareto")
  expect_equal(fit$convergence, 0)

  pars <- predict(fit, type = "response")
  qhat <- predict(fit, type = "quantile", prob = c(0.5, 0.9))
  qref <- qcomppareto(
    c(0.5, 0.9),
    spec = "lnorm",
    meanlog = pars$meanlog[1],
    sdlog = pars$sdlog[1],
    alpha = pars$alpha[1],
    theta = pars$theta[1]
  )

  expect_equal(unname(as.numeric(qhat[1, ])), unname(qref), tolerance = 1e-6)
  expect_true(all(is.finite(rqresid(fit, seed = 1))))
})

test_that("CompPareto GAM supports smooth predictors for the exponential body", {
  set.seed(22)
  n <- 120
  x <- sort(stats::runif(n, -1, 1))
  rate <- exp(-0.1 + 0.6 * sin(pi * x))
  y <- vapply(
    seq_len(n),
    function(i) rcomppareto(1, spec = "exp", rate = rate[i], alpha = 1.3, theta = 1.5),
    numeric(1)
  )
  df <- data.frame(y = y, x = x)

  fit <- suppressMessages(
    egpd(
      list(lograte = y ~ s(x, k = 5), logalpha = ~1, logtheta = ~1),
      data = df,
      family = "comppareto",
      comppareto.args = list(spec = "exp")
    )
  )

  pred <- predict(
    fit,
    newdata = data.frame(x = c(-0.75, 0, 0.75)),
    type = "response"
  )
  qs <- predict(
    fit,
    newdata = data.frame(x = c(-0.75, 0, 0.75)),
    type = "quantile",
    prob = c(0.5, 0.9)
  )

  expect_equal(fit$convergence, 0)
  expect_equal(nrow(pred), 3L)
  expect_true(all(is.finite(as.matrix(pred))))
  expect_true(all(qs[[2]] >= qs[[1]]))
})
