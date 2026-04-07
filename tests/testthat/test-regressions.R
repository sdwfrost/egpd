test_that("model 2 inverse-G helpers solve their defining equations", {
  v <- c(0.2, 0.45, 0.8)
  kappa1 <- c(0.35, 0.6, 1.2)
  kappa2 <- c(1.8, 1.6, 2.4)
  prob <- c(0.25, 0.5, 0.75)

  roots_egpd <- mapply(egpd:::.iG2i_egpd, v, kappa1, kappa2, prob)
  roots_degpd <- mapply(egpd:::.iG2i_degpd, v, kappa1, kappa2, prob)
  roots_zidegpd <- mapply(egpd:::.iG2i_zidegpd, v, kappa1, kappa2, prob)

  expect_equal(
    mapply(egpd:::.G2_egpd, roots_egpd, kappa1, kappa2, prob),
    v,
    tolerance = 1e-6
  )
  expect_equal(
    mapply(egpd:::.G2_degpd, roots_degpd, kappa1, kappa2, prob),
    v,
    tolerance = 1e-6
  )
  expect_equal(
    mapply(egpd:::.G2_zidegpd, roots_zidegpd, kappa1, kappa2, prob),
    v,
    tolerance = 1e-6
  )
})

test_that("type 3 G transformations are internally consistent", {
  probs <- c(0, 0.1, 0.25, 0.5, 0.9, 1)
  quants <- q.G(probs, type = 3, kappa = 2)

  expect_equal(quants[c(1, length(quants))], c(0, 1))
  expect_equal(p.G(quants, type = 3, kappa = 2), probs, tolerance = 1e-10)
})

test_that("ZIDEGPD model 4 quantiles use the fitted kappa-delta ordering", {
  set.seed(11)
  y <- rzidiscegpd(
    1500,
    pi = 0.25,
    sigma = 2,
    xi = 0.1,
    kappa = 2.2,
    delta = 0.7,
    type = 5
  )
  df <- data.frame(y = y)

  fit <- suppressMessages(
    egpd(
      list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1, ldelta = ~1, logitpi = ~1),
      data = df,
      family = "zidegpd",
      zidegpd.args = list(m = 4)
    )
  )

  pars <- predict(fit, type = "response")[1, ]
  probs <- c(0.5, 0.9)
  qhat <- unlist(predict(fit, type = "quantile", prob = probs)[1, ])
  qref <- qzidiscegpd(
    probs,
    pi = pars$pi,
    sigma = pars$scale,
    xi = pars$shape,
    kappa = pars$kappa,
    delta = pars$delta,
    type = 5
  )

  expect_equal(unname(qhat), unname(qref))
})

test_that("compact smooth fits keep expanded fitted values and convergence metadata", {
  data(ny_complaints)
  df <- ny_complaints[rep(seq_len(nrow(ny_complaints)), each = 2), ]

  fit <- suppressMessages(
    egpd(
      list(upheld ~ s(year), ~1, ~1),
      data = df,
      family = "degpd",
      degpd.args = list(m = 1),
      compact = TRUE
    )
  )

  eta <- predict(fit, type = "link")

  expect_true(fit$compacted)
  expect_equal(fit$convergence, 0)
  expect_true(is.finite(fit$negREML))
  expect_equal(length(fit[[1]]$fitted), nrow(df))
  expect_equal(unname(fit[[1]]$fitted), unname(eta[[1]]), tolerance = 1e-8)
  expect_equal(nrow(predict(fit, newdata = df[1:5, ], type = "response")), 5L)
})

test_that("rqresid errors cleanly when model data were removed", {
  set.seed(13)
  y <- rdiscegpd(400, sigma = 3, xi = 0.2, kappa = 2, type = 1)
  df <- data.frame(y = y)

  fit <- suppressMessages(
    egpd(
      list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1),
      data = df,
      family = "degpd",
      removeData = TRUE
    )
  )

  expect_error(rqresid(fit), "requires the original data")
})

test_that("probit inverse links use the CDF", {
  expect_equal(egpd:::.detect_invlink("probitpi")(0), 0.5)
})

test_that("rbegpd supports the xi = 0 limit", {
  set.seed(14)
  Y <- rbegpd(200, kappa = 2, sigma = 1, xi = 0, thL = 5, thU = 5, thw = 0.2)

  expect_equal(dim(Y), c(200L, 2L))
  expect_true(all(is.finite(Y)))
  expect_true(all(Y >= 0))
})
