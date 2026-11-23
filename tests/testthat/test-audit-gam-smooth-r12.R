test_that("a predictor's fx smooth survives a penalty on another predictor", {
  set.seed(495)
  d <- data.frame(v = runif(560), w = runif(560))
  xi <- .18
  k <- exp(.2 + 1.6 * d$w)
  sig <- exp(.7 + .4 * cos(2 * pi * d$v))
  d$y <- sig * expm1(-xi * log1p(-runif(nrow(d))^(1/k))) / xi
  fit <- function(first, ...) egpd(list(first, ~ 1, ~ s(w, k = 5, fx = TRUE)),
                                   d, fixed = list(xi = xi), restarts = FALSE, ...)
  f <- fit(y ~ s(v, k = 6), sp = .7)
  nd <- data.frame(w = c(.09, .31, .72, .95), v = c(.23, .67, .91, .11))
  G <- mgcv::gam(y ~ s(w, k = 5, fx = TRUE), data = d, fit = FALSE)
  X <- cbind(1, mgcv::PredictMat(G$smooth[[1]], nd))
  eta <- drop(X %*% coef(f$logkappa))
  se <- sqrt(rowSums((X %*% f$logkappa$Vc) * X))
  expect_identical(f$convergence, 0L)
  expect_lt(max(abs(predict(f, nd, type = "lpmatrix")$logkappa - X)), 1e-10)
  expect_equal(predict(f, nd, type = "response")$kappa, exp(eta))
  expect_equal(predict(f, nd, se.fit = TRUE)$se.fit$logkappa, se)
  expect_equal(unname(confint(f, "kappa", newdata = nd)),
               exp(cbind(eta - qnorm(.975)*se, eta + qnorm(.975)*se)))

  rows <- c(2L, 19L, 57L)
  expect_equal(unname(simulate(f, 5, seed = 42, newdata = d[rows, ])$kappa),
               unname(simulate(f, 5, seed = 42)$kappa[rows, ]))
  st <- summary(f)
  expect_identical(rownames(st[[1]]$logkappa), "(Intercept)")
  expect_identical(rownames(st[[2]]$logkappa), "s(w)")
  drawn <- character()
  local_mocked_bindings(plot.gam = function(x, select, ...)
    drawn <<- c(drawn, x$smooth[[select]]$label), .package = "mgcv")
  plot(f, onepage = FALSE, ask = FALSE)
  expect_identical(drawn, c("s(v)", "s(w)"))

  ## No penalties anywhere: still a smooth, but no REML search or multipliers.
  no_pen <- fit(y ~ v)
  expect_length(no_pen$sp, 0L)
  expect_null(no_pen$Sdata)
  expect_equal(predict(no_pen, nd)$logkappa, drop(X %*% coef(no_pen$logkappa)))
  expect_identical(rownames(summary(no_pen)[[2]]$logkappa), "s(w)")
  ## A term-fixed multiplier is not a free sp; a top-level zero is still exact.
  term_fixed <- fit(y ~ s(v, k = 6, sp = .7))
  expect_length(term_fixed$sp, 0L)
  expect_equal(unname(term_fixed$full.sp), .7)
  expect_equal(coef(term_fixed), coef(f))
  zero <- fit(y ~ s(v, k = 6), sp = 0)
  expect_equal(unname(zero$full.sp), 0)
  expect_equal(predict(zero, nd)$logkappa, drop(X %*% coef(zero$logkappa)))
})

test_that("a one-column random effect has a scaled one-degree-of-freedom test", {
  set.seed(31)
  d <- data.frame(z = rnorm(180), x = runif(180))
  d$y <- 1 + .21 * d$z - .5 * d$x + rnorm(nrow(d))
  gauss <- list(
    d0 = function(b, ld) -sum(dnorm(ld$y, ld$X[[1L]] %*% b, 1, log = TRUE)),
    d120 = function(b, ld) cbind(drop(ld$X[[1L]] %*% b - ld$y), 1))
  form <- y ~ x + s(z, bs = "re")
  f <- egpd(list(mean = form), d, family = "custom", custom.fns = gauss, sp = .7)
  m <- mgcv::gam(form, data = d, method = "REML", scale = 1, sp = .7)
  j <- m$smooth[[1L]]$first.para
  chisq <- unname(coef(m)[j]^2 / m$Vp[j, j])
  tab <- summary(f)[[2]]$mean
  expect_equal(unname(coef(f)), unname(coef(m)))
  expect_equal(tab["s(z)", "Chi.sq"], chisq)
  ## Its null sampling variance (not its posterior variance) calibrates the test.
  expect_equal(tab["s(z)", "p-value"],
               pchisq(unname(coef(m)[j]^2 / m$Ve[j, j]), 1, lower.tail = FALSE))
})
