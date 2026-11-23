test_that("all native offsets enter the likelihood and new-data prediction", {
  d <- expand.grid(a = c(1, 3), b = c(1, 4), i = 1:30)
  v <- qexp(((1:30) - .5) / 30)
  d$y <- 2.2 * d$a * d$b * v[d$i] / mean(v)
  exposure <- d$a * d$b
  nd <- data.frame(a = 2, b = c(1, 4))
  for (fam in c("egpd", "degpd")) {
    if (fam == "degpd") d$y <- floor(d$y)
    truth <- if (fam == "egpd") mean(d$y / exposure) else {
      score <- function(sc) {
        r <- 1 / (sc * exposure)
        sum(r * (d$y - 1 / expm1(r)))
      }
      uniroot(score, c(.5, 6), tol = 1e-10)$root
    }
    for (f in list(y ~ offset(log(a)) + offset(log(b)),
                   y ~ offset(log(b)) + offset(log(a)))) {
      args <- list(formula = list(f, ~1, ~1), data = d, family = fam,
                   fixed = list(xi = 0, kappa = 1), restarts = FALSE, compact = TRUE)
      args[[paste0(fam, ".args")]] <- list(link = "identity")
      fit <- do.call(egpd, args)
      expect_equal(exp(unname(coef(fit)[1])), truth, tolerance = 5e-6)
      expect_equal(predict(fit, nd, type = "response")$scale, nd$a * nd$b * truth,
                   tolerance = 5e-6)
      rate <- 1 / (truth * exposure)
      ll <- if (fam == "egpd") log(rate) - rate*d$y else log(-expm1(-rate)) - rate*d$y
      expect_equal(as.numeric(logLik(fit)), sum(ll), tolerance = 1e-9)
    }
  }
})

test_that("secondary offset-only blocks retain row masks, names and environments", {
  d <- data.frame(y = 1:32, a = rep(c(0, .1), 16), check.names = FALSE)
  d[["extra offset"]] <- rep(c(0, .04), each = 16)
  d[[".egpd_offset"]] <- 100
  d[["extra offset"]][7] <- NA
  fo <- local({ mult <- 2; ~0 + offset(a) + offset(mult * `extra offset`) })
  expect_warning(f <- egpd(list(y ~ 1, fo, ~1), d, family = "egpd",
                           egpd.args = list(link = "identity"), fixed = list(kappa = 1),
                           compact = TRUE, restarts = FALSE), "Omitting 1")
  expect_equal(f$nobs, 31)
  nd <- d[c(3, 22), c("a", "extra offset")]
  expect_equal(predict(f, nd, type = "response")$shape, nd$a + 2 * nd[["extra offset"]])
})

test_that("hurdle prediction leaves estimable factor rows in their original positions", {
  d <- data.frame(g = factor(rep(c("empty", "large", "small"), c(40, 20, 20))),
                  y = c(rep(0, 40), rep(0, 4), rep(c(10, 14, 18, 22), 4),
                        rep(0, 4), rep(c(1, 2, 2, 3), 4)))
  f <- fit_hurdle(list(y ~ g, ~1, ~1), d, zero.formula = ~1, family = "degpd",
                  degpd.args = list(link = "identity"), fixed = list(xi = 0, kappa = 1))
  means <- tapply(d$y[d$y > 0], d$g[d$y > 0], mean)
  scales <- -1 / log1p(-1 / means)
  prob <- .91
  qs <- ceiling(scales * (log(mean(d$y > 0)) - log1p(-prob)))
  nd <- d[c(41, 1, 61, 43), "g", drop = FALSE]
  expect_warning(par <- predict(f, nd), "empty.*NA")
  expect_equal(par$scale, unname(as.numeric(scales[nd$g])), tolerance = 1e-5)
  expect_equal(par$pzero, rep(mean(d$y == 0), 4), tolerance = 1e-7)
  expect_warning(q <- predict(f, prob = prob), "empty.*NA")
  expect_equal(q[[1]], as.numeric(qs[d$g]))
  expect_silent(known <- predict(f, nd[c(1, 3, 4), , drop = FALSE], prob = prob))
  expect_equal(known[[1]], as.numeric(qs[nd$g[c(1, 3, 4)]]))
  expect_warning(no_positive <- predict(f, nd[2, , drop = FALSE]), "empty.*NA")
  expect_true(is.na(no_positive$scale))
})
