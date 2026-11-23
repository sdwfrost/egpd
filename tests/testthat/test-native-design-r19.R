test_that("an incomplete prediction row does not strip offsets from complete rows", {
  d <- expand.grid(x = 0:1, e = c(.5,2), i = 1:24)
  u <- qexp(((1:24)-.5)/24)
  d$y <- exp(.4 + .65*d$x)*d$e*u[d$i]/mean(u)
  f <- egpd(list(y ~ x + offset(log(e)), ~1, ~1), d, family = "egpd",
            egpd.args = list(link = "identity"), fixed = list(xi = 0, kappa = 1),
            restarts = FALSE)
  nd <- data.frame(x = c(1,.4,0), e = c(3,.8,.5))
  ii <- c(1,3)
  mu <- exp(.4 + .65*nd$x[ii])*nd$e[ii]
  ci <- confint(f, "scale", newdata = nd)[ii,,drop = FALSE]
  draws <- simulate(f, nsim = 12L, seed = 193, type = "link", newdata = nd)$scale[ii,]
  for (v in c("x","e")) {
    missing <- nd
    missing[[v]][2] <- NA
    scale <- predict(f, missing, type = "response")$scale
    expect_true(is.na(scale[2]))
    expect_equal(scale[ii], mu, tolerance = 1e-5)
    X <- predict(f, missing, type = "lpmatrix")
    expect_equal(attr(X, "offsets")[[1L]], c(log(3), NA, log(.5)))
    expect_equal(confint(f, "scale", newdata = missing)[ii,,drop = FALSE], ci)
    expect_equal(simulate(f, nsim = 12L, seed = 193, type = "link", newdata = missing)$scale[ii,],
                 draws)
  }
})

test_that("compaction groups already evaluated terms and offsets", {
  d <- data.frame(x = rep(c(-1,0,2), c(70,30,10)))
  d$e <- .5*(d$x-mean(d$x))^2
  i <- ave(seq_len(nrow(d)), d$x, FUN = seq_along)
  n <- ave(i, d$x, FUN = length)
  u <- qexp((i-.5)/n)
  d$y <- 2.3*exp(d$e)*u/ave(u,d$x,FUN = mean)
  f <- function(form, compact, maxspline = Inf) egpd(list(form,~1,~1), d, family = "egpd",
    egpd.args = list(link = "identity"), fixed = list(xi = 0, kappa = 1), restarts = FALSE,
    compact = compact, maxspline = maxspline)
  for (ms in c(Inf, 50)) {
    set.seed(19)
    a <- f(y ~ offset(I(.5*(x-mean(x))^2)), TRUE, ms)
    expect_equal(exp(unname(coef(a)[1])), 2.3, tolerance = 1e-5)
    expect_equal(as.numeric(logLik(a)), sum(dexp(d$y, rate = exp(-d$e)/2.3, log = TRUE)),
                 tolerance = 1e-9)
    expect_equal(as.numeric(fitted(a)[[1]]), as.numeric(predict(a, d, type = "link")[[1]]))
    expect_equal(length(unique(a$compactid)), 3L)
  }
  ## Identical raw x need not imply identical evaluated design rows.
  d$y <- exp(.3 + .002*cumsum(d$x + 2))*rep(c(.85,1.15), length.out = nrow(d))
  fo <- y ~ I(cumsum(x + 2))
  a <- f(fo, TRUE); b <- f(fo, FALSE)
  expect_equal(coef(a), coef(b), tolerance = 1e-6)
  expect_equal(a$logLik, b$logLik, tolerance = 1e-9)
})

test_that("native fits select rows using evaluated predictor terms across blocks", {
  d <- expand.grid(g = c("a","b","other"), i = 1:36)
  d$valid <- factor(d$g, levels = c("a","b"))
  d$y <- floor(ifelse(d$g == "b",14,2)*qexp((d$i-.5)/36))
  ok <- !is.na(d$valid)
  mu <- tapply(d$y[ok],d$valid[ok],mean)
  truth <- sum(dgeom(d$y[ok],1/(1+mu[d$valid[ok]]),log = TRUE))
  for (comp in c(FALSE, TRUE)) {
    expect_warning(f <- egpd(list(y ~ factor(g, levels = c("a","b")), ~1, ~1),d,
       family = "degpd", degpd.args = list(link = "identity"), fixed = list(xi = 0, kappa = 1),
       restarts = FALSE, compact = comp), "Omitting 36")
    expect_equal(nobs(f), sum(ok))
    expect_equal(as.numeric(logLik(f)), truth, tolerance = 1e-9)
    expect_equal(predict(f, data.frame(g = names(mu)), type = "response")$scale,
                 as.numeric(1/log1p(1/mu)), tolerance = 1e-5)
  }
})
