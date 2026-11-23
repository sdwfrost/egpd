test_that("a source predictor keeps its values in and alongside backfitted terms", {
  skip_if_not_installed("gamlss")
  withr::local_package("gamlss")
  set.seed(87)
  d <- data.frame(t = rep(seq(0, 1, length.out = 40), 4), x = runif(160))
  d$source <- d$t
  d$.src <- seq(0, 1, length.out = nrow(d))
  d$y <- sin(6*d$t) + 2*d$t*d$x + abs(d$.src-.5) + rnorm(nrow(d), sd = .2)
  fit <- function(form) gamlss::gamlss(form, data = d, family = gamlss.dist::NO(),
                                        sigma.start = .2, sigma.fix = TRUE, trace = FALSE)
  nd <- data.frame(t = c(.2, .75, .2, .05), x = c(.95, .4, .1, .6), .src = .25)
  nd$source <- nd$t
  nsim <- 200L
  units <- list(source = 2)
  ## The first control spells out the fitted knot and the sum of both offsets
  ## since predict.gamlss() by itself does not preserve them.
  for (forms in list(list(y ~ pb(x) * source + splines::bs(.src, df=2, degree=1) +
                            offset(source) + offset(.src),
                          y ~ pb(x) * t + splines::bs(.src, knots=.5, degree=1,
                                                      Boundary.knots=c(0,1)) + offset(t + .src)),
                     list(y ~ pb(source * units$source) * x, y ~ pb(t * units$source) * x))) {
    f <- fit(forms[[1]])
    control <- fit(forms[[2]])
    ## Include the upstream prediction in the seeded reference: pb() itself
    ## uses the RNG. The response simulation must use the fixed fitted sigma.
    set.seed(93)
    m <- predict(control, newdata = nd[c("x", "t", ".src")], data = d, type = "link")
    bounds <- vapply(m, function(mu)
      quantile(rnorm(nsim, mu, .2), c(.025, .975), names = FALSE), numeric(2))
    set.seed(93)
    r <- predictive_coverage(f, numeric(nrow(nd)), nd[c("x", "source", ".src")],
                             method = "plug-in", nsim = nsim)
    expect_equal(r$L, unname(bounds[1, ]))
    expect_equal(r$U, unname(bounds[2, ]))
  }
})

test_that("an unused source cannot exchange training and prediction rows", {
  skip_if_not_installed("gamlss")
  withr::local_package("gamlss")
  d <- expand.grid(t = seq(0, 1, length.out = 36), site = 0:1, k = 1:3)
  set.seed(718)
  d$y <- rpois(nrow(d), exp(1 + .85*sin(5*d$t) + .9*d$site))
  d$source <- "data"
  d$source[c(1, 21, 33, 72)] <- "other"
  f <- gamlss::gamlss(y ~ pb(t) + site, data = d, family = gamlss.dist::PO(), trace = FALSE)
  nd <- d[c(45, 28, 45, 59), , drop = FALSE]
  cover <- function(nd) {
    set.seed(109)
    predictive_coverage(f, rep(12, nrow(nd)), nd, data_fit = d, method = "plug-in", nsim = 500)
  }
  ## Exactly as many training rows are labelled "other" as there are
  ## new rows labelled "data": even the wrong result has the right length.
  clean <- cover(nd[c("t", "site")])
  expect_identical(cover(nd[c("source", "t", "site")]), clean)
})

test_that("ridge smoother column selectors keep a predictor named source", {
  skip_if_not_installed("gamlss")
  withr::local_package("gamlss")
  set.seed(193)
  d <- data.frame(a = runif(160), b = rnorm(160))
  d$y <- 3 + d$a - d$b + rnorm(nrow(d), sd = .2)
  ds <- d; names(ds)[1] <- "source"
  selected <- c("source", "b")
  reg <- function(f, data) gamlss::gamlss(f, data = data, family = gamlss.dist::NO(),
                                          sigma.start = .2, sigma.fix = TRUE, trace = FALSE)
  a <- reg(y ~ ri(x.vars = c("a", "b")), d)
  b <- reg(y ~ ri(x.vars = selected), ds)
  at <- c(5, 23, 55)
  set.seed(881)
  mu <- predict(a, newdata = d[at, 1:2], data = d, type = "link")
  ref <- vapply(mu, function(m) quantile(rnorm(256, m, .2), c(.025,.975),
                                          names = FALSE), numeric(2))
  set.seed(881)
  r <- predictive_coverage(b, ds$y[at], ds[at, 1:2], data_fit = ds,
                           method = "plug-in", nsim = 256)
  expect_equal(unname(rbind(r$L, r$U)), unname(ref))
})
