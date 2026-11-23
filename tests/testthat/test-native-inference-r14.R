test_that("a random intercept tests zero variance, not an ordinary smooth null", {
  m <- 32L
  g <- factor(rep(1:2, each = m))
  v <- qexp((seq_len(m) - .5) / m)
  for (share in c(.35, .37, .42)) {
    d <- data.frame(y = c(2 * share * v, 2 * (1 - share) * v) / mean(v), g)
    f <- egpd(list(y ~ s(g, bs = "re"), ~1, ~1), d, fixed = list(xi = 0, kappa = 1),
              restarts = FALSE)
    ref <- mgcv::gam(y ~ s(g, bs = "re"), data = d, family = stats::Gamma("log"),
                     scale = 1, method = "REML")
    p <- summary(f)[[2]]$logscale["s(g)", "p-value"]
    ## Under a common exponential scale, the two group totals are independent
    ## Gammas, so their first share has the exact null law Beta(m,m).
    expect_lt(abs(p - 2 * pbeta(share, m, m)), .004)
    expect_lt(abs(p - summary(ref)$s.table[1, "p-value"]), .003)
  }
})

test_that("boundary tests integrate other random effects across predictor blocks", {
  set.seed(881)
  n <- 180L
  d <- data.frame(g = factor(sample(1:5, n, TRUE)), h = factor(sample(1:4, n, TRUE)),
                  x = runif(n, -1, 1))
  z <- d$z <- sample(c(-1, 1), n, TRUE, prob = c(.2, .8))
  d$y <- .4 * d$x + .12 * as.integer(d$g) + z * (.2 - .15 * as.integer(d$h)) + rnorm(n)
  ## Two predictors in the native likelihood, but a single equivalent
  ## Gaussian regression for mgcv. The cross-block information is nonzero.
  mu <- function(b, ld) drop(ld$X[[1L]] %*% b[ld$idpars == 1L] +
                               z * (ld$X[[2L]] %*% b[ld$idpars == 2L]))
  gaussian2 <- list(
    d0 = function(b, ld) -sum(dnorm(ld$y, mu(b, ld), 1, log = TRUE)),
    d120 = function(b, ld) {
      e <- mu(b, ld) - ld$y[, 1L]
      cbind(e, e * z, 1, z, z^2)
    })
  sp <- c(.7, 5, 8)
  f <- egpd(list(a = y ~ s(x, k = 5) + s(g, bs = "re"), b = ~1 + s(h, bs = "re")),
            d, family = "custom", custom.fns = gaussian2, sp = sp, gamma = 1.5)
  ref <- mgcv::gam(y ~ z + s(x, k = 5) + s(g, bs = "re") + s(h, by = z, bs = "re"),
                   data = d, scale = 1, method = "REML", sp = sp)
  tab <- summary(f)[[2]]
  ref.tab <- summary(ref)$s.table
  for (nm in c("Chi.sq", "p-value"))
    expect_equal(c(tab$a["s(g)", nm], tab$b["s(h)", nm]),
                 unname(ref.tab[c("s(g)", "s(h):z"), nm]), tolerance = 1e-6)
})

test_that("a confounded random effect has no p-value; ordinary smooths keep their test", {
  set.seed(53014)
  x <- runif(180, -1, 1)
  d <- data.frame(y = rexp(length(x), exp(-.2 * cos(pi * x))), x,
                  g = factor(rep(1:3, length.out = length(x))))
  fixed <- list(xi = 0, kappa = 1)
  form <- y ~ g + s(g, bs = "re") + s(x, k = 5)
  f <- egpd(list(form, ~1, ~1), d, fixed = fixed, restarts = FALSE, sp = c(200, .7))
  ## The unpenalised group intercepts already span every random intercept.
  expect_warning(tab <- summary(f)[[2]]$logscale, "s\\(g\\).*no estimable")
  expect_true(all(is.na(tab["s(g)", c("Chi.sq", "p-value")])))
  ordinary <- egpd(list(y ~ g + s(x, k = 5), ~1, ~1), d, fixed = fixed,
                   restarts = FALSE, sp = .7)
  expect_equal(tab["s(x)", ], summary(ordinary)[[2]]$logscale["s(x)", ], tolerance = 1e-6)
})

test_that("magnitude diagnoses cluster fixed effects without changing its scalar factor", {
  set.seed(51463)
  G <- 9L; m <- 28L
  d <- data.frame(g = factor(rep(seq_len(G), each = m)), period = rep(seq_len(m), G),
                  x = rnorm(G * m))
  d$y <- rexp(nrow(d), exp(-.7 - .35 * d$x - rep(seq(-.4, .4, length.out = G), each = m)))
  X <- model.matrix(~x + g, d)
  ## Profile out all group intercepts in the exponential model.
  rho <- function(b) {
    z <- d$y * exp(-b * d$x)
    z / ave(z, d$g, FUN = mean)
  }
  b <- uniroot(function(b) sum(d$x * (rho(b) - 1)), c(-2, 2), tol = 1e-12)$root
  w <- rho(b)
  Hi <- solve(crossprod(X, w * X))
  fit <- function(id, method = "magnitude", force = FALSE)
    egpd(list(y ~ x + g, ~1, ~1), d, fixed = list(xi = 0, kappa = 1), restarts = FALSE,
         sandwich.args = list(id = id, method = method, force = force))
  for (id in c("g", "period")) {
    if (id == "g")
      expect_warning(f <- fit(id), "cluster score covariance is rank deficient \\(rank 1/10\\)")
    else
      expect_silent(f <- fit(id))
    J <- crossprod(rowsum((w - 1) * X, d[[id]]))
    k <- ncol(X) / sum(diag(Hi %*% J))
    expect_equal(f$likdata$k, k, tolerance = 1e-5)
    expect_equal(unname(f$Vp[2, 2]), unname(Hi["x", "x"] / k), tolerance = 1e-5)
  }
  expect_error(fit("g", "curvature"), "cluster score covariance.*rank 1/10")
  ## force is still available for the analytic backend, but does not restore
  ## the nine absent score directions.
  expect_warning(suppressMessages(fc <- fit("g", "curvature", TRUE)),
                 "forced pseudo-inverse inference is unsupported")
  expect_equal(unname(fc$hessian.rank), c(1L, 10L))
  J <- crossprod(rowsum((w - 1) * X, d$g))
  expect_equal(unname(fc$Vp[2, 2]), unname((Hi %*% J %*% Hi)["x", "x"]), tolerance = 1e-5)
})
