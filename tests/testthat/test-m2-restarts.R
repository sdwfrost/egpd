## Independent log-likelihood, not the package d0 or the restart screen.
## G(H) = p H^a + (1-p) H^b; counts are floors, not point densities.
m2_ll <- function(y, eta, family, xi = NULL) {
  s <- exp(eta[[1]])
  if (is.null(xi)) xi <- if (family == "egpd") eta[[2]] else exp(eta[[2]])
  a <- exp(eta[[3]])
  b <- a + exp(eta[[4]])
  p <- plogis(eta[[5]])
  A <- function(x) if (all(xi == 0)) x / s else log1p(xi * (x / s)) / xi
  lh <- function(x) {
    z <- A(x)
    ifelse(z > log(2), log1p(-exp(-z)), log(-expm1(-z)))
  }
  if (family == "egpd") {
    l <- lh(y)
    return(sum(log(p * a * exp((a - 1) * l) + (1 - p) * b * exp((b - 1) * l)) -
                 log(s) - (1 + xi) * A(y)))
  }
  cdf <- function(x) {
    l <- lh(x)
    p * exp(a * l) + (1 - p) * exp(b * l)
  }
  mass <- cdf(y + 1) - cdf(y)
  if (family == "zidegpd") {
    pi <- plogis(eta[[6]])
    mass <- pi * (y == 0) + (1 - pi) * mass
  }
  sum(log(mass))
}

m2_fit <- function(y, family, ...) {
  fl <- c(list(y ~ 1), rep(list(~1), 4L + (family == "zidegpd")))
  egpd(fl, data.frame(y = y), family = family, egpd.args = list(m = 2),
       degpd.args = list(m = 2), zidegpd.args = list(m = 2), ...)
}

## The intentionally stranded single starts and exponential-tail limit can
## have a singular Hessian. Do not muffle nonconvergence or numerical warnings.
m2_allow_singular <- function(expr) withCallingHandlers(expr,
  egpd_singular_hessian = function(w) invokeRestart("muffleWarning"))

m2_profile_data <- function() {
  set.seed(11)
  list(egpd = regpd(400, sigma = 2, xi = 0.15, kappa = 1.5),
       degpd = rdiscegpd(400, sigma = 3, xi = 0.2, kappa = 1.5, type = 1))
}

test_that("model 2 escapes both Round-7 better-fit profile misses", {
  ys <- m2_profile_data()
  ## Feasible points from independent multistart optimization, *not* refits
  ## from the package coefficients. The discrete point also improves on the
  ## published challenger (-1034.323431) by a further 0.21.
  refs <- list(egpd = c(0.6460399, 0.00407764, 0.4830601, 8.269654, 3.846758),
               degpd = c(1.242989, -2.864637, 0.3210772, 5.646950, 3.822108))
  for (fam in names(refs)) {
    y <- ys[[fam]]
    old <- if (fam == "egpd") -802.259371732 else -1034.71237459
    ## The discrete single start is nearly a one-power fit, with singular Hessian.
    one <- m2_allow_singular(m2_fit(y, fam, restarts = FALSE))
    expect_lt(abs(one$logLik - old), 1e-5)
    f <- m2_fit(y, fam)
    expect_equal(f$convergence, 0L)
    ll <- m2_ll(y, as.list(coef(f)), fam)
    expect_lt(abs(f$logLik - ll), 1e-7)
    expect_gte(ll, m2_ll(y, as.list(refs[[fam]]), fam) - 1e-5)
  }
  explicit <- m2_fit(ys$egpd, "egpd", inits = c(0.31, 0.24, -0.05, 0.2, -2.77))
  expect_lt(abs(explicit$logLik + 802.259371732), 1e-5)
})

test_that("a held-out floor mixture is not mistaken for one power", {
  ## Floor independent continuous component draws. A 30-start direct search
  ## found this point; the old five starts all missed it by over 2 logLik units.
  set.seed(137)
  k <- ifelse(runif(320) < 0.9, 1.3, 25)
  u <- runif(320)^(1 / k)
  y <- floor(3 * expm1(-0.18 * log1p(-u)) / 0.18)
  ref <- c(0.8526931, -1.347632, 0.2781144, 2.961276, 1.461595)
  f <- m2_fit(y, "degpd")
  ll <- m2_ll(y, as.list(coef(f)), "degpd")
  expect_lt(abs(f$logLik - ll), 1e-7)
  expect_gte(ll, m2_ll(y, as.list(ref), "degpd") - 1e-5)
})

test_that("an interior ZIDEGPD model-2 warm start can still miss a better fit", {
  set.seed(1101)
  n <- 280
  k <- sample(c(1, 8, 90), n, prob = c(0.65, 0.3, 0.05), replace = TRUE)
  y <- floor(2.5 * expm1(-0.12 * log1p(-runif(n)^(1 / k))) / 0.12)
  y[runif(n) < 0.18] <- 0
  warm <- m2_allow_singular(m2_fit(y, "zidegpd", restarts = FALSE))
  expect_equal(warm$convergence, 0L)
  expect_lt(max(abs(coef(warm)[-1])), 3)
  ## Independent search: the better basin is on the exponential-tail limit.
  ref <- c(1.438154, -22.36789, 0.2460968, 6.684536, 3.986946, -1.264761)
  fit <- m2_allow_singular(m2_fit(y, "zidegpd"))
  expect_equal(fit$convergence, 0L)
  ll <- m2_ll(y, as.list(coef(fit)), "zidegpd")
  expect_lt(abs(fit$logLik - ll), 1e-7)
  expect_gt(ll - warm$logLik, 1.6)
  expect_gte(ll, m2_ll(y, as.list(ref), "zidegpd") - 1e-3)
})

test_that("model-2 restarts fit the same subsample with fixed values and offsets", {
  set.seed(777)
  x <- runif(240)
  k <- ifelse(runif(240) < 0.8, 1, 12)
  y <- floor(exp(1 + x) * expm1(-0.15 * log1p(-runif(240)^(1 / k))) / 0.15)
  d <- data.frame(y = y, x = x, off = log(3) + x / 5)
  fl <- list(y ~ x + offset(off), ~1, ~1, ~1, ~1)
  run <- function(restarts) suppressMessages(egpd(
    fl, d, family = "degpd", degpd.args = list(m = 2, xi.max = 0.25),
    fixed = list(xi = 0.15, p = 0.8), maxdata = 180, maxspline = 120,
    removeData = TRUE, restarts = restarts))
  seen <- list()
  setup <- egpd:::.setup.data
  local_mocked_bindings(.setup.data = function(...) {
    result <- setup(...)
    seen[[length(seen) + 1L]] <<- result$lik.data$y
    result
  }, .package = "egpd")
  set.seed(314)
  one <- run(FALSE)
  state <- .Random.seed
  set.seed(314)
  fit <- run(TRUE)
  expect_identical(.Random.seed, state)
  expect_identical(unique(seen), list(one$likdata$y))
  eta <- lapply(fit[seq_len(fit$ngam)], `[[`, "fitted")
  xi <- -0.25 * expm1(-exp(eta[[2]]) / 0.25)
  expect_lt(max(abs(xi - 0.15)), 1e-12)
  expect_lt(max(abs(plogis(eta[[5]]) - 0.8)), 1e-12)
  expect_lt(abs(fit$logLik - m2_ll(fit$likdata$y[, 1], eta, "degpd", xi)), 1e-7)
})

test_that("a good model-2 screen score is not a converged fit", {
  yd <- m2_profile_data()$degpd
  ws <- list()
  f <- withCallingHandlers(m2_fit(yd, "degpd", control = list(inner = list(itlim = 1))),
    warning = function(w) {
      ws[[length(ws) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  expect_true(any(vapply(ws, inherits, logical(1), "egpd_nonconvergence")))
  expect_true(any(grepl("no restart converged", vapply(ws, conditionMessage, ""))))
  expect_false(f$convergence == 0L)
  expect_true(is.na(f$logLik))
})
