test_that("R2-BAMLSS: .parse_links matches named overrides and validates links", {
  f1 <- egpd_bamlss(1, kappa = "identity")
  expect_identical(f1$links, c(sigma = "log", xi = "log", kappa = "identity"))

  f2 <- egpd_bamlss(2, prob = "probit", delta = "identity")
  expect_identical(
    f2$links,
    c(sigma = "log", xi = "log", kappa = "log", delta = "identity", prob = "probit")
  )

  f_pos <- egpd_bamlss(1, c("identity", "log", "identity"))
  expect_identical(f_pos$links, c(sigma = "identity", xi = "log", kappa = "identity"))

  f_named_vec <- zidegpd_bamlss(1, c(pi = "probit", xi = "identity"))
  expect_identical(
    f_named_vec$links,
    c(sigma = "log", xi = "identity", kappa = "log", pi = "probit")
  )

  expect_error(egpd_bamlss(1, unknown_par = "identity"), "Unknown link parameter")
  expect_error(degpd_bamlss(3, kappa = "identity"), "Unknown link parameter")
  expect_error(egpd_bamlss(1, c("log", "log")), "length 3")
})

test_that("R2-BAMLSS: models 1..6, ZI initializers, and stabilized PIG/ZIPIG work", {
  init_ok <- TRUE
  for (m in 1:6) {
    for (ctor in list(egpd_bamlss, degpd_bamlss, ziegpd_bamlss, zidegpd_bamlss)) {
      fam <- ctor(m)
      if (!inherits(fam, "family.bamlss")) init_ok <- FALSE
      for (y_edge in list(c(0, 0, 0), c(0, 0, 5), c(0, NA_real_, 4), 3)) {
        inits <- lapply(fam$initialize, function(fn) fn(y_edge))
        for (nm in names(inits)) {
          v <- inits[[nm]]
          if (length(v) != length(y_edge) || any(!is.finite(v)) || any(v <= 0)) {
            init_ok <- FALSE
          }
        }
      }
    }
  }
  expect_true(init_ok)

  ## Poisson-limit stabilization in pig_bamlss / zipig_bamlss (mu * sigma < 1e-5)
  fp <- pig_bamlss()
  fz <- zipig_bamlss()
  lp <- fp$d(0:5, list(mu = rep(2, 6), sigma = rep(1e-9, 6)), log = TRUE)
  expect_equal(lp, dpois(0:5, lambda = 2, log = TRUE), tolerance = 1e-6)

  lz <- fz$d(0:5, list(mu = rep(2, 6), sigma = rep(1e-9, 6), pi = rep(0.25, 6)), log = TRUE)
  want_z <- c(
    log(0.25 + 0.75 * dpois(0, 2)),
    log(0.75) + dpois(1:5, 2, log = TRUE)
  )
  expect_equal(lz, want_z, tolerance = 1e-6)
})

test_that("R2-BAMLSS: end-to-end bamlss() fits succeed across univariate families", {
  skip_on_cran()
  skip_if_not_installed("bamlss")

  was_attached <- "package:bamlss" %in% search()
  if (!was_attached) {
    suppressPackageStartupMessages(library(bamlss))
    on.exit(detach("package:bamlss", unload = FALSE, character.only = TRUE), add = TRUE)
  }

  set.seed(2026)
  n <- 80
  df_c <- data.frame(y = regpd(n, sigma = 2, xi = 0.2, kappa = 1.4, type = 1))
  b_egpd <- suppressWarnings(
    bamlss::bamlss(y ~ 1, data = df_c, family = egpd_bamlss(1),
                   optimizer = bamlss::opt_bfit, sampler = FALSE, verbose = FALSE)
  )
  expect_s3_class(b_egpd, "bamlss")

  df_d <- data.frame(y = rdiscegpd(n, sigma = 3, xi = 0.15, kappa = 1.3, type = 1))
  b_degpd <- suppressWarnings(
    bamlss::bamlss(y ~ 1, data = df_d, family = degpd_bamlss(1),
                   optimizer = bamlss::opt_bfit, sampler = FALSE, verbose = FALSE)
  )
  expect_s3_class(b_degpd, "bamlss")

  df_zic <- data.frame(y = rziegpd(n, pi = 0.25, sigma = 2, xi = 0.2, kappa = 1.3, type = 1))
  b_ziegpd <- suppressWarnings(
    bamlss::bamlss(y ~ 1, data = df_zic, family = ziegpd_bamlss(1),
                   optimizer = bamlss::opt_bfit, sampler = FALSE, verbose = FALSE)
  )
  expect_s3_class(b_ziegpd, "bamlss")

  df_zid <- data.frame(y = rzidiscegpd(n, pi = 0.25, sigma = 3, xi = 0.15, kappa = 1.3, type = 1))
  b_zidegpd <- suppressWarnings(
    bamlss::bamlss(y ~ 1, data = df_zid, family = zidegpd_bamlss(1),
                   optimizer = bamlss::opt_bfit, sampler = FALSE, verbose = FALSE)
  )
  expect_s3_class(b_zidegpd, "bamlss")
})

test_that("R2-HURDLE: .ztdegpd1fns handles bounded.xi normalization and derivatives", {
  fns <- egpd:::.ztdegpd1fns
  a <- 0.5
  pars <- c(log(2.5), 0.7, log(1.4))

  ## Check exact PMF normalization over y >= 1 under bounded.xi
  ys <- 1:500
  probs <- vapply(ys, function(yk) {
    ld1 <- list(
      X = replicate(3, matrix(1, 1, 1), simplify = FALSE),
      y = matrix(yk, 1, 1),
      dupid = 0L, duplicate = 0L,
      offsets = replicate(3, numeric(0), simplify = FALSE),
      idpars = 1:3, censored = FALSE,
      xi.max = a, bounded.xi = TRUE
    )
    exp(-fns$d0(pars, ld1))
  }, numeric(1))
  xi_eff <- a * (-expm1(-exp(pars[2] - log(a))))
  f0_eff <- dDEGPD1(0, mu = exp(pars[1]), sigma = xi_eff, nu = exp(pars[3]))
  tail_mass <- pDEGPD1(500, mu = exp(pars[1]), sigma = xi_eff, nu = exp(pars[3]),
                       lower.tail = FALSE) / (1 - f0_eff)
  expect_equal(sum(probs) + tail_mass, 1, tolerance = 1e-12)

  ## Check gradient against central finite differences with covariates + bounded.xi
  y <- c(1, 2, 4, 7, 3)
  x <- seq(-1, 1, length.out = length(y))
  X <- list(cbind(1, x), cbind(1, x), matrix(1, length(y), 1))
  nbk <- vapply(X, ncol, integer(1))
  ld <- list(
    X = X, y = cbind(y), dupid = 0L, duplicate = 0L,
    offsets = replicate(3, numeric(0), simplify = FALSE),
    idpars = rep(1:3, nbk), censored = FALSE,
    xi.max = a, bounded.xi = TRUE
  )
  p_vec <- c(log(3), 0.2, 0.4, -0.3, log(1.3))
  M <- fns$d120(p_vec, ld)
  an_grad <- unname(unlist(lapply(1:3, function(j) colSums(M[, j] * X[[j]]))))
  h <- 1e-6
  fd_grad <- vapply(seq_along(p_vec), function(k) {
    pp <- pm <- p_vec
    pp[k] <- pp[k] + h
    pm[k] <- pm[k] - h
    (fns$d0(pp, ld) - fns$d0(pm, ld)) / (2 * h)
  }, numeric(1))
  expect_equal(an_grad, fd_grad, tolerance = 1e-7)
})

test_that("R2-HURDLE: fit_hurdle validates inputs, omits incomplete cases, and supports BIC", {
  set.seed(404)
  n <- 120
  x1 <- runif(n, -1, 1)
  x2 <- rnorm(n)
  y <- ifelse(runif(n) < 0.35, 0L,
              rdiscegpd(n, sigma = exp(1 + 0.3 * x1), xi = 0.15, kappa = 1.2) + 1L)
  df <- data.frame(y = y, x1 = x1, x2 = x2)

  expect_error(fit_hurdle(y ~ x1, data = df, family = "degpd", degpd.args = list(m = 2)),
               "m == 1")
  df_neg <- df; df_neg$y[1] <- -1
  expect_error(fit_hurdle(y ~ x1, data = df_neg, family = "degpd"),
               "non-negative integers")
  df_frac <- df; df_frac$y[1] <- 1.5
  expect_error(fit_hurdle(y ~ x1, data = df_frac, family = "degpd"),
               "non-negative integers")

  ## Inject NAs in response, positive formula predictor, and zero.formula predictor
  df_na <- df
  df_na$y[2] <- NA
  df_na$x1[5] <- NA
  df_na$x2[9] <- NA

  expect_warning(
    fit <- fit_hurdle(y ~ x1, data = df_na, family = "degpd", zero.formula = ~ x2),
    "omitted 3 observation"
  )
  expect_s3_class(fit, "hurdle.egpd")
  expect_identical(fit$n, as.integer(n - 3L))
  expect_true(is.integer(fit$n) && is.integer(fit$n_zero))
  expect_equal(fit$n_zero, sum(df_na$y[-c(2, 5, 9)] == 0))

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_identical(attr(ll, "nobs"), fit$n)
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
  expect_output(print(fit), "Hurdle model")
})

test_that("R2-SPLICED: Eq. 3 splice eliminates spike, supports negative xi, discrete quantiles, and factor newdata", {
  skip_if_not_installed("qgam")

  set.seed(303)
  n <- 240
  x <- runif(n, -1, 1)
  grp_chr <- sample(c("north", "south", "west"), n, replace = TRUE)
  grp_fac <- factor(sample(c("low", "mid", "high"), n, replace = TRUE))
  ## Bounded tail (negative xi) above bulk counts
  base_y <- rpois(n, lambda = exp(1.1 + 0.4 * x))
  exc_y  <- rdiscegpd(n, sigma = 6, xi = -0.35, kappa = 1, type = 1)
  y <- base_y + rbinom(n, 1, 0.25) * exc_y
  df <- data.frame(y = y, x = x, grp_chr = grp_chr, grp_fac = grp_fac,
                   stringsAsFactors = FALSE)

  fit <- suppressWarnings(
    fit_spliced_degpd(
      y ~ s(x, k = 5),
      data = df,
      bulk_probs = c(0.1, 0.25, 0.5, 0.75),
      threshold_prob = 0.75,
      tail_formula = list(lsigma = ~ grp_chr + grp_fac, xi = ~ 1, lkappa = ~ 1),
      keep_data = TRUE
    )
  )

  expect_s3_class(fit, "spliced_degpd")
  expect_true(is.finite(fit$tail_fit$xi))

  ## Single-row newdata with character and factor predictors
  nd1 <- df[1, , drop = FALSE]
  nd5 <- df[1:5, , drop = FALSE]
  p1 <- predict(fit, newdata = nd1, type = "parameter")
  p5 <- predict(fit, newdata = nd5, type = "parameter")
  expect_equal(as.numeric(p1[1, ]), as.numeric(p5[1, ]), tolerance = 1e-12)

  ## Verify Eq. 3 splice formula and absence of artificial spike at s_i
  s_i <- p1$threshold[1]
  alpha_t <- p1$alpha_tilde[1]
  at_grid <- 0:(s_i + 25)
  cdf1 <- as.numeric(predict(fit, newdata = nd1, type = "cdf", at = at_grid)[1, ])
  pmf1 <- as.numeric(predict(fit, newdata = nd1, type = "pmf", at = at_grid)[1, ])

  if (s_i >= 1) {
    expect_equal(cdf1[s_i], alpha_t, tolerance = 1e-12)
  }
  want_s <- alpha_t + (1 - alpha_t) *
    pdiscegpd(0, sigma = p1$scale[1], xi = p1$shape[1], kappa = p1$kappa[1], type = 1)
  expect_equal(cdf1[s_i + 1L], want_s, tolerance = 1e-12)
  if (s_i >= 1) {
    expect_equal(
      pmf1[s_i + 1L],
      (1 - alpha_t) * ddiscegpd(0, sigma = p1$scale[1], xi = p1$shape[1], kappa = p1$kappa[1], type = 1),
      tolerance = 1e-12
    )
  }

  ## Verify discrete integer quantiles match inf{k : F(k) >= p} across [0, 1]
  probs <- c(0, 0.02, 0.1, 0.3, 0.5, alpha_t, min(alpha_t + 0.05, 0.99), 0.95, 1)
  q1 <- as.numeric(predict(fit, newdata = nd1, type = "quantile", prob = probs)[1, ])
  expect_equal(q1[1], 0)
  fin_q <- which(is.finite(q1))
  expect_equal(q1[fin_q], floor(q1[fin_q]))
  for (j in fin_q) {
    kj <- as.integer(q1[j])
    pj <- probs[j]
    fk <- as.numeric(predict(fit, newdata = nd1, type = "cdf", at = kj)[1, 1])
    expect_gte(fk + 1e-14, pj)
    if (kj > 0 && pj > 0) {
      fkm1 <- as.numeric(predict(fit, newdata = nd1, type = "cdf", at = kj - 1L)[1, 1])
      expect_lt(fkm1, pj + 1e-14)
    }
  }

  ## Also check direct negative-xi tail recovery on synthetic bounded exceedances
  set.seed(304)
  df_b <- data.frame(
    x = runif(200),
    y = rdiscegpd(200, sigma = 5, xi = -0.4, kappa = 1, type = 1)
  )
  fit_b <- suppressWarnings(
    fit_spliced_degpd(y ~ x, data = df_b, bulk_probs = c(0.1, 0.25),
                      threshold_prob = 0.25)
  )
  expect_lt(fit_b$tail_fit$xi, 0)
})

test_that("R2-CARRIERS: Types 4/5 closed form at large delta, Type 2 small-kappa Taylor, and bounds", {
  u <- seq(0.05, 0.95, by = 0.1)

  ## Types 4 & 5 smooth convergence to u as delta -> Inf (e.g. delta = 2000)
  g4_2000 <- p.G(u, type = 4, delta = 2000)
  expect_true(all(g4_2000 < 1))
  expect_lt(max(abs(g4_2000 - u)), 6e-4)
  expect_equal(q.G(g4_2000, type = 4, delta = 2000), u, tolerance = 1e-12)

  g5_2000 <- p.G(u, type = 5, kappa = 1.6, delta = 2000)
  expect_lt(max(abs(g5_2000 - u^0.8)), 1e-3)
  expect_equal(q.G(g5_2000, type = 5, kappa = 1.6, delta = 2000), u, tolerance = 1e-12)
  expect_lt(max(abs(d.G(u, type = 5, kappa = 2, delta = 2000) - 1)), 1e-3)

  ## Type 2 exact kappa -> 0+ limit and first-order Taylor slopes
  expect_equal(p.G(u, type = 2, kappa = 0), u, tolerance = 1e-15)
  expect_equal(d.G(u, type = 2, kappa = 0), rep(1, length(u)), tolerance = 1e-15)
  expect_equal(d.G(u, type = 2, kappa = 0, log = TRUE), rep(0, length(u)), tolerance = 1e-15)
  expect_equal(q.G(u, type = 2, kappa = 0), u, tolerance = 1e-15)

  kap <- 1e-7
  w <- 1 - u
  expect_equal((p.G(u, type = 2, kappa = kap) - u) / kap,
               -w * (1 - w^2) / 6, tolerance = 1e-6)
  expect_equal((d.G(u, type = 2, kappa = kap) - 1) / kap,
               -(3 * w^2 - 1) / 6, tolerance = 1e-6)
  expect_equal((q.G(u, type = 2, kappa = kap) - u) / kap,
               w * (1 - w^2) / 6, tolerance = 1e-6)
  expect_equal(q.G(p.G(u, type = 2, kappa = 1e-5), type = 2, kappa = 1e-5),
               u, tolerance = 1e-14)

  ## Bounds outside [0, 1] for p.G and d.G across all types 1..6
  for (tp in 1:6) {
    expect_equal(
      p.G(c(-Inf, -2, 0, 1, 2, Inf), type = tp, prob = 0.4, kappa = 1.5, delta = 2),
      c(0, 0, 0, 1, 1, 1)
    )
    expect_equal(
      d.G(c(-Inf, -0.5, 1.5, Inf), type = tp, prob = 0.4, kappa = 1.5, delta = 2),
      c(0, 0, 0, 0)
    )
    expect_equal(
      d.G(c(-Inf, -0.5, 1.5, Inf), type = tp, prob = 0.4, kappa = 1.5, delta = 2, log = TRUE),
      rep(-Inf, 4)
    )
  }
})

test_that("R2-DIST: r.G type 6 vectorization, rzidiscegpd unifsamp, zero-length, NA/NaN, and log/tail args", {
  ## Vectorized r.G(..., type = 6)
  set.seed(42)
  rg6 <- r.G(6, prob = c(0, 1, 0.3, 0.7, 0, 1),
             kappa = c(1, 2, 1.5, 2.5, 3, 4),
             delta = c(2, 1, 3.0, 0.5, 0.5, 2),
             type = 6, unifsamp = rep(0.25, 6), direct = FALSE)
  expect_length(rg6, 6)
  expect_equal(rg6[1], 0.25^(1 / 2))
  expect_equal(rg6[2], 0.25^(1 / 2))
  expect_equal(rg6[5], 0.25^(1 / 0.5))
  expect_equal(rg6[6], 0.25^(1 / 4))

  ## rzidiscegpd forwards unifsamp
  u_samp <- c(0.1, 0.4, 0.7, 0.9)
  expect_equal(
    rzidiscegpd(4, pi = 0, sigma = 3, xi = 0.2, kappa = 1.4, type = 1, unifsamp = u_samp),
    rdiscegpd(4, sigma = 3, xi = 0.2, kappa = 1.4, type = 1, unifsamp = u_samp)
  )

  ## Zero-length contract across d/p/q
  expect_identical(p.G(numeric(0), type = 1, kappa = 1.5), numeric(0))
  expect_identical(d.G(numeric(0), type = 1, kappa = 1.5), numeric(0))
  expect_identical(q.G(numeric(0), type = 1, kappa = 1.5), numeric(0))
  for (fn in list(pegpd, degpd_density, qegpd, ddiscegpd, pdiscegpd, qdiscegpd)) {
    expect_identical(fn(numeric(0), sigma = 2, xi = 0.2, kappa = 1.2, type = 1), numeric(0))
    expect_identical(fn(0.5, sigma = numeric(0), xi = 0.2, kappa = 1.2, type = 1), numeric(0))
  }
  for (fn in list(dzidiscegpd, pzidiscegpd, qzidiscegpd, dziegpd, pziegpd, qziegpd)) {
    expect_identical(fn(numeric(0), pi = 0.2, sigma = 2, xi = 0.2, kappa = 1.2, type = 1), numeric(0))
    expect_identical(fn(0.5, pi = numeric(0), sigma = 2, xi = 0.2, kappa = 1.2, type = 1), numeric(0))
  }

  ## .qG6_logsurv NA/NaN guard + NA_real_ vs NaN preservation in q*
  ls_na <- egpd:::.qG6_logsurv(c(0.5, 0.5, NA_real_, NaN),
                              prob = c(NA_real_, 0.5, 0.5, 0.5),
                              kappa = c(1, NaN, 1, 1),
                              delta = 2)
  expect_true(is.na(ls_na[1]) && !is.nan(ls_na[1]))
  expect_true(is.nan(ls_na[2]))
  expect_true(is.na(ls_na[3]) && !is.nan(ls_na[3]))
  expect_true(is.nan(ls_na[4]))

  p_mix <- c(NA_real_, NaN, -0.1, 0.5)
  for (qfn in list(
    function(p) q.G(p, type = 4, delta = 2),
    function(p) qegpd(p, sigma = 2, xi = 0.2, kappa = 1.3, type = 1),
    function(p) qdiscegpd(p, sigma = 2, xi = 0.2, kappa = 1.3, type = 1),
    function(p) qzidiscegpd(p, pi = 0.2, sigma = 2, xi = 0.2, kappa = 1.3, type = 1),
    function(p) qziegpd(p, pi = 0.2, sigma = 2, xi = 0.2, kappa = 1.3, type = 1)
  )) {
    res <- qfn(p_mix)
    expect_true(is.na(res[1]) && !is.nan(res[1]))
    expect_true(is.nan(res[2]))
    expect_true(is.nan(res[3]))
    expect_true(is.finite(res[4]))
  }

  ## log, lower.tail, and log.p contracts including extreme upper tail
  for (tp in 1:6) {
    args <- list(prob = 0.35, kappa = 1.4, delta = 2.2, sigma = 2.5, xi = 0.25, type = tp)
    ## Continuous roundtrip in upper tail on log scale
    lp_u <- do.call(pegpd, c(list(q = c(1, 20, 500), lower.tail = FALSE, log.p = TRUE), args))
    expect_true(all(is.finite(lp_u) & lp_u < 0))
    q_back <- do.call(qegpd, c(list(p = lp_u, lower.tail = FALSE, log.p = TRUE), args))
    expect_equal(q_back, c(1, 20, 500), tolerance = 1e-8)

    ## Discrete log PMF and upper-tail CDF/quantile consistency
    ld <- do.call(ddiscegpd, c(list(x = c(-1, 0, 5, 200), log = TRUE), args))
    expect_equal(ld[1], -Inf)
    expect_true(all(is.finite(ld[2:4])))
    expect_equal(exp(ld[2:3]),
                 do.call(ddiscegpd, c(list(x = c(0, 5), log = FALSE), args)),
                 tolerance = 1e-12)

    s_disc <- do.call(pdiscegpd, c(list(q = 0:6, lower.tail = FALSE), args))
    c_disc <- do.call(pdiscegpd, c(list(q = 0:6, lower.tail = TRUE), args))
    expect_equal(s_disc + c_disc, rep(1, 7), tolerance = 1e-14)
    expect_equal(
      do.call(qdiscegpd, c(list(p = s_disc, lower.tail = FALSE), args)),
      0:6
    )
    ls_mid_d <- do.call(pdiscegpd, c(list(q = 0:5, lower.tail = FALSE, log.p = TRUE), args)) +
                0.5 * (do.call(pdiscegpd, c(list(q = 1:6, lower.tail = FALSE, log.p = TRUE), args)) -
                       do.call(pdiscegpd, c(list(q = 0:5, lower.tail = FALSE, log.p = TRUE), args)))
    expect_equal(
      do.call(qdiscegpd, c(list(p = ls_mid_d, lower.tail = FALSE, log.p = TRUE), args)),
      1:6
    )

    ## ZI discrete & ZI continuous log / lower.tail / log.p
    zi_args <- c(list(pi = 0.25), args)
    lzd <- do.call(dzidiscegpd, c(list(x = c(-1, 0, 3, 100), log = TRUE), zi_args))
    expect_equal(lzd[1], -Inf)
    expect_equal(exp(lzd[2:3]),
                 do.call(dzidiscegpd, c(list(x = c(0, 3), log = FALSE), zi_args)),
                 tolerance = 1e-12)

    szd <- do.call(pzidiscegpd, c(list(q = 0:5, lower.tail = FALSE), zi_args))
    czd <- do.call(pzidiscegpd, c(list(q = 0:5, lower.tail = TRUE), zi_args))
    expect_equal(szd + czd, rep(1, 6), tolerance = 1e-14)
    expect_equal(
      do.call(qzidiscegpd, c(list(p = szd, lower.tail = FALSE), zi_args)),
      0:5
    )
    ls_mid_z <- do.call(pzidiscegpd, c(list(q = 0:4, lower.tail = FALSE, log.p = TRUE), zi_args)) +
                0.5 * (do.call(pzidiscegpd, c(list(q = 1:5, lower.tail = FALSE, log.p = TRUE), zi_args)) -
                       do.call(pzidiscegpd, c(list(q = 0:4, lower.tail = FALSE, log.p = TRUE), zi_args)))
    expect_equal(
      do.call(qzidiscegpd, c(list(p = ls_mid_z, lower.tail = FALSE, log.p = TRUE), zi_args)),
      1:5
    )

    lpc_zi <- do.call(pziegpd, c(list(q = c(2, 50), lower.tail = FALSE, log.p = TRUE), zi_args))
    expect_equal(
      do.call(qziegpd, c(list(p = lpc_zi, lower.tail = FALSE, log.p = TRUE), zi_args)),
      c(2, 50), tolerance = 1e-8
    )
  }
})

test_that("R2-GAMLSS: d/p/q/r wrappers, score recycling, identity sigma, boundary nd fallback, and gamlss() fits", {
  skip_if_not_installed("gamlss")

  ## Check all 3-param and 4-param GAMLSS d/p/q/r wrappers
  for (nm in c("EGPD1", "EGPD3", "EGPD5", "EGPD6")) {
    d_fn <- get(paste0("d", nm)); p_fn <- get(paste0("p", nm))
    q_fn <- get(paste0("q", nm)); r_fn <- get(paste0("r", nm))
    set.seed(1)
    rs <- r_fn(5, mu = 2, sigma = 0.2, nu = 1.4)
    expect_length(rs, 5)
    lp <- p_fn(c(1, 5), mu = 2, sigma = 0.2, nu = 1.4, lower.tail = FALSE, log.p = TRUE)
    expect_equal(q_fn(lp, mu = 2, sigma = 0.2, nu = 1.4, lower.tail = FALSE, log.p = TRUE),
                 c(1, 5), tolerance = 1e-8)
    expect_equal(exp(d_fn(1.5, mu = 2, sigma = 0.2, nu = 1.4, log = TRUE)),
                 d_fn(1.5, mu = 2, sigma = 0.2, nu = 1.4, log = FALSE))
  }
  expect_equal(
    qEGPD4(pEGPD4(2.5, mu = 2, sigma = 0.2, nu = 1.3, tau = 1.6, lower.tail = FALSE, log.p = TRUE),
           mu = 2, sigma = 0.2, nu = 1.3, tau = 1.6, lower.tail = FALSE, log.p = TRUE),
    2.5, tolerance = 1e-8
  )

  for (nm in c("DEGPD1", "DEGPD3", "DEGPD5", "DEGPD6")) {
    d_fn <- get(paste0("d", nm)); p_fn <- get(paste0("p", nm))
    q_fn <- get(paste0("q", nm)); r_fn <- get(paste0("r", nm))
    set.seed(1)
    expect_length(r_fn(5, mu = 3, sigma = 0.2, nu = 1.3), 5)
    s <- p_fn(0:4, mu = 3, sigma = 0.2, nu = 1.3, lower.tail = FALSE, log.p = TRUE)
    expect_equal(q_fn(s, mu = 3, sigma = 0.2, nu = 1.3, lower.tail = FALSE, log.p = TRUE), 0:4)
    expect_equal(exp(d_fn(2, mu = 3, sigma = 0.2, nu = 1.3, log = TRUE)),
                 d_fn(2, mu = 3, sigma = 0.2, nu = 1.3, log = FALSE))
  }
  expect_equal(
    qDEGPD4(pDEGPD4(0:4, mu = 3, sigma = 0.2, nu = 1.3, tau = 1.5, lower.tail = FALSE, log.p = TRUE),
            mu = 3, sigma = 0.2, nu = 1.3, tau = 1.5, lower.tail = FALSE, log.p = TRUE),
    0:4
  )

  for (nm in c("ZIEGPD1", "ZIEGPD3", "ZIEGPD5", "ZIEGPD6",
               "ZIDEGPD1", "ZIDEGPD3", "ZIDEGPD5", "ZIDEGPD6")) {
    d_fn <- get(paste0("d", nm)); p_fn <- get(paste0("p", nm))
    q_fn <- get(paste0("q", nm)); r_fn <- get(paste0("r", nm))
    set.seed(1)
    expect_length(r_fn(5, mu = 2.5, sigma = 0.2, nu = 1.3, tau = 0.2), 5)
    expect_equal(exp(d_fn(0, mu = 2.5, sigma = 0.2, nu = 1.3, tau = 0.2, log = TRUE)),
                 d_fn(0, mu = 2.5, sigma = 0.2, nu = 1.3, tau = 0.2, log = FALSE))
    lp <- p_fn(3, mu = 2.5, sigma = 0.2, nu = 1.3, tau = 0.2, lower.tail = FALSE, log.p = TRUE)
    if (startsWith(nm, "ZIDEGPD")) {
      expect_equal(q_fn(lp, mu = 2.5, sigma = 0.2, nu = 1.3, tau = 0.2, lower.tail = FALSE, log.p = TRUE), 3)
    } else {
      expect_equal(q_fn(lp, mu = 2.5, sigma = 0.2, nu = 1.3, tau = 0.2, lower.tail = FALSE, log.p = TRUE),
                   3, tolerance = 1e-8)
    }
  }

  ## Recycling: scalar y with length-4 parameter vectors in 4-parameter families
  mu_v <- c(2, 2.5, 3, 3.5)
  sig_v <- c(0.1, 0.2, 0.25, 0.3)
  nu_v <- c(1.1, 1.3, 1.5, 1.7)
  tau_v <- c(0.15, 0.2, 0.25, 0.3)
  for (f4 in list(EGPD4(), DEGPD4(), ZIEGPD1(), ZIDEGPD1())) {
    ## y = 0 has infinite density for the continuous EGPD4 (kappa < 1), where
    ## the score is undefined (NaN), so evaluate that family at y = 0.5
    y0 <- if (identical(f4$type, "Continuous")) 0.5 else 0
    for (fn_nm in c("dldm", "dldd", "dldv", "dldt")) {
      got <- f4[[fn_nm]](y0, mu_v, sig_v, nu_v, tau_v)
      expect_length(got, 4)
      want <- vapply(1:4, function(i) f4[[fn_nm]](y0, mu_v[i], sig_v[i], nu_v[i], tau_v[i]), 0)
      expect_equal(got, want, tolerance = 1e-12)
    }
  }

  ## Identity sigma.link permits sigma <= 0 and fallback matches FD at sigma = 0 & boundary
  fam_id <- DEGPD1(sigma.link = "identity")
  expect_true(fam_id$sigma.valid(c(-0.4, 0, 0.3)))
  expect_false(fam_id$sigma.valid(c(0, NA_real_)))
  expect_false(DEGPD1(sigma.link = "log")$sigma.valid(-0.1))

  ## Near-zero sigma score check
  h <- 1e-5
  fd0 <- (dDEGPD1(3, mu = 2.5, sigma = h, nu = 1.4, log = TRUE) -
            dDEGPD1(3, mu = 2.5, sigma = -h, nu = 1.4, log = TRUE)) / (2 * h)
  expect_equal(fam_id$dldd(3, mu = 2.5, sigma = 0, nu = 1.4), fd0, tolerance = 1e-7)

  ## Bounded-tail top-integer score check: y = 9, mu = 2, sigma = -0.21 => 1 + sigma*(y+1)/mu = -0.05 <= 0
  expect_silent(sc_top <- fam_id$dldd(9, mu = 2, sigma = -0.21, nu = 1.5))
  fd_top <- (dDEGPD1(9, mu = 2, sigma = -0.21 + 1e-5, nu = 1.5, log = TRUE) -
               dDEGPD1(9, mu = 2, sigma = -0.21 - 1e-5, nu = 1.5, log = TRUE)) / (2e-5)
  expect_equal(sc_top, fd_top, tolerance = 1e-3)
  expect_gt(abs(sc_top), 1)

  ## ZI mu.initial safe when sum(y > 0) <= 1
  for (zfam in list(ZIEGPD1(), ZIEGPD3(), ZIDEGPD1(), ZIDEGPD3())) {
    for (y_edge in list(c(0, 0, 0), c(0, 0, 4))) {
      e <- list2env(list(y = y_edge), parent = baseenv())
      eval(zfam$mu.initial, envir = e)
      expect_length(e$mu, length(y_edge))
      expect_true(all(is.finite(e$mu) & e$mu > 0))
    }
  }

  ## End-to-end gamlss::gamlss() fits
  skip_on_cran()
  set.seed(505)
  ctrl <- gamlss::gamlss.control(n.cyc = 30, trace = FALSE)

  y_c <- rEGPD1(120, mu = 2.5, sigma = 0.2, nu = 1.4)
  g1 <- suppressWarnings(gamlss::gamlss(y_c ~ 1, family = EGPD1(), control = ctrl))
  expect_s3_class(g1, "gamlss")

  y_d <- rDEGPD1(120, mu = 3, sigma = 0.15, nu = 1.3)
  g2 <- suppressWarnings(gamlss::gamlss(y_d ~ 1, family = DEGPD1(), control = ctrl))
  expect_s3_class(g2, "gamlss")

  y_zid <- rZIDEGPD1(150, mu = 3, sigma = 0.15, nu = 1.3, tau = 0.25)
  g3 <- suppressWarnings(gamlss::gamlss(y_zid ~ 1, family = ZIDEGPD1(), control = ctrl))
  expect_s3_class(g3, "gamlss")
})

test_that("R2-POT: fit_degpd_pot and degpd_pot_stability validate thresholds and fix.kappa", {
  set.seed(606)
  x <- rdiscegpd(250, sigma = 4, xi = 0.15, kappa = 1, type = 1)

  expect_error(fit_degpd_pot(x, threshold = 2.5), "single finite integer")
  expect_error(fit_degpd_pot(x, threshold = c(2L, 4L)), "single finite integer")
  expect_error(fit_degpd_pot(x, threshold = NA_integer_), "single finite integer")
  expect_error(fit_degpd_pot(x, threshold = Inf), "single finite integer")
  expect_error(fit_degpd_pot(x, threshold = 2L, type = 4, fix.kappa = TRUE),
               "fix.kappa = FALSE")

  pot1 <- fit_degpd_pot(x, threshold = 2L, type = 1, fix.kappa = TRUE)
  expect_s3_class(pot1, "degpd_pot")
  expect_identical(pot1$threshold, 2L)
  expect_true(all(is.finite(pot1$estimate)))

  pot4 <- fit_degpd_pot(x, threshold = 2L, type = 4, fix.kappa = FALSE, hessian = FALSE)
  expect_s3_class(pot4, "degpd_pot")
  expect_true("delta" %in% names(pot4$estimate))

  ## Bounded tail POT check (upper.endpoint < Inf when xi < 0)
  x_b <- rdiscegpd(250, sigma = 6, xi = -0.4, kappa = 1, type = 1)
  pot_b <- fit_degpd_pot(x_b, threshold = 2L, type = 1, fix.kappa = TRUE, hessian = FALSE)
  expect_lt(pot_b$estimate[["xi"]], 0)
  expect_true(is.finite(pot_b$upper.endpoint))

  stab <- degpd_pot_stability(x, thresholds = c(1L, 3L, 5L))
  expect_s3_class(stab, "data.frame")
  expect_identical(
    names(stab),
    c("threshold", "n_exceed", "prop_exceed", "sigma", "xi", "xi_lo", "xi_hi", "upper_endpoint")
  )
  expect_equal(nrow(stab), 3L)
})
