## Round-3 univariate audit regressions (U3-xx findings).
##
## Reference values labelled "mpmath" were computed at 40 significant digits
## with an independent mpmath implementation of the GPD and of the carriers
## (written from the source papers, not from the package code) and are
## hard-coded here so the tests need no extra packages.

r3_fd <- function(f, x, h = 1e-3 * max(abs(x), 1e-2)) {
  ## Richardson-extrapolated central difference (independent of the package)
  d1 <- (f(x + h) - f(x - h)) / (2 * h)
  d2 <- (f(x + h / 2) - f(x - h / 2)) / h
  (4 * d2 - d1) / 3
}

## ---------------------------------------------------------------------------
## U3-08: GAMLSS scores
## ---------------------------------------------------------------------------

## mpmath scores: family, y, mu, sigma (xi), nu, tau, log-density, dldm, dldd, dldv, dldt
r3_score_ref <- read.table(header = TRUE, text = "
fam      y      mu  sigma   nu   tau  ld                 dldm                dldd                 dldv                  dldt
DEGPD1   0      2    0.2    0.8  NA  -0.776009203656618  -0.297813601444992  -0.144177454758610   -0.970011504570773    NA
DEGPD1   0      0.7 -0.25   0.4  NA  -0.0749118735641784 -0.261542256205574  -0.149906325496852   -0.187279683910446    NA
DEGPD1   0      3    3e-05  2.5  NA  -3.15164454429103   -0.702143349403276  -0.351070504468574   -1.26065781771641     NA
DEGPD3   0      2    0.2    0.4  NA  -2.21027335977081   -0.781619507963820  -0.378397463047757    0.566743383874941    NA
DEGPD4   0      2    0.2    1.5  1.5 -1.33933526440969   -0.537321009755496  -0.260127728213189    0.198713670215997   -0.892890176273124
DEGPD4   0      0.7  0.9    0.4  0.6 -0.367143007762153  -0.391168143077581  -0.142887445608725    0.134939112797850   -0.611905012936921
DEGPD4   0      3   -0.25   2.5  3   -3.06245350941988   -0.804507651987291  -0.413917768897230    0.282742250011988   -1.02081783647329
DEGPD5   0      2    0.2    0.5  NA  -1.05563194647262   -0.399193653355620  -0.193257543031442   -0.174650944176903    NA
DEGPD6   0      2    0.2    0.5  NA  -0.715845803309183  -0.273895448917240  -0.132598204055326   -0.484920148326313    NA
ZIDEGPD1 0      2    0.2    0.5  0.3 -0.313361575137575  -0.109743414100694  -0.0531289573235986  -0.571914103114379    0.525737260309490
ZIDEGPD3 0      2    0.2    0.5  0.3 -0.964816499093262  -0.164885707016990  -0.0798244319548801   0.110673065428332    2.32043721392658
ZIDEGPD5 0      2    0.2    0.5  0.3 -0.609577063689452  -0.178880291434133  -0.0865994870627672  -0.0782617948231243   1.19950454013788
ZIDEGPD6 0      2    0.2    0.5  0.3 -0.442941193159888  -0.145935815646447  -0.0706504147424337  -0.258373104223355    0.796115361269648
DEGPD1   7      2    0.2    0.8  NA  -4.25553525274000    0.796447781321814   1.16424136445859     1.18631983668001     NA
DEGPD1   1e5    2    0.2    0.8  NA  -56.1789629304836    2.49970003149668    200.264133845281     1.25                 NA
DEGPD5   1000   0.7  0.2    1    NA  -33.4413055168581    7.11297664747782    111.579176593758     0.145562547386397    NA
DEGPD6   1000   0.7  0.2    1    NA  -33.5972293105155    7.11297664747782    111.579176593758     0.448270555210260    NA
DEGPD4   1e5    3    0.2    0.4  1.5 -52.9637129260832    1.66636669816634    190.130253842052    -1.78571329545536     0.666666666666667
EGPD4    1e5    0.7  0.2    0.4  0.6 -61.1557198682910    7.14255714985714    226.505987246381    -1.78571422287497     1.66666666666667
EGPD4    0.01   3    0.2    0.4  0.6  0.562406560647404  -0.199111533067360  -0.00266126280011446  0.213397103061102   -4.21711916921305
EGPD1    1000   3    3e-05  0.4  NA  -333.702548456998    109.680968096810    54493.0271724573     2.5                  NA
DEGPD3   1000   3    3e-05  0.4  NA  -331.693981141383    109.732439902921    54544.3596984292    -1.78571428571429     NA
DEGPD1   11     3   -0.25   1.5  NA  -9.53417354754058    14.6664898355605    136.239850985163     0.666642553721631    NA
ZIEGPD5  0      2    0.2    1.5  0.3 -1.20397280432594    0                   0                    0                    3.33333333333333
ZIEGPD6  2.5    2    0.2    1.5  0.3 -2.24246619123586    0.0478970265235569 -0.481703609158505    0.242086576740791   -1.42857142857143
EGPD5    0.5    2    0.2    0.005 NA -0.986589893110731  -0.357508281698536  -0.208907194832094   -0.140401044750986    NA
EGPD6    40     0.7  0.9    2.5  NA  -7.53584359391098    1.52915846155028    2.58580999170346     0.215341968133613    NA
")

r3_scores <- function(fam, y, mu, sigma, nu, tau = NA) {
  f <- get(fam)(sigma.link = "identity")
  args <- list(y = y, mu = mu, sigma = sigma, nu = nu)
  if (f$nopar == 4) args$tau <- tau
  nm <- c("dldm", "dldd", "dldv", "dldt")[seq_len(f$nopar)]
  c(ld = -do.call(f$G.dev.incr, args) / 2,
    vapply(nm, function(n) do.call(f[[n]], args), 0))
}

test_that("U3-08: GAMLSS scores match high-precision references (y = 0, far tail, small xi)", {
  skip_if_not_installed("gamlss.dist")
  for (i in seq_len(nrow(r3_score_ref))) {
    r <- r3_score_ref[i, ]
    got <- r3_scores(r$fam, r$y, r$mu, r$sigma, r$nu, r$tau)
    want <- unlist(r[c("ld", "dldm", "dldd", "dldv", "dldt")])[seq_along(got)]
    info <- paste(r$fam, "y =", r$y, "mu =", r$mu, "xi =", r$sigma, "nu =", r$nu, "tau =", r$tau)
    ## the log-density is checked elsewhere (U3-01/U3-10); scores to 1e-8
    expect_equal(unname(got[-1]), unname(want[-1]), tolerance = 1e-8, info = info)
  }
})

test_that("U3-08: zero-count scores equal the derivative of log P(0) = nu log F(1)", {
  lp0 <- function(s, xi, k) k * log1p(-(1 + xi / s)^(-1 / xi))
  for (p in list(c(2, .2, .8), c(.7, -.25, .4), c(3, .5, 2.5), c(1.5, .25, .6))) {
    f <- DEGPD1(sigma.link = "identity")
    got <- c(f$dldm(0, p[1], p[2], p[3]), f$dldd(0, p[1], p[2], p[3]),
             f$dldv(0, p[1], p[2], p[3]))
    want <- c(r3_fd(function(s) lp0(s, p[2], p[3]), p[1]),
              r3_fd(function(x) lp0(p[1], x, p[3]), p[2]),
              log1p(-(1 + p[2] / p[1])^(-1 / p[2])))
    expect_equal(got, want, tolerance = 1e-8, info = paste(p, collapse = ","))
    ## the RS working weights use -score^2, so zeros must not get the 1e-15 floor
    expect_equal(f$d2ldm2(0, p[1], p[2], p[3]), -want[1]^2, tolerance = 1e-8)
    expect_equal(f$d2ldd2(0, p[1], p[2], p[3]), -want[2]^2, tolerance = 1e-8)
    expect_equal(f$d2ldmdd(0, p[1], p[2], p[3]), -want[1] * want[2], tolerance = 1e-8)
  }
  ## zero-inflated zeros: d log(pi + (1 - pi) P0) = (1 - pi) P0 d log P0 / f(0)
  fz <- ZIDEGPD1()
  P0 <- exp(lp0(2, .2, .5))
  w0 <- .7 * P0 / (.3 + .7 * P0)
  expect_equal(fz$dldm(0, 2, .2, .5, .3), w0 * r3_fd(function(s) lp0(s, .2, .5), 2), tolerance = 1e-8)
  expect_equal(fz$dldt(0, 2, .2, .5, .3), (1 - P0) / (.3 + .7 * P0), tolerance = 1e-12)
  ## carrier 5 (DEGPD4): log P(0) = (kappa / 2) log G4(F(1))
  G4 <- function(u, d) 1 - (1 - u) * (1 + (1 - (1 - u)^d) / d)
  lp05 <- function(s, xi, d, k) (k / 2) * log(G4(1 - (1 + xi / s)^(-1 / xi), d))
  f4 <- DEGPD4()
  expect_equal(f4$dldm(0, 2, .2, 1.5, 1.2), r3_fd(function(s) lp05(s, .2, 1.5, 1.2), 2), tolerance = 1e-8)
  expect_equal(f4$dldd(0, 2, .2, 1.5, 1.2), r3_fd(function(x) lp05(2, x, 1.5, 1.2), .2), tolerance = 1e-8)
  expect_equal(f4$dldt(0, 2, .2, 1.5, 1.2), lp05(2, .2, 1.5, 1.2) / 1.2, tolerance = 1e-8)
})

test_that("U3-08: vectorised scores recycle against element-wise references", {
  r <- r3_score_ref[r3_score_ref$fam == "DEGPD1", ]
  f <- DEGPD1(sigma.link = "identity")
  got <- f$dldd(r$y, r$mu, r$sigma, r$nu)
  expect_equal(got, r$dldd, tolerance = 1e-8)
  r4 <- r3_score_ref[r3_score_ref$fam == "DEGPD4", ]
  f4 <- DEGPD4(sigma.link = "identity")
  expect_equal(f4$dldm(r4$y, r4$mu, r4$sigma, r4$nu, r4$tau), r4$dldm, tolerance = 1e-8)
  expect_equal(f4$dldt(r4$y, r4$mu, r4$sigma, r4$nu, r4$tau), r4$dldt, tolerance = 1e-8)
})

test_that("U3-08: undefined scores are NaN with a warning, never a silent zero", {
  f <- DEGPD1(sigma.link = "identity")
  ## y = 30 lies beyond the upper endpoint 2 / 0.2 = 10: zero likelihood
  expect_warning(s <- f$dldm(c(30, 3), mu = 2, sigma = -0.2, nu = 1.5), "zero or infinite likelihood")
  expect_true(is.nan(s[1]))
  expect_true(is.finite(s[2]) && s[2] != 0)
  expect_warning(s <- f$dldd(1, mu = -1, sigma = 0.2, nu = 1.5), "NaNs produced")
  expect_true(is.nan(s))
  s <- f$dldv(c(NA, 1), mu = 2, sigma = 0.2, nu = 1.5)
  expect_true(is.na(s[1]) && !is.nan(s[1]))
  ## continuous density infinite at y = 0 when kappa < 1
  expect_warning(s <- EGPD1()$dldm(0, 2, 0.2, 0.5), "zero or infinite likelihood")
  expect_true(is.nan(s))
  ## ...but finite where the zero is an atom of the zero-inflated model
  expect_identical(ZIEGPD1()$dldm(0, 2, 0.2, 0.5, 0.3), 0)
})

test_that("U3-08: gamlss DEGPD1 reaches the closed-form MLE on zero-heavy counts", {
  skip_if_not_installed("gamlss")
  set.seed(2)
  y <- rdiscegpd(400, kappa = .6, sigma = 1.5, xi = .25)
  pmf <- function(y, s, xi, k) {
    F <- function(t) (1 - (1 + xi * t / s)^(-1 / xi))^k
    F(y + 1) - F(y)
  }
  nll <- function(p) -sum(log(pmf(y, exp(p[1]), exp(p[2]), exp(p[3]))))
  o <- optim(log(c(1.5, .25, .6)), nll, method = "BFGS", control = list(reltol = 1e-12))
  m <- gamlss::gamlss(y ~ 1, sigma.formula = ~1, nu.formula = ~1, family = DEGPD1(),
                      data = data.frame(y = y), trace = FALSE)
  expect_true(m$converged)
  expect_equal(as.numeric(logLik(m)), -o$value, tolerance = 0.01 / 528)
  expect_lt(abs(as.numeric(logLik(m)) + o$value), 0.01)
  expect_equal(unname(fitted(m, "sigma")[1]), exp(o$par[2]), tolerance = 0.02)
})

test_that("U3-08: gamlss DEGPD4 (kappa < 2) and EGPD3 fits reach independent MLEs", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  G5 <- function(u, d, k) (1 - (1 - u) * (1 + (1 - (1 - u)^d) / d))^(k / 2)
  H <- function(t, s, xi, d, k) G5(ifelse(t <= 0, 0, 1 - (1 + xi * t / s)^(-1 / xi)), d, k)
  set.seed(2026)
  y4 <- rdiscegpd(400, delta = 1.5, kappa = 1.2, sigma = 2, xi = 0.2, type = 5)
  nll4 <- function(p) {
    s <- exp(p[1]); xi <- exp(p[2]); d <- exp(p[3]); k <- exp(p[4])
    -sum(log(H(y4 + 1, s, xi, d, k) - H(y4, s, xi, d, k)))
  }
  o4 <- optim(log(c(2, .2, 1.5, 1.2)), nll4, method = "BFGS", control = list(reltol = 1e-12, maxit = 2000))
  m4 <- gamlss::gamlss(y4 ~ 1, sigma.formula = ~1, nu.formula = ~1, tau.formula = ~1,
                       family = DEGPD4(), data = data.frame(y4 = y4), trace = FALSE)
  expect_lt(abs(as.numeric(logLik(m4)) + o4$value), 0.01)
  expect_gt(fitted(m4, "sigma")[1], 0.01)
  ## EGPD3 (carrier 4) has a bimodal delta profile on these data
  set.seed(2026)
  invisible(regpd(400, kappa = 1.6, sigma = 2, xi = 0.2))
  y3 <- regpd(400, delta = 2, sigma = 2, xi = 0.2, type = 4)
  nll3 <- function(p) {
    s <- exp(p[1]); xi <- exp(p[2]); d <- exp(p[3])
    u <- 1 - (1 + xi * y3 / s)^(-1 / xi)
    -sum(log((1 + 1 / d) * (1 - (1 - u)^d)) - log(s) - (1 / xi + 1) * log1p(xi * y3 / s))
  }
  o3 <- optim(log(c(2, .1, 2)), nll3, method = "BFGS", control = list(reltol = 1e-12))
  m3 <- gamlss::gamlss(y3 ~ 1, sigma.formula = ~1, nu.formula = ~1, family = EGPD3(),
                       data = data.frame(y3 = y3), trace = FALSE)
  expect_lt(abs(as.numeric(logLik(m3)) + o3$value), 0.01)
  ## zero-inflated DEGPD1 with kappa < 1
  set.seed(2026)
  yz <- rzidiscegpd(400, pi = 0.25, kappa = 0.7, sigma = 2, xi = 0.2)
  Hz <- function(t, s, xi, k) ifelse(t <= 0, 0, (1 - (1 + xi * t / s)^(-1 / xi))^k)
  nllz <- function(p) {
    s <- exp(p[1]); xi <- exp(p[2]); k <- exp(p[3]); pz <- plogis(p[4])
    pm <- Hz(yz + 1, s, xi, k) - Hz(yz, s, xi, k)
    -sum(log(ifelse(yz == 0, pz + (1 - pz) * pm, (1 - pz) * pm)))
  }
  oz <- optim(c(log(2), log(.2), log(.7), qlogis(.25)), nllz, method = "BFGS",
              control = list(reltol = 1e-12, maxit = 2000))
  mz <- gamlss::gamlss(yz ~ 1, sigma.formula = ~1, nu.formula = ~1, tau.formula = ~1,
                       family = ZIDEGPD1(), data = data.frame(yz = yz), trace = FALSE)
  expect_lt(as.numeric(logLik(mz)) + oz$value, 0.01)
  expect_gt(as.numeric(logLik(mz)) + oz$value, -0.01)
})

test_that("numeric score helpers handle parameters at a domain boundary", {
  ## ZIBell zero-inflation probability at 1e-9 (count reviewer): analytic scores
  theta <- uniroot(function(t) t * exp(t) - 5, c(0, 5), tol = 1e-14)$root
  p0 <- exp(1 - exp(theta))
  expect_silent(g <- .egpd_nd2(c(0, 3), 5, 1e-9, dZIBELL, idx = 2L))
  expect_equal(g, c((1 - p0) / (1e-9 + (1 - 1e-9) * p0), -1 / (1 - 1e-9)), tolerance = 1e-6)
  ## a log-density with a singularity just beyond the upper parameter bound
  ld <- function(y, mu, log = TRUE) {
    out <- suppressWarnings(y * log1p(-mu) + log(mu))
    out[!(mu > 0 & mu < 1)] <- NaN
    out
  }
  mu <- 1 - 5e-5
  expect_equal(.egpd_nd1(c(1, 3), mu, ld), -c(1, 3) / (1 - mu) + 1 / mu, tolerance = 1e-6)
  expect_warning(g <- .egpd_nd1(1, 1.5, ld), "derivative undefined")
  expect_true(is.nan(g))
})

## ---------------------------------------------------------------------------
## U3-04 / U3-05 / integer threshold: spliced DEGPD predictions
## ---------------------------------------------------------------------------

## Independent implementation of the documented splice (Maia et al. Eq. 3):
## bulk CDF = linear interpolation of (0, 0) and the fitted quantile knots
## (upper probability at ties), alpha = F_bulk(floor(q)), tail from
## s = max(floor(q) + 1, 0) with a kappa = 1 discrete GPD.
r3_splice_cdf <- function(y, q, p, thr, scale, shape) {
  Fb <- function(t) stats::approx(c(0, pmax(q[seq_len(thr)], 0)), c(0, p[seq_len(thr)]),
                                  xout = t, ties = max, rule = 2)$y
  qt <- q[thr]
  alpha <- Fb(floor(qt))
  s <- max(floor(qt) + 1, 0)
  H <- function(k) {
    z <- 1 + shape * (k + 1) / scale
    ifelse(z <= 0, 1, 1 - z^(-1 / shape))
  }
  ifelse(y < 0, 0, ifelse(y < s, pmin(Fb(pmax(y, 0)), alpha), alpha + (1 - alpha) * H(y - s)))
}

test_that("U3-04/U3-05: spliced predict handles any number of rows and at values", {
  skip_if_not_installed("qgam")
  set.seed(1)
  n <- 300
  x <- runif(n)
  y <- rnbinom(n, mu = exp(1 + x), size = 2)
  hit <- runif(n) < 0.1
  y[hit] <- y[hit] + rdiscegpd(sum(hit), kappa = 1, sigma = 4, xi = 0.2)
  d <- data.frame(y = y, x = x)
  fit <- suppressWarnings(fit_spliced_degpd(y ~ x, data = d, bulk_probs = c(0.1, 0.5, 0.85),
                                            threshold_prob = 0.85))
  nd <- d[c(5, 9, 2), ]
  par <- predict(fit, newdata = nd, type = "parameter")
  bulk <- as.matrix(predict(fit, newdata = nd, type = "bulk"))
  expect_identical(rownames(par), c("5", "9", "2"))
  ## U3-05: a single row keeps its row name (was the threshold probability "0.85")
  expect_identical(rownames(predict(fit, newdata = nd[1, ], type = "parameter")), "5")
  ref <- function(i, at) r3_splice_cdf(at, bulk[i, ], fit$bulk_probs, fit$threshold_index,
                                       par$scale[i], par$shape[i])
  for (rows in list(1L, 1:3)) {
    for (at in list(10, c(5, 10), 0, c(3, 3), c(10, 0, 3), 0:30)) {
      cdf <- predict(fit, newdata = nd[rows, , drop = FALSE], type = "cdf", at = at)
      pmf <- predict(fit, newdata = nd[rows, , drop = FALSE], type = "pmf", at = at)
      expect_identical(dim(cdf), c(length(rows), length(at)))
      expect_identical(dim(pmf), c(length(rows), length(at)))
      expect_identical(names(cdf), paste0("F:", at))
      expect_identical(names(pmf), paste0("P:", at))
      expect_identical(rownames(cdf), rownames(nd)[rows])
      want_cdf <- t(vapply(rows, function(i) ref(i, at), numeric(length(at))))
      want_pmf <- t(vapply(rows, function(i) ref(i, at) - ref(i, at - 1), numeric(length(at))))
      if (length(at) == 1L) {
        want_cdf <- matrix(want_cdf, ncol = 1)
        want_pmf <- matrix(want_pmf, ncol = 1)
      }
      expect_equal(unname(as.matrix(cdf)), unname(want_cdf), tolerance = 1e-12,
                   info = paste("rows", length(rows), "at", paste(at, collapse = ",")))
      expect_equal(unname(as.matrix(pmf)), unname(want_pmf), tolerance = 1e-12,
                   info = paste("rows", length(rows), "at", paste(at, collapse = ",")))
    }
  }
  expect_error(predict(fit, newdata = nd, type = "cdf", at = numeric(0)), "non-negative integers")
})

test_that("spliced splice gives the tail exactly the mass of the fitted exceedances for integer thresholds", {
  probs <- c(.1, .25, .5, .75, .9)
  ob <- list(tail_fit = list(), bulk_probs = probs, threshold_index = 5, threshold_prob = .9)
  tp <- data.frame(scale = 2, shape = 0.2, kappa = 1)
  for (q in list(c(0, 1, 2, 4, 5), c(0, 1, 2, 4, 5.5), c(0.3, 0.6, 1.1, 1.7, 2))) {
    cdf <- drop(egpd:::.spliced_degpd_cdf_matrix(ob, data.frame(x = 0), 0:20,
                                                 bulk_mat = matrix(q, nrow = 1), tail_pars = tp))
    expect_equal(unname(cdf), r3_splice_cdf(0:20, q, probs, 5, 2, 0.2), tolerance = 1e-12,
                 info = paste(q, collapse = ","))
    s <- floor(q[5]) + 1
    ## the tail (fitted to Y - s | Y >= s) carries P(Y >= s) = 1 - alpha_tilde
    alpha <- stats::approx(c(0, q), c(0, probs), xout = floor(q[5]), ties = max)$y
    expect_equal(unname(1 - cdf[s]), 1 - alpha, tolerance = 1e-12)
  }
  expect_equal(egpd:::.spliced_degpd_shift(c(5, 5.5, 0, -0.3, -2, 4.999)), c(6, 6, 1, 0, 0, 5))
})

## ---------------------------------------------------------------------------
## U3-01: log / log.p evaluated on the log scale (no underflow)
## ---------------------------------------------------------------------------

test_that("U3-01: log-scale p/d functions do not underflow (closed forms)", {
  ## exponential carrier: S(y) = exp(-y), P(X = k) = exp(-k) (1 - exp(-1))
  expect_equal(pegpd(c(700, 800, 2000), kappa = 1, sigma = 1, xi = 0,
                     lower.tail = FALSE, log.p = TRUE), -c(700, 800, 2000))
  expect_equal(pdiscegpd(800, kappa = 1, sigma = 1, xi = 0, lower.tail = FALSE, log.p = TRUE), -801)
  expect_equal(ddiscegpd(800, kappa = 1, sigma = 1, xi = 0, log = TRUE), -800 + log1p(-exp(-1)))
  ## kappa = 2, xi = 0.2: S_G = 2 s - s^2 with s = (1 + 0.2 y)^(-5)
  ls <- -5 * log1p(0.2 * 1e70)
  expect_equal(pegpd(1e70, kappa = 2, sigma = 1, xi = 0.2, lower.tail = FALSE, log.p = TRUE),
               log(2) + ls, tolerance = 1e-14)
  ## large kappa: log P(0) = kappa log F(1) far below the double range
  expect_equal(ddiscegpd(0, kappa = 500, sigma = 1, xi = 0.2, log = TRUE),
               500 * log1p(-1.2^-5), tolerance = 1e-14)
  ## zero-inflated variants
  expect_equal(pzidiscegpd(800, pi = 0.2, kappa = 1, sigma = 1, xi = 0,
                           lower.tail = FALSE, log.p = TRUE), log(0.8) - 801)
  expect_equal(dzidiscegpd(0, pi = 0, kappa = 500, sigma = 1, xi = 0.2, log = TRUE),
               500 * log1p(-1.2^-5), tolerance = 1e-14)
  expect_equal(pziegpd(1e-200, pi = 0, kappa = 2, sigma = 1, xi = 0.2, log.p = TRUE),
               2 * log(1e-200), tolerance = 1e-14)
})

test_that("U3-01: log-scale quantiles invert far-tail log probabilities", {
  expect_equal(qegpd(-800, kappa = 1, sigma = 1, xi = 0, lower.tail = FALSE, log.p = TRUE), 800)
  expect_equal(qdiscegpd(-800.5, kappa = 1, sigma = 1, xi = 0, lower.tail = FALSE, log.p = TRUE), 800)
  expect_equal(qziegpd(log(0.8) - 800, pi = 0.2, kappa = 1, sigma = 1, xi = 0,
                       lower.tail = FALSE, log.p = TRUE), 800, tolerance = 1e-12)
  expect_equal(qzidiscegpd(log(0.8) - 800.5, pi = 0.2, kappa = 1, sigma = 1, xi = 0,
                           lower.tail = FALSE, log.p = TRUE), 800)
  ## mpmath: upper-tail quantiles at log p = log(1e-300) and round trips at -800
  want <- c(2.04576856589892689e+31, 2.03954190956501102e+31, 2.08275948798482117e+31,
            2.02369588183561776e+31, 2.12108096322803745e+31)
  for (t in 2:6) {
    args <- list(prob = 0.4, kappa = 1.5, delta = 2, sigma = 2, xi = 0.1, type = t)
    expect_equal(do.call(qegpd, c(list(p = log(1e-300), lower.tail = FALSE, log.p = TRUE), args)),
                 want[t - 1], tolerance = 1e-10, info = paste("type", t))
    q <- do.call(qegpd, c(list(p = -800, lower.tail = FALSE, log.p = TRUE), args))
    expect_true(is.finite(q))
    expect_equal(do.call(pegpd, c(list(q = q, lower.tail = FALSE, log.p = TRUE), args)), -800,
                 tolerance = 1e-12, info = paste("type", t))
  }
})

test_that("U3-01: lower-tail log CDF at 1e-200 for every carrier (mpmath)", {
  want <- c(-461.225093045932851, -462.281528425589780,
            log(2) + 2 * log(1e-200), log(2) + 2 * log(1e-200), -921.950327929492429)
  for (t in 2:6) {
    got <- pegpd(1e-200, prob = 0.4, kappa = 2, delta = 3, sigma = 1, xi = 0.2, type = t, log.p = TRUE)
    expect_equal(got, want[t - 1], tolerance = 1e-13, info = paste("type", t))
  }
  ## and the discrete log-PMF far beyond the double range (mpmath)
  expect_equal(ddiscegpd(1000, kappa = 0.4, sigma = 0.7, xi = 3e-5, log = TRUE),
               -1400.01615669576879, tolerance = 1e-13)
  expect_equal(ddiscegpd(1e5, prob = 0.4, kappa = 0.7, delta = 2.5, sigma = 0.7, xi = 0.2,
                         type = 6, log = TRUE), -60.6279266656315784, tolerance = 1e-13)
})

## ---------------------------------------------------------------------------
## U3-09: lower-tail quantiles of carriers 2, 3 and 6 for tiny p
## ---------------------------------------------------------------------------

test_that("U3-09: small-p quantiles of carriers 2, 3 and 6 are relatively accurate", {
  ## type-6 carrier quantile: root of 0.4 u^2 + 0.6 u^3 = 1e-15
  root <- uniroot(function(lu) log(0.4 * exp(2 * lu) + 0.6 * exp(3 * lu)) - log(1e-15),
                  c(-30, -5), tol = 1e-15)$root
  expect_equal(q.G(1e-15, type = 6, prob = .4, kappa = 2, delta = 3), exp(root), tolerance = 1e-12)
  ## mpmath quantiles
  expect_equal(qegpd(1e-15, kappa = 2, sigma = 1, xi = .2, type = 2),
               2.0300784692787033e-15, tolerance = 1e-11)
  expect_equal(qegpd(1e-15, prob = .3, kappa = .7, delta = 1.5, sigma = 2.5, xi = 0, type = 6),
               5.2040111238179969e-21, tolerance = 1e-11)
  expect_equal(qegpd(1e-12, kappa = 0.5, sigma = 1, xi = .2, type = 3),
               4.5112815393605258e-13, tolerance = 1e-11)
  expect_equal(qegpd(1e-6, kappa = 1e-4, sigma = 1, xi = .3, type = 3),
               2.2180800897815211e-07, tolerance = 1e-11)
  expect_equal(qegpd(1e-12, kappa = 50, sigma = 1, xi = .2, type = 2),
               9.9371696560705192e-03, tolerance = 1e-11)
  ## p(q(p)) / p - 1 stays tiny down to p = 1e-300
  p <- 10^-c(9, 12, 15, 50, 300)
  for (t in c(2, 3, 6)) {
    args <- list(prob = 0.4, kappa = 2, delta = 3, sigma = 1, xi = 0.2, type = t)
    q <- do.call(qegpd, c(list(p = p), args))
    expect_lt(max(abs(do.call(pegpd, c(list(q = q), args)) / p - 1)), 1e-10)
  }
})

## ---------------------------------------------------------------------------
## U3-10 / U3-11: far-tail density and cancellation near y = 0
## ---------------------------------------------------------------------------

test_that("U3-10/U3-11: far-tail densities and near-zero tails are accurate (mpmath)", {
  expect_equal(degpd_density(1e10, delta = 1e-3, sigma = 1, xi = .3, type = 4),
               6.0114612135609698e-40, tolerance = 1e-11)
  expect_equal(degpd_density(1e10, delta = 1e-3, sigma = 1, xi = .3, type = 4, log = TRUE),
               -90.3097358703899100, tolerance = 1e-13)
  expect_equal(pegpd(1e-8, kappa = .3, delta = .5, sigma = 1, xi = .8, type = 5, lower.tail = FALSE),
               0.99618706687332044, tolerance = 1e-13)
  expect_equal(pegpd(1e-8, kappa = 1e-4, sigma = 1, xi = .3, type = 3),
               4.5084110327822204e-08, tolerance = 1e-12)
})

## ---------------------------------------------------------------------------
## U3-02 / U3-03 / U3-12: contracts
## ---------------------------------------------------------------------------

test_that("U3-02: random generators follow base R for zero-length and vector n", {
  a <- list(kappa = 1.2, sigma = 2, xi = 0.1)
  expect_identical(do.call(regpd, c(list(numeric(0)), a)), numeric(0))
  expect_identical(do.call(rdiscegpd, c(list(numeric(0)), a)), numeric(0))
  expect_identical(do.call(rzidiscegpd, c(list(numeric(0), pi = 0.2), a)), numeric(0))
  expect_identical(do.call(rziegpd, c(list(numeric(0), pi = 0.2), a)), numeric(0))
  expect_identical(r.G(numeric(0), kappa = 1.2), numeric(0))
  set.seed(1)
  expect_length(do.call(rzidiscegpd, c(list(c(9, 9, 9), pi = 0.2), a)), 3)
  expect_length(do.call(regpd, c(list(2.7), a)), 2)
  expect_error(do.call(regpd, c(list(-1), a)), "invalid arguments")
  expect_error(do.call(rziegpd, c(list(NA, pi = 0.2), a)), "invalid arguments")
})

test_that("U3-03: missing required parameters are an error; supplied NA gives NA", {
  for (f in list(pegpd, degpd_density, qegpd, ddiscegpd, pdiscegpd, qdiscegpd)) {
    expect_error(f(c(.5, 1), sigma = 1, xi = .1, type = 1), "`kappa' missing")
    expect_error(f(.5, kappa = 1, type = 1), "`sigma', `xi' missing")
    expect_error(f(.5, sigma = 1, xi = .1, kappa = 1, type = 4), "`delta' missing")
    expect_true(all(is.na(f(c(.5, .7), kappa = NA, sigma = 1, xi = .1, type = 1))))
  }
  for (f in list(dzidiscegpd, pzidiscegpd, qzidiscegpd, dziegpd, pziegpd, qziegpd)) {
    expect_error(f(.5, sigma = 1, xi = .1, kappa = 1), "`pi' missing")
    expect_error(f(.5, pi = .1, sigma = 1, xi = .1, kappa = 1, type = 6), "`delta', `prob' missing")
  }
  expect_error(regpd(2, sigma = 1, xi = .1, kappa = 1, type = 4), "`delta' missing")
  expect_error(rdiscegpd(2, sigma = 1, xi = .1, kappa = 1, type = 4), "`delta' missing")
  expect_error(rzidiscegpd(2, sigma = 1, xi = .1, kappa = 1), "`pi' missing")
  expect_error(rziegpd(2, pi = .1, kappa = 1), "`sigma', `xi' missing")
  expect_error(d.G(.5, type = 6, kappa = 1, delta = 2), "`prob' missing")
  expect_error(q.G(.5, type = 1), "`kappa' missing")
})

test_that("U3-12: rzidiscegpd(unifsamp =) is reproducible inverse-CDF sampling for pi > 0", {
  u <- c(.05, .3, .5, .9, .95, .2)
  a <- rzidiscegpd(6, pi = .3, kappa = 1.4, sigma = 3, xi = .2, unifsamp = u)
  set.seed(9)
  b <- rzidiscegpd(6, pi = .3, kappa = 1.4, sigma = 3, xi = .2, unifsamp = u)
  expect_identical(a, b)
  ## smallest k with F(k) >= u, F(k) = pi + (1 - pi) H(k + 1) (closed form)
  F <- function(k) 0.3 + 0.7 * (1 - (1 + 0.2 * (k + 1) / 3)^(-5))^1.4
  want <- vapply(u, function(ui) { k <- 0; while (F(k) < ui) k <- k + 1; k }, 0)
  expect_equal(a, want)
})

test_that("U3-12: degpd_pot_stability passes hessian through and reports failed thresholds", {
  set.seed(42)
  x1 <- rdiscegpd(3000, kappa = 2, sigma = 2, xi = 0.25)
  st <- degpd_pot_stability(x1, thresholds = c(4, 8), hessian = FALSE)
  expect_identical(nrow(st), 2L)
  expect_true(all(is.na(st$xi_lo)))
  for (u in c(4, 8)) {
    k <- x1[x1 >= u] - u
    expect_equal(st$n_exceed[st$threshold == u], length(k))
  }
  expect_warning(st2 <- degpd_pot_stability(x1, thresholds = c(4, 1e6)), "no exceedances")
  expect_identical(st2$threshold, 4L)
  expect_error(degpd_pot_stability(x1, thresholds = 1e6), "no threshold could be fitted")
})

test_that("U3-12: upper.endpoint is the largest count with positive probability", {
  set.seed(42)
  invisible(rdiscegpd(3000, kappa = 2, sigma = 2, xi = 0.25))
  x2 <- c(rpois(1500, 3), 5 + rdiscegpd(800, kappa = 1, sigma = 6, xi = -0.3))
  f2 <- fit_degpd_pot(x2, 5)
  s <- f2$estimate[["sigma"]]
  xi <- f2$estimate[["xi"]]
  expect_lt(xi, 0)
  ## Z = X - u has P(Z >= k) = (1 + xi k / s)^(-1/xi) > 0 exactly when k < -s / xi
  k <- 0:500
  expect_equal(f2$upper.endpoint, 5 + max(k[1 + xi * k / s > 0]))
  expect_equal(f2$upper.endpoint, floor(f2$upper.endpoint))
})

## ---------------------------------------------------------------------------
## Test gap: value-level checks of the model 5/6 carriers (types 2 and 3)
## against closed forms written from arXiv:2409.18719 (independent of the
## package code), including a far lower-tail point.
## ---------------------------------------------------------------------------

test_that("types 2 and 3 d/p/q/r match closed-form carriers", {
  om <- 1 / 32
  G <- list(
    `2` = function(u, k) (pnorm(sqrt(k) * (u - 1)) - pnorm(-sqrt(k))) / (0.5 - pnorm(-sqrt(k))),
    `3` = function(u, k) (pbeta(om + (0.5 - om) * u, k, k) - pbeta(om, k, k)) /
      (pbeta(0.5, k, k) - pbeta(om, k, k)))
  g <- list(
    `2` = function(u, k) sqrt(k) * dnorm(sqrt(k) * (u - 1)) / (0.5 - pnorm(-sqrt(k))),
    `3` = function(u, k) (0.5 - om) * dbeta(om + (0.5 - om) * u, k, k) /
      (pbeta(0.5, k, k) - pbeta(om, k, k)))
  Fg <- function(y, s, xi) 1 - (1 + xi * y / s)^(-1 / xi)
  fg <- function(y, s, xi) (1 + xi * y / s)^(-1 / xi - 1) / s
  s <- 3; xi <- 0.2; k <- 2; pz <- 0.3
  y <- c(0.1, 0.5, 2, 7, 20)
  for (t in c("2", "3")) {
    tp <- as.integer(t)
    H <- function(v) ifelse(v <= 0, 0, G[[t]](Fg(v, s, xi), k))
    expect_equal(degpd_density(y, sigma = s, xi = xi, kappa = k, type = tp),
                 g[[t]](Fg(y, s, xi), k) * fg(y, s, xi), tolerance = 1e-10)
    expect_equal(pegpd(y, sigma = s, xi = xi, kappa = k, type = tp), H(y), tolerance = 1e-10)
    expect_equal(pegpd(y, sigma = s, xi = xi, kappa = k, type = tp, lower.tail = FALSE),
                 1 - H(y), tolerance = 1e-10)
    pp <- c(0.1, 0.5, 0.9)
    want_q <- vapply(pp, function(p) uniroot(function(v) H(v) - p, c(1e-9, 1e4), tol = 1e-13)$root, 0)
    expect_equal(qegpd(pp, sigma = s, xi = xi, kappa = k, type = tp), want_q, tolerance = 1e-9)
    kk <- 0:15
    expect_equal(ddiscegpd(kk, sigma = s, xi = xi, kappa = k, type = tp), H(kk + 1) - H(kk),
                 tolerance = 1e-10)
    expect_equal(pdiscegpd(kk, sigma = s, xi = xi, kappa = k, type = tp), H(kk + 1),
                 tolerance = 1e-10)
    qd <- vapply(pp, function(p) { j <- 0; while (H(j + 1) < p) j <- j + 1; j }, 0)
    expect_equal(qdiscegpd(pp, sigma = s, xi = xi, kappa = k, type = tp), qd)
    u <- c(0.05, 0.4, 0.8, 0.97)
    want_r <- floor(vapply(u, function(p) uniroot(function(v) H(v) - p, c(1e-9, 1e4), tol = 1e-13)$root, 0))
    expect_equal(rdiscegpd(4, sigma = s, xi = xi, kappa = k, type = tp, unifsamp = u), want_r)
    expect_equal(dzidiscegpd(kk, pi = pz, sigma = s, xi = xi, kappa = k, type = tp),
                 (kk == 0) * pz + (1 - pz) * (H(kk + 1) - H(kk)), tolerance = 1e-10)
    expect_equal(pziegpd(y, pi = pz, sigma = s, xi = xi, kappa = k, type = tp),
                 pz + (1 - pz) * H(y), tolerance = 1e-10)
    ## far lower tail: G(u) ~ g(0) u as u -> 0 (relative error O(u))
    expect_equal(pegpd(1e-12, sigma = s, xi = xi, kappa = k, type = tp),
                 g[[t]](0, k) * 1e-12 / s, tolerance = 1e-9)
  }
})

## ---------------------------------------------------------------------------
## Follow-up review: far-tail PMF at large sigma, small-p quantiles of
## carriers 2/3, the subnormal-G4 band, gamlss starting values and
## rzidiscegpd recycling
## ---------------------------------------------------------------------------

test_that("far-tail PMF and scores at large sigma are accurate (mpmath)", {
  skip_if_not_installed("gamlss.dist")
  ## mpmath log P(y) and gamlss scores; mu = GPD scale, sigma = xi
  ref <- read.table(header = TRUE, text = "
fam     y      mu    sigma nu    ld                    dldm                     dldd                   dldv
DEGPD1  1e5    2000  0.2   1.5   -21.5828393652905001  0.00222726691111258086   32.6746402335028994    0.666660457575274341
DEGPD1  5000   3000  0.2   1.5   -9.46257377301773781  0.000101889182582007609  -0.454440766008198019  0.395808896763231125
DEGPD1  1e7    3000  0.2   1.5   -46.6236370410789866  0.00166367116007632122   132.639660021924925    0.666666666666659130
DEGPD6  1e5    2000  0.9   0.3   -16.0697351114042514  0.000532743799798319985  2.43291163795887338    0.669991188898337419
DEGPD5  1e5    2000  0.9   1e-5  -15.6835985125092615  0.000532608806794248724  2.43203756518986694    0.166565531267322940
DEGPD6  400    0.7   0.9   4     -12.1109762669238810  1.58144822041908660      5.36990582137839534    0.132614693061821987
DEGPD3  3e5    3000  0.2   1.5   -25.7626860941302801  0.00157142872231088119   47.5416697769615681   -0.266666664822332850
")
  for (i in seq_len(nrow(ref))) {
    r <- ref[i, ]
    got <- r3_scores(r$fam, r$y, r$mu, r$sigma, r$nu)
    info <- paste(r$fam, "y =", r$y, "mu =", r$mu, "xi =", r$sigma, "nu =", r$nu)
    expect_lt(abs(got[["ld"]] - r$ld), 1e-9, label = info)
    expect_lt(max(abs(got[-1] / unlist(r[c("dldm", "dldd", "dldv")]) - 1)), 1e-7, label = info)
  }
})

test_that("lower-tail quantiles of carriers 2 and 3 are exact for 1e-10 <= p < 0.01 (mpmath)", {
  ## type, kappa, p, mpmath quantile (sigma = 1, xi = 0.2)
  ref <- read.table(text = "
2 0.5    1e-3  0.0011850847242766683
2 0.5    1e-7  1.1845931220527844e-7
2 2.5    1e-10 2.4516986765418805e-10
2 50     1e-5  0.49336069134068257
2 0.0001 1e-9  1.0000333345999995e-9
3 0.0001 1e-8  2.2180762381323476e-9
3 0.5    1e-5  4.5113676075124966e-6
3 2.5    1e-10 1.4883637140599578e-9
3 20     1e-3  0.69378760604662252")
  for (i in seq_len(nrow(ref))) {
    expect_equal(qegpd(ref$V3[i], kappa = ref$V2[i], sigma = 1, xi = 0.2, type = ref$V1[i]),
                 ref$V4[i], tolerance = 1e-12, info = paste(ref[i, 1:3], collapse = " "))
  }
  p <- 10^-c(2.5, 3, 5, 7, 9, 10)
  for (k in c(1e-4, 0.5, 2.5, 50)) for (t in 2:3) {
    q <- qegpd(p, kappa = k, sigma = 1, xi = 0.2, type = t)
    lp <- pegpd(q, kappa = k, sigma = 1, xi = 0.2, type = t, log.p = TRUE)
    expect_lt(max(abs(lp - log(p))), 1e-12, label = paste("type", t, "kappa", k))
  }
  expect_silent(qegpd(c(1e-6, 1e-3), kappa = 20, sigma = 1, xi = 0.2, type = 2))
})

test_that("log CDF and density of carriers 4 and 5 are exact where G4 is subnormal (mpmath)", {
  ## sigma = 1, xi = 0: u = 1 - exp(-x) ~ x and G4(u) ~ (1 + delta) / 2 u^2 < 2.2e-308
  expect_equal(pegpd(2.2e-162, delta = 2, sigma = 1, xi = 0, type = 4, log.p = TRUE),
               -744.0551903012341, tolerance = 1e-14)
  expect_equal(pegpd(1e-157, delta = 0.5, sigma = 1, xi = 0, type = 4, log.p = TRUE),
               -723.29940127258213, tolerance = 1e-14)
  expect_equal(pegpd(1.4e-154, delta = 7, sigma = 1, xi = 0, type = 4, log.p = TRUE),
               -707.13696980780375, tolerance = 1e-14)
  expect_equal(pegpd(2.2e-162, delta = 2, kappa = 3, sigma = 1, xi = 0, type = 5, log.p = TRUE),
               -1116.0827854518511, tolerance = 1e-14)
  expect_equal(pegpd(3e-155, delta = 5, kappa = 0.7, sigma = 1, xi = 0, type = 5, log.p = TRUE),
               -248.67693968675243, tolerance = 1e-14)
  expect_equal(pegpd(1e-158, delta = 0.25, kappa = 6, sigma = 1, xi = 0, type = 5, log.p = TRUE),
               -2184.2606790460925, tolerance = 1e-14)
  expect_equal(degpd_density(1e-158, delta = 2, kappa = 3, sigma = 1, xi = 0, type = 5, log = TRUE),
               -725.91007943528808, tolerance = 1e-14)
  expect_equal(degpd_density(5e-160, delta = 0.5, kappa = 1.2, sigma = 1, xi = 0, type = 5, log = TRUE),
               -73.351123079999739, tolerance = 1e-13)
  expect_equal(d.G(1e-158, type = 5, delta = 2, kappa = 3, log = TRUE),
               -725.91007943528808, tolerance = 1e-14)
})

test_that("gamlss starting values are computed once per data set and reused by gamlssML", {
  skip_if_not_installed("gamlss")
  ## same length, sum, first and last value, different data
  y1 <- c(0, 1, 2, 5, 9, 3, 0, 1, 7, 4, 2, 2, 6, 1, 0)
  y2 <- y1
  y2[2] <- 3
  y2[4] <- 3
  s1 <- .egpd_gamlss_start(y1, "DEGPD1")
  s2 <- .egpd_gamlss_start(y2, "DEGPD1")
  expect_false(isTRUE(all.equal(s1, s2)))
  assign("last", NULL, envir = .egpd_start_cache)
  expect_identical(.egpd_gamlss_start(y2, "DEGPD1"), s2)
  ## gamlssML evaluates each *.initial in a fresh environment (U3-08 data)
  set.seed(2)
  y <- rdiscegpd(400, kappa = .6, sigma = 1.5, xi = .25)
  m <- gamlss::gamlssML(y, family = DEGPD1())
  expect_lt(abs(as.numeric(logLik(m)) - (-528.059827722)), 0.01)
})

test_that("gamlss identity-link fits of a bounded tail start inside the support", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  set.seed(4)
  y <- regpd(5000, kappa = 1.5, sigma = 2, xi = -0.7)
  st <- .egpd_gamlss_start(y, "EGPD1", sigma_link = "identity")
  expect_lt(st[["sigma"]], 0)
  expect_lt(max(y), -st[["mu"]] / st[["sigma"]])
  ## closed-form EGPD1 log-likelihood, maximised independently
  nll <- function(p) {
    s <- exp(p[1]); xi <- p[2]; k <- exp(p[3])
    z <- 1 + xi * y / s
    if (any(z <= 0)) return(1e100)
    -sum(log(k) + (k - 1) * log1p(-z^(-1 / xi)) - log(s) - (1 / xi + 1) * log(z))
  }
  o <- optim(c(log(2), -0.7, log(1.5)), nll, control = list(reltol = 1e-14, maxit = 5000))
  m <- gamlss::gamlss(y ~ 1, family = EGPD1(sigma.link = "identity"),
                      data = data.frame(y = y), trace = FALSE)
  expect_lt(abs(as.numeric(logLik(m)) + o$value), 0.01)
})

test_that("rzidiscegpd(unifsamp =) recycles the parameters to n like base R", {
  u <- c(.1, .95)
  x <- rzidiscegpd(2, pi = c(.1, .2, .3, .4), kappa = 1.4, sigma = 2, xi = .2, unifsamp = u)
  ## smallest k with F(k) >= u, F(k) = pi + (1 - pi) H(k + 1), using pi[1:2]
  F <- function(k, pz) pz + (1 - pz) * (1 - (1 + 0.2 * (k + 1) / 2)^(-5))^1.4
  want <- vapply(1:2, function(i) { k <- 0; while (F(k, c(.1, .2)[i]) < u[i]) k <- k + 1; k }, 0)
  expect_equal(x, want)
  expect_length(rzidiscegpd(3, pi = .2, kappa = c(1, 2), sigma = 2, xi = .2, unifsamp = .5), 3)
})
