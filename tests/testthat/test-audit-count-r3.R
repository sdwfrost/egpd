## Round 3 audit regression tests for the count families (findings C-1 to C-11 of
## files/egpd-review-r3/count/FINDINGS.md). References are independent of the package:
## closed forms, brute-force sums in R, or constants computed with Rmpfr / integrate().

## ---- independent references ------------------------------------------------------

## GPIG is compound Poisson: log G(s) = L0 + sum_j lam_j s^j with
## lam_j = b a Gamma(j - a) c^j / (Gamma(1 - a) j!), so (De Pril) k p_k = sum_j j lam_j p_{k-j}.
## A different recursion from the package's (Zhu & Joe), and every term is positive.
gpig_ref_lpmf <- function(Y, a, b, c) {
  j <- seq_len(Y)
  jl <- j * exp(log(b) + log(a) + lgamma(j - a) - lgamma(1 - a) - lgamma(j + 1) + j * log(c))
  p <- numeric(Y + 1); p[1] <- 1                     # scaled by p_0
  for (k in j) p[k + 1] <- sum(jl[seq_len(k)] * p[k:1]) / k
  b * ((1 - c)^a - 1) + log(p)
}

## log(B_n / n!) by Dobinski's formula, B_n = e^-1 sum_k k^n / k!.
bell_ref_lratio <- function(n) vapply(n, function(m) {
  if (m == 0) return(0)
  k <- 1:max(50, 6 * m)
  t <- m * log(k) - lgamma(k + 1)
  -1 + max(t) + log(sum(exp(t - max(t)))) - lgamma(m + 1)
}, 0)
bell_ref_lpmf <- function(y, theta) y * log(theta) - expm1(theta) + bell_ref_lratio(y)

lse <- function(v) max(v) + log(sum(exp(v - max(v))))

pig_bessel_log <- function(y, mu, s) {
  al <- sqrt(1 / s^2 + 2 * mu / s)
  0.5 * log(2 * al / pi) + y * log(mu) + 1 / s +
    log(besselK(al, y - 0.5, expon.scaled = TRUE)) - al - y * log(al * s) - lgamma(y + 1)
}

gpois_ref_lpmf <- function(y, mu, lam) {
  th <- mu * (1 - lam)
  log(th) + (y - 1) * log(th + lam * y) - th - lam * y - lgamma(y + 1)
}

## Generalized Waring = beta-negative-binomial; for integer k its survivor is the finite sum
## P(Y > n) = sum_{j<k} choose(n+k, j) B(rho + j, a + n + k - j) / B(rho, a).
gw_ref_lsurv <- function(n, mu, k, rho) {
  a <- mu * (rho - 1) / k
  vapply(n, function(m) lse(lchoose(m + k, 0:(k - 1)) +
                              lbeta(rho + 0:(k - 1), a + m + k - 0:(k - 1)) - lbeta(rho, a)), 0)
}

lp_gpig <- function(y, a, b, c) egpd:::gpig_logp_cpp(as.integer(y), a, b, c)

## ---- C-1: GPIG hybrid --------------------------------------------------------------

test_that("GPIG hybrid with ymax = Inf is exact at any count (independent recursion)", {
  old <- gpig_control()
  on.exit(do.call(gpig_control, old))
  gpig_control("hybrid", ymax = Inf)
  expect_equal(lp_gpig(1500, 0.2, 5, 0.999)[["logp"]],
               gpig_ref_lpmf(1500, 0.2, 5, 0.999)[1501], tolerance = 1e-11)
  ## hard-coded from gpig_ref_lpmf (the old hybrid used the saddlepoint here: errors ~0.2)
  refs <- rbind(c(0.2, 5, 0.999, 3000, -12.334036327520),
                c(0.2, 5, 0.999, 10000, -20.598516878903),
                c(0.1, 2, 0.9995, 1500, -10.427344007762),
                c(0.1, 2, 0.9995, 3000, -11.879741607670),
                c(0.1, 2, 0.9995, 10000, -16.609911028466))
  for (i in seq_len(nrow(refs)))
    expect_equal(lp_gpig(refs[i, 4], refs[i, 1], refs[i, 2], refs[i, 3])[["logp"]],
                 refs[i, 5], tolerance = 1e-10)
})

test_that("GPIG log pmf is continuous across the old saddlepoint switches", {
  old <- gpig_control()
  on.exit(do.call(gpig_control, old))
  bs <- seq(4.3175, 4.3190, by = 5e-5)            # |corr| crosses 0.5 near b = 4.3182
  for (m in list(list("hybrid"), list("hybrid", ymax = Inf), list("saddlepoint"))) {
    do.call(gpig_control, m)
    v <- vapply(bs, function(b) lp_gpig(3000, 0.2, b, 0.999)[["logp"]], 0)
    expect_true(all(is.finite(v)))
    expect_lt(max(abs(diff(v))), 1e-3)            # was a 0.69-nat jump
  }
  ## the default (exact branch up to ymax = 1000, saddlepoint above) in the heavy-tailed
  ## regime, across the blend thresholds sptol/2, sptol and the taper (|corr| 0.25-0.5):
  ## log p is steep in b here, so jumps are sought in second differences (smooth: ~1e-3)
  gpig_control("hybrid")
  for (y in c(600, 1500)) {
    bs <- exp(seq(log(0.5), log(400), length.out = 3001))
    v <- vapply(bs, function(b) lp_gpig(y, 0.3, b, 0.9995)[["logp"]], 0)
    expect_true(all(is.finite(v)))
    expect_lt(max(abs(diff(v, differences = 2))), 0.01)
  }
  ## in the hybrid's blend region (sptol/2 < |corr| < sptol) the score is the derivative
  ## of the blended log pmf
  gpig_control("hybrid")
  a <- 0.5; b <- 20; cc <- 0.6; y <- 400
  v <- lp_gpig(y, a, b, cc)
  expect_true(v[["corr"]] > 5e-4 && v[["corr"]] < 1e-3)
  h <- 1e-6
  f <- function(a, b, cc) lp_gpig(y, a, b, cc)[["logp"]]
  fd <- c((f(a * (1 + h), b, cc) - f(a * (1 - h), b, cc)) / (2 * h * a),
          (f(a, b * (1 + h), cc) - f(a, b * (1 - h), cc)) / (2 * h * b),
          (f(a, b, cc * (1 + h)) - f(a, b, cc * (1 - h))) / (2 * h * cc))
  expect_equal(unname(v[c("da", "db", "dc")]), fd, tolerance = 1e-6)
})

test_that("a GPIG fit with ymax = Inf reports the exact log-likelihood at its optimum", {
  skip_on_cran()
  old <- gpig_control()
  on.exit(do.call(gpig_control, old))
  set.seed(7)
  a <- 0.2; c <- 0.999; mu <- 250
  y <- rgpig(500, a, mu * (1 - c)^(1 - a) / (a * c), c)
  expect_gt(sum(y > 1000), 10)
  d <- data.frame(y = y)
  fit <- egpd(list(lmu = y ~ 1, logita = ~1, logitc = ~1), data = d, family = "gpig",
              gpig.args = list(method = "hybrid", ymax = Inf))
  pr <- predict(fit, type = "response")[1, ]
  bb <- pr$mu * (1 - pr$c)^(1 - pr$a) / (pr$a * pr$c)
  exact <- sum(dgpig(y, pr$a, bb, pr$c, log = TRUE))
  expect_equal(as.numeric(logLik(fit)), exact, tolerance = 1e-9)   # was -3110.61 vs -3107.15
  rec <- egpd(list(lmu = y ~ 1, logita = ~1, logitc = ~1), data = d, family = "gpig",
              gpig.args = list(method = "recursion"))
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(rec)), tolerance = 1e-7)
})

test_that("default GPIG hybrid stays fast with covariates and large counts", {
  ## With covariates every observation is its own parameter row; an unbounded exact branch
  ## then costs O(sum y^2) per evaluation (448 s for a Tycho measles fit, 11 s for 150 of
  ## its weeks). The default bounds it at ymax = 1000, also when egpd() is given
  ## gpig.args = list(method = "hybrid"), which resets the other options to the defaults.
  skip_on_cran()
  old <- gpig_control()
  on.exit(do.call(gpig_control, old))
  egpd:::.gpig_apply_control("gpig", list(method = "hybrid"))
  set.seed(5)
  n <- 150; tt <- seq(0, 3, length.out = n)
  X1 <- cbind(1, tt, sin(2 * pi * tt), cos(2 * pi * tt)); X <- matrix(1, n, 1)
  beta <- c(8.5, 0.3, 1.2, -0.4)
  y <- stats::rnbinom(n, size = 3, mu = exp(drop(X1 %*% beta)))  # Tycho-like weekly counts
  expect_gt(max(y), 2e4)
  ld <- list(X = list(X1, X, X), y = cbind(y), dupid = 0L, duplicate = 0L,
             offsets = list(numeric(), numeric(), numeric()), idpars = c(1, 1, 1, 1, 2, 3),
             censored = FALSE)
  p <- c(beta, qlogis(0.6092), qlogis(0.999936))
  el <- system.time(v <- egpd:::.gpig1fns$d0(p, ld))[["elapsed"]]
  expect_true(is.finite(v))
  expect_lt(el, 2)                                   # 11-30 s with the exact branch unbounded
  el <- system.time(M <- egpd:::.gpig1fns$d120(p, ld))[["elapsed"]]
  expect_true(all(is.finite(M)))
  expect_lt(el, 2)
})

## ---- C-10: normalisation ------------------------------------------------------------

test_that("normalised GPIG score is analytic, consistent and fast with covariates", {
  skip_on_cran()
  old <- gpig_control()
  on.exit(do.call(gpig_control, old))
  gpig_control("hybrid", normalize = TRUE)
  set.seed(3)
  n <- 200; x <- seq(-1, 1, length.out = n); a <- 0.5; cc <- 0.6
  mu <- exp(2 + 0.5 * x)
  y <- rgpig(n, a, mu * (1 - cc)^(1 - a) / (a * cc), cc)
  X1 <- cbind(1, x); X <- matrix(1, n, 1); off <- list(numeric(), numeric(), numeric())
  ld <- list(X = list(X1, X, X), y = cbind(y), dupid = 0L, duplicate = 0L, offsets = off,
             idpars = c(1, 1, 2, 3), censored = FALSE)
  fns <- egpd:::.gpig1fns
  p <- c(2, 0.5, qlogis(a), qlogis(cc))
  t12 <- system.time(M <- fns$d120(p, ld))[["elapsed"]]
  expect_lt(t12, 3)                                  # was 11.3 s
  ag <- c(colSums(M[, 1] * X1), sum(M[, 2]), sum(M[, 3]))
  fd <- vapply(seq_along(p), function(i) {
    h <- 1e-5; pp <- p; pm <- p; pp[i] <- pp[i] + h; pm[i] <- pm[i] - h
    (fns$d0(pp, ld) - fns$d0(pm, ld)) / (2 * h)
  }, 0)
  expect_equal(unname(ag), fd, tolerance = 1e-6)
})

## ---- C-7: Bell table ---------------------------------------------------------------

test_that("Bell table is cached, stays exact, and no longer costs O(n^2) per call", {
  skip_on_cran()
  ## the cached, truncated recursion reproduces Dobinski's formula
  y <- c(0, 1, 7, 100, 1000, 5000)
  for (theta in c(0.3, 2, 6.5))
    expect_equal(dbell(y, theta, log = TRUE), bell_ref_lpmf(y, theta), tolerance = 1e-11)
  X <- matrix(1, 100, 1)
  y <- c(rep(3, 99), 16000)
  t1 <- system.time(egpd:::bell1d0(list(log(10)), X, y, integer(), 0L, list(numeric())))
  expect_lt(t1[["elapsed"]], 1)                      # was 3.7 s
  t2 <- system.time(rBELL(20, mu = 1000))
  expect_lt(t2[["elapsed"]], 0.5)                    # was 1.3 s
})

## ---- C-8: CMP score ----------------------------------------------------------------

fd_vs_score <- function(fns, npar, y, p, h = 1e-5) {
  n <- length(y)
  X <- c(list(cbind(1, seq(-1, 1, length.out = n))),
         replicate(npar - 1, matrix(1, n, 1), simplify = FALSE))
  ld <- list(X = X, y = cbind(y), dupid = 0L, duplicate = 0L,
             offsets = replicate(npar, numeric(0), simplify = FALSE),
             idpars = rep(seq_len(npar), vapply(X, ncol, 1L)), censored = FALSE)
  M <- fns$d120(p, ld)
  ag <- unlist(lapply(seq_len(npar), function(i) colSums(M[, i] * X[[i]])))
  fd <- vapply(seq_along(p), function(i) {
    pp <- p; pm <- p; pp[i] <- pp[i] + h; pm[i] <- pm[i] - h
    (fns$d0(pp, ld) - fns$d0(pm, ld)) / (2 * h)
  }, 0)
  max(abs(ag - fd)) / max(1, abs(fd))
}

ycount <- c(0, 0, 0, 1, 2, 3, 5, 8, 12, 4, 1, 25, 60, 0, 7)

test_that("CMP and Poisson-lognormal GAM scores are derivatives of their objectives", {
  cmp <- egpd:::.cmp1fns
  expect_lt(fd_vs_score(cmp, 2, ycount, c(1.5, 0.2, -1)), 1e-8)    # was 1.6e-5
  expect_lt(fd_vs_score(cmp, 2, ycount, c(1.0, 0.2, 0.7)), 1e-8)   # was 7.6e-6
  pln <- egpd:::.plnorm1fns
  expect_lt(fd_vs_score(pln, 2, ycount, c(1.5, 0.2, 1.0)), 1e-8)   # was 1.1e-4
  ## the gridded objective still equals the exact per-observation likelihood
  set.seed(5)
  n <- 300; y <- as.numeric(rpois(n, 20))
  X1 <- cbind(1, as.numeric(scale(seq_len(n)))); X2 <- matrix(1, n, 1)
  pars <- list(c(log(20), 0.15), log(0.6))
  grid <- egpd:::cmp1d0(pars, X1, X2, y, as.integer(0), 0L, list(numeric(0), numeric(0)))
  e1 <- as.vector(X1 %*% pars[[1]])
  slow <- -sum(vapply(seq_len(n), function(i)
    egpd:::cmp_logpmf_cpp(y[i], exp(e1[i]), 0.6), 0))
  expect_equal(grid, slow, tolerance = 1e-9)
})

## ---- C-9: Poisson-lognormal quadrature ---------------------------------------------

pln_ref_lpmf <- function(y, mu, s) {
  muz <- log(mu) - s^2 / 2
  h <- function(z) -exp(z) + y * z - (z - muz)^2 / (2 * s^2)
  zs <- optimize(function(z) -h(z), c(min(muz, log(y + 0.5)) - 10 * s - 5,
                                      max(muz, log(y + 0.5)) + 5), tol = 1e-12)$minimum
  ss <- 1 / sqrt(exp(zs) + 1 / s^2)
  I <- integrate(function(z) exp(h(z) - h(zs)), zs - 60 * ss, zs + 60 * ss,
                 rel.tol = 1e-13, subdivisions = 5000L)$value
  h(zs) + log(I) - lgamma(y + 1) - 0.5 * log(2 * pi) - log(s)
}

test_that("Poisson-lognormal is accurate at y = 0 for large sigma, in d* and the GAM", {
  X <- matrix(1, 1, 1)
  for (pr in list(c(50, 5), c(1, 4), c(50, 3))) {
    for (y in c(0, 1, 2)) {
      ref <- pln_ref_lpmf(y, pr[1], pr[2])
      expect_equal(dpoislnorm(y, pr[1], pr[2], log = TRUE), ref, tolerance = 1e-10)
      expect_equal(-egpd:::plnorm1d0(list(log(pr[1]), log(pr[2])), X, X, y, integer(), 0L,
                                     list(numeric(), numeric())), ref, tolerance = 1e-10)
    }
  }
  ## total mass: P(Y > 2e5) = 2.43730627483872e-06 from the independent representation
  ## P(Y > q) = P(G <= Lambda), G ~ Gamma(q + 1), Lambda lognormal (integrate())
  expect_equal(1 - sum(dpoislnorm(0:2e5, 20, 3)), 2.43730627483872e-06, tolerance = 1e-6)
  expect_equal(ppoislnorm(2e5, 20, 3, lower.tail = FALSE), 2.43730627483872e-06,
               tolerance = 1e-9)
})

## ---- C-4: upper tails --------------------------------------------------------------

test_that("upper tails are summed from the survivor and keep relative accuracy", {
  ## generalized Waring (power law): closed-form survivors for k = 1 and k = 2
  n <- c(10, 1e3, 1e5, 1e6)
  for (cs in list(c(1, 1, 2.5), c(2, 2, 3), c(5, 1, 1.5))) {
    ref <- gw_ref_lsurv(n, cs[1], cs[2], cs[3])
    expect_equal(pgwaring(n, cs[1], cs[2], cs[3], lower.tail = FALSE, log.p = TRUE), ref,
                 tolerance = 1e-11)
    expect_equal(pgwaring(n, cs[1], cs[2], cs[3], lower.tail = FALSE), exp(ref),
                 tolerance = 1e-11)
  }
  expect_equal(pgwaring(1e5, 1, 1, 2.5, lower.tail = FALSE), 2.140774997e-12, tolerance = 1e-9)
  ## Bell: Dobinski tail sum (was -36.737)
  ref <- lse(bell_ref_lpmf(61:300, 1.3))
  expect_equal(pbell(60, 1.3, lower.tail = FALSE, log.p = TRUE), ref, tolerance = 1e-11)
  expect_equal(pbell(60, 1.3, log.p = TRUE), -exp(ref), tolerance = 1e-10)
  expect_equal(pzibell(60, 1.3, 0.3, lower.tail = FALSE, log.p = TRUE), log(0.7) + ref,
               tolerance = 1e-11)
  ## generalized Poisson: closed-form pmf summed in R
  ref <- lse(gpois_ref_lpmf(201:2000, 3, 0.5))
  expect_equal(pgpois(200, 3, 0.5, lower.tail = FALSE, log.p = TRUE), ref, tolerance = 1e-11)
  ## CMP at nu = 1 is Poisson
  expect_equal(pcmp(60, 3, 1, lower.tail = FALSE, log.p = TRUE),
               ppois(60, 3, lower.tail = FALSE, log.p = TRUE), tolerance = 1e-9)
  ## PIG: Bessel closed-form pmf summed (the tail ratio here is 0.44)
  ref <- lse(pig_bessel_log(41:100, 2, 0.2))
  expect_equal(ppig(40, 2, 0.2, lower.tail = FALSE, log.p = TRUE), ref, tolerance = 1e-11)
  ## GPIG: independent recursion summed
  r <- gpig_ref_lpmf(900, 0.5, 20, 0.6)
  expect_equal(pgpig(400, 0.5, 20, 0.6, lower.tail = FALSE, log.p = TRUE), lse(r[402:901]),
               tolerance = 1e-10)
  ## Poisson-lognormal: independent Gamma representation (see above)
  expect_equal(ppoislnorm(1e5, 20, 3, lower.tail = FALSE), 7.15455381755115e-06,
               tolerance = 1e-9)
})

## ---- C-5: generalised inverse near p = 1 --------------------------------------------

test_that("upper-tail quantiles invert the survivor exactly", {
  S <- function(k) lse(gpois_ref_lpmf((k + 1):(k + 3000), 3, 0.5))
  ks <- 0:200; Sk <- vapply(ks, S, 0)
  for (pu in c(1e-9, 1e-12, 1e-13))
    expect_equal(qgpois(pu, 3, 0.5, lower.tail = FALSE), ks[which(Sk <= log(pu))[1]])
  expect_equal(qgpois(c(1e-12, 1e-13), 3, 0.5, lower.tail = FALSE), c(122, 134))
  Sb <- vapply(0:80, function(k) lse(bell_ref_lpmf((k + 1):(k + 300), 1.3)), 0)
  expect_equal(qbell(1e-13, 1.3, lower.tail = FALSE), which(Sb <= log(1e-13))[1] - 1)
  expect_equal(qbell(1e-13, 1.3, lower.tail = FALSE), 51)
  ## round trips on the survivor scale, and on the log scale
  expect_equal(qbell(pbell(0:60, 1.3, lower.tail = FALSE), 1.3, lower.tail = FALSE), 0:60)
  expect_equal(qbell(pbell(0:60, 1.3, lower.tail = FALSE, log.p = TRUE), 1.3,
                     lower.tail = FALSE, log.p = TRUE), 0:60)
  expect_equal(qzibell(pzibell(0:60, 1.3, 0.3, lower.tail = FALSE), 1.3, 0.3,
                       lower.tail = FALSE), 0:60)
  ## heavy tail: Waring closed-form survivor S(n) = B(a + n + 1, rho) / B(a, rho)
  a <- 1.5; rho <- 2.5
  sw <- function(n) lbeta(a + n + 1, rho) - lbeta(a, rho)
  q <- qgwaring(1e-12, 1, 1, rho, lower.tail = FALSE)
  expect_true(sw(q) <= log(1e-12) && sw(q - 1) > log(1e-12))
})

## ---- review of the integrated tree: tail engine ------------------------------------

test_that("an exactly zero or unresolved survivor never spoils the lower tail", {
  ## c = 0 is a documented point mass at zero: its survivor is exactly zero
  expect_silent(v <- pgpig(0:3, a = .5, b = 1, c = 0))
  expect_identical(v, rep(1, 4))                                 # was NA NA NA NA
  expect_silent(v <- pgpig(0:3, a = .5, b = 1, c = 0, lower.tail = FALSE))
  expect_identical(v, rep(0, 4))
  expect_identical(pgpig(0:3, .5, 1, 0, lower.tail = FALSE, log.p = TRUE), rep(-Inf, 4))
  expect_identical(pgpig(0:3, .5, 1, 0, log.p = TRUE), rep(0, 4))
  expect_silent(v <- pzigpig(0:3, .5, 1, 0, pi = .2))
  expect_identical(v, rep(1, 4))
  expect_identical(pzigpig(0:3, .5, 1, 0, pi = .2, lower.tail = FALSE), rep(0, 4))
})

test_that("the lower tail keeps the forward sum when the survivor cannot be certified", {
  skip_on_cran()
  ## c -> 1: the geometric bound on the remainder never certifies the survivor within the
  ## ladder limit, so only the upper tail is unresolved (NaN, with a warning)
  ref <- sum(exp(gpig_ref_lpmf(12000, .5, 1, .999)))             # De Pril recursion
  expect_silent(v <- pgpig(12000, .5, 1, .999))
  expect_equal(v, ref, tolerance = 1e-14)                        # was NA
  expect_warning(v <- pgpig(12000, .5, 1, .999, lower.tail = FALSE), "not resolved")
  expect_true(is.nan(v))
  el <- system.time(expect_silent(q <- qgpig(1 - 1e-9, .5, 1, .9998)))[["elapsed"]]
  expect_lt(el, 30)                                              # was 45 s, with warnings
  Fq <- pgpig(c(q - 1, q), .5, 1, .9998)
  expect_true(Fq[1] < 1 - 1e-9 && Fq[2] >= (1 - 1e-9) * (1 - 64 * .Machine$double.eps))
  ## q beyond the default summation limit (1e7) of the families without an exact survivor
  expect_silent(v <- pgpois(2e7, 5, .5))
  expect_identical(v, 1)                                         # was NA
  expect_identical(pcmp(2e7, 5, 1), 1)                           # was NA
  qq <- c(9.9e6, 1.01e7, 1.2e7)
  expect_equal(ppig(qq, 1e7, 1e-4), cumsum(dpig(0:1.2e7, 1e7, 1e-4))[qq + 1],
               tolerance = 1e-10)                                # was an error
})

test_that("PRL probabilities and quantiles use the closed-form survivor pointwise", {
  ## Alkhairy (2023) Eq. (3)/(4): S(x) = (1+1/tau)^-x tau (x+tau^2) / ((tau-1)(1+tau)^2)
  lS <- function(x, tau)
    -x * log1p(1 / tau) + log(tau) + log(x + tau^2) - log(tau - 1) - 2 * log1p(tau)
  expect_identical(pprl(3e9, 2.5), 1)                            # was an error
  expect_equal(pprl(3e9, 2.5, lower.tail = FALSE, log.p = TRUE), lS(3e9, 2.5),
               tolerance = 1e-14)
  el <- system.time(v <- pprl(3e7, tau = 1e6))[["elapsed"]]
  expect_lt(el, 0.5)                                             # was 3.5 s, O(q) memory
  expect_equal(v, -expm1(lS(3e7, 1e6)), tolerance = 1e-15)
  ## the generalised inverse of the same values, against a brute-force table
  Fr <- cumsum(dprl(0:2000, 2.5))
  p <- c(1e-6, 0.01, 0.1, 0.5, 0.9, 0.999, 1 - 1e-9)
  expect_identical(qprl(p, 2.5), vapply(p, function(pp) min(which(Fr >= pp)) - 1, 0))
  el <- system.time(q <- qprl(1e-9, tau = 1e5, lower.tail = FALSE))[["elapsed"]]
  expect_lt(el, 0.5)
  expect_true(lS(q, 1e5) <= log(1e-9) && lS(q - 1, 1e5) > log(1e-9))
})

## ---- C-2 / C-3: PIG contracts ------------------------------------------------------

test_that("dpig/dzipig follow the dpois conventions for invalid x", {
  expect_identical(dpig(-1, 2, 1), 0)                                   # was 1
  expect_warning(expect_identical(dpig(1.5, 2, 1), 0), "non-integer x") # was 0.195
  v <- dpig(c(NA, 1, NaN, Inf), 2, 1)                                   # NA used to error
  expect_true(is.na(v[1]) && !is.nan(v[1]))
  expect_equal(v[2], exp(pig_bessel_log(1, 2, 1)), tolerance = 1e-12)
  expect_true(is.nan(v[3]))
  expect_identical(v[4], 0)
  expect_identical(dzipig(-1, 2, 1, 0.2), 0)
  expect_warning(expect_identical(dzipig(1.5, 2, 1, 0.2), 0), "non-integer x")
  expect_true(is.na(dzipig(NA, 2, 1, 0.2)))
  expect_equal(dzipig(0, 2, 1, 0.2), 0.2 + 0.8 * exp(pig_bessel_log(0, 2, 1)), tolerance = 1e-12)
})

test_that("qpig/qzipig honour the quantile contract", {
  F <- cumsum(exp(pig_bessel_log(0:150, 100, 10)))
  expect_warning(v <- qpig(c(0.5, 1 - 1e-12), 100, 10, max.value = 200), "not attained")
  expect_equal(v, c(which(F >= 0.5)[1] - 1, NA))                         # was 18, Inf
  F3 <- cumsum(exp(pig_bessel_log(0:50, 3, 1)))
  expect_equal(qpig(c(NA, 0, 1, 0.5), 3, 1), c(NA, 0, Inf, which(F3 >= 0.5)[1] - 1))
  expect_warning(v <- qpig(c(-0.1, 1.2), 3, 1), "NaNs produced")
  expect_true(all(is.nan(v)))
  expect_equal(qpig(1 - ppig(0, 4, 0.7), 4, 0.7, lower.tail = FALSE), 0)  # was 1
  expect_equal(qzipig(c(NA, 0, 1), 3, 1, 0.2), c(NA, 0, Inf))
  expect_warning(expect_true(is.nan(qzipig(1.2, 3, 1, 0.2))), "NaNs produced")
  expect_equal(qzipig(c(0.1, 0.2), 3, 1, 0.2), c(0, 0))
  expect_equal(qzipig(1 - pzipig(0, 4, 0.7, 0.2), 4, 0.7, 0.2, lower.tail = FALSE), 0)
})

## ---- C-6: CompPareto outside its support -------------------------------------------

test_that("CompPareto likelihood and score agree when y <= 0 is supplied", {
  y <- c(0, 0, 1, 2, 5, 9, 3, 14)
  n <- length(y)
  mk <- function(yy, k) list(X = replicate(k, matrix(1, length(yy), 1), simplify = FALSE),
                             y = cbind(yy), dupid = 0L, duplicate = 0L,
                             offsets = replicate(k, numeric(), simplify = FALSE),
                             idpars = seq_len(k), censored = FALSE)
  ## theta = 3.3 is not an observed count, so the FD steps do not cross the splice point
  for (cs in list(list("exp", c(log(0.5), log(1.5), log(3.3))),
                  list("gamma", c(log(1.2), log(2), log(1.5), log(3.3))))) {
    fns <- egpd:::.comppareto_family_setup(list(spec = cs[[1]]))$lik.fns
    p <- cs[[2]]; k <- length(p)
    expect_identical(fns$d0(p, mk(y, k)), Inf)          # density 0 at y = 0 (was finite)
    M <- fns$d120(p, mk(y, k))
    expect_true(all(M[y <= 0, ] == 0))
    yp <- y[y > 0]
    expect_equal(M[y > 0, ], fns$d120(p, mk(yp, k)))
    if (cs[[1]] == "exp") {
      fd <- vapply(seq_len(k), function(i) {
        h <- 1e-6; pp <- p; pm <- p; pp[i] <- pp[i] + h; pm[i] <- pm[i] - h
        (fns$d0(pp, mk(yp, k)) - fns$d0(pm, mk(yp, k))) / (2 * h)
      }, 0)
      expect_equal(colSums(M[, seq_len(k)]), fd, tolerance = 1e-6)
    }
  }
  expect_equal(dcomppareto(c(-1, 0), "exp", rate = 1), c(0, 0))
  expect_equal(dcomppareto(c(-1, 0), "exp", rate = 1, log = TRUE), c(-Inf, -Inf))
  expect_equal(pcomppareto(c(-1, 0), "exp", rate = 1), c(0, 0))
  expect_equal(qcomppareto(0, "exp", rate = 1), 0)
  set.seed(2); expect_true(all(rcomppareto(200, "exp", alpha = 1.5, theta = 2, rate = 1) > 0))
  ## invalid parameters and probabilities give NaN with a warning (were a silent NA)
  expect_warning(expect_true(is.nan(dcomppareto(1, "exp", alpha = -1, rate = 1))), "NaNs produced")
  expect_warning(expect_true(is.nan(qcomppareto(1.1, "exp", rate = 1))), "NaNs produced")
  expect_true(is.na(dcomppareto(1, "exp", alpha = NA_real_, rate = 1)))
})

## ---- C-11: uniform contracts -------------------------------------------------------

test_that("count families share the stats::dpois argument contract", {
  fams <- list(
    list(d = function(x, a) dgpois(x, mu = a, lambda = .3),
         p = function(q, a) pgpois(q, mu = a, lambda = .3),
         q = function(p, a) qgpois(p, mu = a, lambda = .3), r = function(n, a) rgpois(n, mu = a)),
    list(d = function(x, a) dgwaring(x, mu = a), p = function(q, a) pgwaring(q, mu = a),
         q = function(p, a) qgwaring(p, mu = a), r = function(n, a) rgwaring(n, mu = a)),
    list(d = function(x, a) dpoislnorm(x, mu = a), p = function(q, a) ppoislnorm(q, mu = a),
         q = function(p, a) qpoislnorm(p, mu = a), r = function(n, a) rpoislnorm(n, mu = a)),
    list(d = function(x, a) dcmp(x, mu = a), p = function(q, a) pcmp(q, mu = a),
         q = function(p, a) qcmp(p, mu = a), r = function(n, a) rcmp(n, mu = a)),
    list(d = function(x, a) dprl(x, tau = a + 2), p = function(q, a) pprl(q, tau = a + 2),
         q = function(p, a) qprl(p, tau = a + 2), r = function(n, a) rprl(n, tau = a + 2)),
    list(d = function(x, a) dpig(x, mu = a), p = function(q, a) ppig(q, mu = a),
         q = function(p, a) qpig(p, mu = a), r = function(n, a) rpig(n, mu = a)),
    list(d = function(x, a) dzipig(x, mu = a), p = function(q, a) pzipig(q, mu = a),
         q = function(p, a) qzipig(p, mu = a), r = function(n, a) rzipig(n, mu = a)),
    list(d = function(x, a) dbell(x, theta = a), p = function(q, a) pbell(q, theta = a),
         q = function(p, a) qbell(p, theta = a), r = function(n, a) rbell(n, theta = a)),
    list(d = function(x, a) dzibell(x, theta = a), p = function(q, a) pzibell(q, theta = a),
         q = function(p, a) qzibell(p, theta = a), r = function(n, a) rzibell(n, theta = a)),
    list(d = function(x, a) dgpig(x, b = a), p = function(q, a) pgpig(q, b = a),
         q = function(p, a) qgpig(p, b = a), r = function(n, a) rgpig(n, b = a)),
    list(d = function(x, a) dzigpig(x, b = a), p = function(q, a) pzigpig(q, b = a),
         q = function(p, a) qzigpig(p, b = a), r = function(n, a) rzigpig(n, b = a)))
  for (f in fams) {
    ## NA parameter -> NA, silently (dgpois(1, mu = NA) used to stop())
    expect_warning(v <- f$d(1, NA_real_), NA)
    expect_true(is.na(v) && !is.nan(v))
    expect_warning(v <- f$p(1, NA_real_), NA)
    expect_true(is.na(v) && !is.nan(v))
    expect_warning(v <- f$q(0.5, NA_real_), NA)
    expect_true(is.na(v) && !is.nan(v))
    ## invalid parameter -> NaN with a warning
    for (bad in c(-1, NaN, Inf)) {
      expect_warning(expect_true(is.nan(f$d(1, bad))), "NaNs produced")
      expect_warning(expect_true(is.nan(f$p(1, bad))), "NaNs produced")
      expect_warning(expect_true(is.nan(f$q(0.5, bad))), "NaNs produced")
      expect_warning(expect_true(all(is.nan(f$r(2, bad)))), "NaNs produced")
    }
    ## x conventions of stats::dpois
    expect_warning(expect_identical(f$d(2.5, 2), 0), "non-integer x")  # was silent
    expect_identical(f$d(c(-1, Inf), 2), c(0, 0))
    v <- f$d(c(NA, NaN), 2)
    expect_true(is.na(v[1]) && !is.nan(v[1]) && is.nan(v[2]))
    ## probabilities outside [0, 1] -> NaN with a warning
    expect_warning(expect_true(all(is.nan(f$q(c(-0.1, 1.1), 2)))), "NaNs produced")
    expect_equal(f$q(c(0, 1), 2), c(0, Inf))
  }
})

## ---- suspected issues --------------------------------------------------------------

test_that("generalized Waring GAM likelihood is free of lgamma cancellation", {
  X <- matrix(1, 1, 1)
  ## Rmpfr (256-bit) log p(5; mu = 3, k = 0.01, rho = 5e5) = -6.2675556321445223655
  ref <- -6.2675556321445223655
  gam <- -egpd:::gwaring1d0(list(log(3), log(0.01), log(5e5)), X, X, X, 5, integer(), 0L,
                            list(numeric(), numeric(), numeric()))
  expect_equal(gam, ref, tolerance = 1e-13)                               # was 1.4e-6 off
  expect_equal(dgwaring(5, 3, 0.01, 5e5, log = TRUE), ref, tolerance = 1e-13)
})

## ---- cross-area (core-fit): generalized Waring at its negative-binomial / Poisson limits

test_that("generalized Waring log pmf and score are exact across its limits, in dgwaring and the GAM", {
  ## Rmpfr (800-bit) log p and d log p / d(log mu, log k, log rho) from the lgamma/digamma
  ## definition with a = mu (rho - 1)/k. Columns: y, mu, k, rho, log p, g1, g2, g3.
  ref <- rbind(
    ## negative-binomial limit beyond the former clamp rho <= 1e6, and far beyond it
    c(7, 30, 2, 1e7, -3.9175053920127305, -1.437499876328119, -0.67196290118717772, -1.5874998883836859e-07),
    c(60, 30, 2, 1e7, -5.3066149548105281, 1.8749996433594798, -0.027649657939602147, 1.0625001298350159e-07),
    c(0, 30, 2, 1e14, -5.5451774444795348, -1.8749999999999982, -3.6701774444795174, -2.8125000000000063e-14),
    c(60, 30, 2, 1e14, -5.3066148485605318, 1.8749999999999643, -0.027649733330177662, 1.062500000000013e-14),
    ## Poisson limit, rho >> k >> 1
    c(30, 30, 1e8, 1e15, -2.6223150989654784, 5.0000010055550096e-08, 9.9999960444457175e-08, 4.9999995055556149e-08),
    c(45, 30, 1e8, 1e15, -6.0700502643307246, 14.999993700002042, -5.9999938933371582e-07, -2.9999980366674104e-07),
    ## a = 1e12 inside the old box (raw digamma differences lost the score)
    c(3, 1e4, 1e-3, 1e5, -8.0209865830244169, -0.0009996998970264806, 0.98638036028004494, -5.0020798243297199e-09),
    c(12000, 1e4, 1e-3, 1e5, -16.307766296541633, 0.00020001197771913549, 0.99365008954603551, 6.9999564712861216e-09),
    ## a = 2.9e-13, below the former internal clamp a >= 1e-12
    c(0, 1e-6, 1e8, 30, -4.3605104428173692e-12, -4.3605104428173684e-12, 4.0705105283673425e-12, -4.215985933767674e-12),
    c(2, 1e-6, 1e8, 30, -29.562043252494089, 0.99999999999592948, -0.99999939999640552, 1.0344821586169566),
    ## far power-law tail (lgamma(y) ~ 1e8 no longer enters)
    c(1e7, 2, 1, 2.5, -52.232378746920858, 2.064925390450298, -0.38455363490334338, -31.825965438523809),
    c(123456, 0.5, 3, 1.7, -31.672742264385452, 1.0813500081210854, 0.46252427141504293, -14.709779160646852))
  for (i in seq_len(nrow(ref))) {
    r <- ref[i, ]
    expect_equal(dgwaring(r[1], r[2], r[3], r[4], log = TRUE), r[5], tolerance = 5e-13)
    gam <- unname(egpd:::cf_grad_cpp("gwaring", r[1], log(r[2:4])))
    expect_equal(gam[1], r[5], tolerance = 5e-13)
    ## the GAM's parameters are exp(log(.)), within an ulp of the reference ones, so the score
    ## is compared to relative accuracy above an absolute floor of that conditioning
    expect_true(all(abs(gam[2:4] - r[6:8]) <= 1e-9 * abs(r[6:8]) + 1e-14 * (1 + r[1])),
                label = sprintf("score at row %d", i))
  }
})

test_that("generalized Waring probabilities converge to the negative binomial, and the clamp sits there", {
  ## Rmpfr (400-bit) tail sums of the pmf (recursion p(y+1)/p(y) = (a+y)(k+y)/((y+1)(c+y))).
  expect_equal(pgwaring(60, 30, 2, 1e14, lower.tail = FALSE, log.p = TRUE), -2.3656330897779236,
               tolerance = 1e-13)
  expect_equal(pgwaring(1000, 30, 2, 1e14, lower.tail = FALSE, log.p = TRUE), -60.451035984880342,
               tolerance = 1e-13)                                         # was 1.5% off
  expect_equal(pgwaring(7, 30, 2, 1e14, log.p = TRUE), -2.2545496016120739, tolerance = 1e-13)
  expect_equal(pgwaring(60, 30, 1e8, 1e15, lower.tail = FALSE, log.p = TRUE), -14.617430999544336,
               tolerance = 1e-13)
  expect_equal(pgwaring(300, 0.01, 50, 1e12, lower.tail = FALSE, log.p = TRUE), -2424.7988423579709,
               tolerance = 1e-13)
  ## so far out that R's pbeta(log.p = TRUE) underflows (stats::pnbinom gives -Inf with
  ## warnings); P(Y > q) was returned as 1
  expect_silent(lq <- pgwaring(15354, 765.2, 35.77, 2.041e7, lower.tail = FALSE, log.p = TRUE))
  expect_equal(lq, -565.6625713774539, tolerance = 1e-13)
  ## GW(mu, k, rho) -> NB(size = k, mu) as rho -> Inf, with differences O(1/rho)
  y <- c(0, 1, 7, 30, 60, 200, 1000)
  expect_equal(dgwaring(y, 30, 2, 1e14, log = TRUE), stats::dnbinom(y, size = 2, mu = 30, log = TRUE),
               tolerance = 1e-11)
  expect_equal(pgwaring(y, 30, 2, 1e14, lower.tail = FALSE, log.p = TRUE),
               stats::pnbinom(y, size = 2, mu = 30, lower.tail = FALSE, log.p = TRUE), tolerance = 1e-10)
  ## Beyond the bound the GAM likelihood is that AT the bound (as predict() reports), which
  ## is already the negative-binomial limit to double precision: the flat region is the plateau.
  X <- matrix(1, 3, 1); off <- list(numeric(), numeric(), numeric())
  yy <- c(0, 7, 60)
  nll <- function(lrho) egpd:::gwaring1d0(list(log(30), log(2), lrho), X, X, X, yy, integer(), 0L, off)
  rho_hi <- egpd:::cf_bounds_cpp()$gwaring[["rho_hi"]]
  expect_identical(nll(log(1e30)), nll(log(rho_hi) + 1))
  expect_equal(nll(log(1e30)), -sum(dgwaring(yy, exp(log(30)), exp(log(2)), rho_hi, log = TRUE)),
               tolerance = 1e-14)
  expect_equal(nll(log(1e30)), -sum(stats::dnbinom(yy, size = 2, mu = 30, log = TRUE)), tolerance = 1e-12)
})

test_that("generalized Waring fits on Poisson data reach the Poisson-limit optimum", {
  ## On Poisson(30) data the gwaring supremum lies at k, rho -> Inf, where the old log pmf
  ## lost all precision (a positive log-likelihood, +789801 at the point below) and the old
  ## clamps held the fit short. A quasi-Newton run on the package's own objective and
  ## score from egpd()'s default start must now reach at least the Poisson log-likelihood.
  set.seed(11)
  y <- stats::rpois(400, 30)
  X <- matrix(1, length(y), 1); off <- list(numeric(), numeric(), numeric())
  nll <- function(e) egpd:::gwaring1d0(as.list(e), X, X, X, y, integer(), 0L, off)
  gr <- function(e) colSums(egpd:::gwaring1d12(as.list(e), X, X, X, y, integer(), 0L, off)[, 1:3])
  pois <- -sum(stats::dpois(y, mean(y), log = TRUE))
  o <- stats::optim(c(log(mean(y)), 0, log(2.5)), nll, gr, method = "BFGS",
                    control = list(maxit = 1000, reltol = 1e-14))
  expect_equal(o$convergence, 0L)
  expect_lte(o$value, pois)
  expect_gt(o$value, pois - 0.1)
  expect_lt(max(abs(gr(o$par))), 1e-3)
  ## the point where the old code returned a positive log-likelihood: clamped mu and rho
  bd <- egpd:::cf_bounds_cpp()$gwaring
  far <- nll(c(log(3.3e23), log(4.7), log(7e21)))
  expect_gt(far, pois)
  expect_equal(far, -sum(dgwaring(y, bd[["mu_hi"]], exp(log(4.7)), bd[["rho_hi"]], log = TRUE)),
               tolerance = 1e-13)
})

test_that("zero-inflated inverse-quantile helpers use pi", {
  ns <- asNamespace("egpd")
  u <- c(0.2, 0.5, 0.9)
  expect_equal(ns$.zibell1fns$iG(u, mu = 6, pi = 0.5), qZIBELL(u, mu = 6, sigma = 0.5))
  expect_equal(ns$.zibellnat1fns$iG(u, theta = 1.2, pi = 0.5), qzibell(u, 1.2, pi = 0.5))
  expect_equal(ns$.zipig1fns$iG(u, mu = 3, sigma = 0.7, pi = 0.6), qzipig(u, 3, 0.7, 0.6))
  expect_equal(ns$.zigpig1fns$iG(u, mu = 3, a = 0.5, c = 0.6, pi = 0.6),
               qZIGPIG(u, mu = 3, sigma = 0.5, nu = 0.6, tau = 0.6))
  expect_equal(ns$.zigpignat1fns$iG(u, a = 0.5, b = 2, c = 0.6, pi = 0.6),
               qzigpig(u, 0.5, 2, 0.6, 0.6))
  expect_equal(ns$.zibell1fns$iG(0.45, mu = 6, pi = 0.5), 0)               # the zero mass
})
