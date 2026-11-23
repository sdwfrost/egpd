## Round 2 audit regression tests:
##   1. ZIBell / ZIGPIG parameter validation & empty-input quantiles
##   2. Vectorization, unclamped standalone d/p/q/r, PRL tail precision, RNG n conventions
##   3. C++ zero-cell logsumexp stability for zibell, zipig, and zigpig
##   4. Composite Pareto NA/NaN propagation, tail precision, and qcomppareto(lower.tail)
##   5. Direct coverage of previously untested exports (qcmp, rcmp, BELL/ZIBELL/GPIG/ZIGPIG d/p/q/r)

test_that("ZIGPIG and ZIBell validate pi and base parameters across d/p/q/r", {
  for (bad_pi in c(-0.5, -2, 1, 1.2, Inf, -Inf, NaN)) {
    expect_warning(v <- dzigpig(0:2, pi = bad_pi), "NaNs produced")
    expect_true(all(is.nan(v)))
    expect_warning(v <- pzigpig(c(-1, 0, 2), pi = bad_pi), "NaNs produced")
    expect_true(all(is.nan(v)))
    expect_warning(v <- qzigpig(c(0, 0.5, 1), pi = bad_pi), "NaNs produced")
    expect_true(all(is.nan(v)))
    expect_warning(v <- rzigpig(4, pi = bad_pi), "NaNs produced")
    expect_true(all(is.nan(v)))

    expect_warning(v <- dzibell(0:2, pi = bad_pi), "NaNs produced")
    expect_true(all(is.nan(v)))
    expect_warning(v <- pzibell(c(-1, 0, 2), pi = bad_pi), "NaNs produced")
    expect_true(all(is.nan(v)))
    expect_warning(v <- qzibell(c(0, 0.5, 1), pi = bad_pi), "NaNs produced")
    expect_true(all(is.nan(v)))
    expect_warning(v <- rzibell(4, pi = bad_pi), "NaNs produced")
    expect_true(all(is.nan(v)))
  }

  ## NA (not NaN) in parameters propagates NA_real_ without warning in d/p/q
  expect_identical(dzigpig(1, pi = NA_real_), NA_real_)
  expect_identical(pzigpig(1, pi = NA_real_), NA_real_)
  expect_identical(qzigpig(0.5, pi = NA_real_), NA_real_)
  expect_identical(dzibell(1, pi = NA_real_), NA_real_)
  expect_identical(pzibell(1, pi = NA_real_), NA_real_)
  expect_identical(qzibell(0.5, pi = NA_real_), NA_real_)

  ## Invalid base parameters at q < 0 must not be overwritten with 0
  expect_warning(v <- pzigpig(-1, a = 2, pi = 0.2), "NaNs produced")
  expect_true(is.nan(v))
  expect_warning(v <- pzibell(-1, theta = -1, pi = 0.2), "NaNs produced")
  expect_true(is.nan(v))

  ## rzigpig / rzibell validate base parameters even when all uniform draws < pi
  set.seed(1)
  expect_warning(rz1 <- rzigpig(5, a = 2, b = 1, c = 0.5, pi = 0.999999), "NaNs produced")
  expect_true(all(is.nan(rz1)))
  expect_warning(rz2 <- rzigpig(5, a = 0.5, b = -1, c = 0.5, pi = 0.999999), "NaNs produced")
  expect_true(all(is.nan(rz2)))
  expect_warning(rz3 <- rzigpig(5, a = 0.5, b = 1, c = -0.1, pi = 0.999999), "NaNs produced")
  expect_true(all(is.nan(rz3)))
  expect_warning(rb1 <- rzibell(5, theta = -1, pi = 0.999999), "NaNs produced")
  expect_true(all(is.nan(rb1)))
})

test_that("Bell and GPIG quantile functions return numeric(0) on empty inputs", {
  qfns <- list(qgpig, qzigpig, qbell, qzibell, qGPIG, qZIGPIG, qBELL, qZIBELL)
  for (qfn in qfns) {
    res <- qfn(numeric(0))
    expect_identical(res, numeric(0))
  }
})

test_that("CMP, PRL, gpois, gwaring, and poislnorm vectorize over parameters and handle NAs", {
  ## Heterogeneous parameter vectors: compare against scalar elementwise evaluation
  x <- c(0, 2, 5, NA_real_, NaN, -1, 1.5, Inf, -Inf)
  q <- c(-Inf, -1, 0, 2, 6, Inf, NA_real_, NaN)
  p <- c(0, 0.25, 0.75, 1, NA_real_, NaN, -0.1, 1.1)

  ## CMP
  mu_c <- c(0.8, 2.5, 6.0); nu_c <- c(0.6, 1.0, 1.8)
  expect_equal(dcmp(c(0, 2, 5), mu = mu_c, nu = nu_c),
               mapply(dcmp, c(0, 2, 5), mu_c, nu_c))
  expect_equal(pcmp(c(0, 2, 5), mu = mu_c, nu = nu_c),
               mapply(pcmp, c(0, 2, 5), mu_c, nu_c))
  expect_equal(qcmp(c(0.2, 0.5, 0.8), mu = mu_c, nu = nu_c),
               mapply(qcmp, c(0.2, 0.5, 0.8), mu_c, nu_c))

  ## PRL
  tau_v <- c(2.0, 3.5, 8.0)
  expect_equal(dprl(c(0, 3, 9), tau = tau_v),
               mapply(dprl, c(0, 3, 9), tau_v))
  expect_equal(pprl(c(0, 3, 9), tau = tau_v),
               mapply(pprl, c(0, 3, 9), tau_v))
  expect_equal(qprl(c(0.2, 0.5, 0.8), tau = tau_v),
               mapply(qprl, c(0.2, 0.5, 0.8), tau_v))

  ## Generalized Poisson
  mu_gp <- c(1.2, 3.0, 7.5); lam_gp <- c(0, 0.25, 0.6)
  expect_equal(dgpois(c(0, 2, 8), mu = mu_gp, lambda = lam_gp),
               mapply(dgpois, c(0, 2, 8), mu_gp, lam_gp))
  expect_equal(pgpois(c(0, 2, 8), mu = mu_gp, lambda = lam_gp),
               mapply(pgpois, c(0, 2, 8), mu_gp, lam_gp))
  expect_equal(qgpois(c(0.2, 0.5, 0.8), mu = mu_gp, lambda = lam_gp),
               mapply(qgpois, c(0.2, 0.5, 0.8), mu_gp, lam_gp))

  ## Generalized Waring
  mu_gw <- c(1.5, 4.0, 9.0); k_gw <- c(0.7, 1.5, 3.0); rho_gw <- c(1.8, 2.5, 5.0)
  expect_equal(dgwaring(c(0, 2, 7), mu = mu_gw, k = k_gw, rho = rho_gw),
               mapply(dgwaring, c(0, 2, 7), mu_gw, k_gw, rho_gw))
  expect_equal(pgwaring(c(0, 2, 7), mu = mu_gw, k = k_gw, rho = rho_gw),
               mapply(pgwaring, c(0, 2, 7), mu_gw, k_gw, rho_gw))
  expect_equal(qgwaring(c(0.2, 0.5, 0.8), mu = mu_gw, k = k_gw, rho = rho_gw),
               mapply(qgwaring, c(0.2, 0.5, 0.8), mu_gw, k_gw, rho_gw))

  ## Poisson-lognormal
  mu_pl <- c(0.9, 3.5, 10.0); sg_pl <- c(0.3, 0.8, 1.4)
  expect_equal(dpoislnorm(c(0, 2, 9), mu = mu_pl, sigma = sg_pl),
               mapply(dpoislnorm, c(0, 2, 9), mu_pl, sg_pl))
  expect_equal(ppoislnorm(c(0, 2, 9), mu = mu_pl, sigma = sg_pl),
               mapply(ppoislnorm, c(0, 2, 9), mu_pl, sg_pl))
  expect_equal(qpoislnorm(c(0.2, 0.5, 0.8), mu = mu_pl, sigma = sg_pl),
               mapply(qpoislnorm, c(0.2, 0.5, 0.8), mu_pl, sg_pl))

  ## NA / NaN / off-lattice / infinite handling and numeric(0)
  for (dfun in list(dcmp, dprl, dgpois, dgwaring, dpoislnorm)) {
    expect_warning(v <- dfun(x), "non-integer x")      # x = 1.5, as stats::dpois
    expect_true(is.na(v[4]) && !is.nan(v[4]))
    expect_true(is.nan(v[5]))
    expect_equal(v[6:9], c(0, 0, 0, 0))
    expect_identical(dfun(numeric(0)), numeric(0))
  }
  for (pfun in list(pcmp, pprl, pgpois, pgwaring, ppoislnorm)) {
    v <- pfun(q)
    expect_equal(v[c(1, 2, 6)], c(0, 0, 1))
    expect_true(is.na(v[7]) && !is.nan(v[7]))
    expect_true(is.nan(v[8]))
    expect_identical(pfun(numeric(0)), numeric(0))
  }
  for (qfun in list(qcmp, qprl, qgpois, qgwaring, qpoislnorm)) {
    expect_warning(v <- qfun(p), "NaNs produced")      # p = -0.1, 1.1
    expect_equal(v[c(1, 4)], c(0, Inf))
    expect_true(is.na(v[5]) && !is.nan(v[5]))
    expect_true(all(is.nan(v[6:8])))
    expect_identical(qfun(numeric(0)), numeric(0))
  }
})

test_that("standalone dgpois, dgwaring, dpoislnorm, and dcmp are unclamped", {
  ## Generalized Waring near rho -> 1+: a = mu*(rho - 1)/k, P(Y > 0) = a / (a + rho)
  rho <- 1 + 1e-9
  a_ex <- rho - 1
  expect_equal(1 - dgwaring(0, mu = 1, k = 1, rho = rho),
               a_ex / (a_ex + rho), tolerance = 1e-14)
  ## Exact y >= 1 log-PMF at k = 1 (Waring): p(y) = rho * prod_{j=0}^{y-1}(a + j) / prod_{j=0}^y (a + rho + j)
  for (yv in 1:5) {
    poch_a <- log(a_ex) + if (yv == 1L) 0 else lgamma(yv) + sum(log1p(a_ex / seq_len(yv - 1L)))
    ref_lp <- log(rho) + poch_a - sum(log(a_ex + rho + 0:yv))
    expect_equal(dgwaring(yv, mu = 1, k = 1, rho = rho, log = TRUE),
                 ref_lp, tolerance = 1e-13)
  }
  ## Extreme large rho -> NB limit
  expect_equal(dgwaring(0:6, mu = 2, k = 2, rho = 1e8, log = TRUE),
               stats::dnbinom(0:6, size = 2, mu = 2, log = TRUE), tolerance = 1e-6)
  ## Extreme small mu (beyond the GAM's GW_MU_LO = 1e-10), a = 7.5e-16
  mu_sm <- 1e-15; k_sm <- 2; rho_sm <- 2.5; a_sm <- mu_sm * (rho_sm - 1) / k_sm
  ## At y = 0: 1 - P(Y = 0) = a * (digamma(rho + k) - digamma(rho)) + O(a^2)
  expect_equal(1 - dgwaring(0, mu = mu_sm, k = k_sm, rho = rho_sm),
               a_sm * (digamma(rho_sm + k_sm) - digamma(rho_sm)), tolerance = 1e-12)
  for (yv in 1:4) {
    poch_a <- log(a_sm) + if (yv == 1L) 0 else lgamma(yv) + sum(log1p(a_sm / seq_len(yv - 1L)))
    ref_sm <- (a_sm * digamma(rho_sm)) + lgamma(k_sm + rho_sm) + poch_a +
      lgamma(k_sm + yv) - lgamma(k_sm) - lgamma(a_sm + k_sm + rho_sm + yv) - lgamma(yv + 1)
    expect_equal(dgwaring(yv, mu = mu_sm, k = k_sm, rho = rho_sm, log = TRUE),
                 ref_sm, tolerance = 1e-12)
  }

  ## Generalized Poisson at lambda = 0 is exact Poisson, and tiny lambda/mu are unclamped
  expect_equal(dgpois(0:12, mu = 3.7, lambda = 0, log = TRUE),
               stats::dpois(0:12, 3.7, log = TRUE), tolerance = 1e-15)
  expect_equal(dgpois(1, mu = 1e-20, lambda = 0, log = TRUE),
               stats::dpois(1, 1e-20, log = TRUE), tolerance = 1e-14)
  lam_tiny <- 1e-13
  ## Exact Consul-Jain at y = 0: log P(Y = 0) = -mu * (1 - lambda)
  expect_equal(dgpois(0, mu = 2, lambda = lam_tiny, log = TRUE) - (-2),
               2 * lam_tiny, tolerance = 1e-4)

  ## Poisson-lognormal below PLN_SG_LO = 1e-6 and PLN_MU_LO = 1e-10
  for (sg in c(1e-8, 1e-12, 1e-20)) {
    expect_equal(dpoislnorm(0:6, mu = 2.5, sigma = sg, log = TRUE),
                 stats::dpois(0:6, 2.5, log = TRUE) +
                   0.5 * sg^2 * (((0:6) - 2.5)^2 - (0:6)),
                 tolerance = 1e-11)
  }
  expect_equal(dpoislnorm(1, mu = 1e-20, sigma = 1e-10, log = TRUE),
               stats::dpois(1, 1e-20, log = TRUE) - 0.5 * 1e-20, tolerance = 1e-12)

  ## CMP below CMP_NU_LO = 1e-6 and above CMP_NU_HI = 1e3
  expect_equal(dcmp(0:5, mu = 2, nu = 1e-9),
               stats::dgeom(0:5, prob = 1 / 3), tolerance = 1e-6)
  p_ber <- dcmp(0:5, mu = 0.4, nu = 2000)
  expect_equal(p_ber[1:2], c(0.6, 0.4), tolerance = 1e-12)
})

test_that("PRL validates tau >= 2 across d/p/q/r and preserves deep upper-tail precision", {
  ## Round 3 (C-11): invalid parameters give NaN with a warning and NA gives NA, as for
  ## every other count family and stats::dpois, rather than a stopifnot() error.
  for (bad_tau in c(1.999, 1.2, 0, -2, Inf, NaN)) {
    expect_warning(expect_true(is.nan(dprl(0, tau = bad_tau))), "NaNs produced")
    expect_warning(expect_true(is.nan(pprl(0, tau = bad_tau))), "NaNs produced")
    expect_warning(expect_true(is.nan(qprl(0.5, tau = bad_tau))), "NaNs produced")
    expect_warning(expect_true(all(is.nan(rprl(3, tau = bad_tau)))), "NaNs produced")
  }
  expect_identical(dprl(0, tau = NA_real_), NA_real_)
  expect_identical(pprl(0, tau = NA_real_), NA_real_)
  expect_identical(qprl(0.5, tau = NA_real_), NA_real_)

  ## Exact closed-form log-survivor at q = 500, tau = 2.5 (previously underflowed to -Inf)
  x <- 500; tau <- 2.5
  ref_ls <- -x * log1p(1 / tau) + log(tau) + log(x + tau^2) - log(tau - 1) - 2 * log1p(tau)
  expect_equal(pprl(500, 2.5, lower.tail = FALSE, log.p = TRUE), ref_ls, tolerance = 1e-13)
  expect_equal(pprl(500, 2.5, lower.tail = FALSE, log.p = FALSE), exp(ref_ls), tolerance = 1e-13)
  expect_equal(pprl(500.7, 2.5, lower.tail = FALSE, log.p = TRUE), ref_ls, tolerance = 1e-13)
  expect_true(is.finite(pprl(5000, 2.5, lower.tail = FALSE, log.p = TRUE)))
})

test_that("all count and CompPareto RNGs standardize length(n) > 1 and n == 0", {
  rfns <- list(
    function(n) rcmp(n, mu = 2, nu = 1),
    function(n) rprl(n, tau = 3),
    function(n) rgpois(n, mu = 2, lambda = 0.2),
    function(n) rgwaring(n, mu = 2, k = 2, rho = 3),
    function(n) rpoislnorm(n, mu = 2, sigma = 0.5),
    function(n) rpig(n, mu = 2, sigma = 0.5),
    function(n) rpig(n, mu = 2, sigma = 1e-12),
    function(n) rzipig(n, mu = 2, sigma = 0.5, pi = 0.2),
    function(n) rzipig(n, mu = 2, sigma = 1e-12, pi = 0.2),
    function(n) rgpig(n, a = 0.5, b = 2, c = 0.6),
    function(n) rzigpig(n, a = 0.5, b = 2, c = 0.6, pi = 0.2),
    function(n) rbell(n, theta = 0.8),
    function(n) rzibell(n, theta = 0.8, pi = 0.2),
    function(n) rBELL(n, mu = 2),
    function(n) rZIBELL(n, mu = 2, sigma = 0.2),
    function(n) rGPIG(n, mu = 2, sigma = 0.5, nu = 0.6),
    function(n) rZIGPIG(n, mu = 2, sigma = 0.5, nu = 0.6, tau = 0.2),
    function(n) rcomppareto(n, "exp", alpha = 1.5, theta = 2, rate = 1)
  )
  for (rfn in rfns) {
    z0 <- rfn(0)
    expect_length(z0, 0L)
    expect_true(is.numeric(z0))
    zv <- rfn(c(10, 20, 30))
    expect_length(zv, 3L)
    expect_true(all(is.finite(zv)))
    expect_error(rfn(-1))
    expect_error(rfn(NA_real_))
  }
})

test_that("C++ zibell, zipig, and zigpig survive zero-cell underflow in d0 and d12", {
  eval_d0 <- function(fn, e, y) {
    np <- length(e)
    do.call(fn, c(list(as.list(e)),
                  rep(list(matrix(1, 1, 1)), np),
                  list(y, as.integer(0), 0L, as.list(rep(list(numeric(0)), np)))))
  }
  eval_d12 <- function(fn, e, y) {
    np <- length(e)
    as.numeric(do.call(fn, c(list(as.list(e)),
                             rep(list(matrix(1, 1, 1)), np),
                             list(y, as.integer(0), 0L, as.list(rep(list(numeric(0)), np))))))
  }
  fd_g_h <- function(d0fn, e, y) {
    h <- 1e-5
    np <- length(e)
    g <- vapply(seq_len(np), function(i) {
      ep <- em <- e; ep[i] <- ep[i] + h; em[i] <- em[i] - h
      (eval_d0(d0fn, ep, y) - eval_d0(d0fn, em, y)) / (2 * h)
    }, numeric(1))
    H <- numeric(0)
    for (i in seq_len(np)) for (j in i:np) {
      epp <- epm <- emp <- emm <- e
      epp[i] <- epp[i] + h; epp[j] <- epp[j] + h
      epm[i] <- epm[i] + h; epm[j] <- epm[j] - h
      emp[i] <- emp[i] - h; emp[j] <- emp[j] + h
      emm[i] <- emm[i] - h; emm[j] <- emm[j] - h
      H <- c(H, (eval_d0(d0fn, epp, y) - eval_d0(d0fn, epm, y) -
                   eval_d0(d0fn, emp, y) + eval_d0(d0fn, emm, y)) / (4 * h * h))
    }
    list(g = g, H = H)
  }

  ## ZIBell native: theta = exp(3) = 20.0855 => lp0 = 1 - exp(theta) = -5.28e8 (p0 == 0 in double)
  ## With eta_pi = -1000 (pi == 0 in double), w1 = 1, NLL = 1000, g = c(0, -1), H = c(0, 0, 0).
  e_bell_zi <- c(3, -1000)
  expect_equal(eval_d0(egpd:::zibellnat1d0, e_bell_zi, 0), 1000, tolerance = 1e-12)
  expect_equal(eval_d12(egpd:::zibellnat1d12, e_bell_zi, 0), c(0, -1, 0, 0, 0), tolerance = 1e-12)

  ## Tie regime where BOTH pi and p0 underflow to 0 in IEEE double:
  ## theta = 7 => lp0 = 1 - exp(7) = -1095.6331584; set eta_pi = lp0 => w0 = w1 = 0.5!
  lp0_bell <- 1 - exp(7)
  e_bell_tie <- c(log(7), lp0_bell)
  expect_equal(eval_d0(egpd:::zibellnat1d0, e_bell_tie, 0),
               -(lp0_bell + log(2)), tolerance = 1e-11)
  d12_tie <- eval_d12(egpd:::zibellnat1d12, e_bell_tie, 0)
  expect_true(all(is.finite(d12_tie)))
  fd_tie <- fd_g_h(egpd:::zibellnat1d0, e_bell_tie, 0)
  expect_equal(d12_tie[1:2], fd_tie$g, tolerance = 1e-4)
  expect_equal(d12_tie[3:5], fd_tie$H, tolerance = 1e-2)

  ## Also check zibell mean parameterisation and large positive eta_pi at y = 0 and y = 3
  for (eta_pi in c(-1000, 100)) {
    for (yv in c(0, 3)) {
      eb <- c(10, eta_pi)
      expect_true(is.finite(eval_d0(egpd:::zibell1d0, eb, yv)))
      expect_true(all(is.finite(eval_d12(egpd:::zibell1d12, eb, yv))))
    }
  }

  ## ZIPIG: lmu = 16, lsigma = 0 => mu = 8.886e6, b = 4215.71, lp0 = -4214.71 (p0 == 0 in double)
  mu_p <- exp(16); sg_p <- 1; b_p <- sqrt(1 + 2 * sg_p * mu_p)
  lp0_pig <- -2 * mu_p / (1 + b_p)
  ## (a) eta_pi = -1000 > lp0_pig => w1 = 1
  e_pig_w1 <- c(16, 0, -1000)
  expect_equal(eval_d0(egpd:::zipig1d0, e_pig_w1, 0), 1000, tolerance = 1e-12)
  expect_equal(eval_d12(egpd:::zipig1d12, e_pig_w1, 0),
               c(0, 0, -1, 0, 0, 0, 0, 0, 0), tolerance = 1e-12)
  ## (b) Tie point eta_pi == lp0_pig => w0 = w1 = 0.5, exact observed Hessian matches FD
  e_pig_tie <- c(16, 0, lp0_pig)
  expect_equal(eval_d0(egpd:::zipig1d0, e_pig_tie, 0),
               -(lp0_pig + log(2)), tolerance = 1e-11)
  d12_pig_tie <- eval_d12(egpd:::zipig1d12, e_pig_tie, 0)
  expect_true(all(is.finite(d12_pig_tie)))
  fd_pig_tie <- fd_g_h(egpd:::zipig1d0, e_pig_tie, 0)
  expect_equal(d12_pig_tie[1:3], fd_pig_tie$g, tolerance = 1e-4)
  expect_equal(d12_pig_tie[4:9], fd_pig_tie$H, tolerance = 1e-2)
  ## (c) eta_pi = -10000 < lp0_pig => w0 = 1, matches pure pig1d0 / pig1d12
  e_pig_w0 <- c(16, 0, -10000)
  expect_equal(eval_d0(egpd:::zipig1d0, e_pig_w0, 0),
               eval_d0(egpd:::pig1d0, c(16, 0), 0), tolerance = 1e-12)
  p12_ref <- eval_d12(egpd:::pig1d12, c(16, 0), 0)
  expect_equal(eval_d12(egpd:::zipig1d12, e_pig_w0, 0),
               c(p12_ref[1], p12_ref[2], 0, p12_ref[3], p12_ref[4], 0, p12_ref[5], 0, 0),
               tolerance = 1e-12)

  ## ZIGPIG: large b where p0 underflows and eta_pi = -1000
  e_gpig_zi <- c(12, 0, 0, -1000)
  expect_equal(eval_d0(egpd:::zigpig1d0, e_gpig_zi, 0), 1000, tolerance = 1e-12)
  expect_true(all(is.finite(eval_d12(egpd:::zigpig1d12, e_gpig_zi, 0))))
})

test_that("Composite Pareto preserves NA/NaN/invalid, avoids upper-tail cancellation, and supports qcomppareto lower.tail", {
  ## NA vs NaN on x, q, p when parameters are valid
  dx <- dcomppareto(c(NA_real_, NaN, -1, 0, Inf, -Inf), "exp", rate = 1)
  expect_true(is.na(dx[1]) && !is.nan(dx[1]))
  expect_true(is.nan(dx[2]))
  expect_equal(dx[3:6], c(0, 0, 0, 0))

  dx_log <- dcomppareto(c(NA_real_, NaN, -1, 0, Inf), "exp", rate = 1, log = TRUE)
  expect_true(is.na(dx_log[1]) && !is.nan(dx_log[1]))
  expect_true(is.nan(dx_log[2]))
  expect_equal(dx_log[3:5], c(-Inf, -Inf, -Inf))

  pq <- pcomppareto(c(NA_real_, NaN, -1, 0, Inf, -Inf), "exp", rate = 1)
  expect_true(is.na(pq[1]) && !is.nan(pq[1]))
  expect_true(is.nan(pq[2]))
  expect_equal(pq[3:6], c(0, 0, 1, 0))

  qp <- qcomppareto(c(NA_real_, NaN, 0, 1), "exp", rate = 1)
  expect_true(is.na(qp[1]) && !is.nan(qp[1]))
  expect_true(is.nan(qp[2]))
  expect_equal(qp[3:4], c(0, Inf))

  ## Invalid parameters MUST return NA (round 3, C-11: NaN with a warning) even when
  ## x <= 0 or q <= 0
  nan_warn <- function(v) expect_warning(expect_true(all(is.nan(v))), "NaNs produced")
  nan_warn(dcomppareto(c(-2, 0, 1, Inf), "exp", alpha = -1, theta = 1, rate = 1))
  nan_warn(dcomppareto(c(-2, 0, 1), "exp", alpha = 1, theta = 0, rate = 1))
  nan_warn(dcomppareto(c(-2, 0, 1), "exp", alpha = 1, theta = 1, rate = -1))
  nan_warn(pcomppareto(c(-2, 0, 1, Inf, -Inf), "exp", alpha = -1, theta = 1, rate = 1))
  nan_warn(pcomppareto(c(-2, 0, 1), "exp", alpha = 1, theta = -1, rate = 1))
  nan_warn(qcomppareto(c(0, 0.5, 1), "exp", alpha = -1, theta = 1, rate = 1))

  ## Extreme upper-tail precision in pcomppareto and qcomppareto
  phi <- 2 * stats::dexp(1, 1) / stats::pexp(1, 1)
  ls_ref <- log(phi) - log1p(phi) + log(2) - log(1e18) - log1p(1e-18)
  ls_got <- pcomppareto(1e18, "exp", alpha = 1, theta = 1, rate = 1,
                        lower.tail = FALSE, log.p = TRUE)
  expect_true(is.finite(ls_got))
  expect_equal(ls_got, ls_ref, tolerance = 1e-14)
  expect_equal(pcomppareto(1e18, "exp", alpha = 1, theta = 1, rate = 1,
                           lower.tail = FALSE, log.p = FALSE),
               exp(ls_ref), tolerance = 1e-14)

  ## Round-trip qcomppareto <-> pcomppareto across body, splice point, and deep Pareto tail
  qs <- c(0.1, 1.0, 5.0, 1e8, 1e18, 1e80)
  for (lt in c(TRUE, FALSE)) {
    for (lp in c(FALSE, TRUE)) {
      qs_sub <- if (lt) qs[1:3] else qs
      pr <- pcomppareto(qs_sub, "exp", alpha = 1.3, theta = 1.0, rate = 0.9,
                        lower.tail = lt, log.p = lp)
      qr <- qcomppareto(pr, "exp", alpha = 1.3, theta = 1.0, rate = 0.9,
                        lower.tail = lt, log.p = lp)
      expect_equal(qr, qs_sub, tolerance = 1e-10)
    }
  }
  expect_equal(qcomppareto(c(0, 1), "exp", rate = 1, lower.tail = FALSE), c(Inf, 0))
  expect_equal(qcomppareto(c(-Inf, 0), "exp", rate = 1, lower.tail = FALSE, log.p = TRUE),
               c(Inf, 0))
})

test_that("previously untested exports (qcmp, rcmp, BELL/ZIBELL/GPIG/ZIGPIG d/p/q/r) satisfy distribution contracts", {
  ## qcmp & rcmp
  ks <- 0:15
  expect_equal(qcmp(pcmp(ks, mu = 3.2, nu = 0.8), mu = 3.2, nu = 0.8), ks)
  expect_equal(qcmp(pcmp(ks, mu = 3.2, nu = 1.7, lower.tail = FALSE, log.p = TRUE),
                    mu = 3.2, nu = 1.7, lower.tail = FALSE, log.p = TRUE), ks)
  set.seed(101)
  yc <- rcmp(2000, mu = 3.5, nu = 1.2)
  expect_equal(mean(yc), 3.5, tolerance = 0.15)

  ## pBELL, qBELL, rBELL
  pmf_b <- dBELL(0:200, mu = 2.5)
  expect_equal(sum(pmf_b), 1, tolerance = 1e-10)
  expect_equal(pBELL(0:25, mu = 2.5), cumsum(pmf_b[1:26]), tolerance = 1e-10)
  expect_equal(qBELL(pBELL(0:12, mu = 2.5), mu = 2.5), 0:12)
  set.seed(102)
  expect_equal(mean(rBELL(3000, mu = 2.5)), 2.5, tolerance = 0.15)

  ## dZIBELL, pZIBELL, qZIBELL, rZIBELL
  pmf_zb <- dZIBELL(0:200, mu = 2.5, sigma = 0.3)
  expect_equal(sum(pmf_zb), 1, tolerance = 1e-10)
  expect_equal(pmf_zb[1], 0.3 + 0.7 * dBELL(0, mu = 2.5), tolerance = 1e-12)
  expect_equal(pZIBELL(0:25, mu = 2.5, sigma = 0.3), cumsum(pmf_zb[1:26]), tolerance = 1e-10)
  expect_equal(qZIBELL(pZIBELL(0:12, mu = 2.5, sigma = 0.3), mu = 2.5, sigma = 0.3), 0:12)
  set.seed(103)
  expect_equal(mean(rZIBELL(3000, mu = 2.5, sigma = 0.3)), 0.7 * 2.5, tolerance = 0.15)

  ## pGPIG, qGPIG, rGPIG
  pmf_g <- dGPIG(0:500, mu = 2.8, sigma = 0.55, nu = 0.65)
  expect_equal(sum(pmf_g), 1, tolerance = 1e-8)
  expect_equal(pGPIG(0:25, mu = 2.8, sigma = 0.55, nu = 0.65),
               cumsum(pmf_g[1:26]), tolerance = 1e-10)
  expect_equal(qGPIG(pGPIG(0:15, mu = 2.8, sigma = 0.55, nu = 0.65),
                     mu = 2.8, sigma = 0.55, nu = 0.65), 0:15)
  set.seed(104)
  expect_equal(mean(rGPIG(3000, mu = 2.8, sigma = 0.55, nu = 0.65)), 2.8, tolerance = 0.25)

  ## dZIGPIG, pZIGPIG, qZIGPIG, rZIGPIG
  pmf_zg <- dZIGPIG(0:500, mu = 2.8, sigma = 0.55, nu = 0.65, tau = 0.25)
  expect_equal(sum(pmf_zg), 1, tolerance = 1e-8)
  expect_equal(pmf_zg[1], 0.25 + 0.75 * dGPIG(0, mu = 2.8, sigma = 0.55, nu = 0.65),
               tolerance = 1e-12)
  expect_equal(pZIGPIG(0:25, mu = 2.8, sigma = 0.55, nu = 0.65, tau = 0.25),
               cumsum(pmf_zg[1:26]), tolerance = 1e-10)
  expect_equal(qZIGPIG(pZIGPIG(0:15, mu = 2.8, sigma = 0.55, nu = 0.65, tau = 0.25),
                       mu = 2.8, sigma = 0.55, nu = 0.65, tau = 0.25), 0:15)
  set.seed(105)
  expect_equal(mean(rZIGPIG(3000, mu = 2.8, sigma = 0.55, nu = 0.65, tau = 0.25)),
               0.75 * 2.8, tolerance = 0.25)
})
