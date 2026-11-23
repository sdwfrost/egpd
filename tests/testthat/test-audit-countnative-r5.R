## Round-5 regression tests for the native count-family numerics (fix agent `countnative`):
## VB5-01 (CMP mean inversion at integer means), VB5-03 (GPIG zero mass for small a or c) and
## VB5-04 (PIG recursion under/overflow). Every reference is independent of the package: a closed
## form, a limit law, or a high-precision Rmpfr value hard-coded below.

rel_err <- function(x, ref) max(abs(x / ref - 1))

## ---- VB5-01: CMP mean inversion at integer means --------------------------------------------

test_that("the CMP mean inversion finds the centred root at an integer mean (VB5-01)", {
  ## At mu = 1 the j = 0 and j = 2 terms of sum_j (j - 1) lambda^j / (j!)^nu balance exactly at
  ## log lambda = nu log(2) / 2, and the j >= 3 terms are below e^-140 of the j = 1 term, so
  ## log p(0) = log p(2) = -nu log(2) / 2 - log1p(2 exp(-nu log(2) / 2)).
  expect_lt(rel_err(egpd:::cmp_ll_from_mu_cpp(1, 128), 64 * log(2)), 1e-14)
  expect_lt(rel_err(egpd:::cmp_ll_from_mu_cpp(1, 1000), 500 * log(2)), 1e-14)
  lp <- dcmp(0:2, mu = 1, nu = 128, log = TRUE)
  ref <- -64 * log(2) - log1p(2 * 2^-64)
  expect_lt(rel_err(lp[c(1, 3)], ref), 1e-13)
  expect_lt(abs(lp[2]), 1e-18)                          # log p(1) = -2^-63 to first order
  expect_lt(rel_err(dcmp(0, mu = 1, nu = 1000, log = TRUE), -500 * log(2)), 1e-13)
  ## log F(0) = -44.36 < -40 <= log F(1), and -346.57 < -100 <= log F(1)
  expect_identical(qcmp(-40, mu = 1, nu = 128, log.p = TRUE), 1)
  expect_identical(qcmp(-100, mu = 1, nu = 1000, log.p = TRUE), 1)
})

test_that("CMP log-probabilities near integer means match a 600-bit centred-root solve (VB5-01)", {
  ## References: 600-bit Rmpfr bisection on sum_j (j - mu) exp(j L - nu lgamma(j + 1)) = 0 over
  ## j = 0..80, then log p(y) = y L - nu lgamma(y + 1) - log Z.
  refs <- list(
    list(mu = 2, nu = 128, y = c(0, 1, 3, 4),
         lp = c(-140.62237294952879, -25.949766918933265, -25.949766918933265, -88.722839111683744)),
    list(mu = 5, nu = 1000, y = c(3, 4, 6, 7),
         lp = c(-405.46510810816438, -91.160778396977313, -91.160778396977313, -336.47223662121293)),
    list(mu = 1 + 1e-12, nu = 128, y = c(0, 2, 3, 4),
         lp = c(-61.091906892377373, -27.630932219297627, -107.16139827643929, -223.51516960740892)),
    list(mu = 3, nu = 100, y = c(0, 2, 4, 5),
         lp = c(-193.56005167803722, -14.384104755231718, -14.384104755231718, -51.082563509241740)))
  for (r in refs) expect_lt(rel_err(dcmp(r$y, r$mu, r$nu, log = TRUE), r$lp), 1e-12)
  ## tails and quantiles from the same law: log F(4) = -91.16 < -50 <= log F(5)
  expect_lt(rel_err(pcmp(4, 5, 1000, log.p = TRUE), -91.160778396977313), 1e-12)
  expect_lt(rel_err(pcmp(5, 5, 1000, lower.tail = FALSE, log.p = TRUE), -91.160778396977313), 1e-12)
  expect_identical(qcmp(-50, mu = 5, nu = 1000, log.p = TRUE), 5)
})

test_that("the native CMP likelihood uses the centred root at an integer mean (VB5-01)", {
  X <- matrix(1, 1, 1)
  off <- list(numeric(0), numeric(0))
  nu <- exp(log(128))
  L0 <- nu * log(2) / 2                                 # the centred root at mu = exp(0) = 1
  lp0 <- -egpd:::cmp1d0(list(0, log(128)), X, X, 0, integer(0), 0L, off)
  expect_lt(rel_err(lp0, -L0 - log1p(2 * exp(-L0))), 1e-13)
  ## NLL score in log mu at y = 0: (mu - y) mu / Var = 1 / Var, Var = 2 e^-L0 / (1 + 2 e^-L0)
  g <- egpd:::cmp1d12(list(0, log(128)), X, X, 0, integer(0), 0L, off)[1, 1]
  expect_lt(rel_err(g, exp(L0) / 2 + 1), 1e-10)
})

## ---- VB5-03: GPIG zero mass ------------------------------------------------------------------

test_that("the GPIG pmf is exact for tiny a (geometric limit) and tiny c (Poisson limit) (VB5-03)", {
  ## G(s) = exp{b[(1-c)^a - (1-cs)^a]}: with a -> 0, a b = 1 it tends to exp{log((1-c)/(1-cs))},
  ## the geometric law (1-c) c^y; the deviation at a = 1e-20 is O(a).
  p <- dgpig(0:3, a = 1e-20, b = 1e20, c = 0.5)
  expect_lt(rel_err(p, 0.5^(1:4)), 1e-12)
  expect_lt(abs(sum(dgpig(0:200, a = 1e-20, b = 1e20, c = 0.5)) - 1), 1e-12)
  expect_lt(rel_err(dzigpig(0:1, a = 1e-20, b = 1e20, c = 0.5, pi = 0.3), c(0.65, 0.175)), 1e-12)
  expect_lt(abs(pgpig(0, a = 1e-20, b = 1e20, c = 0.5) - 0.5), 1e-12)
  expect_identical(qgpig(0.6, a = 1e-20, b = 1e20, c = 0.5), 1)   # F(0) = 0.5, F(1) = 0.75
  ## c -> 0 with a b c = 3: Poisson(3) up to relative O(c) = 1e-12
  expect_lt(rel_err(dgpig(0:3, a = 0.5, b = 6e12, c = 1e-12), dpois(0:3, 3)), 1e-10)
})

test_that("the GPIG likelihood ladder has the exact zero mass and b-score for tiny a (VB5-03)", {
  ## log p0 = b expm1(a log(1 - c)) and d log p0 / db = expm1(a log(1 - c)) = a log(1 - c) (1 + O(a))
  old <- gpig_control()                                 # default (hybrid) evaluation
  r <- egpd:::gpig_logp_cpp(0L, 1e-20, 1e20, 0.5)
  do.call(gpig_control, old)
  expect_lt(rel_err(r[["logp"]], log(0.5)), 1e-13)
  expect_lt(rel_err(r[["db"]], 1e-20 * log(0.5)), 1e-12)
})

## ---- VB5-04: PIG recursion under/overflow ----------------------------------------------------

test_that("the PIG pmf is finite and exact at tiny mu and at huge sigma mu (VB5-04)", {
  ## mu -> 0: P(y) = mu^y E[Z^y] / y! (1 + O(mu)), Z inverse Gaussian with mean 1 and variance
  ## sigma = 1, whose moments are 1, 1, 2, 7, 37; log P(0) = -2 mu / (1 + sqrt(1 + 2 mu)) = -mu.
  mu <- 1e-200
  ref <- c(-mu, (1:4) * log(mu) + log(c(1, 2, 7, 37)) - lgamma(2:5))
  lp <- dpig(0:4, mu = mu, sigma = 1, log = TRUE)
  expect_true(all(is.finite(lp)))
  expect_lt(rel_err(lp, ref), 1e-13)
  expect_identical(ppig(0, mu = mu, sigma = 1), 1)
  expect_identical(qpig(0.5, mu = mu, sigma = 1), 0)
  expect_lt(rel_err(ppig(3, mu = mu, sigma = 1, lower.tail = FALSE, log.p = TRUE), ref[5]), 1e-13)
  ## sigma mu -> Inf with mu = sigma: b = sqrt(2) sigma, t0 = 1/sqrt(2), t1 = (1 + sqrt(2)) / 2, so
  ## p0 = e^-sqrt(2), p1 = p0 t0, p2 = p0 t0 t1 / 2
  p0 <- exp(-sqrt(2))
  expect_lt(rel_err(dpig(0:2, mu = 1e160, sigma = 1e160),
                    p0 * c(1, 1 / sqrt(2), (1 + sqrt(2)) / (4 * sqrt(2)))), 1e-12)
  expect_identical(qpig(0.5, mu = 1e160, sigma = 1e160), 2)          # F(1) = 0.415 < 0.5 <= F(2)
  ## no NaN anywhere on an extreme grid
  for (m in c(1e-300, 1e-170, 1e-120, 1e120, 1e300)) for (s in c(1e-300, 1e-100, 1, 1e100, 1e300))
    expect_true(all(is.finite(dpig(0:10, mu = m, sigma = s, log = TRUE))))
})

test_that("the PIG upper tail stays finite where 2 sigma mu overflows (VB5-04)", {
  ## For sigma mu -> Inf (mu = sigma) the pmf decays like (2 pi)^-1/2 y^-3/2 (Bessel-K closed form of
  ## gamlss.dist::PIG), so P(Y > k) = sqrt(2 / pi) k^-1/2 (1 + O(1/k)).
  s <- ppig(1e6, mu = 1e160, sigma = 1e160, lower.tail = FALSE)
  expect_lt(abs(s / (sqrt(2 / pi) / 1e3) - 1), 1e-5)
})

test_that("the PIG native score and Hessian are finite at tiny mu (VB5-04)", {
  ## mu -> 0: log P(2) = 2 log mu + log(1 + sigma) - log 2, so the NLL score in (log mu, log sigma)
  ## is (-2, -sigma/(1+sigma)) and the Hessian (0, 0, -sigma/(1+sigma)^2).
  X <- matrix(1, 1, 1)
  off2 <- list(numeric(0), numeric(0))
  r <- egpd:::pig1d12(list(-240, 0), X, X, 2, integer(0), 0L, off2)
  expect_equal(as.numeric(r), c(-2, -0.5, 0, 0, -0.25), tolerance = 1e-8)
  expect_equal(egpd:::pig1d0(list(-400, 0), X, X, 2, integer(0), 0L, off2), 800, tolerance = 1e-13)
  ## ZIPIG at y = 2: the PIG block plus pi = plogis(-1) in the zero-inflation score and
  ## pi (1 - pi) in its Hessian
  off3 <- list(numeric(0), numeric(0), numeric(0))
  rz <- egpd:::zipig1d12(list(-240, 0, -1), X, X, X, 2, integer(0), 0L, off3)
  expect_equal(as.numeric(rz), c(-2, -0.5, plogis(-1), 0, 0, 0, -0.25, 0, plogis(-1) * plogis(1)),
               tolerance = 1e-8)
})
