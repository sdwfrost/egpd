## Round-4 regression (gwtail): the generalized Waring tails when a = mu (rho - 1) / k is small.
## The mixture integral behind pgwaring()/qgwaring() returned P(Y > q) = 1 once a/rho < 2^-53
## (mu <= 1e-16 at k = 1.5, rho = 3; also mu = 1e-6 at k = 1e10, or mu = 1e-10 with rho near 1,
## inside the fitting clamps), and its lower tail lost the mixing mass beyond the double range
## for every a below ~0.056.
##
## Expected values are independent: log S(q) and log F(q) from exact summation of the pmf
##   p(0) = Gamma(a+rho) Gamma(k+rho) / (Gamma(rho) Gamma(a+k+rho)),
##   p(y+1) / p(y) = (a+y) (k+y) / ((a+k+rho+y) (y+1)),   S(q) = 1 - sum_{y <= q} p(y),
## in 512-bit arithmetic (Rmpfr, not needed here), and the small-a expansion of log p(0).

expect_rel <- function(x, ref, tol) {
  err <- abs(x / ref - 1)
  expect_true(all(err < tol), label = sprintf("max relative error %.3g", max(err)))
}

gw_refs <- list(
  list(mu = 1e-12, k = 1.5, rho = 3,
       lS = c(-28.106722887439048, -29.362602120809061, -33.404389024566342, -46.597813413575714),
       lF = c(-6.2144878834795237e-13, -1.770043439037384e-13, -3.1092681294448487e-15,
              -5.7919761614271523e-21)),
  list(mu = 1e-16, k = 1.5, rho = 3,
       lS = c(-37.317063259414709, -38.572942492785721, -42.614729396544156, -55.808153785553891),
       lF = c(-6.2144878834808236e-17, -1.7700434390363792e-17, -3.1092681294397729e-19,
              -5.7919761614155953e-25)),
  list(mu = 1e-20, k = 1.5, rho = 3,
       lS = c(-46.527403631390897, -47.783282864761908, -51.825069768520336, -65.018494157530071),
       lF = c(-6.2144878834808229e-21, -1.770043439036379e-21, -3.1092681294397723e-23,
              -5.7919761614155946e-29)),
  ## a = 1.5e-16 through a large k at a moderate mean
  list(mu = 1e-6, k = 1e10, rho = 2.5,
       lS = c(-33.330292537673607, -33.376124414667984, -33.470946766134077, -33.738757243188715),
       lF = c(-3.3484041434242821e-15, -3.1984041434617821e-15, -2.9090589057040439e-15,
              -2.225583551841727e-15)),
  ## a = 1e-16 inside the fitting clamps (mu >= 1e-10, rho >= 1 + 1e-6)
  list(mu = 1e-10, k = 1, rho = 1 + 1e-6,
       lS = c(-36.841362487986494, -37.534510168546319, -39.239259780661932, -43.750123753771256),
       lF = c(-9.9999899991873356e-17, -4.9999924995974178e-17, -9.0908816367793448e-18,
              -9.9899351995676177e-20)))
gw_q <- c(0, 1, 10, 1000)

test_that("the generalized Waring upper tail is right when a = mu (rho - 1) / k is tiny (gwtail)", {
  for (r in gw_refs) {
    expect_rel(egpd:::gwaring_lcdf_cpp(gw_q, r$mu, r$k, r$rho, FALSE), r$lS, 1e-13)
    expect_rel(pgwaring(gw_q, r$mu, r$k, r$rho, lower.tail = FALSE, log.p = TRUE), r$lS, 1e-13)
    expect_rel(pgwaring(gw_q, r$mu, r$k, r$rho, lower.tail = FALSE), exp(r$lS), 1e-12)
  }
})

test_that("log P(Y > 0) and log P(Y = 0) follow the small-a expansion of p(0) (gwtail)", {
  ## -log p(0) = L = a [psi(k+rho) - psi(rho)] + a^2/2 [psi'(k+rho) - psi'(rho)] + O(a^3),
  ## so log S(0) = log(1 - e^-L) and log F(0) = -L, to relative O(a^2) <= 1e-24 here
  for (kr in list(c(1.5, 3), c(1, 2.5), c(40, 7), c(0.2, 1.1))) {
    k <- kr[1]; rho <- kr[2]
    for (mu in c(1e-12, 1e-16, 1e-20, 1e-50, 1e-100, 1e-300)) {
      a <- mu * ((rho - 1) / k)
      L <- a * (digamma(k + rho) - digamma(rho)) + a^2 / 2 * (trigamma(k + rho) - trigamma(rho))
      expect_rel(egpd:::gwaring_lcdf_cpp(0, mu, k, rho, FALSE), log(-expm1(-L)), 1e-13)
      expect_rel(egpd:::gwaring_lcdf_cpp(0, mu, k, rho, TRUE), -L, 1e-12)
    }
  }
})

test_that("the generalized Waring lower tail is log(1 - P(Y > q)) for small a (gwtail)", {
  for (r in gw_refs)
    expect_rel(egpd:::gwaring_lcdf_cpp(gw_q, r$mu, r$k, r$rho, TRUE), r$lF, 1e-12)
  ## small but ordinary a = 1.3e-4, 1.3e-3, 0.013 (log F(0) was 2.4, 0.49 and 7e-5 too low);
  ## k = 1.5, rho = 3
  lF <- list(
    c(-6.214357923613893e-05, -1.7701439191640855e-05, -1.4972796669197255e-06, -4.0631925941343331e-09),
    c(-0.0006213188613787532, -0.00017710474287256441, -1.4991352167674765e-05, -4.0701930713633166e-08),
    c(-0.0062015279805663253, -0.0017800027226889351, -0.0001517678450051124, -4.1404587566480926e-07))
  mus <- c(1e-4, 1e-3, 0.01)
  for (i in seq_along(mus))
    expect_rel(egpd:::gwaring_lcdf_cpp(c(0, 1, 5, 50), mus[i], 1.5, 3, TRUE), lF[[i]], 1e-12)
  ## a = 0.06 and 0.07, either side of the switch at a = 1/16 (above it the integral of F is
  ## kept): F to 1e-14 relative, i.e. log F to 1e-14 absolute
  lF <- list(
    c(-0.027705304265887694, -0.0081607081929174641, -0.00071530537849654526, -1.9881922723050087e-06),
    c(-0.03227305351405934, -0.0095570956795223531, -0.00084258478559601886, -2.3512897754805499e-06))
  mus <- c(0.045, 0.0525)
  for (i in seq_along(mus))
    expect_lt(max(abs(egpd:::gwaring_lcdf_cpp(c(0, 1, 5, 50), mus[i], 1.5, 3, TRUE) - lF[[i]])), 1e-14)
})

test_that("pgwaring, qgwaring and hurdle quantiles are finite and right at tiny a (gwtail)", {
  ## mu = 1e-16, k = 1.5, rho = 3: log S(0:3) and log F(0:3)
  lS <- c(-37.317063259414709, -38.572942492785721, -39.418459180870649, -40.066904322134206)
  lF <- c(-6.2144878834808236e-17, -1.7700434390363792e-17, -7.5994242893536905e-18,
          -3.973420663350064e-18)
  expect_rel(pgwaring(0:3, 1e-16, 1.5, 3, lower.tail = FALSE, log.p = TRUE), lS, 1e-13)
  expect_rel(pgwaring(0:3, 1e-16, 1.5, 3, log.p = TRUE), lF, 1e-12)
  ## S(1) = 1.770e-17 > 1e-17 >= S(2) = 7.599e-18; S(36) = 1.0333e-20 > 1e-20 >= S(37) = 9.565e-21
  expect_identical(expect_silent(qgwaring(c(1e-20, 1e-17), 1e-16, 1.5, 3, lower.tail = FALSE)),
                   c(37, 2))
  expect_identical(qgwaring(log(c(1e-20, 1e-17)), 1e-16, 1.5, 3, lower.tail = FALSE, log.p = TRUE),
                   c(37, 2))
  ## zero-truncated quantiles at (mu, k, rho) = (1e-10, 1, 1 + 1e-6), inside the fitting clamps:
  ## the smallest q >= 1 with S(q) <= (1 - u) S(0); S(q)/S(0) at q - 1 and q is
  ## (1, 0.49999975), (0.1111109, 0.0999998), (0.0101010, 0.00999996), (0.00100099, 0.000999994)
  pr <- list(mu = rep(1e-10, 4), k = rep(1, 4), rho = rep(1 + 1e-6, 4))
  expect_identical(expect_silent(egpd:::.zt_quantile("gwaring", pr, c(0.5, 0.9, 0.99, 0.999))),
                   c(1, 9, 99, 999))
})

test_that("a near-one upper input does not erase the Waring lower tail", {
  ## Beta-negative-binomial mixture, a=256, k=16, rho=5 (mu=1024).
  i <- 0:15
  j <- 0:30
  p0 <- prod((5+i)/(261+i))
  cdf <- cumsum(cumprod(c(p0, (256+j)*(16+j)/((j+1)*(277+j)))))
  small <- 2^-c(48, 50)
  ref <- vapply(small, function(p) which(cdf >= p)[1L]-1, 0)
  expect_identical(ref, c(13, 12))
  expect_identical(qgwaring(1-small, 1024, 16, 5, lower.tail=FALSE, max.value=32), ref)
  expect_warning(q <- qgwaring(1-small[1], 1024, 16, 5, lower.tail=FALSE, max.value=12),
                 "not attained")
  expect_true(is.na(q))
  ## The same search engine must keep a mixed vector in its original order.
  upper <- c(.2, 1-2^-48, .1, 1-2^-50)
  expect_identical(qgpois(upper, 50, 0, lower.tail=FALSE), qpois(1-upper, 50))
})
