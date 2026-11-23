## Round-5 regression (gwlargek): the generalized Waring tails for a negative-binomial shape k > 1e12.
## gw_lcdf() integrates P(NB(k, w) > q) or P(NB(k, w) <= q) over the Beta(rho, a) mixing law. For a
## large k the integrand matters where 1 - w ~ (q + 40)/k. The continued fraction that cf_lpbeta()
## used for log P(NB(k, w) <= q) there sees only w, which no longer resolves 1 - w. Relative errors
## were 1e-12 at k = 1e13, 1e-6 at k = 1e15 (2.7e-5 in P(Y > q); 7e-3 in log P(Y <= q) on the
## direct lower-tail path), and NaN from k ~ 1e18.
##
## Expected values are independent: log S(q) and log F(q) from exact summation of the pmf
##   p(0) = Gamma(a+rho) Gamma(k+rho) / (Gamma(rho) Gamma(a+k+rho)),
##   p(y+1) / p(y) = (a+y) (k+y) / ((a+k+rho+y) (y+1)),   S(q) = 1 - sum_{y <= q} p(y),
## with a = mu (rho - 1) / k, in 512-bit arithmetic (Rmpfr, not needed here), and the k -> Inf
## expansion of the pmf.

expect_rel <- function(x, ref, tol) {
  err <- abs(x / ref - 1)
  expect_true(all(err < tol), label = sprintf("max relative error %.3g", max(err)))
}
gw_q <- c(0, 1, 10, 1000)

## k = 1e15 (the fitting bound for k), q = 0, 1, 10, 1000
gw_k15 <- list(
  list(mu = 1e-3, rho = 1.5,
       lS = c(-38.598653259540917, -38.62806516657966, -38.687366452339063, -38.843219383774247),
       lF = c(-1.7251143210466054e-17, -1.6751143210466055e-17, -1.5786659083481935e-17, -1.3508407780191633e-17)),
  list(mu = 1, rho = 1.5,
       lS = c(-31.690897980558791, -31.72030988759753, -31.77961117335693, -31.935464104792114),
       lF = c(-1.7251143210466056e-14, -1.6751143210466055e-14, -1.5786659083481934e-14, -1.3508407780191632e-14)),
  list(mu = 1e3, rho = 1.5,
       lS = c(-24.78314270158528, -24.81255460862376, -24.871855894382673, -25.02770882581672),
       lF = c(-1.7251143210465938e-11, -1.6751143210466065e-11, -1.5786659083482013e-11, -1.3508407780191721e-11)),
  list(mu = 1e-3, rho = 3,
       lS = c(-37.238382585432895, -37.268581771260514, -37.329544605427039, -37.490280467605679),
       lF = c(-6.7231984119624448e-17, -6.5231984119624452e-17, -6.1374047611687997e-17, -5.2261042398529753e-17)),
  list(mu = 1, rho = 3,
       lS = c(-30.330627306450793, -30.360826492278409, -30.421789326444934, -30.582525188623567),
       lF = c(-6.7231984119624436e-14, -6.5231984119624445e-14, -6.1374047611688001e-14, -5.2261042398529757e-14)),
  list(mu = 1e3, rho = 3,
       lS = c(-23.422872027502248, -23.453071213328837, -23.514034047493418, -23.674769909667486),
       lF = c(-6.7231984119623647e-11, -6.5231984119625657e-11, -6.1374047611690299e-11, -5.2261042398532249e-11)),
  list(mu = 1e-3, rho = 30,
       lS = c(-34.640282338314677, -34.672907041789855, -34.739014624718031, -34.915075271051904),
       lF = c(-9.0347580960453062e-16, -8.7447580960453143e-16, -8.1853573023945988e-16, -6.8639715464944056e-16)),
  list(mu = 1, rho = 30,
       lS = c(-27.732527059332991, -27.765151762808156, -27.831259345736299, -28.007319992070109),
       lF = c(-9.0347580960453051e-13, -8.7447580960453175e-13, -8.1853573023946047e-13, -6.8639715464944122e-13)),
  list(mu = 1e3, rho = 30,
       lS = c(-20.824771780802156, -20.857396484262356, -20.923504067162241, -21.099564713429842),
       lF = c(-9.0347580960451628e-10, -8.7447580960493762e-10, -8.1853573024009729e-10, -6.8639715465011756e-10)))

## mu = 1, rho = 3 at k = 1e13 ... 1e20
gw_kser <- list(
  list(k = 1e13,
       lS = c(-25.872790099926053, -25.907868097366649, -25.979219219451362, -26.171227670933135),
       lF = c(-5.8021643747648675e-12, -5.6021643747649476e-12, -5.216370723971848e-12, -4.3050702027148311e-12)),
  list(k = 1e14,
       lS = c(-28.098997777923234, -28.131454032552078, -28.197203059005886, -28.372185443048767),
       lF = c(-6.2626813933636383e-13, -6.0626813933636462e-13, -5.6768877425700509e-13, -4.7655872212595734e-13)),
  list(k = 1e16,
       lS = c(-32.566959676254193, -32.595195323113998, -32.6520216773174, -32.800660126579167),
       lF = c(-7.1837154305612532e-15, -6.9837154305612531e-15, -6.5979217797676026e-15, -5.6866212584512442e-15)),
  list(k = 1e18,
       lS = c(-37.051496337014441, -37.076482803622703, -37.12651930670512, -37.255717710742786),
       lF = c(-8.1047494677588715e-17, -7.9047494677588707e-17, -7.5189558169652203e-17, -6.607655295648803e-17)),
  list(k = 1e20,
       lS = c(-41.549031451482499, -41.571439385577101, -41.616135660393134, -41.730395799608104),
       lF = c(-9.0257835049564902e-19, -8.8257835049564894e-19, -8.4399898541628386e-19, -7.5286893328464202e-19)))

test_that("the generalized Waring upper tail is right for k > 1e12 (gwlargek)", {
  for (r in gw_k15) {
    expect_rel(egpd:::gwaring_lcdf_cpp(gw_q, r$mu, 1e15, r$rho, FALSE), r$lS, 1e-13)
    expect_rel(pgwaring(gw_q, r$mu, 1e15, r$rho, lower.tail = FALSE, log.p = TRUE), r$lS, 1e-13)
  }
  for (r in gw_kser) {
    expect_rel(egpd:::gwaring_lcdf_cpp(gw_q, 1, r$k, 3, FALSE), r$lS, 1e-13)
    expect_rel(expect_silent(pgwaring(gw_q, 1, r$k, 3, lower.tail = FALSE, log.p = TRUE)), r$lS, 1e-13)
  }
})

test_that("the generalized Waring lower tail is right for k > 1e12 (gwlargek)", {
  for (r in gw_k15) {
    expect_rel(egpd:::gwaring_lcdf_cpp(gw_q, r$mu, 1e15, r$rho, TRUE), r$lF, 1e-12)
    expect_rel(pgwaring(gw_q, r$mu, 1e15, r$rho, log.p = TRUE), r$lF, 1e-12)
  }
  for (r in gw_kser)
    expect_rel(egpd:::gwaring_lcdf_cpp(gw_q, 1, r$k, 3, TRUE), r$lF, 1e-12)
  ## the direct integral of the lower tail (a >= 1/16, so mu ~ k); there S ~ 1, so log S is
  ## checked to 1e-14 absolute, i.e. S to 1e-14 relative
  direct <- list(
    list(mu = 5e12, k = 1e13, rho = 3,
         lF = c(-28.834993920254785, -28.14184673969504, -26.437098647458413, -21.926239141139565),
         lS = c(-2.99999999999955e-13, -5.9999999999988001e-13, -3.299999999997855e-12, -3.0029999998493997e-10)),
    list(mu = 2e13, k = 1e14, rho = 1.5,
         lF = c(-3.2154286582531757, -3.1201184784488523, -2.929906627588414, -2.4747257000955676),
         lS = c(-0.040965884660563105, -0.045156308537154771, -0.054880799889833416, -0.087942078311243205)),
    list(mu = 10 * 1e15 / 29, k = 1e15, rho = 30,
         lF = c(-310.01304265563175, -307.61514738283341, -297.88625134102966, -255.98509393687834),
         lS = c(-2.3069928930039818e-135, -2.5376921823042879e-134, -4.2623077893968875e-130, -6.7156221633611305e-112)))
  for (r in direct) {
    expect_rel(egpd:::gwaring_lcdf_cpp(gw_q, r$mu, r$k, r$rho, TRUE), r$lF, 1e-12)
    expect_lt(max(abs(egpd:::gwaring_lcdf_cpp(gw_q, r$mu, r$k, r$rho, FALSE) - r$lS)), 1e-14)
  }
})

test_that("both tails follow the k -> Inf expansion at fixed mu and rho (gwlargek)", {
  ## As k -> Inf with mu and rho fixed, a = mu (rho - 1) / k -> 0, and from the pmf
  ##   S(q) / a = psi(k + rho) - psi(rho) - H_q + rho q / k + O(a log(k)^2) + O(q^2 / k^2),
  ## so the law collapses onto 0 (the negative-binomial limit is rho -> Inf instead). With
  ## mu = 1e-3 the remainder is below 1e-15 relative for these k and q.
  q <- c(0, 1, 10, 1000, 1e5)
  for (k in c(1e15, 1e18, 1e20)) for (rho in c(1.5, 3, 30)) {
    a <- 1e-3 * ((rho - 1) / k)
    S <- a * (digamma(k + rho) - digamma(rho) - (digamma(q + 1) - digamma(1)) + rho * q / k)
    expect_rel(egpd:::gwaring_lcdf_cpp(q, 1e-3, k, rho, FALSE), log(S), 1e-13)
    expect_rel(egpd:::gwaring_lcdf_cpp(q, 1e-3, k, rho, TRUE), log1p(-S), 1e-12)
  }
})

test_that("the lower tail is right when t* = a/rho exceeds 2^53 (gwlargek)", {
  ## log F(0) = log p(0); 512-bit Rmpfr. The mixing law's centre t* = a/rho = 6.7e16 rounds
  ## 1 - w to 1, so the second logarithm of the mixing drop was -Inf there: F(0) came out as 1.
  expect_rel(egpd:::gwaring_lcdf_cpp(0, 1e34, 1e17, 3, TRUE), -1.9095425048844374e+17, 1e-12)
})
