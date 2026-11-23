## Round-4 native audit: regression tests for VA-06, VA-07, VA-09 and VA-20
## (round-4 verification, area A).
##
## Independent references: own 512-bit Rmpfr implementations of the continuous
## EGPD, DEGPD and ZIDEGPD likelihoods (no package code). The carrier-5
## normaliser is erf(sqrt(kappa/2))/2 and the carrier-6 normaliser
## Z(kappa) = int_0^1 (1 - C t^2)^(kappa-1) dt, C = (15/16)^2, is summed as its
## hypergeometric series sum_n (1-kappa)_n/n! C^n/(2n+1). DEGPD cells are
## differences of G(H(y+1)) and G(H(y)) at 512 bits. Derivatives are mpfr
## central differences (step 2^-70, 2^-90 at y = 1e12). Values are hard-coded.

r4n_ld <- function(eta, y, xi.max = Inf) {
  n <- length(y)
  list(X = rep(list(matrix(1, n, 1)), length(eta)), y = cbind(y),
       idpars = seq_along(eta), dupid = 0L, duplicate = 0L,
       offsets = lapply(eta, rep, n), censored = FALSE, xi.max = xi.max)
}
r4n_fns <- function(name) get(name, envir = asNamespace("egpd"))
r4n_d0 <- function(fns, eta, y, xi.max = Inf)
  fns$d0(rep(0, length(eta)), r4n_ld(eta, y, xi.max))
r4n_d12 <- function(fns, eta, y, xi.max = Inf)
  fns$d120(rep(0, length(eta)), r4n_ld(eta, y, xi.max))[1, ]
## largest error relative to the row's scale (max(1, |reference|))
r4n_rowerr <- function(got, want) max(abs(got - want)) / max(1, abs(want))

test_that("VA-20: bounded shape link uses t/(e^t - 1) = 1 - t/2 + t^2/12 for small t", {
  ## xi = xi.max (1 - exp(-exp(eta)/xi.max)) with t = exp(eta)/xi.max = 9e-5 < 1e-4;
  ## the series sign error made the log-shape score 1.35e-9 too small (relative).
  eta1 <- c(log(2), log(9e-5), log(1.5))
  ref1 <- c(-0.54448493173336132, 3.4446356751319602e-05, -0.70368369478668891,
            1.9099647829206754, 0.00011573841907891278, 0.55317568116746108,
            3.4446696351782377e-05, 4.2898414195979206e-05, 0.29416958108564745)
  got1 <- r4n_d12(r4n_fns(".degpd1fns"), eta1, 3, xi.max = 1)   # chain rule in C++
  expect_lt(abs(got1[2] / ref1[2] - 1), 1e-12)
  expect_lt(r4n_rowerr(got1, ref1), 1e-13)
  eta3 <- c(log(2), log(9e-5), log(0.7))
  ref3 <- c(-0.22184102506377551, 5.9625998487999983e-05, 0.074804650549536139,
            2.0761636203237219, 0.00010435605459480561, -0.36626485529897074,
            5.9621888962093184e-05, -2.8834387842527261e-05, 0.1285052847485281)
  got3 <- r4n_d12(r4n_fns(".degpd3fns"), eta3, 3, xi.max = 1)   # chain rule in R
  expect_lt(abs(got3[2] / ref3[2] - 1), 1e-12)
  expect_lt(r4n_rowerr(got3, ref3), 1e-13)
})

test_that("VA-06: continuous carriers 5/6 keep their normaliser as kappa -> 0", {
  f5 <- r4n_fns(".egpd5fns"); f6 <- r4n_fns(".egpd6fns")
  e <- function(lk) c(0.1, 0.2, lk)
  ## carrier 5: 1/2 - pnorm(-sqrt(kappa)) lost every digit below kappa ~ 1e-32 (1e20 penalty)
  expect_lt(abs(r4n_d0(f5, e(-37), 1.9) - 1.873179171872944), 1e-12)
  expect_lt(abs(r4n_d0(f5, e(-80), 1.9) - 1.873179171872944), 1e-12)
  expect_lt(abs(r4n_d0(f5, e(-80), 0) - 0.1), 1e-12)
  ## carrier 6: pbeta(1/2) - pbeta(1/32) cancelled (0.28 too small at log kappa = -37)
  expect_lt(abs(r4n_d0(f6, e(-12), 1.9) - 2.4314493828196837), 1e-12)
  expect_lt(abs(r4n_d0(f6, e(-30), 1.9) - 2.431453955569653), 1e-12)
  expect_lt(abs(r4n_d0(f6, e(-37), 1.9) - 2.4314539555697228), 1e-12)
  expect_lt(abs(r4n_d0(f6, e(-45), 1.9) - 2.4314539555697228), 1e-12)
  expect_lt(abs(r4n_d0(f6, e(-37), 0) - -1.4060768630739156), 1e-12)
  ## both sides of the switch to the stable forms (kappa = 0.01), and below it
  expect_lt(abs(r4n_d0(f5, e(log(0.0099)), 1.9) - 1.8717879729308546), 1e-12)
  expect_lt(abs(r4n_d0(f5, e(log(0.0101)), 1.9) - 1.8717598903319048), 1e-12)
  expect_lt(abs(r4n_d0(f5, e(-5), 1.9) - 1.8722320842840219), 1e-12)
  expect_lt(abs(r4n_d0(f6, e(log(0.0099)), 1.9) - 2.4241067536388123), 1e-12)
  expect_lt(abs(r4n_d0(f6, e(log(0.0101)), 1.9) - 2.4239587527050168), 1e-12)
  expect_lt(abs(r4n_d0(f6, e(-5), 1.9) - 2.4264489354382617), 1e-12)
  ## the carrier-6 density integrates to one (it was 1.32 at log kappa = -37)
  dens <- function(y) vapply(y, function(yy) exp(-r4n_d0(f6, c(0, 0.2, -37), yy)), 0)
  expect_lt(abs(integrate(dens, 0, Inf, rel.tol = 1e-11)$value - 1), 1e-8)
  ## carrier-5 gradient/Hessian use the same normaliser (log-kappa score was off by 1e-8)
  ref <- c(-0.53517684210664529, 0.28763766106262639, -5.9746777258773149e-19,
           1.1423821860165566, 0.68465924536157197, 2.8296192763847794e-19,
           0.096897338010231385, 2.1934149918749958e-19, -5.9746777258773149e-19)
  expect_lt(r4n_rowerr(r4n_d12(f5, e(-40), 1.9), ref), 1e-12)
  ref0 <- c(1, 0, 1.4161180850971963e-18, 0, 0, 0, 0, 0, 1.4161180850971963e-18)
  expect_lt(r4n_rowerr(r4n_d12(f5, e(-40), 0), ref0), 1e-12)
})

test_that("VA-07: continuous carriers 3/4 derivatives stay exact as delta -> 0", {
  f3 <- r4n_fns(".egpd3fns"); f4 <- r4n_fns(".egpd4fns")
  e3 <- function(ld) c(log(3), 0.25, ld)
  e4 <- function(ld) c(log(3), 0.25, ld, log(1.7))
  ref3 <- list(
    `-12` = c(1.2124495580305623, 0.68378331357031275, -4.2499067757625968e-06,
              0.67674189753367453, -0.421246282472305, -1.7554870268954639e-06,
              -0.3273595378309867, -5.551184138864551e-07, -4.2498702209660975e-06),
    `-20` = c(1.2124513129297965, 0.68378386850285544, -1.4256971540701656e-09,
              0.67674039333535807, -0.42124728527154265, -5.8890103473296445e-10,
              -0.32735996609297563, -1.8622171700988727e-10, -1.425697149956413e-09),
    `-36` = c(1.2124513135186974, 0.68378386868907715, -1.604410787929784e-16,
              0.67674039283058585, -0.42124728560805752, -6.6272080864101975e-17,
              -0.32735996623668984, -2.095649346910605e-17, -1.6044107879297835e-16))
  for (ld in names(ref3)) {
    got <- r4n_d12(f3, e3(as.numeric(ld)), 2)
    expect_lt(r4n_rowerr(got, ref3[[ld]]), 1e-12)
  }
  ## the log-delta score at delta = exp(-36) is -1.6e-16 (it was +0.356)
  expect_lt(abs(r4n_d12(f3, e3(-36), 2)[3] - -1.604410787929784e-16), 1e-15)
  ref4 <- list(
    `-9` = c(0.98849115491846706, 0.61296340342150768, -7.0440092521551438e-05,
             0.75128809498741378, 0.61397835448711569, -0.37381985927569233,
             -3.2116027255395676e-05, 1.2689077615858706, -0.29817656034017664,
             -1.0155699151977298e-05, 0.40125278807982101, -7.0427586657127684e-05,
             -8.4476148429496879e-05, 1.7512880949874137),
    `-16` = c(0.98852324206447473, 0.61297354998785936, -6.4244444949686656e-08,
              0.75137250043303427, 0.6139504800059512, -0.37383831217236341,
              -2.9286764953201566e-08, 1.2689255560407913, -0.29818442794223865,
              -9.2610325565537708e-09, 0.40125841502518728, -6.424443454804681e-08,
              -7.70438307465493e-08, 1.7513725004330343),
    `-25` = c(0.98852327134762574, 0.61297355924774921, -7.9283956484043685e-12,
              0.7513725774673623, 0.61395045456699227, -0.37383832901282726,
              -3.6142740084243216e-12, 1.2689255722802548, -0.29818443512234977,
              -1.1429022397594879e-12, 0.4012584201604148, -7.9283956482459521e-12,
              -9.5079653599936436e-12, 1.7513725774673623))
  for (ld in names(ref4)) {
    got <- r4n_d12(f4, e4(as.numeric(ld)), 2)
    expect_lt(r4n_rowerr(got, ref4[[ld]]), 1e-12)
  }
  ## carrier 4 also loses G3 to cancellation when A = -log(1 - H) is small
  ref_small <- c(1.6975201276637659, 0.0026969040699679693, -0.0056828605502312907,
                 10.149862508326921, 0.0024786811082420442, -0.0026938099397618522,
                 -6.0613189929106885e-06, 1.6984387838450448, -5.1568095162286808e-06,
                 -6.0603090255577528e-09, 0.0016981557814614631, -0.0056447854483741722,
                 -0.005685106714598911, 11.149862508326921)
  expect_lt(r4n_rowerr(r4n_d12(f4, e4(-5), 0.006), ref_small), 1e-12)
})

test_that("carrier 3 keeps third and fourth derivatives and the REML gradient at tiny delta", {
  fns <- r4n_fns(".egpd3fns")
  delta <- 1e-15
  for (xi in c(0, .2)) {
    p <- c(log(1.8), xi, log(delta))
    d <- as.numeric(fns$d340(rep(0, 3), r4n_ld(p, 2)))
    z <- 2/1.8
    A <- if (xi == 0) z else log1p(xi*z)/xi
    ## NLL = lpsi + (1+xi)*A - log(A) + delta*(A/2-1) + O(delta^2).
    ## Divide the pure and mixed log-delta entries by delta: zero is not accurate here.
    expect_equal(d[c(10,25)] / delta, rep(A/2 - 1, 2), tolerance = 1e-12)
    expect_equal(d[c(6,20)] / delta, rep(-z/(2*(1+xi*z)), 2), tolerance = 1e-12)
    ## 90-decimal differentiation of the log-density (derivatives in scale and xi).
    ref <- if (xi == 0) c(-1.11111111111111, 1.48605395518976, 1.11111111111111, -4.20413554844281)
           else c(-.623647033908267, .890044777333816, .115342298880274, -2.07383671239837)
    expect_equal(d[c(1,7,11,21)], ref, tolerance = 1e-11)
  }

  ## Compare the analytic smoothing gradient to differences of REML itself
  ## near the minimum; a pinned delta does not remove its effect on scale curvature.
  set.seed(172)
  d <- data.frame(x = seq_len(180)/180)
  d$y <- exp(.2 + 1.1*cos(2*pi*d$x))*expm1(.2*rgamma(180, 2))/.2
  rho <- -2.9270631
  f <- egpd(list(y ~ s(x, k = 7), ~1, ~1), d, egpd.args = list(m = 3),
            fixed = list(xi = .2, delta = delta), restarts = FALSE, correctV = FALSE,
            outer = "fixed", rho0 = rho, control = list(inner = list(gradtol = 1e-8)))
  reml <- function(r) as.numeric(egpd:::.reml0(r, f$likfns, f$likdata, f$Sdata, beta = f$coefficients))
  grad <- egpd:::.reml12(rho, f$likfns, f$likdata, f$Sdata, beta = f$coefficients)$d1
  expect_lt(abs(grad - (reml(rho + .001) - reml(rho - .001))/.002), 2e-6)
})

test_that("VA-09: DEGPD and ZIDEGPD objectives are exact at very large counts", {
  e1 <- c(log(3), log(0.25), log(1.7))
  f1 <- r4n_fns(".degpd1fns")
  ## the cell difference lost ~eps*y*log(y): 5.4e-3 at 1e12, a 1e20 penalty at 1e16
  expect_lt(abs(r4n_d0(f1, e1, 1e9) - 91.75978003589799), 1e-12)
  expect_lt(abs(r4n_d0(f1, e1, 1e12) - 126.29855636837118), 1e-12)
  expect_lt(abs(r4n_d0(f1, e1, 1e16) - 172.3502582281896), 1e-12)
  expect_lt(abs(r4n_d0(r4n_fns(".degpd1idfns"), c(log(3), 0.25, log(1.7)), 1e10) -
                  103.27270544461823), 1e-12)
  expect_lt(abs(r4n_d0(r4n_fns(".degpd3fns"), c(log(3), log(0.4), log(0.6)), 1e10) -
                  73.656100721350072), 1e-12)
  expect_lt(abs(r4n_d0(r4n_fns(".degpd5idfns"), c(log(3), 0.3, log(2)), 1e11) -
                  100.58537409932214), 1e-12)
  expect_lt(abs(r4n_d0(r4n_fns(".zidegpd1fns"), c(e1, -1), 1e12) -
                  126.6118180558894), 1e-12)
  ## gradient and Hessian (the log-sigma curvature was 0.005 instead of 6e-11)
  ref <- c(-3.99999999994, -95.584457864672189, -1, 5.9999999998529996e-11,
           3.9999999998919997, 1.4109707707415664e-43, 92.584457864768183,
           3.4049902420574944e-42, 3.5535290570962862e-44)
  got <- r4n_d12(f1, e1, 1e12)
  expect_lt(r4n_rowerr(got, ref), 1e-12)
  expect_lt(abs(got[4] - ref[4]), 1e-13)
})

test_that("VA-08: d34 next to the xi < 0 endpoint and in tiny-probability cells", {
  ## largest error of each block (third, fourth derivatives) relative to its scale
  blockerr <- function(got, want) c(max(abs(got[1:10] - want[1:10])) / max(abs(want[1:10])),
                                    max(abs(got[11:25] - want[11:25])) / max(abs(want[11:25])))
  d34 <- function(fns, eta, y) as.numeric(fns$d340(rep(0, length(eta)), r4n_ld(eta, y)))
  ## DEGPD carrier 1, identity shape link, top attainable bin: P = 1 - F^kappa = 1.6e-14
  ## (F^kappa - 1 cancelled: third derivatives 16% and fourth 67% off)
  ref_a <- c(-77505.689576686244, -881201.6385980472, 2.1014959288680425e-10, -9693218.0245791525,
             1.9716426176298013e-09, 1.3659723537642082e-12, -103485418.50424647, 1.8577621825503369e-08,
             1.3016587821577567e-11, 8.0716548176975178e-15, 3459625.3632587888, 39789643.219789378,
             2.9076620978884677e-08, 447379293.44227862, 2.6776585622372424e-07, 2.1014959288680189e-10,
             4921172227.8652859, 2.4765800383486027e-06, 1.971642617629779e-09, 1.3659723537641934e-12,
             53026241809.981354, 2.3010567318217284e-05, 1.8577621825503154e-08, 1.3016587821577427e-11,
             8.071654817697431e-15)
  expect_true(all(blockerr(d34(r4n_fns(".degpd1idfns"), c(0, -0.085, 0.4), 11), ref_a) < 1e-9))
  ## log shape link, far-tail cell P(Y = 160) = 1.6e-13 (7% and 160% off)
  ref_c <- c(1.0515710972519183, -2.3121937861716479, 2.8830722514496878e-10, -7.5045741839092797,
             3.9976944944394018e-10, 2.629668889053244e-11, 1.2332365591833712, 5.9956965449849502e-10,
             3.8192336624201913e-11, 2.3738967468655868e-12, 0.5376312557910139, -1.5113145782801269,
             3.1250913799687181e-09, 3.7456801347220892, 4.1274321314214299e-09, 2.883072251449677e-10,
             3.8364823561806189, 5.9098514692116728e-09, 3.9976944944393873e-10, 2.6296688890532343e-11,
             -5.1484729511684897, 9.1640446771210077e-09, 5.9956965449849295e-10, 3.8192336624201777e-11,
             2.3738967468655779e-12)
  expect_true(all(blockerr(d34(r4n_fns(".degpd1fns"), c(0.5, log(0.08), 0.4), 160), ref_c) < 1e-8))
  ## continuous carrier 3 with 1 + xi y/psi = 2e-7: the generated forms overflow and the
  ## stencil fallback was 0.3% (third) and 160% (fourth derivatives) off; the remaining
  ## error is the input rounding eps/(1 + xi y/psi)
  ref_b <- c(-1.224999631716036e+22, -6.1249975335803045e+23, 0, -3.0624981480404015e+25, 0, 0,
             -1.5312486115205733e+27, 0, 0, 0.066647734775181794, 1.8374992634321087e+29,
             9.1874956921607318e+30, 0, 4.593747227330614e+32, 0, 0, 2.2968731511655862e+34, 0, 0, 0,
             1.1484362677705716e+36, 0, 0, 0, 0.08526913163132023)
  got_b <- d34(r4n_fns(".egpd3fns"), c(0, -0.02, 0.6), 49.99999)
  expect_true(all(blockerr(got_b, ref_b) < 1e-7))
  ## the log-delta-only entries are exact
  expect_lt(abs(got_b[10] - ref_b[10]), 1e-14)
  expect_lt(abs(got_b[25] - ref_b[25]), 1e-14)
})
