## Round-3 native audit: regression tests for findings N1-N7 and the bounded-xi
## (finite xi.max) derivative chain.
##
## Independent references:
## * SymPy derivatives of the closed-form NLLs evaluated at 40 digits
##   (derive_egpd_d34.py, mode "check");
## * a 320-bit Rmpfr implementation of the likelihoods (carrier 6 normaliser by
##   composite Gauss-Legendre quadrature), differentiated by mpfr central
##   differences with step 2^-60 (round-3 review reference library);
## * closed forms and Richardson-extrapolated differences of lower-order
##   analytic derivatives.

r3_ld <- function(eta, y, xi.max = Inf) {
  n <- length(y)
  list(X = rep(list(matrix(1, n, 1)), length(eta)), y = cbind(y),
       idpars = seq_along(eta), dupid = 0L, duplicate = 0L,
       offsets = lapply(eta, rep, n), censored = FALSE,
       xi.max = xi.max, bounded.xi = is.finite(xi.max))
}

r3_fns <- function(family, m, link = "log") {
  nm <- switch(family,
    degpd = if (link == "identity") sprintf(".degpd%didfns", m) else sprintf(".degpd%dfns", m),
    zidegpd = sprintf(".zidegpd%dfns", m),
    egpd = sprintf(".egpd%dfns", m))
  get(nm, envir = asNamespace("egpd"))
}

## NLL, gradient and packed Hessian of one observation.
r3_row <- function(fns, eta, y, xi.max = Inf) {
  ld <- r3_ld(eta, y, xi.max)
  p <- rep(0, length(eta))
  c(fns$d0(p, ld), fns$d120(p, ld)[1, ])
}

r3_unpack <- function(v, n) {
  H <- matrix(0, n, n)
  k <- 1
  for (i in 1:n) for (j in i:n) {
    H[i, j] <- H[j, i] <- v[k]
    k <- k + 1
  }
  H
}

r3_relerr <- function(got, want) max(abs(got - want) / pmax(1, abs(want)))

## Third then fourth derivatives in graded-lexicographic order (as .indices()),
## by Richardson-extrapolated central differences of the analytic Hessian
## (steps h for the third and h4 for the fourth derivatives).
r3_fd34 <- function(fns, eta, y, h = rep(2e-4, length(eta)), h4 = 5 * h) {
  n <- length(eta)
  Hf <- function(e) r3_unpack(fns$d120(rep(0, n), r3_ld(e, y))[1, -(1:n)], n)
  H0 <- Hf(eta)
  sh <- function(e, k, d) { e[k] <- e[k] + d; e }
  d1 <- function(k, s) (Hf(sh(eta, k, s)) - Hf(sh(eta, k, -s))) / (2 * s)
  d2 <- function(k, l, s, t) {
    if (k == l) return((Hf(sh(eta, k, s)) - 2 * H0 + Hf(sh(eta, k, -s))) / s^2)
    (Hf(sh(sh(eta, k, s), l, t)) - Hf(sh(sh(eta, k, s), l, -t)) -
       Hf(sh(sh(eta, k, -s), l, t)) + Hf(sh(sh(eta, k, -s), l, -t))) / (4 * s * t)
  }
  T3 <- lapply(1:n, function(k) (4 * d1(k, h[k] / 2) - d1(k, h[k])) / 3)
  Q <- list()
  for (k in 1:n) for (l in k:n)
    Q[[paste(k, l)]] <- (4 * d2(k, l, h4[k] / 2, h4[l] / 2) - d2(k, l, h4[k], h4[l])) / 3
  t3 <- t4 <- numeric(0)
  for (i in 1:n) for (j in i:n) for (k in j:n) t3 <- c(t3, T3[[k]][i, j])
  for (i in 1:n) for (j in i:n) for (k in j:n) for (l in k:n)
    t4 <- c(t4, Q[[paste(k, l)]][i, j])
  c(t3, t4)
}

## Gradient and packed Hessian by Richardson differences of the objective.
r3_fd12 <- function(f, b, h = 1e-3) {
  one <- function(h) {
    n <- length(b)
    f0 <- f(b)
    sh <- function(i, d) { z <- b; z[i] <- z[i] + d; z }
    g <- vapply(seq_len(n), function(i) (f(sh(i, h)) - f(sh(i, -h))) / (2 * h), 0)
    H <- numeric(0)
    for (i in 1:n) for (j in i:n) {
      if (i == j) {
        H <- c(H, (f(sh(i, h)) - 2 * f0 + f(sh(i, -h))) / h^2)
      } else {
        pp <- pm <- mp <- mm <- b
        pp[c(i, j)] <- pp[c(i, j)] + h
        mm[c(i, j)] <- mm[c(i, j)] - h
        pm[i] <- pm[i] + h; pm[j] <- pm[j] - h
        mp[i] <- mp[i] - h; mp[j] <- mp[j] + h
        H <- c(H, (f(pp) - f(pm) - f(mp) + f(mm)) / (4 * h^2))
      }
    }
    c(g, H)
  }
  (4 * one(h / 2) - one(h)) / 3
}

## ---------------------------------------------------------------------------
## N1: analytic third/fourth derivatives

test_that("N1: every analytic d340 matches Richardson differences of the analytic d120", {
  cases <- list(
    list("degpd", 1, "log", c(log(2), log(0.3), log(1.5)), c(0, 1, 5, 40)),
    list("degpd", 1, "log", c(log(5), log(0.08), log(0.6)), c(0, 3, 12)),
    list("degpd", 1, "identity", c(log(2), 0.3, log(1.5)), c(0, 2, 9)),
    list("degpd", 1, "identity", c(log(2), -0.3, log(1.5)), c(0, 3, 6)),
    list("egpd", 1, "identity", c(log(2), 0.3, log(1.5)), c(0.7, 3.3, 20)),
    list("egpd", 1, "identity", c(log(2), -0.2, log(1.5)), c(0.9, 5.1, 9)),
    list("egpd", 1, "identity", c(log(0.5), 0.15, log(0.7)), c(0.3, 2.2)),
    list("egpd", 3, "identity", c(log(2), -0.1, log(0.5)), c(1.1, 6)),
    list("egpd", 3, "identity", c(log(2), 0.3, log(1.5)), c(0.8, 3.3, 25)),
    list("egpd", 3, "identity", c(log(1), 0.2, log(3)), c(0.6, 4)))
  for (cs in cases) {
    fns <- r3_fns(cs[[1]], cs[[2]], cs[[3]])
    for (yy in cs[[5]]) {
      got <- fns$d340(rep(0, 3), r3_ld(cs[[4]], yy))[1, ]
      want <- r3_fd34(fns, cs[[4]], yy)
      lab <- sprintf("%s m%d %s y=%g", cs[[1]], cs[[2]], cs[[3]], yy)
      expect_true(all(is.finite(got)), label = lab)
      expect_lt(r3_relerr(got[1:10], want[1:10]), 1e-7, label = paste(lab, "third"))
      expect_lt(r3_relerr(got[11:25], want[11:25]), 5e-5, label = paste(lab, "fourth"))
    }
  }
})

test_that("N1: the (unwired) carrier-4 d34 matches differences of its d120", {
  fns <- list(d120 = egpd:::.egpd4.d12)
  for (eta in list(c(log(1.2), 0.15, log(2), log(3)), c(log(2), -0.2, log(0.5), log(0.8))))
    for (yy in c(0.7, 2.4, 6)) {
      got <- egpd:::.egpd4.d34(rep(0, 4), r3_ld(eta, yy))[1, ]
      want <- r3_fd34(fns, eta, yy)
      lab <- sprintf("egpd4 eta=(%s) y=%g", paste(signif(eta, 3), collapse = ","), yy)
      expect_lt(r3_relerr(got[1:20], want[1:20]), 1e-7, label = paste(lab, "third"))
      expect_lt(r3_relerr(got[21:55], want[21:55]), 5e-5, label = paste(lab, "fourth"))
    }
})

test_that("N1: continuous EGPD d34 equal high-precision symbolic values in every entry", {
  ## derive_egpd_d34.py check (SymPy, 40 digits)
  e1 <- c(-0.315207517180254, -1.48236212540321, 0.484842255229684, -1.08032161222525,
    -0.0204061140938435, 0.586927876905482, -0.295612567905477, -0.189876920488404,
    0.419658092747152, 0.455179749469019, -0.446975497767743, 1.95348106157858,
    0.0269373313124114, 3.94336593315959, -0.306484894755516, 0.484842255229684,
    3.60159421929206, -0.0632939323040065, -0.0204061140938435, 0.586927876905482,
    1.92889343952168, 0.175904153685651, -0.189876920488404, 0.419658092747152,
    0.455179749469019)
  got1 <- egpd:::.egpd1.d34(rep(0, 3), r3_ld(c(log(2), 0.3, log(1.5)), 3.3))[1, ]
  expect_lt(r3_relerr(got1, e1), 1e-10)
  e1b <- c(-25.9567254604911, -78.2128494464232, 0.752113675710409, -172.521244587648,
    0.730786894682249, 0.226913379590115, -326.50410015216, 0.672528652255469,
    0.356962676329325, 0.0429812432981471, 115.536949387113, 418.806116593561,
    0.678727561611548, 1160.62528799628, -1.20265003159164, 0.752113675710409,
    2682.68930005314, -3.54958999631641, 0.730786894682249, 0.226913379590115,
    5499.31480201193, -6.52494644223614, 0.672528652255469, 0.356962676329325,
    0.0429812432981471)
  got1b <- egpd:::.egpd1.d34(rep(0, 3), r3_ld(c(log(2), -0.2, log(1.5)), 5.1))[1, ]
  expect_lt(max(abs(got1b - e1b) / pmax(1, abs(e1b))), 1e-10)
  e3 <- c(-0.727203136142096, 0.14172111200384, 0.125491033574669, 0.526187746040974,
    0.0727979172953991, -0.118211241845129, 0.33802511059161, 0.027427393910811,
    -0.0331210241785059, 0.0408250653532041, 0.86416827890871, 0.650834046145556,
    -0.111510434585116, -0.663344312986872, -0.141987330035816, 0.097311701581534,
    -0.869914308190143, -0.0817263672373438, 0.0610167760166825, -0.0912100239798657,
    -0.550181281414325, -0.0347092368133962, 0.0234006557343821, -0.0255556862647384,
    0.162728529427466)
  got3 <- egpd:::.egpd3.d34(rep(0, 3), r3_ld(c(log(2), -0.1, log(0.5)), 1.1))[1, ]
  expect_lt(r3_relerr(got3, e3), 1e-10)
  ## The unwired carrier-4 and direct-parametrisation carrier-2 d34 share the
  ## generated GPD term; carrier 4 also had a wrong (2,2,3,3) entry.
  e4 <- c(-0.619313303200812, -2.10777883086302, 0.0998564280881686, 0.699262620748812,
    -1.51505746881694, 0.286192966195427, -0.0402676402211983, -0.157538393306598,
    -0.151050777572709, 0.828483703919403, -0.366308031802586, 0.323177793841897,
    -0.31851542277002, -0.11931835600958, -0.114404686220635, 0.627487125249554,
    0.144975935728785, -0.0314193755726376, -0.167438528097305, 0.592429602913584,
    -0.363653735477478, 3.39043546537949, 0.151175851140214, -0.0365270068259695,
    6.3758891717209, -0.196522355876419, -0.519688297765704, -0.111297046313664,
    0.000973801143429619, 0.699262620748812, 5.89858980239466, -0.565000153360545,
    -0.109442161311394, 0.0240693956138532, 0.104639869618372, -0.0402676402211983,
    0.106194287995575, -0.0868229206478572, -0.151050777572709, 0.828483703919403,
    3.23366307289504, -0.7681901436787, 0.327471040070281, 0.0729971813435892,
    0.131765254996374, -0.31851542277002, 0.0804307292672458, -0.0657590060315707,
    -0.114404686220635, 0.627487125249554, -0.00728640616574343, 0.0955900237253617,
    -0.0314193755726376, -0.167438528097305, 0.592429602913584)
  got4 <- egpd:::.egpd4.d34(rep(0, 4), r3_ld(c(0.4, 0.2, -0.3, 0.7), 2.5))[1, ]
  expect_lt(r3_relerr(got4, e4), 1e-10)
  e2 <- c(-0.72267925621419, -2.04898759027992, 0.0883112937682377, 0.257018597033236,
    -0.115856319565452, -1.60870600406219, 0.0476551556162599, -0.145992822634559,
    -0.00441912213030283, 0.252065356386732, -0.180121405902435, 0.0580617899207882,
    0.546661639692414, -0.184518825298012, -0.0155332590885064, -0.610956287772284,
    0.0263743049317915, -0.282740413374653, 0.0387674254279208, 0.190912343967482,
    -0.13642255442193, 0.0439755489116817, 0.414037281786016, -0.139753125731895,
    -0.0117647698380699, 0.188229815958732, -0.0131831758366561, 0.0287371060614694,
    -0.11190202971976, 0.011583870175798, -0.0252508887060886, 0.491530887423502,
    -0.139722861073594, 0.0144638259674192, -0.0789502495837408, -0.449276910886037,
    3.26636393640519, 0.187930944446197, -0.513262136140245, -0.00839677981536226,
    6.28939023186097, 0.0367047254388087, -0.461395553637277, 0.0843068332880099,
    0.28817060320552, -0.084886620693836, 0.0834887373159519, 0.0247491507174072,
    -0.202038384805046, 0.0144956038697293, 6.16739786522037, -0.0188231718456491,
    -0.039764766528816, 0.0341229403020706, 0.0448716132276141, 0.0596064775666623,
    0.0232951340259556, -0.357283787409331, -0.0260985346372286, 0.0216636201701714,
    0.330786132652329, -0.0349510831729722, 0.125010263975419, -0.0923947135444609,
    -0.0449841322836407, 0.0551579552872798, 0.317162023099663, -0.165387442401151,
    -0.0509898352639994, 0.0575837040522647, 3.9321674780587, -0.0344873904679379,
    0.299953401576592, -0.0269930481073742, -0.0536434864980104, 0.107763530222277,
    -0.00254125372543662, -0.460647365581182, 0.0443799539029924, 0.0218078952235904,
    0.250534848746548, -0.0264716790454818, 0.0946818033927851, -0.0699790387137294,
    -0.0340706325483179, 0.0417762071047637, 0.240216053944859, -0.125263164856131,
    -0.0386193053591648, 0.0436134503865773, 0.286638610008306, -0.0569243725986736,
    0.12408555822355, 0.032652613219454, -0.0372722531214503, 0.081247243029471,
    -0.0968160682964993, 0.00603304862332365, -0.0429428330643488, 0.0936081535769896,
    0.317069050752117, -0.120886261796685, 0.00753297162489303, -0.0536191837930267,
    0.0107901925988322)
  got2 <- egpd:::.egpd2.d34(rep(0, 5), r3_ld(c(0.4, 0.2, -0.3, 0.7, 0.4), 2.5))[1, ]
  expect_lt(r3_relerr(got2, e2), 1e-10)
})

test_that("N1: the analytic REML gradient is the gradient of the REML criterion (egpd m1, m3)", {
  ns <- asNamespace("egpd")
  set.seed(7); n <- 600; x <- runif(n); x2 <- runif(n)
  set.seed(11)
  y <- regpd(n, sigma = exp(1 + sin(2 * pi * x)), xi = 0.1 + 0.2 * x2, kappa = 1.5, type = 1)
  d <- data.frame(y = y, x = x, x2 = x2)
  for (m in c(1, 3)) {
    f <- suppressMessages(egpd(list(y ~ s(x, k = 8), ~ s(x2, k = 6), ~ 1), data = d,
                               family = "egpd", egpd.args = list(m = m), outer = "fd"))
    ld <- f$likdata; Sd <- f$Sdata
    lf <- get(sprintf(".egpd%dfns", m), envir = ns)
    expect_false(is.null(lf$d340))
    ldt <- ld; ldt$control$inner$gradtol <- 1e-10
    rho <- log(f$sp) + c(0.7, -0.5)
    r0 <- function(rr) as.numeric(ns$.reml0(structure(rr, beta = f$coefficients), lf, ldt, Sd))
    h <- 1e-4
    fdg <- vapply(seq_along(rho), function(k) {
      rp <- rm <- rho; rp[k] <- rp[k] + h; rm[k] <- rm[k] - h
      (r0(rp) - r0(rm)) / (2 * h)
    }, 0)
    bt <- attr(ns$.reml0(structure(rho, beta = f$coefficients), lf, ldt, Sd), "beta")
    an <- ns$.reml1(structure(rho, beta = bt), lf, ld, Sd)
    expect_lt(max(abs(an - fdg) / abs(fdg)), 2e-3, label = paste("REML gradient m", m))
  }
})

test_that("N1: analytic DEGPD d34 still requests finite differences where unsafe", {
  fns <- egpd:::.degpd1fns
  expect_error(fns$d340(rep(0, 3), r3_ld(c(log(3), log(0.3), log(1.5)), 2, xi.max = 1)),
               class = "egpd_d34_unavailable")
  expect_error(fns$d340(rep(0, 3), r3_ld(c(log(3), log(5e-4), log(1.5)), 2)),
               class = "egpd_d34_unavailable")
  expect_error(egpd:::.degpd1idfns$d340(rep(0, 3), r3_ld(c(log(3), 0.01, log(1.5)), 2)),
               class = "egpd_d34_unavailable")
})

## ---------------------------------------------------------------------------
## N3: carrier 6 (truncated beta) kappa derivatives

test_that("N3: continuous EGPD model 6 derivatives are exact for large kappa", {
  ## Rmpfr references (NLL, gradient, packed Hessian) at sigma = 1, xi = 0.3, y = 7
  refs <- list(
    `2` = c(4.5565279767622, -1.93337928529354, -2.78148693655061, -0.53299607472967,
      0.955760392097769, 4.38693717941475, 0.00420917134841157, 11.5652629568819,
      0.00940279802717438, -0.0421286970134673),
    `10` = c(3.71189418881201, -1.91654259989989, -2.74387574444192, -0.507825381096706,
      1.02640128122123, 4.51885501015415, 0.0210458567420579, 11.8208669670005,
      0.0470139901358719, 0.0171125220498461),
    `30` = c(3.16357784821625, -1.87445088641577, -2.64984776417017, -0.490188882849025,
      1.20300350402989, 4.84864958700267, 0.0631375702261736, 12.4598769922972,
      0.141041970407616, 0.0181421375954767),
    `50` = c(2.91581665779862, -1.83235917293166, -2.55581978389843, -0.479204532633221,
      1.37960572683854, 5.17844416385118, 0.105229283710289, 13.0988870175938,
      0.235069950679359, 0.0257949675166937),
    `100` = c(2.59128844636391, -1.72712988922137, -2.32074983321907, -0.454659299592245,
      1.82111128386018, 6.00293060597246, 0.210458567420579, 14.6964120808354,
      0.470139901358719, 0.0478406379124417),
    `200` = c(2.29068054542461, -1.51667132180079, -1.85060993186035, -0.407443628479827,
      2.70412239790345, 7.65190349021503, 0.420917134841157, 17.8914622073186,
      0.940279802717437, 0.0938063637078191),
    `500` = c(1.97193223444751, -0.885295619539055, -0.440190227784193, -0.26729657595732,
      5.35315574003327, 12.5988221429427, 1.05229283710289, 27.4766125867683,
      2.35069950679359, 0.233203423542681))
  fns <- egpd:::.egpd6fns
  for (k in names(refs)) {
    eta <- c(0, 0.3, log(as.numeric(k)))
    got <- r3_row(fns, eta, 7)
    expect_lt(max(abs(got - refs[[k]]) / pmax(1e-3, abs(refs[[k]]))), 1e-9,
              label = paste("model 6 vs Rmpfr, kappa =", k))
    fd <- r3_fd12(function(e) fns$d0(rep(0, 3), r3_ld(e, 7)), eta)
    expect_lt(r3_relerr(got[-1], fd), 1e-6, label = paste("model 6 vs FD, kappa =", k))
  }
  ## a cell with a negative shape and a large carrier penalty (y = 0.4, sigma = 2)
  ref2 <- c(172.844388599749, 117.391684814486, 12.0016476810062, 175.069958762286,
    -4.40475906710191, -12.6139684989057, 117.144071840355, -2.10781785372572,
    11.8738001520936, 175.571208754473)
  got2 <- r3_row(fns, c(log(2), -0.2, log(200)), 0.4)
  expect_lt(max(abs(got2 - ref2) / abs(ref2)), 1e-9)
  ## y = 0 boundary branch (Rmpfr: NLL, d/dlkappa, d2/dlkappa2; d/dlpsi = 1, others 0)
  ## and a vector of rows sharing kappa (cached normaliser)
  ref0 <- list(`2` = c(2.457657007340625, 3.688452591563491, 4.179319969279694),
               `50` = c(102.1317136794129, 105.0570121246958, 105.5620116248457),
               `500` = c(1051.013779172023, 1055.094869997333, 1055.595369996833))
  for (k in names(ref0)) {
    eta <- c(log(2), 0.25, log(as.numeric(k)))
    want0 <- c(ref0[[k]][1], 1, 0, ref0[[k]][2], 0, 0, 0, 0, 0, ref0[[k]][3])
    expect_lt(r3_relerr(r3_row(fns, eta, 0), want0), 1e-12,
              label = paste("model 6 y = 0, kappa =", k))
    yy <- c(0.2, 3, 0, 11)
    rows <- fns$d120(rep(0, 3), r3_ld(eta, yy))
    one <- t(vapply(yy, function(v) fns$d120(rep(0, 3), r3_ld(eta, v))[1, ], numeric(9)))
    expect_equal(rows, one, tolerance = 1e-14)
  }
})

test_that("N3: model 6 with kappa = Inf returns the out-of-support row promptly", {
  ## lkappa > 709.78 overflows kappa; the normaliser quadrature used to bisect a
  ## NaN integrand to depth 40 (the unfixed code hangs here as well).
  fns <- egpd:::.egpd6fns
  elapsed <- system.time({
    for (lk in c(709.8, 710, 800)) for (yy in c(0, 1e-6, 7)) for (xi in c(0.3, 1e-5)) {
      ld <- r3_ld(c(0, xi, lk), yy)
      expect_equal(fns$d0(rep(0, 3), ld), 1e20)
      expect_equal(unname(fns$d120(rep(0, 3), ld)[1, ]), rep(0, 9))
    }
  })[["elapsed"]]
  expect_lt(elapsed, 5)
})

## ---------------------------------------------------------------------------
## N4: rows next to the xi < 0 upper endpoint

test_that("N4: derivatives next to the xi < 0 endpoint never difference the support penalty", {
  ## Per-entry tolerances: next to the endpoint 1 + xi*y/sigma = g carries an input
  ## rounding of ~eps/g relative, so the g = 1e-9 cells are checked to 1e-6 and the
  ## g = 1e-6 cell to 1e-8, entry by entry (NLL, gradient and every Hessian entry).
  cells <- list(
    list(3, c(log(0.05), -0.001, log(1e-6)), 50 * (1 - 1e-9),
         c(20689.6181103383, -998952183395.473, -998952162674447, 0.0103248453340837,
           9.98954452917271e+20, 9.98954451918022e+23, -496546.142021046,
           9.98954449918693e+26, -496546131.731112, 0.0102890582914891), 1e-6),
    list(4, c(log(0.05), -0.45, log(1e3), log(200)), 0.05 / 0.45 * (1 - 1e-9),
         c(17.7265340594398, -1222222216.38042, -2716049269.61938, 0.000999000999000999,
           -1, 1.22222221376085e+18, 2.71604935897472e+18, -2.20000000576307e-12,
           2.22444445027155e-09, 6.03566522515644e+18, -4.88888880038197e-12,
           4.94320978705288e-09, -0.000998002996004994, -1.00000000658105e-21,
           1.00100000658763e-18), 1e-6),
    list(1, c(0, -0.3, log(1.5)), (1 / 0.3) * (1 - 1e-6),
         c(31.8307261938105, -2333330.00000198, -7777616.49433375, -1, 2333331000003.97,
           7777758888913.23, 4.99999499999007e-14, 25925800001167.5, 1.66664197414576e-13,
           1.49999999999575e-20), 1e-8))
  for (cs in cells) {
    got <- r3_row(r3_fns("egpd", cs[[1]]), cs[[2]], cs[[3]])
    ref <- cs[[4]]
    bad <- which(!(abs(got - ref) <= cs[[5]] * abs(ref) + 1e-300))
    expect_equal(bad, integer(0),
                 label = sprintf("endpoint model %d entries off by more than %g", cs[[1]], cs[[5]]))
  }
  ## d34 stencils shrink inside the support: compare with differences of d12
  d34_cells <- list(
    list(3, c(log(0.05), -0.001, log(1e-6)), 50 * (1 - 1e-9), c(1e-11, 1e-14, 1e-4)),
    list(1, c(0, -0.3, log(1.5)), (1 / 0.3) * (1 - 1e-6), c(1e-8, 1e-8, 1e-4)),
    list(3, c(0, -0.3, log(2)), (1 / 0.3) * (1 - 1e-7), c(1e-9, 1e-9, 1e-4)))
  for (cs in d34_cells) {
    fns <- r3_fns("egpd", cs[[1]])
    got <- fns$d340(rep(0, 3), r3_ld(cs[[2]], cs[[3]]))[1, ]
    want <- r3_fd34(fns, cs[[2]], cs[[3]], h = cs[[4]])
    expect_true(all(is.finite(got)))
    expect_lt(max(abs(got[1:10] - want[1:10])) / max(abs(want[1:10])), 1e-4,
              label = paste("endpoint d34 model", cs[[1]], "y =", cs[[3]]))
  }
})

## ---------------------------------------------------------------------------
## N5: small y / sigma

test_that("N5: continuous EGPD d12 is exact for y / sigma -> 0 (no stencil noise)", {
  cells <- list(
    list(1, c(0, 0.3, log(5)), c(53.6526082194223, 4.99999610000147, 2.99999770000099e-06,
      68.0775560398207, 3.89999705333479e-06, -2.99999540000297e-06, 4.99999675000135,
      -2.66666316666897e-12, 2.49999812500088e-06, 69.0775560398207)),
    list(3, c(0, 0.3, log(2)), c(12.7169007192956, 1.99999755000106, 1.49999857500069e-06,
      -0.66666566666715, 2.44999787833443e-06, -1.49999715000207e-06, -9.99999033333723e-07,
      -1.4166645083349e-12, -4.99999466666917e-13, -0.222221222223039)),
    list(4, c(0, 0.3, log(1), log(2)), c(13.122365327404, 1.99999805000066,
      1.49999882500047e-06, -0.499999500000158, 26.6310224159283, 1.94999867833396e-06,
      -1.4999976500014e-06, -4.99999683333453e-07, 1.99999870000054, -1.41666484166778e-12,
      -2.49999816666746e-13, 9.99999250000351e-07, -0.249999500000242, -0.4999996666668,
      27.6310224159283)),
    list(5, c(0, 0.3, log(5)), c(1.89539460070142, 1.00000369998889, 1.00000169999429e-06,
      2.07512423754119, -3.6999777800445e-06, -1.00000339998287e-06, 4.99998850001495e-06,
      -1.00000206665597e-12, 2.49999400000809e-12, 2.34409212260222)),
    list(6, c(0, 0.3, log(5)), c(7.60869358342085, 1.00005676354011, 1.00002823176725e-06,
      10.0310757138081, -5.67625641141904e-05, -1.00005646304369e-06, 7.25794246549417e-05,
      -1.00003744210628e-12, 3.62897086985002e-11, 10.5797659118753)),
    list(2, c(0, 0.3, log(0.7), 0.3, 0.2), c(-3.18983826815067, 0.699998920757698,
      8.49999325378931e-07, 8.67085785814035, 3.43286276035458e-07, -0.450165983605966,
      1.07023058390957e-06, -8.49998655263802e-07, 0.699999528038613, -4.37630304038209e-07,
      -2.57573886524955e-08, -8.74998920258375e-13, 3.49999729019335e-07,
      -2.18815130137592e-13, -1.28786930383785e-14, 9.67085783300944, -2.21767555468629e-07,
      -1.2565453434942e-08, -5.82432274487306e-06, -3.43286269485022e-07, 0.247516553630305)))
  for (cs in cells) {
    got <- r3_row(r3_fns("egpd", cs[[1]]), cs[[2]], 1e-6)
    ## every entry, including the O(1e-6) and O(1e-12) curvatures, to 1e-9 absolute
    expect_lt(max(abs(got - cs[[3]])), 1e-9, label = paste("small y model", cs[[1]]))
  }
})

## ---------------------------------------------------------------------------
## N6: value-only objective

test_that("N6: DEGPD/ZIDEGPD d0 does not pay for derivatives it discards", {
  skip_on_cran()
  set.seed(3)
  n <- 20000
  y <- rpois(n, 4)
  car <- list(log(1.5), c(log(0.7), log(1.2), 0.2), log(2), c(log(1.3), log(2.5)), log(3))
  t0 <- t12 <- 0
  for (m in 1:5) {
    for (fam in c("degpd", "zidegpd")) {
      eta <- c(log(3), log(0.2), car[[m]], if (fam == "zidegpd") -1)
      fns <- r3_fns(fam, m)
      ld <- r3_ld(eta, y)
      p <- rep(0, length(eta))
      t0 <- t0 + system.time(for (i in 1:3) fns$d0(p, ld))[["elapsed"]]
      t12 <- t12 + system.time(for (i in 1:3) fns$d120(p, ld))[["elapsed"]]
    }
  }
  ## d0 used to evaluate the full second-order Jets and was as slow as d120
  expect_lt(t0, 0.6 * t12)
})

## ---------------------------------------------------------------------------
## N2: zero-inflated y = 0 with tiny probabilities

test_that("N2: ZIDEGPD y = 0 NLL is non-negative and accurate for tiny probabilities", {
  fns <- egpd:::.zidegpd1fns
  ## closed form, carrier 1: P(Y = 0) = H(1)^kappa, H(1) = 1 - (1 + xi/sigma)^(-1/xi);
  ## NLL = -log(pi + (1 - pi) P0) = -log(1 - (1 - pi) (1 - P0))
  want_nll <- function(sigma, xi, kappa, lp) {
    lP0 <- kappa * log1p(-(1 + xi / sigma)^(-1 / xi))
    la <- plogis(lp, log.p = TRUE)
    lb <- plogis(-lp, log.p = TRUE) + lP0
    lsum <- max(la, lb) + log1p(exp(-abs(la - lb)))
    if (lsum < -0.5) -lsum else -log1p(-plogis(-lp) * -expm1(lP0))
  }
  for (sigma in c(0.05, 0.2, 3, 100)) for (xi in c(0.01, 0.3)) for (kappa in c(0.01, 1.5, 200))
    for (lp in c(qlogis(1e-9, lower.tail = FALSE), 8, 0, -5, -1000)) {
      want <- want_nll(sigma, xi, kappa, lp)
      got <- fns$d0(rep(0, 4), r3_ld(c(log(sigma), log(xi), log(kappa), lp), 0))
      lab <- sprintf("sigma=%g xi=%g kappa=%g logitpi=%g", sigma, xi, kappa, lp)
      expect_gte(got, 0, label = lab)
      expect_lt(abs(got - want), 1e-12 * want + 1e-300, label = lab)
    }
  ## the reviewer's cell: pi = 1 - 1e-9, sigma = 0.05
  lp <- qlogis(1e-9, lower.tail = FALSE)
  got <- fns$d0(rep(0, 4), r3_ld(c(log(0.05), log(1e-3), log(1e-3), lp), 0))
  expect_gte(got, 0)
  expect_equal(got, want_nll(0.05, 1e-3, 1e-3, lp), tolerance = 1e-10)
})

test_that("N2: ZIDEGPD y = 0 stays finite when P(0) rounds to (or above) 1", {
  ## sigma ~ 0.01 puts P(Y = 0) within an ulp of 1; the DEGPD cell NLL
  ## L = -log P(0) can then be -1e-16 (mixture log-sum-exp, quadrature).
  ## Exact reproducers first (DEGPD carrier 6 alone gives d0 = -1.1e-16 there).
  cells <- list(
    list(6, c(log(0.01), log(0.05), -1.63, -3), TRUE),
    list(2, c(-5.4704412458062963, -3.2799074057458482, -0.79928845562972128,
              1.6487807277590036, -1.9945460809394717, -2.181596415117383), TRUE))
  car <- list(`2` = c(log(0.45), log(5.2), -1.99), `5` = log(3), `6` = c(-1.63))
  for (m in names(car)) for (sigma in c(0.005, 0.01, 0.02)) for (xi in c(0.05, 0.3))
    for (lp in c(-3, 0, 2))
      cells[[length(cells) + 1L]] <-
        list(as.integer(m), c(log(sigma), log(xi), car[[m]], lp), FALSE)
  cells[[length(cells) + 1L]] <- list(6L, c(log(0.01), log(0.05), log(40), -1), FALSE)
  for (cs in cells) {
    fns <- r3_fns("zidegpd", cs[[1]])
    eta <- cs[[2]]
    lab <- sprintf("zidegpd%d eta=(%s)", cs[[1]], paste(signif(eta, 4), collapse = ","))
    v <- fns$d0(rep(0, length(eta)), r3_ld(eta, 0))
    expect_true(is.finite(v) && v >= 0, label = paste(lab, "d0"))
    ## the exact reproducers have 1 - P(0) below double precision
    if (cs[[3]]) expect_lt(v, 1e-14, label = paste(lab, "d0 size"))
    g <- fns$d120(rep(0, length(eta)), r3_ld(eta, 0))[1, ]
    expect_true(all(is.finite(g)), label = paste(lab, "d12"))
    fd <- r3_fd12(function(e) fns$d0(rep(0, length(e)), r3_ld(e, 0)), eta)
    expect_lt(r3_relerr(g, fd), 1e-6, label = paste(lab, "d12 vs d0"))
  }
})

## ---------------------------------------------------------------------------
## N7: ginv.egpd relative tolerance

test_that("N7: ginv.egpd uses a MASS::ginv-style relative tolerance", {
  svd_ginv <- function(x, tol = sqrt(.Machine$double.eps)) {
    s <- svd(x)
    keep <- s$d > max(tol * s$d[1], 0)
    if (!any(keep)) return(matrix(0, ncol(x), nrow(x)))
    s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
  }
  set.seed(3)
  A <- matrix(rnorm(20), 5, 4) %*% diag(c(1, 1, 1, 0)) %*% matrix(rnorm(16), 4, 4)
  mats <- list(A = A, tA = t(A), small = 1e-12 * A, tiny_wide = 1e-30 * t(A),
               rect = matrix(rnorm(12), 3, 4), big = 1e12 * A, zero = matrix(0, 3, 2))
  for (nm in names(mats)) {
    x <- mats[[nm]]
    want <- svd_ginv(x)
    got <- ginv.egpd(x)
    expect_equal(dim(got), rev(dim(x)), label = nm)
    expect_lt(max(abs(got - want)), 1e-10 * max(1, max(abs(want))), label = nm)
  }
  ## scale invariance: ginv(c * A) = ginv(A) / c
  expect_equal(ginv.egpd(1e-12 * A) * 1e-12, ginv.egpd(A), tolerance = 1e-10)
  ## an explicit relative tolerance drops the smaller singular values
  D <- diag(c(1, 1e-3, 1e-9))
  expect_equal(ginv.egpd(D, tol = 1e-6), diag(c(1, 1e3, 0)))
  expect_equal(ginv.egpd(D, tol = 0), diag(c(1, 1e3, 1e9)))
  expect_error(ginv.egpd(D, tol = -1), "non-negative")
  expect_error(ginv.egpd(D, tol = c(1e-8, 1e-6)), "non-negative")
})

## ---------------------------------------------------------------------------
## Bounded shape link (finite xi.max) and independent d12 anchors

test_that("finite xi.max: DEGPD d120 is the derivative of d0 through the bounded link", {
  car <- list(log(1.5), c(log(0.7), 0.3, 0.2), log(2), c(log(1.3), log(2.5)), log(3), log(4))
  for (m in 1:6) for (xm in c(0.5, 3)) for (lt in c(log(1e-6), log(0.05), 0, log(5), log(800)))
    for (yy in c(0, 2, 15)) {
      fns <- r3_fns("degpd", m)
      eta <- c(log(3), log(xm) + lt, car[[m]])
      got <- fns$d120(rep(0, length(eta)), r3_ld(eta, yy, xm))[1, ]
      fd <- r3_fd12(function(e) fns$d0(rep(0, length(eta)), r3_ld(e, yy, xm)), eta)
      expect_lt(r3_relerr(got, fd), 1e-6,
                label = sprintf("m%d xi.max=%g t=%g y=%g", m, xm, exp(lt), yy))
    }
  ## Rmpfr anchors with xi = 0.8 (1 - exp(-exp(eta)/0.8))
  a1 <- c(2.22927447261545, -0.321558238061801, 0.0876937141107296, -0.590941825682583,
    1.20798158970313, 0.127117304071914, 0.560091541947341, 0.0555358171587568,
    0.0878705181967553, 0.406869509429512)
  expect_lt(r3_relerr(r3_row(egpd:::.degpd1fns, c(log(2), log(0.3), log(1.5)), 3, 0.8), a1), 1e-10)
  a3 <- c(2.20038834637635, -0.286811690969046, 0.107251907421801, 0.133978858088739,
    1.12676091568276, 0.14915976216513, -0.28941133220068, 0.0425248963717695,
    -0.0559926396164471, 0.147085201002854)
  expect_lt(r3_relerr(r3_row(egpd:::.degpd3fns, c(log(2), log(0.5), log(2)), 3, 0.8), a3), 1e-10)
})

test_that("d12 agree with independent Rmpfr values for all three native families", {
  cells <- list(
    list("egpd", 1, "identity", 3.3, c(log(2), 0.3, log(1.5)),
         c(2.18195555192627, -0.239139983060491, 0.454425761649561, -0.544820250530981,
           1.12133489356864, 0.473058366214756, 0.586927876905482, -0.080828759988928,
           0.419658092747152, 0.455179749469019)),
    list("egpd", 3, "identity", 1.1, c(log(2), -0.1, log(0.5)),
         c(1.5045949510447, 1.36636415415663, 0.668353278295252, -0.198565674544378,
           0.663767491530096, -0.500636572267439, -0.131820566746258, -0.322881189643891,
           -0.0369341537258122, -0.0940950708485458)),
    list("egpd", 5, "identity", 2.5, c(log(1.5), 0.2, log(3)),
         c(1.80557671611696, -0.288824319839477, 0.467098934302735, -0.247345630552943,
           1.49455744028091, 0.758907263792493, 0.211175680160522, 0.0737151542934549,
           0.159150745597258, -0.111999715931122)),
    list("degpd", 1, "log", 3, c(log(2), log(0.3), log(1.5)),
         c(2.25047616608, -0.292628541282343, 0.127611948421185, -0.569933468304204,
           1.12502160874999, 0.163992826696916, 0.560175892853725, 0.126421914007999,
           0.125119542353312, 0.427878713681593)),
    list("degpd", 1, "identity", 5, c(log(2), -0.2, log(1.5)),
         c(3.46777805138977, -3.72199843533333, -4.19852118611377, -0.968555199414352,
           10.4911765821932, 22.570629983993, 0.178490278293703, 40.9021728691909,
           0.29805183578818, 0.0313584883013728)),
    list("degpd", 6, "log", 4, c(log(3), log(0.2), log(40)),
         c(3.33170922082147, 5.56442711812014, 0.896782851505575, 2.04730331547743,
           10.8233094622405, 0.729317090352167, 5.85244553098123, 0.915334794141884,
           0.829819937172583, 2.40598989089438)),
    list("zidegpd", 1, "log", 0, c(log(2), log(0.3), log(1.5), -0.5),
         c(0.655840842037801, 0.299563438700435, 0.0214212497392835, 0.403845687792471,
           -0.349890322220347, -0.150250491882245, -0.0311904642793489, -0.0233005304643854,
           -0.217911729086765, 0.0182159435843899, -0.00166617957218241, -0.0155824809267016,
           -0.0314117730526154, -0.293769868889422, 0.0367285678772481)),
    list("zidegpd", 5, "log", 2, c(log(2), log(0.3), log(4), 0.3),
         c(2.65324468806426, 0.271655234133447, 0.19171407942841, -0.137041554118953,
           0.574442516811659, 1.29178627492103, 0.0592579643765743, 0.433968591466336, 0,
           0.169882228841893, 0.0720637053695402, 0, 0.0646561176145799, 0,
           0.244458311690746)),
    list("egpd", 2, "identity", 1.7, c(log(1.5), 0.2, log(0.8), 0.4, -0.3),
         c(1.52525038336013, 0.303196711421792, 0.653499464953966, -0.209898153383334,
           0.00996077211155615, 0.090764240771412, 1.02877948194947, -0.0368002195052183,
           0.303437986731656, 0.513270690737872, -0.172705361328366, -0.386284022061888,
           0.160254761333342, 0.271073746993775, -0.091210915148675, 0.111585019890639,
           0.148977002969955, -0.144967194211097, 0.291768240874501, -0.003334799192204,
           0.0217515844418574)),
    list("egpd", 4, "identity", 2.4, c(log(1.2), 0.15, log(2), log(3)),
         c(1.56295212106982, -0.410860298539079, 0.461338823591192, 0.173558639895271,
           -0.551853107359685, 1.90600699545773, 1.38121017031451, -0.33567383506367,
           0.787149968864987, 0.575027639370112, -0.306384754879982, 0.718467527320799,
           0.0963848696847654, -0.151937313928674, 0.448146892640315)),
    list("egpd", 6, "identity", 1.3, c(log(0.9), -0.15, log(4)),
         c(0.663530085307231, -0.179430923685937, 0.706969473613397, -0.392608036035172,
           2.98654260827458, 1.48728027374636, 0.517259950451848, 0.241922729578371,
           0.403925073954532, 0.164025840078063)),
    list("degpd", 2, "log", 3, c(log(2), log(0.25), log(0.8), 0.4, -0.3),
         c(2.24713778665865, -0.182325074060285, 0.132930358973089, -0.313317375635261,
           -0.176699513622784, 0.145827419190627, 1.26800770806351, 0.142678707868531,
           0.225291783426153, 0.427287004190713, -0.11220917216835, 0.128846112307506,
           0.0429534596730772, 0.0815641535447621, -0.0213508515250992, -0.0315122385717157,
           0.194941489851992, -0.130916905212334, 0.114806450519444, 0.0494281662540663,
           0.042977156397197)),
    list("degpd", 3, "log", 4, c(log(2.5), log(0.2), log(1.5)),
         c(2.35175828785005, -0.357815738620152, 0.086712863048487, 0.141835323894422,
           1.47755682606705, 0.172493942051156, -0.341162777550066, 0.0936733554531907,
           -0.0548454625635344, 0.157880390915923)),
    list("degpd", 4, "log", 2, c(log(1.5), log(0.3), log(2), log(3)),
         c(1.65759435721486, 0.0693239348483094, 0.197731322920512, 0.0512507793143031,
           -0.25865473434831, 1.50186711293065, 0.171540996908395, -0.39271500052552,
           0.971631437753363, 0.190960153215195, -0.0841290554117299, 0.207449022335023,
           0.183997054625995, -0.236531792904616, 0.727939013184982)),
    list("degpd", 5, "log", 6, c(log(3), log(0.15), log(4)),
         c(2.7862616010445, -0.719237488228538, 0.0276264986836903, -0.338875164544301,
           1.79767300218955, 0.249719575109248, 0.155197470434925, 0.0511107087408832,
           0.0226878152079755, -0.134544500649191)),
    list("degpd", 2, "identity", 2, c(log(2), -0.2, log(0.8), 0.4, -0.3),
         c(1.59410287378895, 0.133862255318096, 0.816128899845733, -0.306745850704363,
           -0.165121396404102, 0.142185631296946, 1.96561465999783, 0.412441758841702,
           0.30740298097125, 0.577551998181712, -0.154938891685498, -0.496210691580051,
           0.210050861868493, 0.396353497377481, -0.105165317278067, -0.0233155818426954,
           0.189883590646539, -0.131275365299013, 0.121783250311847, 0.0467907558859235,
           0.0413860662437097)),
    list("degpd", 3, "identity", 4, c(log(2.5), 0.2, log(1.5)),
         c(2.35175828785005, -0.357815738620152, 0.433564315242435, 0.141835323894422,
           1.47755682606705, 0.86246971025578, -0.341162777550066, 0.174012310117591,
           -0.274227312817672, 0.157880390915923)),
    list("degpd", 4, "identity", 3, c(log(1.5), -0.15, log(2), log(3)),
         c(2.0636846345491, -1.77513353184871, -0.865428731490404, 0.293775834026844,
           -0.851516358918839, 4.93061013712805, 7.23604648004991, -0.177463875025032,
           0.509335152890833, 8.84723573275583, -0.224666679399435, 0.650166539199184,
           -0.107596233940444, -0.0501225631797806, 0.146673923720034)),
    list("degpd", 5, "identity", 1, c(log(3), 0.25, log(4)),
         c(1.96443151264123, 1.12522066814408, 0.498730612059829, 0.401677684753454,
           0.485254364005853, -0.336586763173224, 0.669182677198833, -0.201242082619498,
           0.173639332137034, 0.587958495859349)),
    list("degpd", 6, "identity", 5, c(log(3), -0.1, log(20)),
         c(1.67446360944696, 0.333953855544221, 1.35141596038913, -0.186823266604683,
           6.77924069574513, 5.11537516447256, 1.39541793307795, 3.32797180440476,
           1.3556769127725, 0.320090546722306)),
    list("zidegpd", 2, "log", 3, c(log(2), log(0.25), log(0.8), 0.4, -0.3, -1),
         c(2.56039947417687, -0.182325074060285, 0.132930358973089, -0.313317375635261,
           -0.176699513622784, 0.145827419190627, 0.268941421369995, 1.26800770806351,
           0.142678707868531, 0.225291783426153, 0.427287004190713, -0.11220917216835, 0,
           0.128846112307506, 0.0429534596730772, 0.0815641535447621, -0.0213508515250992, 0,
           -0.0315122385717157, 0.194941489851992, -0.130916905212334, 0, 0.114806450519444,
           0.0494281662540663, 0, 0.042977156397197, 0, 0.196611933241482)),
    list("zidegpd", 3, "log", 0, c(log(2.5), log(0.2), log(1.5), 0.2),
         c(0.50409257309392, 0.137299475900641, 0.00535111737159323, -0.0387009069874884,
           -0.360406623191176, -0.155085470940073, -0.0112599029676616, 0.0412318981575596,
           -0.124975560138626, 0.00477674915980571, 0.00160697428047248, -0.00487080439670691,
           -0.0230532230641761, 0.0352271375903456, 0.165813939423083)),
    list("zidegpd", 4, "log", 2, c(log(1.5), log(0.3), log(3), log(2), -0.5),
         c(2.13167134139497, 0.0693239348483094, 0.197731322920512, -0.25865473434831,
           0.0512507793143031, 0.377540668798145, 1.50186711293065, 0.171540996908395,
           0.971631437753363, -0.39271500052552, 0, 0.190960153215195, 0.207449022335023,
           -0.0841290554117299, 0, 0.727939013184982, -0.236531792904616, 0,
           0.183997054625995, 0, 0.235003712201594)),
    list("zidegpd", 6, "log", 0, c(log(3), log(0.2), log(40), -1.5),
         c(1.70141327797209, 2.31062967767153e-10, 7.53642790245228e-12, 2.68659731060786e-10,
           -0.817574476182983, -4.96059223163387e-09, -1.6917238823181e-10,
           -5.59515998151531e-09, -2.3106296776469e-10, 1.85825049803499e-12,
           -1.82493630246581e-10, -7.53642790237194e-12, -6.50705555481874e-09,
           -2.68659731057922e-10, 0.149146452059673)))
  for (cs in cells) {
    got <- r3_row(r3_fns(cs[[1]], cs[[2]], cs[[3]]), cs[[5]], cs[[4]])
    expect_lt(r3_relerr(got, cs[[6]]), 1e-10,
              label = sprintf("%s m%d %s y=%g", cs[[1]], cs[[2]], cs[[3]], cs[[4]]))
  }
})
