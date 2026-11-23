## Round-5 native audit (corenative): regression tests for VA5-01, VA5-02 and VA5-08
## (round-5 verification, area A).
##
## Independent references, hard-coded (no package code): own Rmpfr implementations
## of the DEGPD carrier-1 cell probability -log(H(y+1)^kappa - H(y)^kappa) (1600 bits)
## and of the continuous carrier-3 density (1200 bits), with nested mpfr central
## differences, plus closed forms where noted.

c5_ld <- function(eta, y) {
  list(X = rep(list(matrix(1, 1, 1)), length(eta)), y = cbind(y),
       idpars = seq_along(eta), dupid = 0L, duplicate = 0L,
       offsets = lapply(eta, rep, 1), censored = FALSE, xi.max = Inf)
}
c5_fns <- function(name) get(name, envir = asNamespace("egpd"))
## largest error relative to the largest reference entry of the block
c5_err <- function(got, want) max(abs(got - want)) / max(abs(want))

test_that("VA5-01: large-count cells next to the xi = -1 endpoint keep exact derivatives", {
  ## Identity shape link, y >= 65536. At xi = -1 the density is constant, so the 2- and
  ## 3-point Gauss-Legendre values agree while the xi-derivatives (log singularity at
  ## the endpoint, here 0.01 or 0.1 beyond y + 1) did not: Hessians were 25% and 3%
  ## (per entry 16%) off.
  f <- c5_fns(".degpd1idfns")
  d12 <- function(eta, y) f$d120(c(0, 0, 0), c5_ld(eta, y))[1, ]
  ## columns: d/d(log sigma, xi, log kappa), then the Hessian (11, 12, 13, 22, 23, 33)
  ref1 <- c(1, 12.034268764487024, -0.99999221809752836, 0, -302460.19987296831, 1,
            -604897.1160851008, 0.99991342070728562, 7.781883069415913e-06)
  got1 <- d12(c(log(65537.01), -1, 0), 65536)
  expect_lt(c5_err(got1, ref1), 1e-8)
  ## closed forms at xi = -1, kappa = 1 (d0 = sigma - y, d1 = d0 - 1):
  ## score 1 + log(sigma) - d0 log(d0) + d1 log(d1), curvature 1 - sigma log(d0/d1)
  expect_lt(abs(got1[2] / 12.0342687644886 - 1), 1e-9)
  expect_lt(abs(got1[5] / -302460.199875122 - 1), 1e-8)
  ref2 <- c(2.4976020822690486, 15.978001537386993, -0.99999850000143142, 9090.9347049833086,
            -2388826.8654009281, 2.4999999652021052, -4786714.57923768, 2.4999790876793164,
            1.4999980476992041e-06)
  got2 <- d12(c(log((1 - 1e-9) * (1e6 + 1.1)), -1 + 1e-9, log(2.5)), 1e6)
  expect_lt(c5_err(got2, ref2), 1e-8)
})

test_that("VA5-02: the DEGPD d34 tiny-cell guard is limited to the upper tail", {
  blockerr <- function(got, want) c(c5_err(got[1:10], want[1:10]), c5_err(got[11:25], want[11:25]))
  d34 <- function(fns, eta, y) as.numeric(fns$d340(c(0, 0, 0), c5_ld(eta, y)))
  ## Lower tail, log link: y = 1, sigma = 1e8, xi = 0.2, kappa = 3 (H(y) = 1e-8, P = 7e-24).
  ## The guard mixed recomputed and plain symbols: third and fourth derivatives were
  ## 4.5e-4 and 2.5e-2 of their scale off (d3333 53.88 for 52.55).
  ref_a <- c(-3.857142772016321e-08, 6.4285711095510151e-09, 3.2460416213369751e-08,
             -6.4285712277142784e-09, -5.4100693006935718e-09, 2.9999999634234036,
             6.4285712867959101e-09, 5.4100693567932728e-09, 6.0960993921513653e-09,
             53.115989647246622, 3.8571426868897923e-08, -6.4285707905306172e-09,
             -3.2460415732900714e-08, 6.4285710268571414e-09, 5.4100691524140115e-09,
             3.6576596310170304e-08, -6.428571145020404e-09, -5.4100692646134143e-09,
             -6.0960993031063964e-09, 2.9999999584861294, 6.4285712041020356e-09,
             5.4100693207131153e-09, 6.096099370483705e-09, 6.9189784041307501e-09,
             52.545613373680908)
  got_a <- d34(c5_fns(".degpd1fns"), c(log(1 / 1e-8), log(0.2), log(3)), 1)
  expect_true(all(blockerr(got_a, ref_a) < 1e-7))
  expect_lt(abs(got_a[10] / ref_a[10] - 1), 1e-7)
  expect_lt(abs(got_a[25] / ref_a[25] - 1), 1e-7)
  ## identity link: y = 20, sigma = 2e5, xi = 0.2, kappa = 10 (H(y) = 1e-4, P = 6e-41;
  ## fourth derivatives were 1.3e-4 off)
  ref_b <- c(-0.00067761633752179049, 0.00056459592430476899, 0.00061737764202818702,
             1.0029123529050659e-07, -0.00051444605379699583, 9.9993802190744887,
             9.5102877628972121e-12, -4.4243039917986247e-08, 0.00051646631121836817,
             91.781971377579467, 0.0006775299348200333, -0.00056443957375130202,
             -0.0006173349715255566, -2.0056477869009708e-07, 0.00051437518323858844,
             0.00061975928743903698, -2.8529378618244496e-11, 8.8480756369983897e-08,
             -0.00051643048704935626, 9.9993755680577028, -2.7807840212975542e-15,
             8.2370182530331661e-12, -4.4678974357003672e-08, 0.00052034189369296068,
             91.706329568918846)
  got_b <- d34(c5_fns(".degpd1idfns"), c(log(20 / 1e-4), 0.2, log(10)), 20)
  expect_true(all(blockerr(got_b, ref_b) < 1e-6))
  ## The guard still applies where F(y) = H(y)^kappa is near 1 although H(y) = 1e-4 is
  ## small (kappa = 1e-8, P = 4.9e-10): without it the error is 1.2e-3.
  ref_c <- c(-6.1489809227930725e-05, 1.0245500187873145e-05, 6.1499484921373641e-13,
             -1.0246550949618664e-05, -1.024921329006834e-13, 9.9993849841242792e-09,
             1.0247076314339763e-05, 1.0249738945780668e-13, 1.0249914154322221e-13,
             9.1859467899257242e-08, 6.1488125359566198e-05, -1.0242418290118455e-05,
             -6.1495279727463054e-13, 1.0244519824296432e-05, 1.0247811614868484e-13,
             6.1499484921617545e-13, -1.024557054293077e-05, -1.0248862893924961e-13,
             -1.0249213290108984e-13, 9.9993849841242743e-09, 1.0246095878911752e-05,
             1.0249388535251875e-13, 1.0249738945821319e-13, 1.0249914154403526e-13,
             9.1859467899177899e-08)
  sig_c <- 0.2 * 20 / ((1 - 1e-4)^(-0.2) - 1)
  got_c <- d34(c5_fns(".degpd1fns"), c(log(sig_c), log(0.2), log(1e-8)), 20)
  expect_true(all(blockerr(got_c, ref_c) < 1e-5))
})

test_that("VA5-08: carrier-3 d34 next to a negative-shape endpoint with small delta", {
  blockerr <- function(got, want) c(c5_err(got[1:10], want[1:10]), c5_err(got[11:25], want[11:25]))
  d34 <- function(eta, y) as.numeric(c5_fns(".egpd3fns")$d340(c(0, 0, 0), c5_ld(eta, y))[1, ])
  ## The generated forms hold powers exp((delta - 4) A) (A = -log of the GPD survival)
  ## that underflow once (4 - delta) A > ~708, and the tail row only covered
  ## delta A > 40. A-opus point: 1 + xi y/psi = 1e-9, delta = 0.05 (delta A = 10.4):
  ## d222 was 30% and d2222 129 times its size off (wrong sign). The remaining error is
  ## the input rounding eps/(1 + xi y/psi) ~ 2e-7.
  ref_a <- c(-1.7999986644081807e+28, -1.7999986634081727e+29, 81916397069030.484,
             -1.7999986615081651e+30, 819163982910955.88, 1222065.1054727822,
             -1.7999986587081586e+31, 8191639936513039, 12220650.801476013,
             -0.066357466662647732, 5.3999964364138736e+37, 5.3999964344138613e+38,
             -1.2682718915517379e+23, 5.3999964306138495e+39, -1.2682718990624817e+24,
             -751074389408320.38, 5.3999964250138393e+40, -1.2682719057540616e+25,
             -7510743979467207, -8538400.2833855357, 5.3999964176138307e+41,
             -1.2682719116264778e+26, -75107440526305568, -85384001.064419925,
             0.14393260685631684)
  got_a <- d34(c(0, -0.1, log(0.05)), 9.99999999)
  expect_true(all(blockerr(got_a, ref_a) < 2e-6))
  expect_lt(abs(got_a[21] / ref_a[21] - 1), 2e-6)
  ## 1 + xi y/psi = 2^-36 exactly (xi = -1/8, y = 8 - 2^-33): no input rounding, so the
  ## row is exact (the stencil fallback was 5.6e-5 off)
  ref_e <- c(-4.5432539868551034e+33, -3.6346031894538591e+34, 5.0581487378774893e+17,
             -2.9076825515177526e+35, 4.0465189910211226e+18, 89891377.416159004,
             -2.326146041161009e+36, 3.2372151933191107e+19, 719131019.06814241,
             -0.073673488684243074, 9.3663026840817698e+44, 7.4930421472238775e+45,
             -5.705551705805668e+28, 5.9944337177167934e+46, -4.5644413649884377e+29,
             -4.2987882697277476e+18, 4.7955469741003272e+47, -3.6515530922335004e+30,
             -3.4390306162559791e+19, -592226374.16475081, 3.8364375792031665e+48,
             -2.9212424739551025e+31, -2.7512244933262777e+20, -4737810991.597621,
             0.18203728819037207)
  got_e <- d34(c(0, -0.125, log(0.05)), 8 - 2^-33)
  expect_true(all(blockerr(got_e, ref_e) < 1e-12))
  expect_lt(abs(got_e[25] - ref_e[25]), 1e-12)
  ## 1 + xi y/psi = 1e-12, delta = 0.02: the row was all zero
  ref_g <- c(-9.800087329468734e+37, -4.9000436647318672e+39, 1000002970354.1348,
             -2.4500218323634586e+41, 50000148551785.758, 681.58041980130565,
             -1.2250109161798791e+43, 2500007429226661.5, 34079.020989123645,
             -0.01846953658983113, 2.9400349318393742e+50, 1.4700174659191872e+52,
             -1.0000059407703586e+24, 7.3500873295909853e+53, -5.0000297041131044e+25,
             -52262203408749.805, 3.6750436647917925e+55, -2.5000148521847077e+27,
             -2613110171272846.5, -16707.130664701283, 1.8375218323934338e+57,
             -1.2500074261551816e+29, -1.3065550860370621e+17, -835356.53321198234,
             -0.017005684968537493)
  expect_true(all(blockerr(d34(c(0, -0.02, log(0.02)), (1 - 1e-12) / 0.02), ref_g) < 1e-3))
  ## positive shape far in the tail (xi = 0.1, A = 200, delta = 1e-3): garbage (1e26) before
  ref_f <- c(0.00025002090125466302, -0.0023344136404456479, -1.6600148746280128e-05,
             -2000.6293608840681, 0.043517074481600734, -0.0043351072994137936,
             -1086617.0673055393, 16.669043049271149, -0.82367038867568521,
             0.085706137182957345, 3.7522638192913139e-05, -0.00037536624441997371,
             -6.6252870127023173e-09, 0.023794497866141514, 0.0003307423289742585,
             -3.3067607610413149e-05, 60020.44075328928, -0.80916042558038515,
             0.037068225747463525, -0.0036737549671359395, 42903045.123518854,
             -473.75081754946757, 14.846156081984846, -0.6980134452702631,
             0.072483072558037473)
  expect_true(all(blockerr(d34(c(0, 0.1, log(1e-3)), expm1(20) / 0.1), ref_f) < 1e-12))
})
