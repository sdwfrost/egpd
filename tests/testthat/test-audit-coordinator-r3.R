## Coordinator integration tests (round 3): behaviour that spans several areas.

test_that("distfit(family = 'cpegpd') warns for off-lattice data with a coarse h", {
  set.seed(12)
  x <- rcpegpd(400, lambda = 3, sigma = 1, xi = 0.1, kappa = 1, h = 0)
  expect_warning(distfit(x, family = "cpegpd", cpegpd.h = 0.5),
                 "not on the cpegpd.h = 0.5 lattice")
  xg <- rcpegpd(400, lambda = 3, sigma = 1, xi = 0.1, kappa = 1, h = 0.2)
  w <- character()
  withCallingHandlers(distfit(xg, family = "cpegpd", cpegpd.h = 0.2),
                      warning = function(cnd) {
                        w <<- c(w, conditionMessage(cnd))
                        invokeRestart("muffleWarning")
                      })
  expect_false(any(grepl("lattice", w)))
})

test_that("distfit() validates cpegpd.h", {
  x <- c(0, 0.2, 0.4, 1.2)
  expect_error(distfit(x, family = "cpegpd", cpegpd.h = 0), "single positive finite")
  expect_error(distfit(x, family = "cpegpd", cpegpd.h = c(0.1, 0.2)), "single positive finite")
})

test_that("distfit() forwards npe.control to the neural fit", {
  captured <- NULL
  local_mocked_bindings(.distfit_neuralbayes = function(x, ...) {
    captured <<- list(...)
    structure(list(), class = "distfit")
  })
  Y <- matrix(rexp(20), 10, 2)
  distfit(Y, family = "mdgpd", method = "neuralbayes", npe.control = list(K = 600L))
  expect_identical(captured$npe.control, list(K = 600L))
  distfit(Y, family = "mdgpd", method = "neuralbayes")
  expect_identical(captured$npe.control, list())
})

test_that("zero-truncated Waring f0 matches dgwaring() inside the C++ bounds, even for tiny a", {
  ## (mu, k, rho) chosen so that a = mu (rho - 1) / k spans 1e-26 .. 1e3; the C++ bounds
  ## (cf_bounds_cpp) no longer floor a, so neither may the truncation term.
  pts <- rbind(c(3, 2, 4), c(1e-10, 1e10, 1 + 1e-6), c(1e-6, 1e8, 1.5),
               c(50, 1e12, 1e9), c(1e3, 0.5, 2.5))
  ## log p(0) = lgamma(a + rho) + lgamma(k + rho) - lgamma(rho) - lgamma(a + k + rho),
  ## evaluated with Rmpfr at 512 bits
  ref <- c(-1.395625688548775, -2.3603064948017382e-25, -9.1920953899868932e-14,
           -0.34543773864404959, -3.59509433491688)
  for (i in seq_len(nrow(pts))) {
    p <- pts[i, ]
    got <- egpd:::.zt_gwaring_logf0(as.list(log(p)), deriv = TRUE)
    expect_equal(got$value, ref[i], tolerance = 1e-10, label = sprintf("log f0 at point %d", i))
    expect_equal(got$value, dgwaring(0, mu = p[1], k = p[2], rho = p[3], log = TRUE),
                 tolerance = 1e-10, label = sprintf("vs dgwaring at point %d", i))
    num <- numDeriv::grad(function(e) egpd:::.zt_gwaring_logf0(as.list(e)),
                          log(p), method.args = list(r = 6))
    expect_equal(as.numeric(got$grad), num, tolerance = 1e-6 * max(1, abs(num)),
                 label = sprintf("gradient at point %d", i))
  }
})
