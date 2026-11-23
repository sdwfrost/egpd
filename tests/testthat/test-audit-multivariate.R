test_that("MDGPD realizes the original-author spectral CDF and conditional margins", {
  # Aka--Kratz--Naveau v1, Theorem 2.1 / Eq. (7), Definition 2.1.
  set.seed(416)
  n <- 80000L
  tmat <- matrix(rpois(2 * n, 0.8), ncol = 2)
  spectral <- tmat - apply(tmat, 1, max)
  e <- rexp(n)
  expected <- ceiling(2 * expm1(0.2 * (spectral + e)) / 0.2)
  set.seed(416)
  y <- rmdgpd(n, 2, 0.2, rate = 0.8)
  expect_equal(unname(y), expected, ignore_attr = TRUE)
  k <- c(2, 4)
  h <- log1p(0.2 * k / 2) / 0.2
  cdf <- 1 - mean(pmin(1, exp(pmax(spectral[, 1] - h[1],
                                    spectral[, 2] - h[2]))))
  expect_equal(mean(y[, 1] <= k[1] & y[, 2] <= k[2]), cdf,
               tolerance = 0.008)
  pos <- y[y[, 1] > 0, 1]
  expect_equal(mean(pos > 3), (1 + 0.2 * 3 / 2)^(-1 / 0.2),
               tolerance = 0.008)
  expect_true(any(y[, 1] == 1)) # no artificial even-only lattice
})

test_that("standard maximum has positive, not zero-based, geometric support", {
  set.seed(245)
  y <- rmdgpd(60000, 1, 0, rate = 0.7, d = 3)
  maximum <- apply(y, 1, max)
  expect_true(all(maximum >= 1))
  expect_equal(mean(maximum > 2), exp(-2), tolerance = 0.008)
  positive <- y[y[, 1] > 0, 1]
  expect_equal(mean(positive > 1), exp(-1), tolerance = 0.008)
  expect_gt(mean(y[, 1] <= 0), 0.1)
})

test_that("positive compound likelihoods normalize using the actual zero atom", {
  for (family in c("cpegpd", "cpdegpd")) {
    h <- if (family == "cpegpd") 0.2 else 1
    args <- list(family = family, lambda = 1, sigma = 1, xi = 0,
                 kappa = 1, K = 400)
    mass <- exp(do.call(egpd:::.aggregate_positive_logdens,
                        c(list(x = (1:400) * h), args)))
    expect_equal(sum(mass), 1, tolerance = 1e-8)
    expect_equal(do.call(egpd:::.aggregate_positive_logdens,
                        c(list(x = 0), args)), -Inf)
    if (family == "cpegpd") {
      expect_equal(do.call(egpd:::.aggregate_positive_logdens,
                           c(list(x = 0.1), args)),
                   do.call(egpd:::.aggregate_positive_logdens,
                           c(list(x = 0.2), args)))
    }
    p <- c(0, 0.1, 0.5, 0.9, 1)
    q <- do.call(egpd:::.aggregate_positive_quantile, c(list(p = p), args))
    expect_equal(q[c(1, 5)], c(h, Inf))
    expect_true(all(q > 0))
    f <- get(paste0("p", family))
    p0 <- f(0, lambda = 1, sigma = 1, xi = 0, kappa = 1, K = 400)
    fq <- f(q[2:4], lambda = 1, sigma = 1, xi = 0, kappa = 1, K = 400)
    expect_true(all((fq - p0) / (1 - p0) >= p[2:4] - 1e-12))
  }
  expect_equal(dcpegpd(0, 1, sigma = 1, xi = 0, kappa = 1),
               exp(-1), tolerance = 1e-12)
})

test_that("compound support and grid exhaustion cannot masquerade as zero mass or quantiles", {
  for (family in c("cpegpd", "cpdegpd")) {
    args <- list(lambda = 1, sigma = 1, xi = 0, kappa = 1, K = 200)
    d <- get(paste0("d", family))
    p <- get(paste0("p", family))
    q <- get(paste0("q", family))
    expect_equal(do.call(d, c(list(x = c(-Inf, -1, Inf)), args)), c(0, 0, 0))
    expect_equal(do.call(p, c(list(q = c(-Inf, -1, Inf)), args)), c(0, 0, 1))
    expect_equal(do.call(d, c(list(x = -1, log = TRUE), args)), -Inf)
    expect_equal(do.call(q, c(list(p = c(0, 1, NA_real_)), args)),
                 c(0, Inf, NA_real_))
    expect_warning(bad <- do.call(q, c(list(p = c(-0.1, 1.1)), args)), "NaNs")
    expect_true(all(is.nan(bad)))
    args$sigma <- 100
    args$K <- 10
    expect_error(do.call(q, c(list(p = 0.9), args)), "increase 'K'")
    expect_error(do.call(p, c(list(q = 100), args)), "increase 'K'")
  }
  expect_equal(dcpdegpd(c(0.5, 1.5), 1, sigma = 1, xi = 0, kappa = 1), c(0, 0))
  expect_equal(pcpdegpd(0.5, 1, sigma = 1, xi = 0, kappa = 1),
               pcpdegpd(0, 1, sigma = 1, xi = 0, kappa = 1))
})

test_that("sorting durations preserves paired unnamed and named aggregates", {
  samples <- list(rep(80, 5), rep(10, 5))
  a <- egpd:::.prepare_aggregate_input(samples, c(8, 1))
  b <- egpd:::.prepare_aggregate_input(setNames(samples, c("8", "1")), c(8, 1))
  expect_equal(a$durations, c(1L, 8L))
  expect_equal(vapply(a$full_samples, mean, 0), c(10, 80))
  expect_equal(unname(b$full_samples), a$full_samples)
  expect_error(egpd:::.prepare_aggregate_input(samples, c(1, 1)), "same length|distinct")
  expect_error(egpd:::.prepare_aggregate_input(setNames(samples, c("8", "2")),
                                               c(8, 1)), "named by duration")
  expect_error(fit_aggregate_egpd(list(c(0, 0.01, 1), c(0, 1, 2)), c(1, 2),
                                  p = 1, q = 1, cpegpd.h = 0), "positive and finite")
})

test_that("blend validation resolves concentrated Beta transitions", {
  expect_error(pbgpd(0.4, xi = -1, pa = 0.5, pb = 0.99,
                     c1 = 100500001, c2 = 99500001), "negative")
  x <- seq(0, 2, length.out = 120)
  expect_true(all(diff(pbgpd(x, xi = 0.2)) >= -1e-12))
})

test_that("large discrete counts retain finite integer-valued double storage", {
  set.seed(45)
  b <- rbdegpd(100, 2, 1e10, 0.1, 5, 5, 0.2)
  m <- rmdgpd(100, 1e10, 0, rate = 1)
  for (y in list(b, m)) {
    expect_type(y, "double")
    expect_true(all(is.finite(y)))
    expect_equal(y, floor(y))
    expect_true(any(y > .Machine$integer.max))
  }
})

test_that("neural simulation is reproducible across repeated and parallel calls", {
  fun <- function(i) c(sample(1:10, 1), runif(3))
  set.seed(501)
  a <- egpd:::.neural_simulate(1:8, fun, 1)
  set.seed(501)
  b <- egpd:::.neural_simulate(1:8, fun, 1)
  expect_identical(a, b)
  if (.Platform$OS.type != "windows") {
    set.seed(501)
    expect_identical(egpd:::.neural_simulate(1:8, fun, 2), a)
  }
})

test_that("neural public constructor and explicit Julia seed are forwarded", {
  calls <- list()
  local_mocked_bindings(
    .julia_let = function(...) { calls[[length(calls) + 1L]] <<- list(...); "mock" },
    .package = "egpd"
  )
  expect_identical(egpd:::.init_npe_begpd(3L, 2L), "mock")
  egpd:::.neural_seed(19)
  expect_match(calls[[1]][[1]], "initialize_ratio_network")
  expect_equal(calls[[1]]$d, 3L)
  expect_match(calls[[2]][[1]], "Random.seed!", fixed = TRUE)
  expect_equal(calls[[2]]$seed, 19L)
  architecture <- readLines(system.file("julia", "begpd_architecture.jl",
                                       package = "egpd", mustWork = TRUE))
  expect_true(any(grepl("RatioEstimator(d, summary; num_summaries = 128)", architecture,
                       fixed = TRUE)))
})

test_that("mocked learned log-ratio reaches the NPE importance sampler intact", {
  ## Round 3 (D3): the fixed prior-grid SIR (sampleposterior over 20000 prior
  ## draws) was replaced by adaptive importance sampling on the log-ratio.
  path <- test_path("audit-neural-mock.bson")
  on.exit(unlink(c(path, paste0(path, ".metadata.rds"))), add = TRUE)
  file.create(path)
  egpd:::.write_neural_metadata(path, "mdgpd", "npe")
  seen <- NULL
  centre <- c(log(2), 0.2, log(0.5))
  local_mocked_bindings(
    .check_julia_deps = function() NULL,
    .init_julia_begpd = function() NULL,
    .init_npe_begpd = function(...) "mock",
    .neural_loadstate = function(...) NULL,
    .neural_logratio = function(estimator, Z, theta) {
      seen <<- list(Z = Z, theta_rows = nrow(theta))
      -0.5 * colSums(((theta - centre) / 0.01)^2)
    },
    .package = "egpd"
  )
  set.seed(1)
  fit <- egpd:::.distfit_neuralbayes(matrix(1:20, ncol = 2), "mdgpd",
                                     model.path = path, nsamples = 7L)
  expect_equal(dim(seen$Z), c(2L, 10L))
  expect_equal(seen$theta_rows, 3L)
  expect_equal(dim(fit$posterior_samples), c(3L, 7L))
  expect_equal(fit$estimate, c(sigma = 2, xi = 0.2, rate = 0.5), tolerance = 0.03)
  expect_gte(fit$importance_sampling$ess, 200)
})

test_that("mocked NBE uses the public estimate call and validates output dimension", {
  path <- test_path("audit-neural-nbe.bson")
  on.exit(unlink(c(path, paste0(path, ".metadata.rds"))), add = TRUE)
  file.create(path)
  egpd:::.write_neural_metadata(path, "mdgpd", "nbe")
  code <- NULL
  local_mocked_bindings(
    .check_julia_deps = function() NULL,
    .init_julia_begpd = function() NULL,
    .init_nbe_begpd = function(d, n) {
      expect_equal(c(d, n), c(3L, 2L))
      "mock"
    },
    .neural_loadstate = function(...) NULL,
    .julia_let = function(...) {
      code <<- list(...)[[1]]
      matrix(c(log(2), 0.2, log(0.5)), ncol = 1)
    },
    .package = "egpd"
  )
  fit <- egpd:::.distfit_neuralbayes(matrix(1:20, ncol = 2), "mdgpd",
                                     model.path = path, estimator = "nbe")
  expect_match(code, "egpd_estimate(", fixed = TRUE)
  expect_equal(fit$estimate, c(sigma = 2, xi = 0.2, rate = 0.5))
  expect_null(fit$posterior_samples)
  local_mocked_bindings(.julia_let = function(...) 1:4, .package = "egpd")
  expect_error(egpd:::.distfit_neuralbayes(matrix(1:20, ncol = 2), "mdgpd",
                   model.path = path, estimator = "nbe"), "incompatible parameter")
})

test_that("metadata rejects stale, mismatched and malformed model contracts", {
  path <- test_path("audit-neural-contract.bson")
  on.exit(unlink(c(path, paste0(path, ".metadata.rds"))), add = TRUE)
  file.create(path)
  expect_error(egpd:::.validate_neural_model(path, "mdgpd", "npe", 2L), "[Rr]etrain")
  egpd:::.write_neural_metadata(path, "mdgpd", "npe", 2L)
  expect_silent(egpd:::.validate_neural_model(path, "mdgpd", "npe", 2L))
  expect_error(egpd:::.validate_neural_model(path, "mdgpd", "npe", 3L), "Incompatible")
  expect_error(egpd:::.validate_neural_model(path, "zimdgpd", "npe", 2L), "Incompatible")
  expect_error(egpd:::.plot_mdgpd(list(estimate = c(sigma = 2, xi = 0.2,
                           lambda = 1, rho = 0.5), data = matrix(1:20, ncol = 2),
                           family = "mdgpd", estimator_type = "npe")), "Legacy")
})

test_that("mocked trainers share parameter dimensions, priors, seeds and metadata", {
  path <- test_path("audit-training-models")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  records <- list()
  julia_seeds <- integer()
  init_calls <- 0L
  local_mocked_bindings(
    .check_julia_deps = function() NULL,
    .init_julia_begpd = function() {
      init_calls <<- init_calls + 1L
      runif(init_calls)
    },
    .julia_let = function(code, seed) {
      expect_match(code, "Random.seed!", fixed = TRUE)
      julia_seeds <<- c(julia_seeds, seed)
    },
    .init_npe_begpd = function(d = 6L, n = 2L) list(d = d, n = n),
    .init_nbe_begpd = function(d = 6L, n = 2L) list(d = d, n = n),
    .neural_train = function(estimator, ...) {
      args <- list(...)
      expect_equal(nrow(args$theta_train), estimator$d)
      expect_true(all(vapply(args$Z_train, nrow, 0L) == estimator$n))
      records[[length(records) + 1L]] <<- args
      estimator
    },
    .neural_savestate = function(estimator, path) file.create(path),
    .package = "egpd"
  )
  for (family in c("begpd", "bdegpd", "bzidegpd", "mdgpd", "zimdgpd")) {
    args <- list(savepath = path, K = 12L, m = c(10L, 12L), epochs = 1L,
                 estimator = "both", seed = 193L, mc.cores = 1L, verbose = FALSE)
    train <- if (family == "begpd") train_begpd else if (family %in%
                   c("bdegpd", "bzidegpd")) train_bdegpd else train_mdgpd
    if (family != "begpd") args$family <- family
    data_dim <- if (family %in% c("mdgpd", "zimdgpd")) 3L else 2L
    if (data_dim == 3L) args$data_dim <- data_dim
    paths <- do.call(train, args)
    cfg <- egpd:::.neuralbayes_config(family, data_dim)
    for (estimator in c("npe", "nbe"))
      expect_silent(egpd:::.validate_neural_model(paths[[estimator]], family,
                                                  estimator, data_dim))
    a <- records[[length(records)]]
    theta <- cfg$backtransform(a$theta_train)
    prior <- cfg$backtransform(egpd:::.sample_theta_prior(family, 500L, data_dim))
    expect_equal(nrow(prior), cfg$d)
    if (data_dim == 3L) {
      for (z in list(theta, prior)) {
        expect_true(all(z[1, ] >= 0.1 & z[1, ] <= 10))
        expect_true(all(z[2, ] >= -0.4 & z[2, ] <= 0.5))
        expect_true(all(z[3, ] >= 0.01 & z[3, ] <= 5))
        if (family == "zimdgpd")
          expect_true(all(z[4, ] >= 0.01 & z[4, ] <= 0.9))
      }
    }
    args$mc.cores <- if (.Platform$OS.type == "windows") 1L else 2L
    do.call(train, args)
    b <- records[[length(records)]]
    expect_identical(a$theta_train, b$theta_train)
    expect_identical(a$theta_val, b$theta_val)
    expect_identical(a$Z_train, b$Z_train)
    expect_identical(a$Z_val, b$Z_val)
  }
  expect_identical(julia_seeds, rep(193L, 10))
})
