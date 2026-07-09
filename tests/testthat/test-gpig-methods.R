## Evaluation routes for the GPIG pmf (see gpig_control / src/gpig.cpp).
## The reference is the "legacy" recursion wherever it is numerically valid; beyond
## that the routes are checked against one another and against finite differences.

lp <- function(y, a, b, c) egpd:::gpig_logp_cpp(y, a, b, c)

PARS <- list(c(0.5, 20, 0.6), c(0.45, 300, 0.7), c(0.7, 50, 0.4),
             c(0.3, 5, 0.85), c(0.2, 2, 0.95))

restore <- function() gpig_control("hybrid")

test_that("log-space recursion reproduces the legacy recursion exactly", {
  on.exit(restore())
  for (p in PARS) for (y in c(0L, 1L, 5L, 50L, 300L)) {
    gpig_control("legacy");    old <- lp(y, p[1], p[2], p[3])[1:4]
    gpig_control("recursion"); new <- lp(y, p[1], p[2], p[3])[1:4]
    expect_equal(unname(new), unname(old), tolerance = 1e-9)
  }
})

test_that("the log-space recursion survives the p_0 underflow that kills legacy", {
  on.exit(restore())
  a <- 0.45; c <- 0.7
  for (mu in c(5000, 50000)) {
    b <- mu * (1 - c)^(1 - a) / (a * c); y <- as.integer(round(mu))
    gpig_control("legacy")
    expect_true(is.infinite(lp(y, a, b, c)[["logp"]]))     # p_0 flushes to zero
    gpig_control("recursion")
    expect_true(is.finite(lp(y, a, b, c)[["logp"]]))
    gpig_control("hybrid")
    expect_true(is.finite(lp(y, a, b, c)[["logp"]]))
  }
})

test_that("relative severity truncation matches the exact recursion", {
  on.exit(restore())
  for (p in PARS) for (y in c(200L, 1000L)) {
    gpig_control("recursion");          ex <- lp(y, p[1], p[2], p[3])[["logp"]]
    gpig_control("trunc", eps = 1e-12); tr <- lp(y, p[1], p[2], p[3])[["logp"]]
    if (is.finite(ex)) expect_equal(tr, ex, tolerance = 1e-8)
  }
})

test_that("the analytic saddlepoint gradient matches central finite differences", {
  on.exit(restore())
  gpig_control("saddlepoint", order = 2)
  for (p in PARS) for (y in c(50L, 500L)) {
    a <- p[1]; b <- p[2]; c <- p[3]
    g <- lp(y, a, b, c)[c("da", "db", "dc")]
    h <- 1e-6
    fd <- c(da = (lp(y, a + h * a, b, c)[["logp"]] - lp(y, a - h * a, b, c)[["logp"]]) / (2 * h * a),
            db = (lp(y, a, b * (1 + h), c)[["logp"]] - lp(y, a, b * (1 - h), c)[["logp"]]) / (2 * h * b),
            dc = (lp(y, a, b, c + h * c)[["logp"]] - lp(y, a, b, c - h * c)[["logp"]]) / (2 * h * c))
    expect_equal(unname(g), unname(fd), tolerance = 1e-5)
  }
})

test_that("d log p / dc has the closed form -a b (1-c)^(a-1) + y/c", {
  on.exit(restore())
  # z^ solves a b z (1-z)^(a-1) = y and so is free of c, forcing dt^/dc = -1/c and
  # cancelling every other c-dependence in the saddlepoint. Holds for the UNNORMALISED
  # pmf: normalising adds a -dlogS/dc term.
  gpig_control("saddlepoint", order = 2, normalize = FALSE)
  for (p in PARS) for (y in c(50L, 500L)) {
    a <- p[1]; b <- p[2]; c <- p[3]
    expect_equal(lp(y, a, b, c)[["dc"]], -a * b * (1 - c)^(a - 1) + y / c, tolerance = 1e-10)
  }
})

test_that("|corr| flags the saddlepoint's own unreliability", {
  on.exit(restore())
  for (p in PARS) for (y in c(300L, 1000L)) {
    gpig_control("recursion");              ex <- lp(y, p[1], p[2], p[3])[["logp"]]
    gpig_control("saddlepoint", order = 2); sp <- lp(y, p[1], p[2], p[3])
    if (is.finite(ex) && sp[["corr"]] <= 1e-3)
      expect_equal(sp[["logp"]], ex, tolerance = 1e-4)     # small corr => accurate
  }
})

test_that("the exact recursion fails loudly rather than returning a saturated value", {
  on.exit(restore())
  # p_y / max_j p_j < ~1e-308: no common scale spans the ladder. The recursion must
  # return -Inf, never a finite number, and hybrid must fall back to the saddlepoint.
  a <- 0.5; b <- 20; c <- 0.6; y <- 2000L
  gpig_control("recursion")
  expect_true(is.infinite(lp(y, a, b, c)[["logp"]]))
  gpig_control("hybrid")
  h <- lp(y, a, b, c)[["logp"]]
  expect_true(is.finite(h))
  gpig_control("saddlepoint")
  expect_equal(h, lp(y, a, b, c)[["logp"]], tolerance = 1e-10)
})

test_that("hybrid agrees with the exact recursion wherever the recursion is valid", {
  on.exit(restore())
  for (p in PARS) for (y in c(250L, 500L)) {
    gpig_control("recursion"); ex <- lp(y, p[1], p[2], p[3])[["logp"]]
    gpig_control("hybrid");    hy <- lp(y, p[1], p[2], p[3])[["logp"]]
    if (is.finite(ex)) expect_equal(hy, ex, tolerance = 1e-3)
  }
})

test_that("gpig_control validates its arguments and round-trips", {
  on.exit(restore())
  expect_error(gpig_control("nope"))
  expect_error(gpig_control("hybrid", order = 3))
  expect_error(gpig_control("hybrid", eps = 0))
  expect_error(gpig_control("hybrid", sptol = -1))
  old <- gpig_control("saddlepoint", yswitch = 50L, order = 1L)
  cur <- egpd:::gpig_get_opts()
  expect_identical(cur$method, 3L)
  expect_identical(cur$yswitch, 50L)
  expect_identical(cur$order, 1L)
  do.call(gpig_control, old)
  expect_identical(egpd:::gpig_get_opts()$method, unname(egpd:::GPIG_METHODS[[old$method]]))
})

test_that("a gpig GAM fits under every method and agrees with legacy", {
  on.exit(restore())
  set.seed(42)
  d <- data.frame(y = rgpig(200, a = 0.5, b = 30, c = 0.6))
  aics <- sapply(c("legacy", "recursion", "trunc", "hybrid"), function(m) {
    f <- try(egpd(list(lmu = y ~ 1, logita = ~1, logitc = ~1), data = d,
                  family = "gpig", gpig.args = list(method = m)), silent = TRUE)
    if (inherits(f, "try-error")) NA_real_ else as.numeric(f$AIC)
  })
  expect_true(all(is.finite(aics)))
  expect_lt(max(aics) - min(aics), 0.05)
})

test_that("egpd() restores the gpig method after a fit and warns off-family", {
  on.exit(restore())
  gpig_control("hybrid")
  set.seed(1)
  d <- data.frame(y = rgpig(120, a = 0.5, b = 20, c = 0.6))
  invisible(egpd(list(lmu = y ~ 1, logita = ~1, logitc = ~1), data = d,
                 family = "gpig", gpig.args = list(method = "recursion")))
  expect_identical(egpd:::gpig_get_opts()$method, unname(egpd:::GPIG_METHODS[["hybrid"]]))
  expect_warning(
    egpd(list(lmu = y ~ 1, lsigma = ~1), data = d, family = "pig",
         gpig.args = list(method = "recursion")),
    "ignored for family")
})

test_that("the saddlepoint never returns a log-pmf above zero", {
  on.exit(restore())
  # As c -> 1 (the discrete-stable boundary) |corr| reaches ~1e8. Applying the
  # second-order Daniels term there adds log1p(corr) ~ +20 to log p, yielding "pmf"
  # values above 1 and an unbounded pseudo-likelihood. Regression: a GPIG GAM on a
  # large-count series once returned AIC = -265478 with c pinned at 1.
  a <- 0.952; c <- 1 - 1e-6
  for (mth in c("saddlepoint", "hybrid")) {
    gpig_control(mth, order = 2)
    for (mu in c(1000, 5000, 77616)) {
      b <- mu * (1 - c)^(1 - a) / (a * c)
      for (y in as.integer(c(mu, 2 * mu))) expect_lte(lp(y, a, b, c)[["logp"]], 0)
    }
  }
})

test_that("the saddlepoint pmf never sums to more than one", {
  on.exit(restore())
  ys <- 1:4000
  tot <- function(mth, a, b, c) {
    gpig_control(mth, order = 2)
    sum(exp(vapply(ys, function(y) lp(as.integer(y), a, b, c)[["logp"]], 0)))
  }
  for (p in list(c(0.5, 20, 0.6), c(0.45, 300, 0.7), c(0.952, 500, 1 - 1e-6))) {
    # With normalize=TRUE the residual is the quadrature error in S, ~1e-3, not 0.
    for (mth in c("saddlepoint", "hybrid"))
      expect_lte(tot(mth, p[1], p[2], p[3]), 1 + 2e-3)
  }
})

test_that("the normalising quadrature reproduces the lattice sum of the same pmf", {
  on.exit(restore())
  # The integrand must be the very pmf the likelihood uses, second-order gate included.
  # Integrating an order-1 integrand against an order-2 pmf was off by 7% at c=0.999936.
  gpig_control("saddlepoint", nquad = 1024L, normalize = FALSE, order = 2L)
  for (p in list(c(0.45, 300, 0.7), c(0.7, 50, 0.4), c(0.609, 187, 0.999936))) {
    a <- p[1]; b <- p[2]; cc <- p[3]
    mu <- a * b * cc * (1 - cc)^(a - 1)
    ys <- unique(round(exp(seq(log(1), log(max(3e4, 60 * mu)), length.out = 20000))))
    brute <- sum(c(diff(ys), 1) * exp(vapply(ys, function(y) lp(as.integer(y), a, b, cc)[["logp"]], 0)))
    lS <- egpd:::gpig_logS_cpp(a, b, cc)
    # The contract is: either S is right, or the quadrature declines (returns exactly 0)
    # because the peak is unresolvable within the node cap. It is never wrong.
    if (lS != 0) expect_equal(exp(lS), brute, tolerance = 5e-3)
  }
})

test_that("the normalising quadrature declines rather than returning a wrong constant", {
  on.exit(restore())
  gpig_control("saddlepoint", nquad = 256L, normalize = FALSE, order = 2L)
  # At b = 4e4 a 256-node uniform rule returns S = 0 -- dividing the likelihood by zero.
  # The refinement must either converge or decline (exactly 0.0), never emit that value.
  for (p in list(c(0.5, 3000, 0.6), c(0.45, 40000, 0.7), c(0.5, 1e6, 0.6))) {
    lS <- egpd:::gpig_logS_cpp(p[1], p[2], p[3])
    expect_true(lS == 0 || (is.finite(lS) && abs(lS) < 1))   # never a huge/negative-inf shift
  }
})

test_that("normalisation sharpens the saddlepoint where its series has converged", {
  on.exit(restore())
  for (p in list(c(0.5, 0.6), c(0.45, 0.7))) {
    a <- p[1]; cc <- p[2]
    for (y in c(300L, 800L)) {
      b <- y * (1 - cc)^(1 - a) / (a * cc)          # GAM regime: b tracks the mean
      gpig_control("recursion"); ex <- lp(y, a, b, cc)[["logp"]]
      gpig_control("saddlepoint", normalize = FALSE, nquad = 1024L)
      un <- abs(lp(y, a, b, cc)[["logp"]] - ex)
      gpig_control("saddlepoint", normalize = TRUE, nquad = 1024L)
      nm <- abs(lp(y, a, b, cc)[["logp"]] - ex)
      expect_lt(nm, un)                              # strictly better, typically 10-100x
    }
  }
})
