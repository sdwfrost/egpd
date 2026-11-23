test_that("MDGPD supports arbitrary dimension and signed integer excesses", {
  for (d in c(2L, 3L, 5L)) {
    set.seed(42)
    y <- rmdgpd(1000, sigma = 2, xi = 0.2, rate = 2, d = d)
    expect_equal(dim(y), c(1000L, d))
    expect_equal(colnames(y), paste0("Y", seq_len(d)))
    expect_true(all(is.finite(y)))
    expect_equal(y, ceiling(y))
    expect_true(any(y < 0))
    expect_true(all(apply(y, 1, max) > 0))
  }
})

test_that("MDGPD validates parameters and refuses ambiguous legacy calls", {
  expect_error(rmdgpd(-1, 2, 0.2, rate = 1), "'n'")
  expect_error(rmdgpd(10, -1, 0.2, rate = 1), "'sigma'")
  expect_error(rmdgpd(10, 2, Inf, rate = 1), "'xi'")
  expect_error(rmdgpd(10, 2, 0.2, rate = -0.1), "'rate'")
  expect_error(rmdgpd(10, 2, 0.2, rate = Inf), "'rate'")
  expect_error(rmdgpd(10, 2, 0.2, rate = 1, d = 2.5), "'d'")
  expect_error(rmdgpd(10, 2, 0.2, 1, 0.5), "Legacy.*retraining")
  expect_error(rmdgpd(10, 2, 0.2, lambda = 1), "Legacy")
  expect_error(rmdgpd(10, 2, 0.2, rho = 0.5, rate = 1), "Legacy")
  expect_no_error(rmdgpd(10, 2, -0.2, rate = 1))
  expect_no_error(rmdgpd(10, 2, 0.2, rate = 0))
})

test_that("smaller identifiable rates yield more nearly equal components", {
  set.seed(42)
  high <- rmdgpd(10000, 2, 0.2, rate = 0.01)
  low <- rmdgpd(10000, 2, 0.2, rate = 5)
  expect_gt(mean(high[, 1] == high[, 2]), mean(low[, 1] == low[, 2]))
})

test_that("zero inflation adds exactly the joint zero atom", {
  set.seed(42)
  y <- rzimdgpd(20000, 2, 0.2, rate = 1, pi0 = 0.3, d = 3)
  expect_equal(mean(apply(y, 1, function(x) all(x == 0))), 0.3,
               tolerance = 0.015)
  expect_error(rzimdgpd(10, 2, 0.2, rate = 1, pi0 = 1), "'pi0'")
  expect_error(rzimdgpd(10, 2, 0.2, rate = 1, pi0 = NA), "'pi0'")
  set.seed(5)
  a <- rzimdgpd(100, 2, 0.2, rate = 1, pi0 = 0.3)
  set.seed(5)
  expect_identical(a, rzimdgpd(100, 2, 0.2, rate = 1, pi0 = 0.3))
})

test_that("multivariate fits still require neural Bayes", {
  y <- matrix(1:20, ncol = 2)
  expect_error(distfit(y, family = "mdgpd", method = "mle"),
               "requires method='neuralbayes'")
  expect_error(distfit(y, family = "zimdgpd", method = "mle"),
               "requires method='neuralbayes'")
})

test_that("legacy bundled MDGPD assets are rejected before Julia is called", {
  for (estimator in c("npe", "nbe")) {
    expect_error(
      distfit(matrix(1:20, ncol = 2), family = "mdgpd",
              method = "neuralbayes", estimator = estimator),
      "[Rr]etrain"
    )
  }
})
