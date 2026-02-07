# Fit a DEGPD-1 model used across multiple test blocks
set.seed(42)
y <- rdiscegpd(500, sigma = 3, xi = 0.2, kappa = 2, type = 1)
df <- data.frame(y = y)
fit <- suppressMessages(
  egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
       data = df, family = "degpd", degpd.args = list(m = 1))
)

# ---- vcov.egpd ----

test_that("vcov returns a named matrix matching Vp", {
  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), length(fit$coefficients))
  expect_equal(ncol(V), length(fit$coefficients))
  expect_equal(unname(V), unname(fit$Vp), tolerance = 1e-12)
  expect_equal(rownames(V), names(fit$coefficients))
  expect_equal(colnames(V), names(fit$coefficients))
})

# ---- confint.egpd — Wald method ----

test_that("Wald CI returns correct structure", {
  ci <- confint(fit, method = "wald")
  expect_true(is.matrix(ci))
  expect_equal(dim(ci), c(3L, 2L))
  expect_equal(colnames(ci), c("lower", "upper"))
  expect_true(all(rownames(ci) %in% c("scale", "shape", "kappa")))
})

test_that("Wald CI values are valid", {
  ci <- confint(fit, method = "wald")
  expect_true(all(ci[, "lower"] < ci[, "upper"]))
  expect_true(all(ci > 0))
})

test_that("Wald CI contains true values", {
  ci <- confint(fit, method = "wald")
  expect_lt(ci["scale", "lower"], 3)
  expect_gt(ci["scale", "upper"], 3)
  expect_lt(ci["shape", "lower"], 0.2)
  expect_gt(ci["shape", "upper"], 0.2)
  expect_lt(ci["kappa", "lower"], 2)
  expect_gt(ci["kappa", "upper"], 2)
})

# ---- parm argument ----

test_that("parm selects parameters by name", {
  ci <- confint(fit, parm = "scale", method = "wald")
  expect_equal(nrow(ci), 1L)
  expect_equal(rownames(ci), "scale")

  ci2 <- confint(fit, parm = c("scale", "kappa"), method = "wald")
  expect_equal(nrow(ci2), 2L)
  expect_equal(rownames(ci2), c("scale", "kappa"))
})

test_that("parm selects parameters by index", {
  ci_name <- confint(fit, parm = "scale", method = "wald")
  ci_idx  <- confint(fit, parm = 1, method = "wald")
  expect_equal(nrow(ci_idx), 1L)
  expect_equal(unname(ci_name), unname(ci_idx), tolerance = 1e-12)
})

test_that("parm errors on invalid input", {
  expect_error(confint(fit, parm = "nonexistent"), "not found")
  expect_error(confint(fit, parm = 99), "out of range")
})

# ---- level argument ----

test_that("level controls CI width", {
  ci_90 <- confint(fit, method = "wald", level = 0.90)
  ci_95 <- confint(fit, method = "wald", level = 0.95)
  ci_99 <- confint(fit, method = "wald", level = 0.99)

  width_90 <- ci_90[, "upper"] - ci_90[, "lower"]
  width_95 <- ci_95[, "upper"] - ci_95[, "lower"]
  width_99 <- ci_99[, "upper"] - ci_99[, "lower"]

  expect_true(all(width_90 < width_95))
  expect_true(all(width_99 > width_95))
})

test_that("level = 0 collapses CI to point estimate", {
  ci <- confint(fit, method = "wald", level = 0)
  expect_equal(ci[, "lower"], ci[, "upper"])
})

# ---- confint.egpd — Profile method ----

test_that("Profile CI returns correct structure", {
  ci <- confint(fit, method = "profile")
  expect_true(is.matrix(ci))
  expect_equal(dim(ci), c(3L, 2L))
  expect_equal(colnames(ci), c("lower", "upper"))
})

test_that("Profile CI values are valid", {
  ci <- confint(fit, method = "profile")
  expect_true(all(ci[, "lower"] < ci[, "upper"]))
  expect_true(all(ci > 0))
})

test_that("Profile CI contains true values", {
  ci <- confint(fit, method = "profile")
  expect_lt(ci["scale", "lower"], 3)
  expect_gt(ci["scale", "upper"], 3)
  expect_lt(ci["shape", "lower"], 0.2)
  expect_gt(ci["shape", "upper"], 0.2)
  expect_lt(ci["kappa", "lower"], 2)
  expect_gt(ci["kappa", "upper"], 2)
})

test_that("Profile and Wald CI widths roughly agree", {
  ci_w <- confint(fit, method = "wald")
  ci_p <- confint(fit, method = "profile")
  for (p in rownames(ci_w)) {
    w_width <- ci_w[p, "upper"] - ci_w[p, "lower"]
    p_width <- ci_p[p, "upper"] - ci_p[p, "lower"]
    expect_lt(abs(p_width - w_width) / w_width, 0.5,
              label = sprintf("Profile/Wald width ratio for %s", p))
  }
})

test_that("Profile CI works with parm subset", {
  ci <- confint(fit, parm = "scale", method = "profile")
  expect_equal(nrow(ci), 1L)
  expect_equal(rownames(ci), "scale")
})

# ---- Profile/Wald convergence at large n ----

test_that("Profile and Wald converge at large sample size", {
  set.seed(99)
  y_big <- rdiscegpd(5000, sigma = 3, xi = 0.2, kappa = 2, type = 1)
  df_big <- data.frame(y = y_big)
  fit_big <- suppressMessages(
    egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
         data = df_big, family = "degpd", degpd.args = list(m = 1))
  )

  ci_w <- confint(fit_big, method = "wald")
  ci_p <- confint(fit_big, method = "profile")

  for (p in rownames(ci_w)) {
    w_width <- ci_w[p, "upper"] - ci_w[p, "lower"]
    p_width <- ci_p[p, "upper"] - ci_p[p, "lower"]
    expect_lt(abs(p_width - w_width) / w_width, 0.2,
              label = sprintf("Large-n convergence for %s", p))
  }
})

# ---- Profile CI restrictions ----

test_that("Profile CI errors for models with covariates", {
  set.seed(50)
  df2 <- data.frame(
    y = rdiscegpd(200, sigma = 3, xi = 0.2, kappa = 2, type = 1),
    x = rnorm(200)
  )
  fit_cov <- suppressMessages(
    egpd(list(lsigma = y ~ x, lxi = ~ 1, lkappa = ~ 1),
         data = df2, family = "degpd", degpd.args = list(m = 1))
  )

  expect_error(confint(fit_cov, method = "profile"), "intercept-only")
  expect_true(is.matrix(confint(fit_cov, method = "wald")))
})

# ---- EGPD (continuous) model ----

test_that("confint and vcov work for continuous EGPD", {
  set.seed(77)
  y_egpd <- regpd(500, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
  df_egpd <- data.frame(y = y_egpd)
  fit_egpd <- suppressMessages(
    egpd(list(lpsi = y ~ 1, xi = ~1, lkappa = ~1),
         data = df_egpd, family = "egpd", egpd.args = list(m = 1))
  )

  ci_w <- confint(fit_egpd, method = "wald")
  expect_equal(dim(ci_w), c(3L, 2L))
  expect_true("scale" %in% rownames(ci_w))
  expect_true("kappa" %in% rownames(ci_w))

  ci_p <- confint(fit_egpd, method = "profile")
  expect_equal(dim(ci_p), c(3L, 2L))

  V <- vcov(fit_egpd)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 3L)
})
