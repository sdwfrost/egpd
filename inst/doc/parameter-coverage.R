## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)


## ----load-packages, message = FALSE, warning = FALSE--------------------------
library(egpd)


## ----sim-helper---------------------------------------------------------------
sim_param_coverage <- function(sigma, xi, kappa, n = 500, n_rep = 200,
                               level = 0.95) {
  true_vals <- c(sigma, xi, kappa)

  # Matrices to store results: n_rep x 3
  covered  <- matrix(NA, nrow = n_rep, ncol = 3)
  estimates <- matrix(NA, nrow = n_rep, ncol = 3)
  colnames(covered)  <- c("sigma", "xi", "kappa")
  colnames(estimates) <- c("sigma", "xi", "kappa")

  for (i in seq_len(n_rep)) {
    y <- rdiscegpd(n, sigma = sigma, xi = xi, kappa = kappa, type = 1)
    df <- data.frame(y = y)

    fit <- tryCatch(
      suppressMessages(egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
           data = df, family = "degpd", degpd.args = list(m = 1))),
      error = function(e) NULL
    )

    if (is.null(fit)) next

    ci <- tryCatch(
      suppressWarnings(confint(fit, method = "wald", level = level)),
      error = function(e) NULL
    )
    if (is.null(ci)) next

    # Response-scale point estimates (level = 0 collapses CI to MLE)
    estimates[i, ] <- confint(fit, method = "wald", level = 0)[, 1]

    for (k in 1:3) {
      if (!any(is.na(ci[k, ])))
        covered[i, k] <- (true_vals[k] >= ci[k, 1]) &
                          (true_vals[k] <= ci[k, 2])
    }
  }

  # Coverage = fraction of successful fits where CI covers truth
  ok <- complete.cases(covered)
  cov_rate <- colMeans(covered[ok, , drop = FALSE])

  # Relative bias = (estimate - true) / true
  ok_est <- complete.cases(estimates)
  rel_bias <- sweep(estimates[ok_est, , drop = FALSE], 2, true_vals, "-")
  rel_bias <- sweep(rel_bias, 2, true_vals, "/")

  list(coverage = cov_rate, n_ok = sum(ok), n_rep = n_rep,
       estimates = estimates[ok_est, , drop = FALSE],
       rel_bias = rel_bias,
       true_vals = true_vals)
}


## ----vary-sigma, cache = TRUE-------------------------------------------------
set.seed(101)
sigma_vals <- c(1, 2, 3, 5, 10)
xi_fix <- 0.2
kappa_fix <- 2

res_sigma <- lapply(sigma_vals, function(s) {
  sim_param_coverage(sigma = s, xi = xi_fix, kappa = kappa_fix)
})


## ----sigma-table--------------------------------------------------------------
tab_sigma <- do.call(rbind, lapply(seq_along(sigma_vals), function(j) {
  r <- res_sigma[[j]]
  data.frame(sigma_true = sigma_vals[j],
             cov_sigma = round(r$coverage["sigma"], 3),
             cov_xi    = round(r$coverage["xi"], 3),
             cov_kappa = round(r$coverage["kappa"], 3),
             n_ok      = r$n_ok)
}))
tab_sigma


## ----sigma-coverage-plot, fig.width = 6, fig.height = 5-----------------------
se_cov <- sqrt(0.95 * 0.05 / 200)

plot(sigma_vals, tab_sigma$cov_sigma, type = "b", pch = 16,
     col = "steelblue", lwd = 2, ylim = c(0.85, 1),
     xlab = expression(sigma[true]), ylab = "Coverage",
     main = expression("Coverage vs " * sigma[true]))
lines(sigma_vals, tab_sigma$cov_xi, type = "b", pch = 17,
      col = "firebrick", lwd = 2)
lines(sigma_vals, tab_sigma$cov_kappa, type = "b", pch = 15,
      col = "forestgreen", lwd = 2)
abline(h = 0.95, lty = 2, col = "grey40")
rect(min(sigma_vals) - 1, 0.95 - 2 * se_cov,
     max(sigma_vals) + 1, 0.95 + 2 * se_cov,
     col = adjustcolor("grey70", alpha.f = 0.3), border = NA)
legend("bottomright",
       legend = c(expression(sigma), expression(xi), expression(kappa)),
       col = c("steelblue", "firebrick", "forestgreen"),
       pch = c(16, 17, 15), lwd = 2)


## ----vary-xi, cache = TRUE----------------------------------------------------
set.seed(102)
xi_vals <- c(0.05, 0.1, 0.2, 0.5, 1.0)
sigma_fix <- 3

res_xi <- lapply(xi_vals, function(x) {
  sim_param_coverage(sigma = sigma_fix, xi = x, kappa = kappa_fix)
})


## ----xi-table-----------------------------------------------------------------
tab_xi <- do.call(rbind, lapply(seq_along(xi_vals), function(j) {
  r <- res_xi[[j]]
  data.frame(xi_true   = xi_vals[j],
             cov_sigma = round(r$coverage["sigma"], 3),
             cov_xi    = round(r$coverage["xi"], 3),
             cov_kappa = round(r$coverage["kappa"], 3),
             n_ok      = r$n_ok)
}))
tab_xi


## ----xi-coverage-plot, fig.width = 6, fig.height = 5--------------------------
plot(xi_vals, tab_xi$cov_sigma, type = "b", pch = 16,
     col = "steelblue", lwd = 2, ylim = c(0.85, 1),
     xlab = expression(xi[true]), ylab = "Coverage",
     main = expression("Coverage vs " * xi[true]))
lines(xi_vals, tab_xi$cov_xi, type = "b", pch = 17,
      col = "firebrick", lwd = 2)
lines(xi_vals, tab_xi$cov_kappa, type = "b", pch = 15,
      col = "forestgreen", lwd = 2)
abline(h = 0.95, lty = 2, col = "grey40")
rect(min(xi_vals) - 0.1, 0.95 - 2 * se_cov,
     max(xi_vals) + 0.1, 0.95 + 2 * se_cov,
     col = adjustcolor("grey70", alpha.f = 0.3), border = NA)
legend("bottomright",
       legend = c(expression(sigma), expression(xi), expression(kappa)),
       col = c("steelblue", "firebrick", "forestgreen"),
       pch = c(16, 17, 15), lwd = 2)


## ----vary-kappa, cache = TRUE-------------------------------------------------
set.seed(103)
kappa_vals <- c(0.5, 1, 2, 3, 5)
xi_fix2 <- 0.2

res_kappa <- lapply(kappa_vals, function(k) {
  sim_param_coverage(sigma = sigma_fix, xi = xi_fix2, kappa = k)
})


## ----kappa-table--------------------------------------------------------------
tab_kappa <- do.call(rbind, lapply(seq_along(kappa_vals), function(j) {
  r <- res_kappa[[j]]
  data.frame(kappa_true = kappa_vals[j],
             cov_sigma  = round(r$coverage["sigma"], 3),
             cov_xi     = round(r$coverage["xi"], 3),
             cov_kappa  = round(r$coverage["kappa"], 3),
             n_ok       = r$n_ok)
}))
tab_kappa


## ----kappa-coverage-plot, fig.width = 6, fig.height = 5-----------------------
plot(kappa_vals, tab_kappa$cov_sigma, type = "b", pch = 16,
     col = "steelblue", lwd = 2, ylim = c(0.85, 1),
     xlab = expression(kappa[true]), ylab = "Coverage",
     main = expression("Coverage vs " * kappa[true]))
lines(kappa_vals, tab_kappa$cov_xi, type = "b", pch = 17,
      col = "firebrick", lwd = 2)
lines(kappa_vals, tab_kappa$cov_kappa, type = "b", pch = 15,
      col = "forestgreen", lwd = 2)
abline(h = 0.95, lty = 2, col = "grey40")
rect(min(kappa_vals) - 0.5, 0.95 - 2 * se_cov,
     max(kappa_vals) + 0.5, 0.95 + 2 * se_cov,
     col = adjustcolor("grey70", alpha.f = 0.3), border = NA)
legend("bottomright",
       legend = c(expression(sigma), expression(xi), expression(kappa)),
       col = c("steelblue", "firebrick", "forestgreen"),
       pch = c(16, 17, 15), lwd = 2)


## ----bias-data----------------------------------------------------------------
# Combine all relative bias matrices with scenario labels
make_bias_df <- function(res_list, varied_name, varied_vals, fixed_label) {
  do.call(rbind, lapply(seq_along(varied_vals), function(j) {
    rb <- as.data.frame(res_list[[j]]$rel_bias)
    rb$scenario <- paste0(varied_name, " = ", varied_vals[j])
    rb$group <- fixed_label
    rb
  }))
}

bias_sigma <- make_bias_df(res_sigma, "sigma", sigma_vals, "Vary sigma")
bias_xi    <- make_bias_df(res_xi, "xi", xi_vals, "Vary xi")
bias_kappa <- make_bias_df(res_kappa, "kappa", kappa_vals, "Vary kappa")
bias_all   <- rbind(bias_sigma, bias_xi, bias_kappa)


## ----bias-sigma-plot, fig.width = 7, fig.height = 5---------------------------
par(mfrow = c(1, 3), mar = c(7, 4, 3, 1))

boxplot(sigma ~ scenario, data = bias_sigma,
        main = expression("Relative bias: " * sigma),
        ylab = "(est - true) / true", las = 2, col = "steelblue",
        outline = FALSE)
abline(h = 0, lty = 2, col = "grey40")

boxplot(xi ~ scenario, data = bias_sigma,
        main = expression("Relative bias: " * xi),
        ylab = "", las = 2, col = "firebrick",
        outline = FALSE)
abline(h = 0, lty = 2, col = "grey40")

boxplot(kappa ~ scenario, data = bias_sigma,
        main = expression("Relative bias: " * kappa),
        ylab = "", las = 2, col = "forestgreen",
        outline = FALSE)
abline(h = 0, lty = 2, col = "grey40")


## ----bias-xi-plot, fig.width = 7, fig.height = 5------------------------------
par(mfrow = c(1, 3), mar = c(7, 4, 3, 1))

boxplot(sigma ~ scenario, data = bias_xi,
        main = expression("Relative bias: " * sigma),
        ylab = "(est - true) / true", las = 2, col = "steelblue",
        outline = FALSE)
abline(h = 0, lty = 2, col = "grey40")

boxplot(xi ~ scenario, data = bias_xi,
        main = expression("Relative bias: " * xi),
        ylab = "", las = 2, col = "firebrick",
        outline = FALSE)
abline(h = 0, lty = 2, col = "grey40")

boxplot(kappa ~ scenario, data = bias_xi,
        main = expression("Relative bias: " * kappa),
        ylab = "", las = 2, col = "forestgreen",
        outline = FALSE)
abline(h = 0, lty = 2, col = "grey40")


## ----bias-kappa-plot, fig.width = 7, fig.height = 5---------------------------
par(mfrow = c(1, 3), mar = c(7, 4, 3, 1))

boxplot(sigma ~ scenario, data = bias_kappa,
        main = expression("Relative bias: " * sigma),
        ylab = "(est - true) / true", las = 2, col = "steelblue",
        outline = FALSE)
abline(h = 0, lty = 2, col = "grey40")

boxplot(xi ~ scenario, data = bias_kappa,
        main = expression("Relative bias: " * xi),
        ylab = "", las = 2, col = "firebrick",
        outline = FALSE)
abline(h = 0, lty = 2, col = "grey40")

boxplot(kappa ~ scenario, data = bias_kappa,
        main = expression("Relative bias: " * kappa),
        ylab = "", las = 2, col = "forestgreen",
        outline = FALSE)
abline(h = 0, lty = 2, col = "grey40")


## ----profile-example----------------------------------------------------------
set.seed(200)
y <- rdiscegpd(500, sigma = 3, xi = 0.2, kappa = 2, type = 1)
df <- data.frame(y = y)
fit <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
            data = df, family = "degpd", degpd.args = list(m = 1))

ci_wald    <- confint(fit, method = "wald")
ci_profile <- confint(fit, method = "profile")

cat("Wald CIs (response scale):\n")
print(ci_wald)
cat("\nProfile CIs (response scale):\n")
print(ci_profile)


## ----sim-both-helper----------------------------------------------------------
sim_param_coverage_both <- function(sigma, xi, kappa, n = 500, n_rep = 100,
                                    level = 0.95) {
  true_vals <- c(sigma, xi, kappa)
  nms <- c("sigma", "xi", "kappa")

  wald_covered    <- matrix(NA, nrow = n_rep, ncol = 3)
  profile_covered <- matrix(NA, nrow = n_rep, ncol = 3)
  colnames(wald_covered) <- colnames(profile_covered) <- nms

  for (i in seq_len(n_rep)) {
    y <- rdiscegpd(n, sigma = sigma, xi = xi, kappa = kappa, type = 1)
    df <- data.frame(y = y)

    fit <- tryCatch(
      suppressMessages(egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
           data = df, family = "degpd", degpd.args = list(m = 1))),
      error = function(e) NULL
    )

    if (is.null(fit)) next

    ci_w <- tryCatch(suppressWarnings(confint(fit, method = "wald", level = level)),
                     error = function(e) NULL)
    ci_p <- tryCatch(suppressWarnings(confint(fit, method = "profile", level = level)),
                     error = function(e) NULL)

    if (!is.null(ci_w)) {
      for (k in 1:3) {
        if (!any(is.na(ci_w[k, ])))
          wald_covered[i, k] <- (true_vals[k] >= ci_w[k, 1]) &
                                 (true_vals[k] <= ci_w[k, 2])
      }
    }
    if (!is.null(ci_p)) {
      for (k in 1:3) {
        if (!any(is.na(ci_p[k, ])))
          profile_covered[i, k] <- (true_vals[k] >= ci_p[k, 1]) &
                                    (true_vals[k] <= ci_p[k, 2])
      }
    }
  }

  ok_w <- complete.cases(wald_covered)
  ok_p <- complete.cases(profile_covered)

  list(
    wald_coverage    = if (sum(ok_w) > 0) colMeans(wald_covered[ok_w, , drop = FALSE]) else rep(NA, 3),
    profile_coverage = if (sum(ok_p) > 0) colMeans(profile_covered[ok_p, , drop = FALSE]) else rep(NA, 3),
    n_ok_wald    = sum(ok_w),
    n_ok_profile = sum(ok_p),
    true_vals    = true_vals
  )
}


## ----both-sigma, cache = TRUE-------------------------------------------------
set.seed(301)
sigma_vals_b <- c(1, 3, 10)
xi_fix_b <- 0.2
kappa_fix_b <- 2

res_both_sigma <- lapply(sigma_vals_b, function(s) {
  sim_param_coverage_both(sigma = s, xi = xi_fix_b, kappa = kappa_fix_b)
})


## ----both-xi, cache = TRUE----------------------------------------------------
set.seed(302)
xi_vals_b <- c(0.1, 0.2, 0.5)

res_both_xi <- lapply(xi_vals_b, function(x) {
  sim_param_coverage_both(sigma = 3, xi = x, kappa = kappa_fix_b)
})


## ----both-kappa, cache = TRUE-------------------------------------------------
set.seed(303)
kappa_vals_b <- c(0.5, 2, 5)

res_both_kappa <- lapply(kappa_vals_b, function(k) {
  sim_param_coverage_both(sigma = 3, xi = xi_fix_b, kappa = k)
})


## ----both-tables--------------------------------------------------------------
make_both_table <- function(res_list, varied_name, varied_vals) {
  do.call(rbind, lapply(seq_along(varied_vals), function(j) {
    r <- res_list[[j]]
    data.frame(
      param_value = varied_vals[j],
      wald_sigma = round(r$wald_coverage["sigma"], 3),
      prof_sigma = round(r$profile_coverage["sigma"], 3),
      wald_xi    = round(r$wald_coverage["xi"], 3),
      prof_xi    = round(r$profile_coverage["xi"], 3),
      wald_kappa = round(r$wald_coverage["kappa"], 3),
      prof_kappa = round(r$profile_coverage["kappa"], 3)
    )
  }))
}

cat("Varying sigma:\n")
make_both_table(res_both_sigma, "sigma", sigma_vals_b)

cat("\nVarying xi:\n")
make_both_table(res_both_xi, "xi", xi_vals_b)

cat("\nVarying kappa:\n")
make_both_table(res_both_kappa, "kappa", kappa_vals_b)


## ----both-coverage-plot, fig.width = 9, fig.height = 5------------------------
se_cov_b <- sqrt(0.95 * 0.05 / 100)
par(mfrow = c(1, 3), mar = c(5, 4, 3, 1))

# --- Vary sigma ---
tab_bs <- make_both_table(res_both_sigma, "sigma", sigma_vals_b)
plot(sigma_vals_b, tab_bs$wald_sigma, type = "b", pch = 16,
     col = "steelblue", lwd = 2, ylim = c(0.80, 1),
     xlab = expression(sigma[true]), ylab = "Coverage",
     main = expression("Vary " * sigma))
lines(sigma_vals_b, tab_bs$prof_sigma, type = "b", pch = 16,
      col = "steelblue", lwd = 2, lty = 2)
lines(sigma_vals_b, tab_bs$wald_xi, type = "b", pch = 17,
      col = "firebrick", lwd = 2)
lines(sigma_vals_b, tab_bs$prof_xi, type = "b", pch = 17,
      col = "firebrick", lwd = 2, lty = 2)
lines(sigma_vals_b, tab_bs$wald_kappa, type = "b", pch = 15,
      col = "forestgreen", lwd = 2)
lines(sigma_vals_b, tab_bs$prof_kappa, type = "b", pch = 15,
      col = "forestgreen", lwd = 2, lty = 2)
abline(h = 0.95, lty = 3, col = "grey40")
rect(min(sigma_vals_b) - 1, 0.95 - 2 * se_cov_b,
     max(sigma_vals_b) + 1, 0.95 + 2 * se_cov_b,
     col = adjustcolor("grey70", alpha.f = 0.3), border = NA)

# --- Vary xi ---
tab_bx <- make_both_table(res_both_xi, "xi", xi_vals_b)
plot(xi_vals_b, tab_bx$wald_sigma, type = "b", pch = 16,
     col = "steelblue", lwd = 2, ylim = c(0.80, 1),
     xlab = expression(xi[true]), ylab = "Coverage",
     main = expression("Vary " * xi))
lines(xi_vals_b, tab_bx$prof_sigma, type = "b", pch = 16,
      col = "steelblue", lwd = 2, lty = 2)
lines(xi_vals_b, tab_bx$wald_xi, type = "b", pch = 17,
      col = "firebrick", lwd = 2)
lines(xi_vals_b, tab_bx$prof_xi, type = "b", pch = 17,
      col = "firebrick", lwd = 2, lty = 2)
lines(xi_vals_b, tab_bx$wald_kappa, type = "b", pch = 15,
      col = "forestgreen", lwd = 2)
lines(xi_vals_b, tab_bx$prof_kappa, type = "b", pch = 15,
      col = "forestgreen", lwd = 2, lty = 2)
abline(h = 0.95, lty = 3, col = "grey40")
rect(min(xi_vals_b) - 0.05, 0.95 - 2 * se_cov_b,
     max(xi_vals_b) + 0.05, 0.95 + 2 * se_cov_b,
     col = adjustcolor("grey70", alpha.f = 0.3), border = NA)

# --- Vary kappa ---
tab_bk <- make_both_table(res_both_kappa, "kappa", kappa_vals_b)
plot(kappa_vals_b, tab_bk$wald_sigma, type = "b", pch = 16,
     col = "steelblue", lwd = 2, ylim = c(0.80, 1),
     xlab = expression(kappa[true]), ylab = "Coverage",
     main = expression("Vary " * kappa))
lines(kappa_vals_b, tab_bk$prof_sigma, type = "b", pch = 16,
      col = "steelblue", lwd = 2, lty = 2)
lines(kappa_vals_b, tab_bk$wald_xi, type = "b", pch = 17,
      col = "firebrick", lwd = 2)
lines(kappa_vals_b, tab_bk$prof_xi, type = "b", pch = 17,
      col = "firebrick", lwd = 2, lty = 2)
lines(kappa_vals_b, tab_bk$wald_kappa, type = "b", pch = 15,
      col = "forestgreen", lwd = 2)
lines(kappa_vals_b, tab_bk$prof_kappa, type = "b", pch = 15,
      col = "forestgreen", lwd = 2, lty = 2)
abline(h = 0.95, lty = 3, col = "grey40")
rect(min(kappa_vals_b) - 0.5, 0.95 - 2 * se_cov_b,
     max(kappa_vals_b) + 0.5, 0.95 + 2 * se_cov_b,
     col = adjustcolor("grey70", alpha.f = 0.3), border = NA)

legend("bottomright",
       legend = c(expression(sigma * " Wald"), expression(sigma * " Profile"),
                  expression(xi * " Wald"), expression(xi * " Profile"),
                  expression(kappa * " Wald"), expression(kappa * " Profile")),
       col = rep(c("steelblue", "firebrick", "forestgreen"), each = 2),
       pch = rep(c(16, 17, 15), each = 2), lwd = 2,
       lty = rep(c(1, 2), 3), cex = 0.8)

