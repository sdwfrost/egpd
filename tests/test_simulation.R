#!/usr/bin/env Rscript
#
# Simulation-based tests for the egpd package
#
# For each model class (EGPD, DEGPD, ZIDEGPD) and model variants (1--6),
# this script:
#   1. Simulates data from known parameters
#   2. Fits the model
#   3. Computes bias (estimate - truth on link scale)
#   4. Computes coverage (does 95% Wald CI contain the true value?)
#   5. Checks SE calibration (mean model SE vs empirical SD)
#
# Distribution type mapping:
#   Fitting model 1 -> distribution type 1  (G(u) = u^kappa)
#   Fitting model 2 -> distribution type 6  (G(u) = p*u^kappa1 + (1-p)*u^kappa2)
#   Fitting model 3 -> distribution type 4  (G(u) = 1 - barB(1-u; 1/delta, 2)^(1/delta))
#   Fitting model 4 -> distribution type 5  (G(u) = [1 - barB(1-u; 1/delta, 2)^(1/delta)]^(kappa/2))
#   Fitting model 5 -> distribution type 2  (truncated normal G-function)
#   Fitting model 6 -> distribution type 3  (truncated beta G-function)
#
# Usage:
#   Rscript tests/test_simulation.R              # default: nsim=100, n=2000
#   Rscript tests/test_simulation.R 50 1000      # nsim=50, n=1000

library(egpd)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args) >= 1) as.integer(args[1]) else 100L
n    <- if (length(args) >= 2) as.integer(args[2]) else 2000L

logit <- function(p) log(p / (1 - p))

cat(sprintf("Running simulation tests: nsim=%d, n=%d\n\n", nsim, n))

# ---------------------------------------------------------------------------
# Helper: run one simulation study
# ---------------------------------------------------------------------------

run_simulation <- function(model_name, family, model_num, true_link,
                           param_names, sim_fn, fmla_fn,
                           family_args_name, nsim, n) {

  npar    <- length(true_link)
  estimates <- matrix(NA, nsim, npar, dimnames = list(NULL, param_names))
  ses       <- matrix(NA, nsim, npar, dimnames = list(NULL, param_names))
  converged <- logical(nsim)

  cat(sprintf("  Fitting %s ...", model_name))

  for (i in seq_len(nsim)) {
    tryCatch({
      set.seed(i)

      # Simulate data
      df <- data.frame(x = rnorm(n))
      df$y <- sim_fn(n)

      # Build formula and family args
      fmla <- fmla_fn()
      fargs <- list()
      fargs[[family_args_name]] <- list(m = model_num)

      # Fit model
      fit_args <- c(
        list(formula = fmla, data = df, family = family, trace = 0),
        fargs
      )
      m <- do.call(egpd, fit_args)

      # Extract estimates and SEs (handle non-PD Hessian)
      estimates[i, ] <- m$coefficients
      se_vals <- suppressWarnings(sqrt(diag(m$Vp)))
      se_vals[is.nan(se_vals)] <- NA
      ses[i, ] <- se_vals
      converged[i] <- TRUE
    }, error = function(e) {
      converged[i] <<- FALSE
    })
  }

  ok <- converged & complete.cases(estimates) & complete.cases(ses)
  nconv <- sum(converged)
  nok   <- sum(ok)

  cat(sprintf(" %d/%d converged, %d usable\n", nconv, nsim, nok))

  if (nok < 10) {
    cat(sprintf("    *** Too few usable fits -- skipping diagnostics\n\n"))
    return(invisible(NULL))
  }

  est <- estimates[ok, , drop = FALSE]
  se  <- ses[ok, , drop = FALSE]

  # Bias
  bias     <- colMeans(est) - true_link
  rel_bias <- ifelse(abs(true_link) > 0.01, bias / abs(true_link) * 100, NA)

  # RMSE
  rmse <- sqrt(colMeans(sweep(est, 2, true_link)^2))

  # Coverage: does 95% Wald CI contain true value?
  lo <- est - 1.96 * se
  hi <- est + 1.96 * se
  coverage <- colMeans(sweep(lo, 2, true_link, "<=") &
                       sweep(hi, 2, true_link, ">=")) * 100

  # Mean SE vs empirical SD
  mean_se <- colMeans(se)
  emp_sd  <- apply(est, 2, sd)
  se_ratio <- mean_se / emp_sd

  # Print results
  cat(sprintf("    %-12s %10s %10s %10s %10s %10s %10s\n",
              "Parameter", "True", "Mean.Est", "Bias", "RMSE",
              "Cover95%", "SE/EmpSD"))
  for (j in seq_len(npar)) {
    cat(sprintf("    %-12s %10.4f %10.4f %10.4f %10.4f %10.1f %10.2f\n",
                param_names[j], true_link[j], colMeans(est)[j],
                bias[j], rmse[j], coverage[j], se_ratio[j]))
  }
  cat("\n")

  invisible(data.frame(
    model       = model_name,
    parameter   = param_names,
    true        = true_link,
    mean_est    = colMeans(est),
    bias        = bias,
    rel_bias_pct = rel_bias,
    rmse        = rmse,
    coverage    = coverage,
    se_ratio    = se_ratio,
    n_converged = nconv,
    n_usable    = nok,
    stringsAsFactors = FALSE
  ))
}

# ===========================================================================
# EGPD (Continuous Extended GPD)
# ===========================================================================

cat("============================================================\n")
cat("EGPD (Continuous Extended GPD)\n")
cat("============================================================\n\n")

# True parameters (natural scale)
egpd_sigma <- 2; egpd_xi <- 0.1; egpd_kappa <- 1.5
egpd4_delta <- 2.0; egpd4_kappa <- 2.0

# EGPD model 1: G(u) = u^kappa -> distribution type 1
res_egpd1 <- run_simulation(
  model_name = "EGPD-1", family = "egpd", model_num = 1,
  true_link   = c(log(egpd_sigma), egpd_xi, log(egpd_kappa)),
  param_names = c("lpsi", "xi", "lkappa"),
  sim_fn   = function(nn) regpd(nn, sigma = egpd_sigma, xi = egpd_xi,
                                kappa = egpd_kappa, type = 1),
  fmla_fn  = function() list(lpsi = y ~ 1, xi = ~1, lkappa = ~1),
  family_args_name = "egpd.args", nsim = nsim, n = n
)

# EGPD model 2: G(u) = p*u^kappa1 + (1-p)*u^kappa2 -> distribution type 6
# Reparameterized: ldkappa = log(kappa2 - kappa1)
egpd2_kappa1 <- 1.5; egpd2_kappa2 <- 3.0; egpd2_p <- 0.6
res_egpd2 <- run_simulation(
  model_name = "EGPD-2", family = "egpd", model_num = 2,
  true_link   = c(log(egpd_sigma), egpd_xi, log(egpd2_kappa1), log(egpd2_kappa2 - egpd2_kappa1), logit(egpd2_p)),
  param_names = c("lpsi", "xi", "lkappa1", "ldkappa", "logitp"),
  sim_fn   = function(nn) regpd(nn, sigma = egpd_sigma, xi = egpd_xi,
                                kappa = egpd2_kappa1, delta = egpd2_kappa2,
                                prob = egpd2_p, type = 6),
  fmla_fn  = function() list(lpsi = y ~ 1, xi = ~1, lkappa1 = ~1, ldkappa = ~1, logitp = ~1),
  family_args_name = "egpd.args", nsim = nsim, n = n
)

# EGPD model 3: G via incomplete beta -> distribution type 4
egpd3_delta <- 0.5
res_egpd3 <- run_simulation(
  model_name = "EGPD-3", family = "egpd", model_num = 3,
  true_link   = c(log(egpd_sigma), egpd_xi, log(egpd3_delta)),
  param_names = c("lpsi", "xi", "ldelta"),
  sim_fn   = function(nn) regpd(nn, sigma = egpd_sigma, xi = egpd_xi,
                                delta = egpd3_delta, type = 4),
  fmla_fn  = function() list(lpsi = y ~ 1, xi = ~1, ldelta = ~1),
  family_args_name = "egpd.args", nsim = nsim, n = n
)

# EGPD model 4: G via beta + power -> distribution type 5
res_egpd4 <- run_simulation(
  model_name = "EGPD-4", family = "egpd", model_num = 4,
  true_link   = c(log(egpd_sigma), egpd_xi, log(egpd4_delta), log(egpd4_kappa)),
  param_names = c("lpsi", "xi", "ldelta", "lkappa"),
  sim_fn   = function(nn) regpd(nn, sigma = egpd_sigma, xi = egpd_xi,
                                delta = egpd4_delta, kappa = egpd4_kappa, type = 5),
  fmla_fn  = function() list(lpsi = y ~ 1, xi = ~1, ldelta = ~1, lkappa = ~1),
  family_args_name = "egpd.args", nsim = nsim, n = n
)

# EGPD model 5: truncated normal G -> distribution type 2
egpd5_kappa <- 2.0
res_egpd5 <- run_simulation(
  model_name = "EGPD-5", family = "egpd", model_num = 5,
  true_link   = c(log(egpd_sigma), egpd_xi, log(egpd5_kappa)),
  param_names = c("lpsi", "xi", "lkappa"),
  sim_fn   = function(nn) regpd(nn, sigma = egpd_sigma, xi = egpd_xi,
                                kappa = egpd5_kappa, type = 2),
  fmla_fn  = function() list(lpsi = y ~ 1, xi = ~1, lkappa = ~1),
  family_args_name = "egpd.args", nsim = nsim, n = n
)

# EGPD model 6: truncated beta G -> distribution type 3
egpd6_kappa <- 2.0
res_egpd6 <- run_simulation(
  model_name = "EGPD-6", family = "egpd", model_num = 6,
  true_link   = c(log(egpd_sigma), egpd_xi, log(egpd6_kappa)),
  param_names = c("lpsi", "xi", "lkappa"),
  sim_fn   = function(nn) regpd(nn, sigma = egpd_sigma, xi = egpd_xi,
                                kappa = egpd6_kappa, type = 3),
  fmla_fn  = function() list(lpsi = y ~ 1, xi = ~1, lkappa = ~1),
  family_args_name = "egpd.args", nsim = nsim, n = n
)

# ===========================================================================
# DEGPD (Discrete Extended GPD)
# ===========================================================================

cat("============================================================\n")
cat("DEGPD (Discrete Extended GPD)\n")
cat("============================================================\n\n")

degpd_sigma <- 3; degpd_xi <- 0.2; degpd_kappa <- 2
degpd3_delta <- 0.5
degpd4_delta <- 2.0; degpd4_kappa <- 2.0

# DEGPD model 1 -> distribution type 1
res_degpd1 <- run_simulation(
  model_name = "DEGPD-1", family = "degpd", model_num = 1,
  true_link   = c(log(degpd_sigma), log(degpd_xi), log(degpd_kappa)),
  param_names = c("lsigma", "lxi", "lkappa"),
  sim_fn   = function(nn) rdiscegpd(nn, sigma = degpd_sigma, xi = degpd_xi,
                                     kappa = degpd_kappa, type = 1),
  fmla_fn  = function() list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1),
  family_args_name = "degpd.args", nsim = nsim, n = n
)

# DEGPD model 2: mixture -> distribution type 6
# Reparameterized: ldkappa = log(kappa2 - kappa1)
degpd2_kappa1 <- 1.5; degpd2_kappa2 <- 3.0; degpd2_p <- 0.6
res_degpd2 <- run_simulation(
  model_name = "DEGPD-2", family = "degpd", model_num = 2,
  true_link   = c(log(degpd_sigma), log(degpd_xi), log(degpd2_kappa1), log(degpd2_kappa2 - degpd2_kappa1), logit(degpd2_p)),
  param_names = c("lsigma", "lxi", "lkappa1", "ldkappa", "logitp"),
  sim_fn   = function(nn) rdiscegpd(nn, sigma = degpd_sigma, xi = degpd_xi,
                                     kappa = degpd2_kappa1, delta = degpd2_kappa2,
                                     prob = degpd2_p, type = 6),
  fmla_fn  = function() list(lsigma = y ~ 1, lxi = ~1, lkappa1 = ~1, ldkappa = ~1, logitp = ~1),
  family_args_name = "degpd.args", nsim = nsim, n = n
)

# DEGPD model 3 -> distribution type 4
res_degpd3 <- run_simulation(
  model_name = "DEGPD-3", family = "degpd", model_num = 3,
  true_link   = c(log(degpd_sigma), log(degpd_xi), log(degpd3_delta)),
  param_names = c("lsigma", "lxi", "ldelta"),
  sim_fn   = function(nn) rdiscegpd(nn, sigma = degpd_sigma, xi = degpd_xi,
                                     delta = degpd3_delta, type = 4),
  fmla_fn  = function() list(lsigma = y ~ 1, lxi = ~1, ldelta = ~1),
  family_args_name = "degpd.args", nsim = nsim, n = n
)

# DEGPD model 4 -> distribution type 5
res_degpd4 <- run_simulation(
  model_name = "DEGPD-4", family = "degpd", model_num = 4,
  true_link   = c(log(degpd_sigma), log(degpd_xi), log(degpd4_delta), log(degpd4_kappa)),
  param_names = c("lsigma", "lxi", "ldelta", "lkappa"),
  sim_fn   = function(nn) rdiscegpd(nn, sigma = degpd_sigma, xi = degpd_xi,
                                     delta = degpd4_delta, kappa = degpd4_kappa, type = 5),
  fmla_fn  = function() list(lsigma = y ~ 1, lxi = ~1, ldelta = ~1, lkappa = ~1),
  family_args_name = "degpd.args", nsim = nsim, n = n
)

# DEGPD model 5: truncated normal -> distribution type 2
degpd5_kappa <- 2.0
res_degpd5 <- run_simulation(
  model_name = "DEGPD-5", family = "degpd", model_num = 5,
  true_link   = c(log(degpd_sigma), log(degpd_xi), log(degpd5_kappa)),
  param_names = c("lsigma", "lxi", "lkappa"),
  sim_fn   = function(nn) rdiscegpd(nn, sigma = degpd_sigma, xi = degpd_xi,
                                     kappa = degpd5_kappa, type = 2),
  fmla_fn  = function() list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1),
  family_args_name = "degpd.args", nsim = nsim, n = n
)

# DEGPD model 6: truncated beta -> distribution type 3
degpd6_kappa <- 2.0
res_degpd6 <- run_simulation(
  model_name = "DEGPD-6", family = "degpd", model_num = 6,
  true_link   = c(log(degpd_sigma), log(degpd_xi), log(degpd6_kappa)),
  param_names = c("lsigma", "lxi", "lkappa"),
  sim_fn   = function(nn) rdiscegpd(nn, sigma = degpd_sigma, xi = degpd_xi,
                                     kappa = degpd6_kappa, type = 3),
  fmla_fn  = function() list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1),
  family_args_name = "degpd.args", nsim = nsim, n = n
)

# ===========================================================================
# ZIDEGPD (Zero-Inflated Discrete Extended GPD)
# ===========================================================================

cat("============================================================\n")
cat("ZIDEGPD (Zero-Inflated Discrete Extended GPD)\n")
cat("============================================================\n\n")

zi_sigma <- 3; zi_xi <- 0.2; zi_kappa <- 2; zi_pi <- 0.3
zi3_delta <- 0.5
zi4_delta <- 2.0; zi4_kappa <- 2.0

# ZIDEGPD model 1 -> distribution type 1
res_zi1 <- run_simulation(
  model_name = "ZIDEGPD-1", family = "zidegpd", model_num = 1,
  true_link   = c(log(zi_sigma), log(zi_xi), log(zi_kappa), logit(zi_pi)),
  param_names = c("lsigma", "lxi", "lkappa", "logitpi"),
  sim_fn   = function(nn) rzidiscegpd(nn, pi = zi_pi, sigma = zi_sigma, xi = zi_xi,
                                       kappa = zi_kappa, type = 1),
  fmla_fn  = function() list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1, lpi = ~1),
  family_args_name = "zidegpd.args", nsim = nsim, n = n
)

# ZIDEGPD model 2: mixture with ZI -> distribution type 6
# Reparameterized: ldkappa = log(kappa2 - kappa1)
zi2_kappa1 <- 1.5; zi2_kappa2 <- 3.0; zi2_p <- 0.6
res_zi2 <- run_simulation(
  model_name = "ZIDEGPD-2", family = "zidegpd", model_num = 2,
  true_link   = c(log(zi_sigma), log(zi_xi), log(zi2_kappa1), log(zi2_kappa2 - zi2_kappa1), logit(zi2_p), logit(zi_pi)),
  param_names = c("lsigma", "lxi", "lkappa1", "ldkappa", "logitp", "logitpi"),
  sim_fn   = function(nn) rzidiscegpd(nn, pi = zi_pi, sigma = zi_sigma, xi = zi_xi,
                                       kappa = zi2_kappa1, delta = zi2_kappa2,
                                       prob = zi2_p, type = 6),
  fmla_fn  = function() list(lsigma = y ~ 1, lxi = ~1, lkappa1 = ~1, ldkappa = ~1, logitp = ~1, lpi = ~1),
  family_args_name = "zidegpd.args", nsim = nsim, n = n
)

# ZIDEGPD model 3 -> distribution type 4
res_zi3 <- run_simulation(
  model_name = "ZIDEGPD-3", family = "zidegpd", model_num = 3,
  true_link   = c(log(zi_sigma), log(zi_xi), log(zi3_delta), logit(zi_pi)),
  param_names = c("lsigma", "lxi", "ldelta", "logitpi"),
  sim_fn   = function(nn) rzidiscegpd(nn, pi = zi_pi, sigma = zi_sigma, xi = zi_xi,
                                       delta = zi3_delta, type = 4),
  fmla_fn  = function() list(lsigma = y ~ 1, lxi = ~1, ldelta = ~1, lpi = ~1),
  family_args_name = "zidegpd.args", nsim = nsim, n = n
)

# ZIDEGPD model 4 -> distribution type 5
res_zi4 <- run_simulation(
  model_name = "ZIDEGPD-4", family = "zidegpd", model_num = 4,
  true_link   = c(log(zi_sigma), log(zi_xi), log(zi4_kappa), log(zi4_delta), logit(zi_pi)),
  param_names = c("lsigma", "lxi", "lkappa", "ldelta", "logitpi"),
  sim_fn   = function(nn) rzidiscegpd(nn, pi = zi_pi, sigma = zi_sigma, xi = zi_xi,
                                       delta = zi4_delta, kappa = zi4_kappa, type = 5),
  fmla_fn  = function() list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1, ldelta = ~1, lpi = ~1),
  family_args_name = "zidegpd.args", nsim = nsim, n = n
)

# ZIDEGPD model 5: truncated normal with ZI -> distribution type 2
res_zi5 <- run_simulation(
  model_name = "ZIDEGPD-5", family = "zidegpd", model_num = 5,
  true_link   = c(log(zi_sigma), log(zi_xi), log(zi_kappa), logit(zi_pi)),
  param_names = c("lsigma", "lxi", "lkappa", "logitpi"),
  sim_fn   = function(nn) rzidiscegpd(nn, pi = zi_pi, sigma = zi_sigma, xi = zi_xi,
                                       kappa = zi_kappa, type = 2),
  fmla_fn  = function() list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1, lpi = ~1),
  family_args_name = "zidegpd.args", nsim = nsim, n = n
)

# ZIDEGPD model 6: truncated beta with ZI -> distribution type 3
res_zi6 <- run_simulation(
  model_name = "ZIDEGPD-6", family = "zidegpd", model_num = 6,
  true_link   = c(log(zi_sigma), log(zi_xi), log(zi_kappa), logit(zi_pi)),
  param_names = c("lsigma", "lxi", "lkappa", "logitpi"),
  sim_fn   = function(nn) rzidiscegpd(nn, pi = zi_pi, sigma = zi_sigma, xi = zi_xi,
                                       kappa = zi_kappa, type = 3),
  fmla_fn  = function() list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1, lpi = ~1),
  family_args_name = "zidegpd.args", nsim = nsim, n = n
)

# ===========================================================================
# Combined summary
# ===========================================================================

cat("============================================================\n")
cat("OVERALL SUMMARY\n")
cat("============================================================\n\n")

all_results <- do.call(rbind, Filter(Negate(is.null),
  list(res_egpd1, res_egpd2, res_egpd3, res_egpd4, res_egpd5, res_egpd6,
       res_degpd1, res_degpd2, res_degpd3, res_degpd4, res_degpd5, res_degpd6,
       res_zi1, res_zi2, res_zi3, res_zi4, res_zi5, res_zi6)))

if (!is.null(all_results) && nrow(all_results) > 0) {
  # Per-model summary
  models <- unique(all_results$model)
  cat(sprintf("%-12s %8s %12s %12s %12s\n",
              "Model", "Usable", "MaxAbsBias", "MinCover%", "SE.ratio"))
  cat(paste(rep("-", 60), collapse = ""), "\n")
  for (mod in models) {
    sub <- all_results[all_results$model == mod, ]
    cat(sprintf("%-12s %7d%% %12.4f %11.1f%% %8.2f-%.2f\n",
                mod,
                round(sub$n_usable[1] / nsim * 100),
                max(abs(sub$bias)),
                min(sub$coverage),
                min(sub$se_ratio),
                max(sub$se_ratio)))
  }

  cat("\nDiagnostics:\n")

  # Flag large bias (>0.1 on link scale)
  bad_bias <- all_results[abs(all_results$bias) > 0.2, ]
  if (nrow(bad_bias) > 0) {
    cat("  NOTE: Parameters with |bias| > 0.2 on link scale:\n")
    for (i in seq_len(nrow(bad_bias))) {
      cat(sprintf("    %s / %s: bias=%.4f (true=%.4f, est=%.4f)\n",
                  bad_bias$model[i], bad_bias$parameter[i],
                  bad_bias$bias[i], bad_bias$true[i], bad_bias$mean_est[i]))
    }
  } else {
    cat("  All parameters have |bias| < 0.2 on link scale.\n")
  }

  # Flag poor coverage
  bad_cov <- all_results[all_results$coverage < 85, ]
  if (nrow(bad_cov) > 0) {
    cat("  NOTE: Parameters with <85% coverage:\n")
    for (i in seq_len(nrow(bad_cov))) {
      cat(sprintf("    %s / %s: %.1f%%\n",
                  bad_cov$model[i], bad_cov$parameter[i], bad_cov$coverage[i]))
    }
  } else {
    cat("  All parameters have >=85% coverage.\n")
  }

  # Flag bad SE calibration
  bad_se <- all_results[all_results$se_ratio < 0.7 | all_results$se_ratio > 1.5, ]
  if (nrow(bad_se) > 0) {
    cat("  NOTE: Parameters with SE ratio outside [0.7, 1.5]:\n")
    for (i in seq_len(nrow(bad_se))) {
      cat(sprintf("    %s / %s: SE/EmpSD = %.2f\n",
                  bad_se$model[i], bad_se$parameter[i], bad_se$se_ratio[i]))
    }
  } else {
    cat("  All SE ratios within [0.7, 1.5].\n")
  }
}

cat("\nDone.\n")
