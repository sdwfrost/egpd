## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)


## ----multiscale-sim-----------------------------------------------------------
library(egpd)
set.seed(2)

durations <- c(1L, 2L, 4L, 8L)
sigma_coef_true <- c(log(0.6), 0.12)
lambda_coef_true <- c(log(0.8), 0.25)
kappa_true <- 0.35
xi_true <- 0.2

sigma_true <- exp(sigma_coef_true[1] + sigma_coef_true[2] * log(durations))
lambda_true <- exp(lambda_coef_true[1] + lambda_coef_true[2] * log(durations))

agg_samples <- setNames(lapply(seq_along(durations), function(i) {
  rcpegpd(
    300,
    sigma = sigma_true[i],
    xi = xi_true,
    kappa = kappa_true,
    lambda = lambda_true[i],
    type = 1
  )
}), durations)

data.frame(
  duration = durations,
  n = vapply(agg_samples, length, integer(1)),
  prop_zero = round(vapply(agg_samples, function(x) mean(x == 0), numeric(1)), 3),
  mean = round(vapply(agg_samples, mean, numeric(1)), 2)
)


## ----multiscale-fit-----------------------------------------------------------
fit_agg <- fit_aggregate_egpd(
  agg_samples,
  durations = durations,
  family = "cpegpd",
  p = 1,
  q = 1,
  start = list(
    sigma_coef = sigma_coef_true,
    lambda_coef = lambda_coef_true,
    kappa = kappa_true,
    xi = xi_true
  ),
  optim.method = "BFGS",
  control = list(maxit = 200)
)

fit_agg


## ----multiscale-summary-------------------------------------------------------
summary(fit_agg)


## ----multiscale-parameters----------------------------------------------------
fit_pars <- predict(fit_agg, type = "parameter")

comparison <- data.frame(
  duration = durations,
  sigma_true = sigma_true,
  sigma_hat = fit_pars$sigma,
  lambda_true = lambda_true,
  lambda_hat = fit_pars$lambda,
  kappa_true = kappa_true,
  kappa_hat = fit_pars$kappa,
  xi_true = xi_true,
  xi_hat = fit_pars$xi
)

round(comparison, 3)


## ----multiscale-predict-------------------------------------------------------
new_durations <- c(1L, 2L, 4L, 8L, 16L)

round(
  predict(fit_agg, durations = new_durations, type = "parameter"),
  3
)


## ----multiscale-quantiles-----------------------------------------------------
round(
  predict(
    fit_agg,
    durations = new_durations,
    type = "quantile",
    prob = c(0.5, 0.9, 0.99)
  ),
  2
)


## ----multiscale-plot, fig.height = 4------------------------------------------
plot(fit_agg)


## ----series-fit, eval = FALSE-------------------------------------------------
# set.seed(4)
# x <- rcpegpd(300, sigma = 0.8, xi = 0.2, kappa = 0.35, lambda = 0.7, type = 1)
# 
# fit_from_series <- fit_aggregate_egpd(
#   x,
#   durations = c(1L, 2L, 4L, 8L),
#   family = "cpegpd",
#   p = 1,
#   q = 1,
#   overlap = TRUE
# )
# 
# predict(fit_from_series, type = "parameter")


## ----discrete-fit, eval = FALSE-----------------------------------------------
# set.seed(3)
# count_samples <- setNames(lapply(durations, function(d) {
#   rcpdegpd(
#     250,
#     sigma = exp(log(1.2) + 0.08 * log(d)),
#     xi = 0.18,
#     kappa = 0.4,
#     lambda = exp(log(0.9) + 0.22 * log(d)),
#     type = 1
#   )
# }), durations)
# 
# fit_counts <- fit_aggregate_egpd(
#   count_samples,
#   durations = durations,
#   family = "cpdegpd",
#   p = 1,
#   q = 1
# )
# 
# predict(fit_counts, type = "quantile", prob = c(0.5, 0.9))

