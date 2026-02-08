## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)


## ----degpd1-sim---------------------------------------------------------------
library(egpd)
library(gamlss)

set.seed(99)
sigma_true <- 3
xi_true    <- 0.15
kappa_true <- 2
n <- 2000

y <- rdiscegpd(n, sigma = sigma_true, xi = xi_true, kappa = kappa_true, type = 1)
df <- data.frame(y = y)
cat("Range of y:", range(y), "\n")


## ----degpd1-egpd--------------------------------------------------------------
fit_egpd <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
                 data = df, family = "degpd", degpd.args = list(m = 1))
pars_egpd <- predict(fit_egpd, type = "response")[1, ]


## ----degpd1-gamlss, warning = FALSE-------------------------------------------
fit_gamlss <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                     data = df, family = DEGPD1(),
                     control = gamlss.control(n.cyc = 200, trace = FALSE))


## ----degpd1-compare-----------------------------------------------------------
mu_g    <- exp(coef(fit_gamlss, what = "mu"))
sigma_g <- exp(coef(fit_gamlss, what = "sigma"))
nu_g    <- exp(coef(fit_gamlss, what = "nu"))

data.frame(
  parameter = c("sigma", "xi", "kappa"),
  true      = c(sigma_true, xi_true, kappa_true),
  egpd      = round(as.numeric(unlist(pars_egpd)), 4),
  gamlss    = round(c(mu_g, sigma_g, nu_g), 4),
  row.names = NULL
)


## ----zidegpd1-sim-------------------------------------------------------------
set.seed(77)
sigma_true <- 3
xi_true    <- 0.15
kappa_true <- 2
pi_true    <- 0.25
n <- 2000

y <- rzidiscegpd(n, pi = pi_true, sigma = sigma_true, xi = xi_true,
                 kappa = kappa_true, type = 1)
df <- data.frame(y = y)
cat("Proportion of zeros:", mean(y == 0), "\n")


## ----zidegpd1-egpd------------------------------------------------------------
fit_zidegpd <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1, logitpi = ~ 1),
                    data = df, family = "zidegpd", zidegpd.args = list(m = 1))
pars_zidegpd <- predict(fit_zidegpd, type = "response")[1, ]


## ----zidegpd1-gamlss, warning = FALSE-----------------------------------------
fit_zi_gamlss <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                        tau.formula = ~ 1, data = df, family = ZIDEGPD1(),
                        control = gamlss.control(n.cyc = 200, trace = FALSE))


## ----zidegpd1-compare---------------------------------------------------------
mu_zg    <- exp(coef(fit_zi_gamlss, what = "mu"))
sigma_zg <- exp(coef(fit_zi_gamlss, what = "sigma"))
nu_zg    <- exp(coef(fit_zi_gamlss, what = "nu"))
tau_zg   <- plogis(coef(fit_zi_gamlss, what = "tau"))

data.frame(
  parameter = c("sigma", "xi", "kappa", "pi"),
  true      = c(sigma_true, xi_true, kappa_true, pi_true),
  egpd      = round(as.numeric(unlist(pars_zidegpd)), 4),
  gamlss    = round(c(mu_zg, sigma_zg, nu_zg, tau_zg), 4),
  row.names = NULL
)


## ----egpd1-sim----------------------------------------------------------------
set.seed(42)
sigma_true <- 2
xi_true    <- 0.2
kappa_true <- 1.5
n <- 2000

y <- regpd(n, sigma = sigma_true, xi = xi_true, kappa = kappa_true, type = 1)
df <- data.frame(y = y, x = rep(1, n))


## ----egpd1-egpd---------------------------------------------------------------
fit_egpd_c <- egpd(list(lpsi = y ~ 1, xi = ~ 1, lkappa = ~ 1),
                   data = df, family = "egpd", egpd.args = list(m = 1))
pars_egpd_c <- predict(fit_egpd_c, type = "response")[1, ]


## ----egpd1-gamlss, warning = FALSE--------------------------------------------
fit_gamlss_c <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                       data = df, family = EGPD1(),
                       control = gamlss.control(n.cyc = 200, trace = FALSE))


## ----egpd1-compare------------------------------------------------------------
mu_gc    <- exp(coef(fit_gamlss_c, what = "mu"))
sigma_gc <- exp(coef(fit_gamlss_c, what = "sigma"))
nu_gc    <- exp(coef(fit_gamlss_c, what = "nu"))

data.frame(
  parameter = c("sigma", "xi", "kappa"),
  true      = c(sigma_true, xi_true, kappa_true),
  egpd      = round(as.numeric(unlist(pars_egpd_c)), 4),
  gamlss    = round(c(mu_gc, sigma_gc, nu_gc), 4),
  row.names = NULL
)

