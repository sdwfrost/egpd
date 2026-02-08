## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)


## ----egpd-m1------------------------------------------------------------------
library(egpd)
library(bamlss)

set.seed(42)
sigma_true <- 2
xi_true    <- 0.2
kappa_true <- 1.5
n <- 2000

y <- regpd(n, sigma = sigma_true, xi = xi_true, kappa = kappa_true, type = 1)
df <- data.frame(y = y, x = rep(1, n))


## ----egpd-m1-egpd-------------------------------------------------------------
fit_egpd <- egpd(list(lpsi = y ~ 1, xi = ~ 1, lkappa = ~ 1),
                 data = df, family = "egpd", egpd.args = list(m = 1))
pars_egpd <- predict(fit_egpd, type = "response")[1, ]


## ----egpd-m1-bamlss, results = "hide", warning = FALSE------------------------
fit_bamlss <- bamlss(list(y ~ 1, ~ 1, ~ 1),
                     data = df, family = egpd_bamlss(m = 1),
                     verbose = FALSE)


## ----egpd-m1-compare----------------------------------------------------------
sigma_b <- exp(fit_bamlss$parameters$sigma$p)
xi_b    <- exp(fit_bamlss$parameters$xi$p)
kappa_b <- exp(fit_bamlss$parameters$kappa$p)

data.frame(
  parameter = c("sigma", "xi", "kappa"),
  true      = c(sigma_true, xi_true, kappa_true),
  egpd      = round(as.numeric(unlist(pars_egpd)), 4),
  bamlss    = round(c(sigma_b, xi_b, kappa_b), 4),
  row.names = NULL
)


## ----egpd-smooth-sim----------------------------------------------------------
set.seed(7)
n <- 2000
x <- runif(n, 0, 1)
sigma_x <- exp(0.5 + 1.5 * sin(2 * pi * x))
xi_true  <- 0.1
kappa_true <- 1.5

y <- regpd(n, sigma = sigma_x, xi = xi_true, kappa = kappa_true, type = 1)
df <- data.frame(y = y, x = x)


## ----egpd-smooth-egpd---------------------------------------------------------
fit_egpd_s <- egpd(list(lpsi = y ~ s(x, k = 15), xi = ~ 1, lkappa = ~ 1),
                   data = df, family = "egpd", egpd.args = list(m = 1))


## ----egpd-smooth-bamlss, results = "hide", warning = FALSE--------------------
fit_bamlss_s <- bamlss(list(y ~ s(x, k = 15), ~ 1, ~ 1),
                       data = df, family = egpd_bamlss(m = 1),
                       verbose = FALSE)


## ----egpd-smooth-compare, fig.width = 7, fig.height = 5-----------------------
xgrid <- data.frame(x = seq(0, 1, length = 200))
pred_egpd <- predict(fit_egpd_s, newdata = xgrid, type = "response")

pred_bamlss <- predict(fit_bamlss_s, newdata = xgrid, type = "parameter")

plot(xgrid$x, exp(0.5 + 1.5 * sin(2 * pi * xgrid$x)), type = "l",
     lwd = 2, col = "black", ylim = c(0, 15),
     xlab = "x", ylab = expression(sigma(x)),
     main = "Smooth scale recovery: egpd vs bamlss")
lines(xgrid$x, pred_egpd$scale, col = "steelblue", lwd = 2, lty = 2)
lines(xgrid$x, pred_bamlss$sigma, col = "firebrick", lwd = 2, lty = 3)
legend("topright", legend = c("True", "egpd", "bamlss"),
       col = c("black", "steelblue", "firebrick"), lwd = 2, lty = 1:3)


## ----ziegpd-sim---------------------------------------------------------------
set.seed(123)
sigma_true <- 2
xi_true    <- 0.2
kappa_true <- 1.5
pi_true    <- 0.3
n <- 2000

y <- rziegpd(n, pi = pi_true, sigma = sigma_true, xi = xi_true,
             kappa = kappa_true, type = 1)
df <- data.frame(y = y)
cat("Proportion of zeros:", mean(y == 0), "\n")


## ----ziegpd-fit, results = "hide", warning = FALSE----------------------------
fit_zi <- bamlss(list(y ~ 1, ~ 1, ~ 1, ~ 1),
                 data = df, family = ziegpd_bamlss(m = 1),
                 verbose = FALSE)


## ----ziegpd-compare-----------------------------------------------------------
sigma_zb <- exp(fit_zi$parameters$sigma$p)
xi_zb    <- exp(fit_zi$parameters$xi$p)
kappa_zb <- exp(fit_zi$parameters$kappa$p)
pi_zb    <- 1 / (1 + exp(-fit_zi$parameters$pi$p))

data.frame(
  parameter = c("sigma", "xi", "kappa", "pi"),
  true      = c(sigma_true, xi_true, kappa_true, pi_true),
  bamlss    = round(c(sigma_zb, xi_zb, kappa_zb, pi_zb), 4),
  row.names = NULL
)


## ----degpd-sim----------------------------------------------------------------
set.seed(99)
sigma_true <- 3
xi_true    <- 0.15
kappa_true <- 2
n <- 2000

y <- rdiscegpd(n, sigma = sigma_true, xi = xi_true, kappa = kappa_true, type = 1)
df <- data.frame(y = y)
cat("Range of y:", range(y), "\n")


## ----degpd-egpd---------------------------------------------------------------
fit_degpd <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
                  data = df, family = "degpd", degpd.args = list(m = 1))
pars_degpd <- predict(fit_degpd, type = "response")[1, ]


## ----degpd-bamlss, results = "hide", warning = FALSE--------------------------
fit_degpd_b <- bamlss(list(y ~ 1, ~ 1, ~ 1),
                      data = df, family = degpd_bamlss(m = 1),
                      verbose = FALSE)


## ----degpd-compare------------------------------------------------------------
sigma_db <- exp(fit_degpd_b$parameters$sigma$p)
xi_db    <- exp(fit_degpd_b$parameters$xi$p)
kappa_db <- exp(fit_degpd_b$parameters$kappa$p)

data.frame(
  parameter = c("sigma", "xi", "kappa"),
  true      = c(sigma_true, xi_true, kappa_true),
  egpd      = round(as.numeric(unlist(pars_degpd)), 4),
  bamlss    = round(c(sigma_db, xi_db, kappa_db), 4),
  row.names = NULL
)


## ----zidegpd-sim--------------------------------------------------------------
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


## ----zidegpd-egpd-------------------------------------------------------------
fit_zidegpd <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1, logitpi = ~ 1),
                    data = df, family = "zidegpd", zidegpd.args = list(m = 1))
pars_zidegpd <- predict(fit_zidegpd, type = "response")[1, ]


## ----zidegpd-bamlss, results = "hide", warning = FALSE------------------------
fit_zidegpd_b <- bamlss(list(y ~ 1, ~ 1, ~ 1, ~ 1),
                        data = df, family = zidegpd_bamlss(m = 1),
                        verbose = FALSE)


## ----zidegpd-compare----------------------------------------------------------
sigma_zdb <- exp(fit_zidegpd_b$parameters$sigma$p)
xi_zdb    <- exp(fit_zidegpd_b$parameters$xi$p)
kappa_zdb <- exp(fit_zidegpd_b$parameters$kappa$p)
pi_zdb    <- 1 / (1 + exp(-fit_zidegpd_b$parameters$pi$p))

data.frame(
  parameter = c("sigma", "xi", "kappa", "pi"),
  true      = c(sigma_true, xi_true, kappa_true, pi_true),
  egpd      = round(as.numeric(unlist(pars_zidegpd)), 4),
  bamlss    = round(c(sigma_zdb, xi_zdb, kappa_zdb, pi_zdb), 4),
  row.names = NULL
)

