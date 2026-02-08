## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)


## ----egpd-intercept-----------------------------------------------------------
library(egpd)
set.seed(1)

# True parameters
sigma_true <- 2
xi_true    <- 0.2
kappa_true <- 1.5

# Simulate
n <- 2000
y <- regpd(n, sigma = sigma_true, xi = xi_true, kappa = kappa_true, type = 1)
df <- data.frame(y = y, x = rep(1, n))

# Fit
fit <- egpd(list(lpsi = y ~ 1, xi = ~ 1, lkappa = ~ 1),
            data = df, family = "egpd", egpd.args = list(m = 1))
summary(fit)


## ----egpd-intercept-check-----------------------------------------------------
pars <- predict(fit, type = "response")[1, ]
truth <- c(scale = sigma_true, shape = xi_true, kappa = kappa_true)
cbind(true = truth, fitted = round(unlist(pars), 4))


## ----egpd-intercept-qq, fig.width = 5, fig.height = 5-------------------------
set.seed(1)
r <- rqresid(fit)
qqnorm(r, main = "Q-Q Plot (EGPD-1)", pch = 20, col = "grey60")
qqline(r, col = "red")


## ----egpd-smooth-sim----------------------------------------------------------
set.seed(42)
n <- 2000
x <- runif(n, 0, 1)

# Scale varies smoothly with x
sigma_x <- exp(0.5 + 1.5 * sin(2 * pi * x))
xi_true  <- 0.1
kappa_true <- 1.5

y <- regpd(n, sigma = sigma_x, xi = xi_true, kappa = kappa_true, type = 1)
df <- data.frame(y = y, x = x)


## ----egpd-smooth-fit----------------------------------------------------------
fit_smooth <- egpd(list(lpsi = y ~ s(x, k = 15), xi = ~ 1, lkappa = ~ 1),
                   data = df, family = "egpd", egpd.args = list(m = 1))
summary(fit_smooth)


## ----egpd-smooth-plot, fig.width = 7, fig.height = 5--------------------------
plot(fit_smooth)


## ----egpd-smooth-compare, fig.width = 7, fig.height = 5-----------------------
xgrid <- data.frame(x = seq(0, 1, length = 200))
pred <- predict(fit_smooth, newdata = xgrid, type = "response")

plot(xgrid$x, exp(0.5 + 1.5 * sin(2 * pi * xgrid$x)), type = "l",
     lwd = 2, col = "black", ylim = c(0, 15),
     xlab = "x", ylab = expression(sigma(x)),
     main = "Recovered smooth scale function")
lines(xgrid$x, pred$scale, col = "steelblue", lwd = 2, lty = 2)
legend("topright", legend = c("True", "Fitted"),
       col = c("black", "steelblue"), lwd = 2, lty = c(1, 2))


## ----degpd-intercept----------------------------------------------------------
set.seed(2)
sigma_true <- 3
xi_true    <- 0.3
kappa_true <- 2.0

n <- 2000
y <- rdiscegpd(n, sigma = sigma_true, xi = xi_true, kappa = kappa_true,
               type = 1)
df <- data.frame(y = y, x = rep(1, n))

fit_d <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
              data = df, family = "degpd", degpd.args = list(m = 1))
summary(fit_d)


## ----degpd-intercept-check----------------------------------------------------
pars_d <- predict(fit_d, type = "response")[1, ]
truth_d <- c(scale = sigma_true, shape = xi_true, kappa = kappa_true)
cbind(true = truth_d, fitted = round(unlist(pars_d), 4))


## ----degpd-pmf, fig.width = 7, fig.height = 5---------------------------------
xvals <- 0:20
emp_pmf <- tabulate(y + 1, nbins = max(xvals) + 1) / n
fit_pmf <- ddiscegpd(xvals, sigma = pars_d$scale[1], xi = pars_d$shape[1],
                     kappa = pars_d$kappa[1], type = 1)

plot(xvals, emp_pmf[seq_along(xvals)], type = "h", lwd = 3, col = "grey60",
     main = "Empirical vs fitted PMF (DEGPD-1)",
     xlab = "Count", ylab = "Probability")
lines(xvals + 0.2, fit_pmf, type = "h", lwd = 3, col = "steelblue")
legend("topright", legend = c("Empirical", "Fitted"),
       col = c("grey60", "steelblue"), lwd = 3)


## ----degpd-intercept-qq, fig.width = 5, fig.height = 5------------------------
set.seed(1)
r_d <- rqresid(fit_d)
qqnorm(r_d, main = "Q-Q Plot (DEGPD-1)", pch = 20, col = "grey60")
qqline(r_d, col = "red")


## ----degpd2-sim---------------------------------------------------------------
set.seed(22)
sigma_true <- 3
xi_true    <- 0.3
kappa1_true <- 1.5
kappa2_true <- 3.0
prob_true   <- 0.6

n <- 2000
y <- rdiscegpd(n, sigma = sigma_true, xi = xi_true, kappa = kappa1_true,
               delta = kappa2_true, prob = prob_true, type = 6)
df <- data.frame(y = y, x = rep(1, n))

cat("Range:", range(y), "\n")
cat("Mean:", mean(y), "\n")


## ----degpd2-fit---------------------------------------------------------------
fit_d2 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa1 = ~ 1, ldkappa = ~ 1,
                    logitp = ~ 1),
               data = df, family = "degpd", degpd.args = list(m = 2))
summary(fit_d2)


## ----degpd2-check-------------------------------------------------------------
pars_d2 <- predict(fit_d2, type = "response")[1, ]
truth_d2 <- c(scale = sigma_true, shape = xi_true, kappa1 = kappa1_true,
              kappa2 = kappa2_true, p = prob_true)
fitted_d2 <- round(unlist(pars_d2), 4)
cbind(true = truth_d2, fitted = fitted_d2)


## ----degpd2-qq, fig.width = 5, fig.height = 5---------------------------------
set.seed(1)
r_d2 <- rqresid(fit_d2)
qqnorm(r_d2, main = "Q-Q Plot (DEGPD-2)", pch = 20, col = "grey60")
qqline(r_d2, col = "red")


## ----degpd5-sim---------------------------------------------------------------
set.seed(55)
sigma_true <- 3
xi_true    <- 0.3
kappa_true <- 2.0

n <- 2000
y <- rdiscegpd(n, sigma = sigma_true, xi = xi_true, kappa = kappa_true,
               type = 2)
df <- data.frame(y = y, x = rep(1, n))

cat("Range:", range(y), "\n")
cat("Mean:", mean(y), "\n")


## ----degpd5-fit---------------------------------------------------------------
fit_d5 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
               data = df, family = "degpd", degpd.args = list(m = 5))
summary(fit_d5)


## ----degpd5-check-------------------------------------------------------------
pars_d5 <- predict(fit_d5, type = "response")[1, ]
truth_d5 <- c(scale = sigma_true, shape = xi_true, kappa = kappa_true)
cbind(true = truth_d5, fitted = round(unlist(pars_d5), 4))


## ----degpd5-qq, fig.width = 5, fig.height = 5---------------------------------
set.seed(1)
r_d5 <- rqresid(fit_d5)
qqnorm(r_d5, main = "Q-Q Plot (DEGPD-5)", pch = 20, col = "grey60")
qqline(r_d5, col = "red")


## ----degpd6-sim---------------------------------------------------------------
set.seed(66)
sigma_true <- 3
xi_true    <- 0.3
kappa_true <- 2.0

n <- 2000
y <- rdiscegpd(n, sigma = sigma_true, xi = xi_true, kappa = kappa_true,
               type = 3)
df <- data.frame(y = y, x = rep(1, n))

cat("Range:", range(y), "\n")
cat("Mean:", mean(y), "\n")


## ----degpd6-fit---------------------------------------------------------------
fit_d6 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
               data = df, family = "degpd", degpd.args = list(m = 6))
summary(fit_d6)


## ----degpd6-check-------------------------------------------------------------
pars_d6 <- predict(fit_d6, type = "response")[1, ]
truth_d6 <- c(scale = sigma_true, shape = xi_true, kappa = kappa_true)
cbind(true = truth_d6, fitted = round(unlist(pars_d6), 4))


## ----degpd6-qq, fig.width = 5, fig.height = 5---------------------------------
set.seed(1)
r_d6 <- rqresid(fit_d6)
qqnorm(r_d6, main = "Q-Q Plot (DEGPD-6)", pch = 20, col = "grey60")
qqline(r_d6, col = "red")


## ----degpd-smooth-sim---------------------------------------------------------
set.seed(7)
n <- 2000
x <- runif(n, 0, 1)
sigma_x <- exp(1 + 2 * x)
xi_true  <- 0.3
kappa_true <- 1.5

y <- rdiscegpd(n, sigma = sigma_x, xi = xi_true, kappa = kappa_true, type = 1)
df <- data.frame(y = y, x = x)


## ----degpd-smooth-fit---------------------------------------------------------
fit_ds <- egpd(list(lsigma = y ~ s(x, k = 10), lxi = ~ 1, lkappa = ~ 1),
               data = df, family = "degpd", degpd.args = list(m = 1))
summary(fit_ds)


## ----degpd-smooth-plot, fig.width = 7, fig.height = 5-------------------------
plot(fit_ds)


## ----degpd-smooth-compare, fig.width = 7, fig.height = 5----------------------
xgrid <- data.frame(x = seq(0, 1, length = 200))
pred_ds <- predict(fit_ds, newdata = xgrid, type = "response")

plot(xgrid$x, exp(1 + 2 * xgrid$x), type = "l", lwd = 2, col = "black",
     xlab = "x", ylab = expression(sigma(x)),
     main = "Recovered smooth scale function (DEGPD-1)")
lines(xgrid$x, pred_ds$scale, col = "steelblue", lwd = 2, lty = 2)
legend("topleft", legend = c("True", "Fitted"),
       col = c("black", "steelblue"), lwd = 2, lty = c(1, 2))


## ----zidegpd-sim--------------------------------------------------------------
set.seed(3)
sigma_true <- 3
xi_true    <- 0.3
kappa_true <- 1.5
pi_true    <- 0.3

n <- 2000
y <- rzidiscegpd(n, pi = pi_true, sigma = sigma_true, xi = xi_true,
                 kappa = kappa_true, type = 1)
df <- data.frame(y = y, x = rep(1, n))

cat("Proportion of zeros:", mean(y == 0), "\n")
cat("Expected proportion: pi + (1-pi)*P(Y=0) =",
    round(pi_true + (1 - pi_true) * ddiscegpd(0, sigma = sigma_true,
          xi = xi_true, kappa = kappa_true, type = 1), 3), "\n")


## ----zidegpd-fit--------------------------------------------------------------
fit_zi <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1, logitpi = ~ 1),
               data = df, family = "zidegpd", zidegpd.args = list(m = 1))
summary(fit_zi)


## ----zidegpd-check------------------------------------------------------------
pars_zi <- predict(fit_zi, type = "response")[1, ]
truth_zi <- c(scale = sigma_true, shape = xi_true, kappa = kappa_true,
              pi = pi_true)
cbind(true = truth_zi, fitted = round(unlist(pars_zi), 4))


## ----zidegpd-qq, fig.width = 5, fig.height = 5--------------------------------
set.seed(1)
r_zi <- rqresid(fit_zi)
qqnorm(r_zi, main = "Q-Q Plot (ZIDEGPD-1)", pch = 20, col = "grey60")
qqline(r_zi, col = "red")


## ----model-selection----------------------------------------------------------
fit_nozi <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
                 data = df, family = "degpd", degpd.args = list(m = 1))

data.frame(
  Model = c("DEGPD-1 (no ZI)", "ZIDEGPD-1"),
  logLik = c(logLik(fit_nozi), logLik(fit_zi)),
  AIC = c(AIC(fit_nozi), AIC(fit_zi))
)


## ----quantile-pred------------------------------------------------------------
probs <- c(0.5, 0.9, 0.95, 0.99)

# Theoretical quantiles from the known distribution
true_q <- qzidiscegpd(probs, pi = pi_true, sigma = sigma_true,
                       xi = xi_true, kappa = kappa_true, type = 1)

# Fitted quantiles
fit_q <- predict(fit_zi, type = "quantile", prob = probs)

data.frame(prob = probs,
           true = true_q,
           fitted = unlist(fit_q[1, ]),
           empirical = quantile(y, probs))

