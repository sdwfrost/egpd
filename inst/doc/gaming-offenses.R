## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)


## ----data---------------------------------------------------------------------
library(egpd)
data(nsw_offenses)
str(nsw_offenses)


## ----eda----------------------------------------------------------------------
d <- nsw_offenses$offenses
plot(table(d), main = "Gaming and betting offenses (NSW)",
     xlab = "Number of offenses", ylab = "Frequency")
cat("n =", length(d), " range:", range(d), "\n")


## ----thresh10-----------------------------------------------------------------
u10 <- floor(quantile(d, 0.10))
cat("Threshold u (10th percentile):", u10, "\n")
y10 <- d[d >= u10] - u10
cat("Exceedances: n =", length(y10), "\n")

hist(y10, breaks = 50, col = "lightblue", border = "white",
     main = paste0("Threshold exceedances (u = ", u10, ")"),
     xlab = "Excess count", ylab = "Frequency")


## ----fit10--------------------------------------------------------------------
df10 <- data.frame(y = y10, x = rep(1, length(y10)))

fit10_m1 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
                 data = df10, family = "degpd", degpd.args = list(m = 1))

fit10_m2 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa1 = ~ 1, ldkappa = ~ 1,
                      logitp = ~ 1),
                 data = df10, family = "degpd", degpd.args = list(m = 2))

fit10_m3 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, ldelta = ~ 1),
                 data = df10, family = "degpd", degpd.args = list(m = 3))

fit10_m4 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, ldelta = ~ 1, lkappa = ~ 1),
                 data = df10, family = "degpd", degpd.args = list(m = 4))

fit10_m5 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
                 data = df10, family = "degpd", degpd.args = list(m = 5))

fit10_m6 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
                 data = df10, family = "degpd", degpd.args = list(m = 6))

aic10 <- data.frame(
  Model = c("DEGPD-1", "DEGPD-2", "DEGPD-3", "DEGPD-4", "DEGPD-5", "DEGPD-6"),
  npar = c(3, 5, 3, 4, 3, 3),
  logLik = c(logLik(fit10_m1), logLik(fit10_m2), logLik(fit10_m3),
             logLik(fit10_m4), logLik(fit10_m5), logLik(fit10_m6)),
  AIC = c(AIC(fit10_m1), AIC(fit10_m2), AIC(fit10_m3),
          AIC(fit10_m4), AIC(fit10_m5), AIC(fit10_m6))
)
aic10


## ----summary10----------------------------------------------------------------
summary(fit10_m1)


## ----thresh20-----------------------------------------------------------------
u20 <- floor(quantile(d, 0.20))
cat("Threshold u (20th percentile):", u20, "\n")
y20 <- d[d >= u20] - u20
cat("Exceedances: n =", length(y20), "\n")

hist(y20, breaks = 50, col = "lightblue", border = "white",
     main = paste0("Threshold exceedances (u = ", u20, ")"),
     xlab = "Excess count", ylab = "Frequency")


## ----fit20--------------------------------------------------------------------
df20 <- data.frame(y = y20, x = rep(1, length(y20)))

fit20_m1 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
                 data = df20, family = "degpd", degpd.args = list(m = 1))

fit20_m2 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa1 = ~ 1, ldkappa = ~ 1,
                      logitp = ~ 1),
                 data = df20, family = "degpd", degpd.args = list(m = 2))

fit20_m3 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, ldelta = ~ 1),
                 data = df20, family = "degpd", degpd.args = list(m = 3))

fit20_m4 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, ldelta = ~ 1, lkappa = ~ 1),
                 data = df20, family = "degpd", degpd.args = list(m = 4))

fit20_m5 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
                 data = df20, family = "degpd", degpd.args = list(m = 5))

fit20_m6 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
                 data = df20, family = "degpd", degpd.args = list(m = 6))

aic20 <- data.frame(
  Model = c("DEGPD-1", "DEGPD-2", "DEGPD-3", "DEGPD-4", "DEGPD-5", "DEGPD-6"),
  npar = c(3, 5, 3, 4, 3, 3),
  logLik = c(logLik(fit20_m1), logLik(fit20_m2), logLik(fit20_m3),
             logLik(fit20_m4), logLik(fit20_m5), logLik(fit20_m6)),
  AIC = c(AIC(fit20_m1), AIC(fit20_m2), AIC(fit20_m3),
          AIC(fit20_m4), AIC(fit20_m5), AIC(fit20_m6))
)
aic20


## ----gof, fig.width = 7, fig.height = 5---------------------------------------
# Fitted parameters for the 10th-percentile model 1
pars <- predict(fit10_m1, type = "response")
sigma <- pars$scale[1]; xi <- pars$shape[1]; kappa <- pars$kappa[1]

xvals <- 0:max(y10)
emp_pmf <- tabulate(y10 + 1, nbins = max(xvals) + 1) / length(y10)
fit_pmf <- ddiscegpd(xvals, sigma = sigma, xi = xi, kappa = kappa, type = 1)

plot(xvals[1:40], emp_pmf[1:40], type = "h", lwd = 2, col = "grey60",
     main = "Empirical vs fitted PMF (DEGPD-1, u = 10th pctl)",
     xlab = "Excess count", ylab = "Probability")
lines(xvals[1:40] + 0.2, fit_pmf[1:40], type = "h", lwd = 2, col = "steelblue")
legend("topright", legend = c("Empirical", "DEGPD-1"),
       col = c("grey60", "steelblue"), lwd = 2)


## ----qq10, fig.width = 7, fig.height = 10-------------------------------------
set.seed(1)
par(mfrow = c(3, 2))

r10_1 <- rqresid(fit10_m1)
qqnorm(r10_1, main = "Q-Q Plot (DEGPD-1, u = 10th pctl)", pch = 20, col = "grey60")
qqline(r10_1, col = "red")

r10_2 <- rqresid(fit10_m2)
qqnorm(r10_2, main = "Q-Q Plot (DEGPD-2, u = 10th pctl)", pch = 20, col = "grey60")
qqline(r10_2, col = "red")

r10_3 <- rqresid(fit10_m3)
qqnorm(r10_3, main = "Q-Q Plot (DEGPD-3, u = 10th pctl)", pch = 20, col = "grey60")
qqline(r10_3, col = "red")

r10_4 <- rqresid(fit10_m4)
qqnorm(r10_4, main = "Q-Q Plot (DEGPD-4, u = 10th pctl)", pch = 20, col = "grey60")
qqline(r10_4, col = "red")

r10_5 <- rqresid(fit10_m5)
qqnorm(r10_5, main = "Q-Q Plot (DEGPD-5, u = 10th pctl)", pch = 20, col = "grey60")
qqline(r10_5, col = "red")

r10_6 <- rqresid(fit10_m6)
qqnorm(r10_6, main = "Q-Q Plot (DEGPD-6, u = 10th pctl)", pch = 20, col = "grey60")
qqline(r10_6, col = "red")

par(mfrow = c(1, 1))


## ----qq20, fig.width = 7, fig.height = 10-------------------------------------
set.seed(1)
par(mfrow = c(3, 2))

r20_1 <- rqresid(fit20_m1)
qqnorm(r20_1, main = "Q-Q Plot (DEGPD-1, u = 20th pctl)", pch = 20, col = "grey60")
qqline(r20_1, col = "red")

r20_2 <- rqresid(fit20_m2)
qqnorm(r20_2, main = "Q-Q Plot (DEGPD-2, u = 20th pctl)", pch = 20, col = "grey60")
qqline(r20_2, col = "red")

r20_3 <- rqresid(fit20_m3)
qqnorm(r20_3, main = "Q-Q Plot (DEGPD-3, u = 20th pctl)", pch = 20, col = "grey60")
qqline(r20_3, col = "red")

r20_4 <- rqresid(fit20_m4)
qqnorm(r20_4, main = "Q-Q Plot (DEGPD-4, u = 20th pctl)", pch = 20, col = "grey60")
qqline(r20_4, col = "red")

r20_5 <- rqresid(fit20_m5)
qqnorm(r20_5, main = "Q-Q Plot (DEGPD-5, u = 20th pctl)", pch = 20, col = "grey60")
qqline(r20_5, col = "red")

r20_6 <- rqresid(fit20_m6)
qqnorm(r20_6, main = "Q-Q Plot (DEGPD-6, u = 20th pctl)", pch = 20, col = "grey60")
qqline(r20_6, col = "red")

par(mfrow = c(1, 1))


## ----sensitivity--------------------------------------------------------------
p10 <- predict(fit10_m1, type = "response")[1, ]
p20 <- predict(fit20_m1, type = "response")[1, ]
rbind("u = 10th pctl" = unlist(p10), "u = 20th pctl" = unlist(p20))

