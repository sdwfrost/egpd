## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)


## ----data---------------------------------------------------------------------
library(egpd)
library(evgam)
data(FCtmax)

FCtmax$month <- as.integer(format(FCtmax$date, "%m"))
FCtmax$year  <- as.integer(format(FCtmax$date, "%Y"))

summer <- FCtmax[FCtmax$month %in% 6:8, ]
thresh <- quantile(summer$tmax, 0.9)
cat("90th percentile threshold:", thresh, "°C\n")

exc_idx <- summer$tmax > thresh
df <- data.frame(
  y    = summer$tmax[exc_idx] - thresh,
  year = summer$year[exc_idx]
)
cat("Number of exceedances:", nrow(df), "out of", nrow(summer), "summer days\n")


## ----eda, fig.width = 7, fig.height = 4---------------------------------------
hist(df$y, breaks = 30, col = "steelblue", border = "white",
     main = "Summer temperature exceedances (Fort Collins, CO)",
     xlab = "Exceedance above threshold (°C)", ylab = "Frequency")


## ----egpd1--------------------------------------------------------------------
fit1 <- egpd(list(lpsi = y ~ 1, xi = ~ 1, lkappa = ~ 1),
             data = df, family = "egpd", egpd.args = list(m = 1))
summary(fit1)


## ----egpd1-pars---------------------------------------------------------------
pars1 <- predict(fit1, type = "response")[1, ]
cat("sigma =", round(pars1$scale, 3),
    " xi =", round(pars1$shape, 3),
    " kappa =", round(pars1$kappa, 3), "\n")


## ----egpd3--------------------------------------------------------------------
fit3 <- egpd(list(lpsi = y ~ 1, xi = ~ 1, ldelta = ~ 1),
             data = df, family = "egpd", egpd.args = list(m = 3))
summary(fit3)


## ----egpd4--------------------------------------------------------------------
fit4 <- egpd(list(lpsi = y ~ 1, xi = ~ 1, ldelta = ~ 1, lkappa = ~ 1),
             data = df, family = "egpd", egpd.args = list(m = 4))
summary(fit4)


## ----egpd5--------------------------------------------------------------------
fit5 <- egpd(list(lpsi = y ~ 1, xi = ~ 1, lkappa = ~ 1),
             data = df, family = "egpd", egpd.args = list(m = 5))
summary(fit5)


## ----egpd6--------------------------------------------------------------------
fit6 <- egpd(list(lpsi = y ~ 1, xi = ~ 1, lkappa = ~ 1),
             data = df, family = "egpd", egpd.args = list(m = 6))
summary(fit6)


## ----compare------------------------------------------------------------------
aic_table <- data.frame(
  Model  = c("EGPD-1", "EGPD-3", "EGPD-4", "EGPD-5", "EGPD-6"),
  npar   = c(3, 3, 4, 3, 3),
  logLik = round(c(logLik(fit1), logLik(fit3), logLik(fit4),
                   logLik(fit5), logLik(fit6)), 2),
  AIC    = round(c(AIC(fit1), AIC(fit3), AIC(fit4),
                   AIC(fit5), AIC(fit6)), 2)
)
aic_table


## ----qq, fig.width = 7, fig.height = 10---------------------------------------
set.seed(1)
par(mfrow = c(3, 2))

r1 <- rqresid(fit1)
qqnorm(r1, main = "Q-Q Plot (EGPD-1)", pch = 20, col = "grey60")
qqline(r1, col = "red")

r3 <- rqresid(fit3)
qqnorm(r3, main = "Q-Q Plot (EGPD-3)", pch = 20, col = "grey60")
qqline(r3, col = "red")

r4 <- rqresid(fit4)
qqnorm(r4, main = "Q-Q Plot (EGPD-4)", pch = 20, col = "grey60")
qqline(r4, col = "red")

r5 <- rqresid(fit5)
qqnorm(r5, main = "Q-Q Plot (EGPD-5)", pch = 20, col = "grey60")
qqline(r5, col = "red")

r6 <- rqresid(fit6)
qqnorm(r6, main = "Q-Q Plot (EGPD-6)", pch = 20, col = "grey60")
qqline(r6, col = "red")

par(mfrow = c(1, 1))


## ----survivor, fig.width = 7, fig.height = 5----------------------------------
y_sorted <- sort(df$y)
n <- length(y_sorted)
emp_surv <- 1 - (1:n) / (n + 1)

p1 <- predict(fit1, type = "response")[1, ]
surv1 <- 1 - pegpd(y_sorted, sigma = p1$scale, xi = p1$shape,
                    kappa = p1$kappa, type = 1)

p4 <- predict(fit4, type = "response")[1, ]
surv4 <- 1 - pegpd(y_sorted, sigma = p4$scale, xi = p4$shape,
                    delta = p4$delta, kappa = p4$kappa, type = 5)

plot(y_sorted, emp_surv, log = "y", pch = 20, col = "grey50",
     xlab = "Exceedance (°C)", ylab = "Survival probability",
     main = "Empirical vs fitted survivor functions")
lines(y_sorted, surv1, col = "steelblue", lwd = 2)
lines(y_sorted, surv4, col = "firebrick", lwd = 2, lty = 2)
p5 <- predict(fit5, type = "response")[1, ]
surv5 <- 1 - pegpd(y_sorted, sigma = p5$scale, xi = p5$shape,
                    kappa = p5$kappa, type = 2)

p6 <- predict(fit6, type = "response")[1, ]
surv6 <- 1 - pegpd(y_sorted, sigma = p6$scale, xi = p6$shape,
                    kappa = p6$kappa, type = 3)

lines(y_sorted, surv5, col = "darkgreen", lwd = 2, lty = 3)
lines(y_sorted, surv6, col = "purple", lwd = 2, lty = 4)
legend("topright", legend = c("Empirical", "EGPD-1", "EGPD-4", "EGPD-5", "EGPD-6"),
       col = c("grey50", "steelblue", "firebrick", "darkgreen", "purple"),
       pch = c(20, NA, NA, NA, NA), lty = c(NA, 1, 2, 3, 4),
       lwd = c(NA, 2, 2, 2, 2))


## ----trend--------------------------------------------------------------------
fit1_yr <- egpd(list(lpsi = y ~ s(year, k = 5), xi = ~ 1, lkappa = ~ 1),
                data = df, family = "egpd", egpd.args = list(m = 1))
summary(fit1_yr)
cat("\nAIC (intercept-only):", round(AIC(fit1), 2),
    "\nAIC (year trend):    ", round(AIC(fit1_yr), 2), "\n")


## ----trend-plot, fig.width = 7, fig.height = 5--------------------------------
year_grid <- data.frame(year = 1970:2019)
pred_yr <- predict(fit1_yr, newdata = year_grid, type = "response")

plot(year_grid$year, pred_yr$scale, type = "l", lwd = 2, col = "steelblue",
     xlab = "Year", ylab = expression(hat(sigma)(year)),
     main = "Estimated scale parameter over time")


## ----quantile-----------------------------------------------------------------
probs <- c(0.5, 0.9, 0.95, 0.99)
qpred <- predict(fit1, type = "quantile", prob = probs)
emp_q <- quantile(df$y, probs)

data.frame(
  probability = probs,
  empirical   = round(as.numeric(emp_q), 3),
  fitted      = round(as.numeric(unlist(qpred[1, ])), 3),
  row.names   = NULL
)

