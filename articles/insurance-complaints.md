Discrete EGPD Models for Insurance Complaint Counts
================

This vignette demonstrates fitting Discrete Extended Generalized Pareto
Distribution (DEGPD) models to insurance complaint count data using the
`egpd` package. The data record the number of upheld complaints against
automobile insurance companies in New York State.

## Data

``` r
library(egpd)
data(ny_complaints)
str(ny_complaints)
```

    'data.frame':   1942 obs. of  4 variables:
     $ year    : int  2020 2020 2020 2020 2020 2020 2020 2020 2020 2020 ...
     $ upheld  : int  0 0 0 0 0 0 0 0 0 0 ...
     $ total   : int  20 2 13 5 0 4 5 5 3 5 ...
     $ premiums: num  107.1 65 63.5 53.2 43.8 ...

The response variable is the number of upheld complaints per insurer per
year. These are non-negative integer counts with a heavy right tail.

``` r
y <- ny_complaints$upheld
plot(table(y[y <= 30]), main = "Upheld insurance complaints (NY)",
     xlab = "Number of upheld complaints", ylab = "Frequency")
```

![](insurance-complaints_files/figure-gfm/eda-1.png)

``` r
cat("n =", length(y), " range:", range(y), "\n")
```

    n = 1942  range: 0 265 

## Fitting DEGPD models

The `egpd()` function fits DEGPD models within a GAM framework. For
intercept-only models (no covariates), each parameter gets a `~ 1`
formula. The DEGPD family uses a log link for the scale (`sigma`) and
shape (`xi`) parameters, so the fitted coefficients are on the log
scale.

### DEGPD Model 1: Power transformation

Model 1 uses G(u) = u^kappa with three parameters: sigma, xi, kappa.

``` r
df <- data.frame(y = y, x = rep(1, length(y)))
fit1 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
             data = df, family = "degpd", degpd.args = list(m = 1))
summary(fit1)
```


    ** Parametric terms **

    logscale
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)    -0.23       0.27   -0.85    0.397

    logshape
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)    -0.32       0.08   -4.06 4.89e-05

    logkappa
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)     0.35       0.22    1.57    0.116

    ** Smooth terms **

``` r
cat("Log-likelihood:", logLik(fit1), "\n")
```

    Log-likelihood: -3642.467 

``` r
cat("AIC:", AIC(fit1), "\n")
```

    AIC: 7290.933 

### DEGPD Model 2: Mixture of power transformations

Model 2 uses G(u) = p·u^kappa1 + (1-p)·u^kappa2 with five parameters:
sigma, xi, kappa1, dkappa, p.

``` r
fit2 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa1 = ~ 1, ldkappa = ~ 1,
                  logitp = ~ 1),
             data = df, family = "degpd", degpd.args = list(m = 2))
summary(fit2)
```


    ** Parametric terms **

    logscale
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)     0.29       0.16    1.82   0.0688

    logshape
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)    -0.78       0.15   -5.36 8.36e-08

    logkappa1
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)    -0.03       0.13   -0.22    0.829

    logdkappa
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)     5.64       0.39   14.63   <2e-16

    logitp
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)     4.03       0.29   13.82   <2e-16

    ** Smooth terms **

``` r
cat("Log-likelihood:", logLik(fit2), "\n")
```

    Log-likelihood: -3638.066 

``` r
cat("AIC:", AIC(fit2), "\n")
```

    AIC: 7286.132 

The two-power mixture gains 4.40 log-likelihood units over model 1. The
second power has weight `1 - p`, not `p`; a small fitted weight does not
by itself mean the components coincide.

### DEGPD Model 3: Incomplete beta transformation

Model 3 uses an incomplete beta transformation with three parameters:
sigma, xi, delta.

``` r
fit3 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, ldelta = ~ 1),
             data = df, family = "degpd", degpd.args = list(m = 3))
summary(fit3)
```


    ** Parametric terms **

    logscale
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)    -0.02       0.15   -0.15    0.878

    logshape
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)    -0.34       0.08   -4.29 1.79e-05

    logdelta
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)     2.02       0.73    2.75  0.00589

    ** Smooth terms **

``` r
cat("Log-likelihood:", logLik(fit3), "\n")
```

    Log-likelihood: -3642.691 

``` r
cat("AIC:", AIC(fit3), "\n")
```

    AIC: 7291.382 

### DEGPD Model 4: Power-beta transformation

Model 4 combines the incomplete beta and power transformations with four
parameters: sigma, xi, delta, kappa.

``` r
fit4 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, ldelta = ~ 1, lkappa = ~ 1),
             data = df, family = "degpd", degpd.args = list(m = 4))
summary(fit4)
```

    Warning in sqrt(obj$Vp[cbind(id2, id2)]): NaNs produced


    ** Parametric terms **

    logscale
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)    -0.23       0.27   -0.85    0.397

    logshape
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)    -0.32       0.08   -4.06 4.89e-05

    logdelta
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)   100.69        NaN     NaN      NaN

    logkappa
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)     1.04       0.22    4.73  2.3e-06

    ** Smooth terms **

    Note: the Hessian behind the standard errors is singular (numerical rank 3 of 4 estimated coefficients; `hessian.rank`, see ?egpd): they do not reflect the uncertainty along the unidentified direction(s).

``` r
cat("Log-likelihood:", logLik(fit4), "\n")
```

    Log-likelihood: -3642.467 

``` r
cat("AIC:", AIC(fit4), "\n")
```

    AIC: 7290.933 

Here the summary notes a singular Hessian: `logdelta` runs to a plateau
of its link (delta grows without bound, see `fit4$plateau`), where model
4 reduces to model 1, so delta is not identified; the log-likelihood and
AIC equal model 1’s.

### DEGPD Model 5: Truncated normal transformation

Model 5 uses a truncated normal G-transformation with three parameters:
sigma, xi, kappa.

``` r
fit5 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
             data = df, family = "degpd", degpd.args = list(m = 5))
summary(fit5)
```


    ** Parametric terms **

    logscale
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)     -0.2       0.21   -0.94    0.348

    logshape
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)    -0.31       0.08   -3.96 7.55e-05

    logkappa
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)     0.65       0.59    1.11    0.266

    ** Smooth terms **

``` r
cat("Log-likelihood:", logLik(fit5), "\n")
```

    Log-likelihood: -3642.309 

``` r
cat("AIC:", AIC(fit5), "\n")
```

    AIC: 7290.617 

### DEGPD Model 6: Truncated beta transformation

Model 6 uses a truncated beta G-transformation with three parameters:
sigma, xi, kappa.

``` r
fit6 <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
             data = df, family = "degpd", degpd.args = list(m = 6))
summary(fit6)
```


    ** Parametric terms **

    logscale
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)    -0.13       0.19   -0.67    0.502

    logshape
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)    -0.32       0.08   -4.07 4.68e-05

    logkappa
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)     0.45       0.24    1.84   0.0651

    ** Smooth terms **

``` r
cat("Log-likelihood:", logLik(fit6), "\n")
```

    Log-likelihood: -3642.442 

``` r
cat("AIC:", AIC(fit6), "\n")
```

    AIC: 7290.884 

## Model comparison

``` r
aic_table <- data.frame(
  Model = c("DEGPD-1", "DEGPD-2", "DEGPD-3", "DEGPD-4", "DEGPD-5", "DEGPD-6"),
  npar = c(3, 5, 3, 4, 3, 3),
  logLik = c(logLik(fit1), logLik(fit2), logLik(fit3), logLik(fit4),
             logLik(fit5), logLik(fit6)),
  AIC = c(AIC(fit1), AIC(fit2), AIC(fit3), AIC(fit4),
          AIC(fit5), AIC(fit6))
)
aic_table
```

        Model npar    logLik      AIC
    1 DEGPD-1    3 -3642.467 7290.933
    2 DEGPD-2    5 -3638.066 7286.132
    3 DEGPD-3    3 -3642.691 7291.382
    4 DEGPD-4    4 -3642.467 7290.933
    5 DEGPD-5    3 -3642.309 7290.617
    6 DEGPD-6    3 -3642.442 7290.884

## Goodness of fit

We can assess model fit by comparing the empirical and fitted PMFs.

``` r
# Extract fitted parameters on the response scale
pars1 <- predict(fit1, type = "response")
sigma1 <- pars1$scale[1]; xi1 <- pars1$shape[1]; kappa1 <- pars1$kappa[1]

xvals <- 0:30
emp_pmf <- tabulate(y + 1, nbins = max(xvals) + 1) / length(y)
emp_pmf <- emp_pmf[seq_along(xvals)]

fit_pmf1 <- ddiscegpd(xvals, sigma = sigma1, xi = xi1, kappa = kappa1, type = 1)

plot(xvals, emp_pmf, type = "h", lwd = 2, col = "grey60",
     main = "Empirical vs fitted PMF (DEGPD-1)",
     xlab = "Count", ylab = "Probability")
lines(xvals + 0.2, fit_pmf1, type = "h", lwd = 2, col = "steelblue")
legend("topright", legend = c("Empirical", "DEGPD-1"),
       col = c("grey60", "steelblue"), lwd = 2)
```

![](insurance-complaints_files/figure-gfm/gof-1.png)

## Q-Q plots

Randomized quantile residuals provide a model-agnostic diagnostic: if
the model is correct, the residuals should follow a standard normal
distribution. The `rqresid()` function computes these by drawing a
uniform variate between the lower and upper CDF bounds at each
observation and transforming to the normal scale.

``` r
set.seed(1)
par(mfrow = c(3, 2))

r1 <- rqresid(fit1)
qqnorm(r1, main = "Q-Q Plot (DEGPD-1)", pch = 20, col = "grey60")
qqline(r1, col = "red")

r2 <- rqresid(fit2)
qqnorm(r2, main = "Q-Q Plot (DEGPD-2)", pch = 20, col = "grey60")
qqline(r2, col = "red")

r3 <- rqresid(fit3)
qqnorm(r3, main = "Q-Q Plot (DEGPD-3)", pch = 20, col = "grey60")
qqline(r3, col = "red")

r4 <- rqresid(fit4)
qqnorm(r4, main = "Q-Q Plot (DEGPD-4)", pch = 20, col = "grey60")
qqline(r4, col = "red")

r5 <- rqresid(fit5)
qqnorm(r5, main = "Q-Q Plot (DEGPD-5)", pch = 20, col = "grey60")
qqline(r5, col = "red")

r6 <- rqresid(fit6)
qqnorm(r6, main = "Q-Q Plot (DEGPD-6)", pch = 20, col = "grey60")
qqline(r6, col = "red")
```

![](insurance-complaints_files/figure-gfm/qq-1.png)

``` r
par(mfrow = c(1, 1))
```

## Quantile predictions

The `predict` method with `type = "quantile"` computes quantile
estimates from the fitted model.

``` r
probs <- c(0.5, 0.9, 0.95, 0.99)
qpred <- predict(fit1, type = "quantile", prob = probs)
qpred[1, ]
```

      q:0.5 q:0.9 q:0.95 q:0.99
    1     1     6     11     38

Compare with empirical quantiles:

``` r
emp_q <- quantile(y, probs)
cbind(prob = probs, empirical = emp_q, fitted = unlist(qpred[1, ]))
```

        prob empirical fitted
    50% 0.50      1.00      1
    90% 0.90      5.00      6
    95% 0.95     10.00     11
    99% 0.99     47.18     38
