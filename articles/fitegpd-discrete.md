# Fitting Discrete Distributions with fitegpd

The [`fitegpd()`](https://sdwfrost.github.io/egpd/reference/fitegpd.md)
function supports two discrete distribution families for fitting
non-negative integer (count) data:

1.  **Discrete EGPD** (`family = "degpd"`) — for counts without excess
    zeros
2.  **Zero-Inflated Discrete EGPD** (`family = "zidegpd"`) — for counts
    with more zeros than the base model can accommodate

Both discretize the continuous EGPD by placing the probability mass of
each integer $`k`$ at
$`P(X = k) = F_{\mathrm{EGPD}}(k+1) - F_{\mathrm{EGPD}}(k)`$, inheriting
the flexible body and Pareto tail of the continuous model. This makes
them natural alternatives to the Poisson, negative binomial, or
Conway–Maxwell–Poisson for heavy-tailed or over-dispersed count data.

## 1. Discrete EGPD

### The model

The discrete EGPD (DEGPD) has PMF

``` math
P(X = k) = G\!\bigl(H(k+1)\bigr) - G\!\bigl(H(k)\bigr), \qquad k = 0, 1, 2, \ldots
```

where $`H`$ is the standard GPD CDF and $`G`$ is the EGPD transformation
function. As with the continuous EGPD, the `type` argument selects the
parametric form of $`G`$.

### Simulating and fitting

We generate data from a DEGPD Type 1 model and fit it back.

``` r

library(egpd)
set.seed(1)

sigma_true <- 3
xi_true    <- 0.1
kappa_true <- 1.5

x <- rdiscegpd(1000, sigma = sigma_true, xi = xi_true,
               kappa = kappa_true, type = 1)

cat("Range:", range(x), "\n")
#> Range: 0 51
cat("Mean:", round(mean(x), 2), "  Var:", round(var(x), 2), "\n")
#> Mean: 3.89   Var: 19
cat("Proportion zeros:", mean(x == 0), "\n")
#> Proportion zeros: 0.145
```

``` r

barplot(table(x) / length(x), main = "Simulated DEGPD data",
        xlab = "x", ylab = "Proportion", col = "lightblue", border = "grey")
```

![](fitegpd-discrete_files/figure-html/degpd-hist-1.png)

``` r

fit_degpd <- fitegpd(x, type = 1, family = "degpd")
summary(fit_degpd)
#> Fitting of the distribution 'degpd' (type 1)
#> Method: mle
#> 
#> Estimated parameters:
#>       Estimate Std. Error z value Pr(>|z|)    
#> sigma  2.79011    0.27440  10.168  < 2e-16 ***
#> xi     0.14362    0.04309   3.333 0.000859 ***
#> kappa  1.59145    0.13492  11.796  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Convergence:  successful 
#> Loglikelihood:  -2457.38   AIC:  4920.77   BIC:  4935.49 
#> Number of observations:  1000
```

### Parameter recovery

``` r

truth <- c(sigma = sigma_true, xi = xi_true, kappa = kappa_true)
est <- fit_degpd$estimate
cbind(true = truth, estimate = round(est, 4),
      SE = round(fit_degpd$sd, 4))
#>       true estimate     SE
#> sigma  3.0   2.7901 0.2744
#> xi     0.1   0.1436 0.0431
#> kappa  1.5   1.5915 0.1349
```

### Diagnostics

The diagnostic plot for discrete families shows the empirical PMF with
fitted probabilities (panel 1), a step-function CDF comparison (panel
2), a normal Q-Q plot of randomized quantile residuals (panel 3), and a
P-P plot using the randomized probability integral transform (panel 4).
The randomized PIT (Dunn & Smyth, 1996) spreads tied discrete values
into continuous uniform variates, producing clean diagonal plots instead
of the banded patterns that arise from integer-valued quantiles.

``` r

plot(fit_degpd)
```

![](fitegpd-discrete_files/figure-html/degpd-plot-1.png)

### Comparing types by AIC

Different G-transformation types suit different data shapes. We can use
AIC to select among them:

``` r

aic_table <- data.frame(
  type = c(1, 4, 5),
  AIC = c(
    AIC(fitegpd(x, type = 1, family = "degpd")),
    AIC(fitegpd(x, type = 4, family = "degpd")),
    AIC(fitegpd(x, type = 5, family = "degpd"))
  )
)
aic_table
#>   type      AIC
#> 1    1 4920.769
#> 2    4 4921.979
#> 3    5 4922.722
```

### Confidence intervals

``` r

confint(fit_degpd)
#>            2.5 %    97.5 %
#> sigma 2.25229832 3.3279173
#> xi    0.05916484 0.2280755
#> kappa 1.32702000 1.8558842
```

### Fixing parameters

As with continuous fits, individual parameters can be fixed at known
values:

``` r

fit_fix <- fitegpd(x, type = 1, family = "degpd",
                    fix.arg = list(xi = 0.1))
summary(fit_fix)
#> Fitting of the distribution 'degpd' (type 1)
#> Method: mle
#> 
#> Estimated parameters:
#>       Estimate Std. Error z value Pr(>|z|)    
#> sigma  3.03816    0.14717   20.64   <2e-16 ***
#> kappa  1.49640    0.08492   17.62   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Fixed parameters:
#>    Value
#> xi   0.1
#> 
#> Convergence:  successful 
#> Loglikelihood:  -2457.94   AIC:  4919.87   BIC:  4929.69 
#> Number of observations:  1000
```

## 2. Zero-Inflated Discrete EGPD

### The model

Many count datasets have more zeros than any standard count model can
explain — e.g. the number of insurance claims (many policyholders file
none), species counts (many sites have zero individuals), or disease
counts (many days with no cases). The zero-inflated discrete EGPD
(ZIDEGPD) handles this by mixing a point mass at zero with the DEGPD:

``` math
P(X = k) = \begin{cases}
\pi + (1-\pi)\,P_{\mathrm{DEGPD}}(X=0) & k = 0 \\
(1-\pi)\,P_{\mathrm{DEGPD}}(X=k) & k \ge 1
\end{cases}
```

where $`\pi \in (0,1)`$ is the zero-inflation probability. The overall
zero probability is $`\pi + (1-\pi)\,P_{\mathrm{DEGPD}}(0)`$, which is
always larger than $`P_{\mathrm{DEGPD}}(0)`$ alone.

### Simulating and fitting

``` r

set.seed(42)

sigma_true <- 2
xi_true    <- 0.1
kappa_true <- 1.5
pi_true    <- 0.3

y <- rzidiscegpd(1000, pi = pi_true, sigma = sigma_true,
                  xi = xi_true, kappa = kappa_true, type = 1)

cat("Range:", range(y), "\n")
#> Range: 0 27
cat("Mean:", round(mean(y), 2), "  Var:", round(var(y), 2), "\n")
#> Mean: 1.72   Var: 7.42
cat("Proportion zeros:", mean(y == 0), "\n")
#> Proportion zeros: 0.477
```

The zero proportion should be notably higher than for the non-inflated
model because of the extra $`\pi = 0.3`$ point mass.

``` r

barplot(table(y) / length(y), main = "Simulated ZIDEGPD data",
        xlab = "y", ylab = "Proportion", col = "lightblue", border = "grey")
```

![](fitegpd-discrete_files/figure-html/zidegpd-hist-1.png)

``` r

fit_zi <- fitegpd(y, type = 1, family = "zidegpd")
summary(fit_zi)
#> Fitting of the distribution 'zidegpd' (type 1)
#> Method: mle
#> 
#> Estimated parameters:
#>       Estimate Std. Error z value Pr(>|z|)    
#> sigma  1.75201    0.45741   3.830 0.000128 ***
#> xi     0.15800    0.07385   2.140 0.032389 *  
#> kappa  1.97969    0.84112   2.354 0.018590 *  
#> pi     0.36187    0.05998   6.034  1.6e-09 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Convergence:  successful 
#> Loglikelihood:  -1745.38   AIC:  3498.76   BIC:  3518.39 
#> Number of observations:  1000
```

### Parameter recovery

``` r

truth <- c(sigma = sigma_true, xi = xi_true,
           kappa = kappa_true, pi = pi_true)
est <- fit_zi$estimate
cbind(true = truth, estimate = round(est, 4),
      SE = round(fit_zi$sd, 4))
#>       true estimate     SE
#> sigma  2.0   1.7520 0.4574
#> xi     0.1   0.1580 0.0738
#> kappa  1.5   1.9797 0.8411
#> pi     0.3   0.3619 0.0600
```

### Diagnostics

``` r

plot(fit_zi)
```

![](fitegpd-discrete_files/figure-html/zidegpd-plot-1.png)

### Confidence intervals

``` r

confint(fit_zi)
#>            2.5 %    97.5 %
#> sigma 0.85549768 2.6485174
#> xi    0.01326359 0.3027374
#> kappa 0.33113076 3.6282584
#> pi    0.24431952 0.4794218
```

## 3. Comparing DEGPD and ZIDEGPD

When faced with zero-heavy count data, a natural question is whether the
zero inflation is needed or whether the base DEGPD already accounts for
the zeros adequately. We can compare the two models by AIC.

### Fitting both models to zero-inflated data

``` r

fit_no_zi <- fitegpd(y, type = 1, family = "degpd")
fit_with_zi <- fit_zi  # already fitted above

cat("DEGPD AIC:   ", AIC(fit_no_zi), "\n")
#> DEGPD AIC:    3503.567
cat("ZIDEGPD AIC: ", AIC(fit_with_zi), "\n")
#> ZIDEGPD AIC:  3498.756
```

Since the data were generated with $`\pi = 0.3`$, the zero-inflated
model should have a substantially lower AIC.

### Fitting both models to non-inflated data

Conversely, when the data have no excess zeros, the extra $`\pi`$
parameter should not help:

``` r

fit_degpd_on_x <- fit_degpd  # fitted to non-inflated x above
fit_zi_on_x <- fitegpd(x, type = 1, family = "zidegpd")

cat("DEGPD AIC:   ", AIC(fit_degpd_on_x), "\n")
#> DEGPD AIC:    4920.769
cat("ZIDEGPD AIC: ", AIC(fit_zi_on_x), "\n")
#> ZIDEGPD AIC:  4922.401
cat("Estimated pi:", round(fit_zi_on_x$estimate["pi"], 4), "\n")
#> Estimated pi: 0.0278
```

The estimated $`\pi`$ should be near zero and the ZIDEGPD AIC slightly
higher (penalised for the extra parameter).

## 4. Higher-order types

Types 5 and 6 provide additional flexibility. Type 5 has both $`\delta`$
and $`\kappa`$ parameters; Type 6 adds a mixing probability $`p`$.

``` r

set.seed(99)
x5 <- rdiscegpd(1000, sigma = 2, xi = 0.1, delta = 1.2,
                 kappa = 1.8, type = 5)
fit5 <- fitegpd(x5, type = 5, family = "degpd")
summary(fit5)
#> Fitting of the distribution 'degpd' (type 5)
#> Method: mle
#> 
#> Estimated parameters:
#>       Estimate Std. Error z value Pr(>|z|)    
#> sigma  1.63595    0.60546   2.702  0.00689 ** 
#> xi     0.15591    0.03919   3.978 6.94e-05 ***
#> delta  0.41132    0.97364   0.422  0.67269    
#> kappa  1.65941    0.15856  10.466  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Convergence:  successful 
#> Loglikelihood:  -2117.78   AIC:  4243.56   BIC:  4263.19 
#> Number of observations:  1000
```

``` r

plot(fit5)
```

![](fitegpd-discrete_files/figure-html/type5-plot-1.png)

## 5. Compound Poisson-Discrete EGPD

### The model

The Compound Poisson-Discrete EGPD (`family = "cpdegpd"`) models the
aggregate sum $`S = X_1 + \cdots + X_N`$ where

- $`N \sim \mathrm{Poisson}(\lambda)`$ is the random number of events,
  and
- $`X_i \sim \mathrm{Discrete\text{-}EGPD}(\sigma, \xi, \kappa, \ldots)`$
  are i.i.d. non-negative integer severities.

Since the individual claims are already integer-valued, the Panjer
recursion computes an **exact** compound distribution on
$`\{0, 1, 2, \ldots\}`$ — no discretization step is needed and there is
no bin-width parameter `h` to tune (unlike `family = "cpegpd"`).

### Simulating and fitting

``` r

set.seed(42)

sigma_true  <- 3
xi_true     <- 0.1
kappa_true  <- 1.5
lambda_true <- 2

z <- rcpdegpd(500, sigma = sigma_true, xi = xi_true,
               kappa = kappa_true, lambda = lambda_true, type = 1)

cat("Range:", range(z), "\n")
#> Range: 0 43
cat("Mean:", round(mean(z), 2), "  Var:", round(var(z), 2), "\n")
#> Mean: 7.03   Var: 61.4
cat("Proportion zeros:", mean(z == 0), "\n")
#> Proportion zeros: 0.214
```

``` r

barplot(table(z) / length(z), main = "Simulated CPDEGPD data",
        xlab = "z", ylab = "Proportion", col = "lightblue", border = "grey")
```

![](fitegpd-discrete_files/figure-html/cpdegpd-hist-1.png)

``` r

fit_cpdegpd <- fitegpd(z, type = 1, family = "cpdegpd")
summary(fit_cpdegpd)
#> Fitting of the distribution 'cpdegpd' (type 1)
#> Method: mle
#> 
#> Estimated parameters:
#>        Estimate Std. Error z value Pr(>|z|)    
#> sigma    2.7440     1.1024   2.489   0.0128 *  
#> xi       0.1482     0.1137   1.303   0.1924    
#> kappa    1.6625     1.0139   1.640   0.1011    
#> lambda   1.7767     0.2376   7.477  7.6e-14 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Convergence:  successful 
#> Loglikelihood:  -1488.58   AIC:  2985.17   BIC:  3002.02 
#> Number of observations:  500
```

### Parameter recovery

``` r

truth <- c(sigma = sigma_true, xi = xi_true,
           kappa = kappa_true, lambda = lambda_true)
est <- fit_cpdegpd$estimate
cbind(true = truth, estimate = round(est, 4),
      SE = round(fit_cpdegpd$sd, 4))
#>        true estimate     SE
#> sigma   3.0   2.7440 1.1024
#> xi      0.1   0.1482 0.1137
#> kappa   1.5   1.6625 1.0139
#> lambda  2.0   1.7767 0.2376
```

### Diagnostics

The diagnostic plots use discrete-style panels: a barplot with fitted
PMF points, a step-function CDF, and randomized PIT-based Q-Q and P-P
plots.

``` r

plot(fit_cpdegpd)
```

![](fitegpd-discrete_files/figure-html/cpdegpd-plot-1.png)

### Comparing with plain DEGPD on the same data

Since `cpdegpd` is a compound sum model while `degpd` is a single-event
model, we can compare their AIC on the same data to see which fits
better:

``` r

fit_degpd_on_z <- fitegpd(z, type = 1, family = "degpd")

cat("DEGPD AIC:   ", AIC(fit_degpd_on_z), "\n")
#> DEGPD AIC:    2995.597
cat("CPDEGPD AIC: ", AIC(fit_cpdegpd), "\n")
#> CPDEGPD AIC:  2985.166
```

The compound model should fit noticeably better when the data-generating
process truly involves random summation.

## 6. Practical guidelines

**Choosing between DEGPD and ZIDEGPD.** Start with `family = "degpd"`.
If the diagnostic plots show the model underestimates the zero
probability, refit with `family = "zidegpd"` and compare by AIC.

**Choosing the type.** Fit several types and compare AIC. Types 1 and 4
have three parameters and are good defaults. Type 5 (four parameters)
and Type 6 (five parameters) offer more flexibility but need larger
samples.

**Starting values.** The automatic starting values work well in most
cases. For difficult datasets (e.g. very heavy tails), supply custom
start values via `start = list(sigma = ..., xi = ...)`.

**Fixed parameters.** Use `fix.arg` when external information is
available — e.g. `fix.arg = list(xi = 0)` for an exponential tail.

## Summary

| Family | Parameters | Zero mechanism | Use case |
|----|----|----|----|
| `"degpd"` | $`\sigma, \xi`$ + type-specific | Natural zeros from DEGPD | Standard count data |
| `"zidegpd"` | $`\sigma, \xi`$ + type-specific + $`\pi`$ | Extra point mass at zero | Counts with excess zeros |
| `"cpdegpd"` | $`\sigma, \xi`$ + type-specific + $`\lambda`$ | Poisson zero ($`N=0`$) | Aggregate integer sums |

All families:

- Return S3 objects of class `"fitegpd"` with the same interface
  (`summary`, `plot`, `AIC`, `confint`, `coef`, `vcov`, `logLik`)
- Support `fix.arg` for fixing parameters at known values
- Estimate standard errors via the delta method on the Hessian
- Provide four-panel diagnostic plots tailored for discrete data (bar
  plots for PMF, step functions for CDF, randomized PIT-based Q-Q and
  P-P plots)
