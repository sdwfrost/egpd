# Comparing egpd and gamlss Fits

The `egpd` package provides gamlss family constructors
([`EGPD1()`](https://sdwfrost.github.io/egpd/reference/EGPD1.md),
[`DEGPD1()`](https://sdwfrost.github.io/egpd/reference/DEGPD1.md),
[`ZIEGPD1()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD1.md),
[`ZIDEGPD1()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD1.md),
etc.) that allow fitting EGPD models using the `gamlss` package’s
distributional regression framework. This vignette compares parameter
estimates from the two fitting approaches on simulated data.

## Parameter Mapping

The gamlss framework uses standard parameter names (mu, sigma, nu, tau).
These map to EGPD parameters as follows:

| gamlss param | Model 1 (EGPD)  | Model 3 (EGPD)  | Model 4 (EGPD)  | ZI Models    |
|:-------------|:----------------|:----------------|:----------------|:-------------|
| mu           | sigma (scale)   | sigma (scale)   | sigma (scale)   | sigma        |
| sigma        | xi (shape)      | xi (shape)      | xi (shape)      | xi           |
| nu           | kappa (G-param) | delta (G-param) | delta (G-param) | kappa/delta  |
| tau          | –               | –               | kappa (G-param) | pi (ZI prob) |

## Fitting notes

- **gamlss convergence**: The EGPD families use numerical derivatives,
  which can require more iterations than the gamlss default
  (`n.cyc = 20`). We recommend setting `n.cyc = 200` via
  [`gamlss.control()`](https://rdrr.io/pkg/gamlss/man/gamlss.control.html)
  to ensure convergence.
- **Penalised vs raw log-likelihood**:
  [`egpd()`](https://sdwfrost.github.io/egpd/reference/egpd.md) uses
  penalised likelihood estimation (as in `mgcv`), so
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) on an egpd object
  returns the penalised log-likelihood, which is slightly lower than the
  raw log-likelihood. The gamlss
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) returns the raw
  log-likelihood. Despite this difference in objective function, both
  approaches recover very similar parameter estimates.

## Discrete EGPD Model 1

We simulate from a discrete EGPD with G(u) = u^kappa and fit using both
[`egpd()`](https://sdwfrost.github.io/egpd/reference/egpd.md) (penalised
likelihood via GAM) and
[`gamlss()`](https://rdrr.io/pkg/gamlss/man/gamlss.html) (with the
[`DEGPD1()`](https://sdwfrost.github.io/egpd/reference/DEGPD1.md)
family).

``` r

library(egpd)
library(gamlss)
#> Loading required package: splines
#> Loading required package: gamlss.data
#> 
#> Attaching package: 'gamlss.data'
#> The following object is masked from 'package:datasets':
#> 
#>     sleep
#> Loading required package: gamlss.dist
#> Loading required package: nlme
#> Loading required package: parallel
#>  **********   GAMLSS Version 5.5-0  **********
#> For more on GAMLSS look at https://www.gamlss.com/
#> Type gamlssNews() to see new features/changes/bug fixes.

set.seed(99)
sigma_true <- 3
xi_true    <- 0.15
kappa_true <- 2
n <- 2000

y <- rdiscegpd(n, sigma = sigma_true, xi = xi_true, kappa = kappa_true, type = 1)
df <- data.frame(y = y)
cat("Range of y:", range(y), "\n")
#> Range of y: 0 147
```

### egpd fit

``` r

fit_egpd <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
                 data = df, family = "degpd", degpd.args = list(m = 1))
pars_egpd <- predict(fit_egpd, type = "response")[1, ]
```

### gamlss fit

``` r

fit_gamlss <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                     data = df, family = DEGPD1(),
                     control = gamlss.control(n.cyc = 200, trace = FALSE))
```

### Comparison

``` r

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
#>   parameter true   egpd gamlss
#> 1     sigma 3.00 2.7332 2.7205
#> 2        xi 0.15 0.1974 0.1989
#> 3     kappa 2.00 2.1228 2.1310
```

Both approaches recover similar parameter estimates for the discrete
model.

## Zero-Inflated Discrete EGPD Model 1

The
[`ZIDEGPD1()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD1.md)
family supports zero-inflated discrete EGPD models via gamlss. Here we
compare [`egpd()`](https://sdwfrost.github.io/egpd/reference/egpd.md)
and [`gamlss()`](https://rdrr.io/pkg/gamlss/man/gamlss.html) on
simulated data with excess zeros.

``` r

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
#> Proportion of zeros: 0.292
```

### egpd fit

``` r

fit_zidegpd <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1, logitpi = ~ 1),
                    data = df, family = "zidegpd", zidegpd.args = list(m = 1))
pars_zidegpd <- predict(fit_zidegpd, type = "response")[1, ]
```

### gamlss fit

``` r

fit_zi_gamlss <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                        tau.formula = ~ 1, data = df, family = ZIDEGPD1(),
                        control = gamlss.control(n.cyc = 200, trace = FALSE))
```

### Comparison

``` r

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
#>   parameter true   egpd gamlss
#> 1     sigma 3.00 3.0878 3.0460
#> 2        xi 0.15 0.1465 0.1505
#> 3     kappa 2.00 1.7892 1.8202
#> 4        pi 0.25 0.2161 0.2178
```

Both fitting approaches recover all four parameters for the
zero-inflated discrete model.

## Continuous EGPD Model 1

We also demonstrate the continuous EGPD family.

``` r

set.seed(42)
sigma_true <- 2
xi_true    <- 0.2
kappa_true <- 1.5
n <- 2000

y <- regpd(n, sigma = sigma_true, xi = xi_true, kappa = kappa_true, type = 1)
df <- data.frame(y = y, x = rep(1, n))
```

### egpd fit

``` r

fit_egpd_c <- egpd(list(lpsi = y ~ 1, xi = ~ 1, lkappa = ~ 1),
                   data = df, family = "egpd", egpd.args = list(m = 1))
pars_egpd_c <- predict(fit_egpd_c, type = "response")[1, ]
```

### gamlss fit

``` r

fit_gamlss_c <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                       data = df, family = EGPD1(),
                       control = gamlss.control(n.cyc = 200, trace = FALSE))
```

### Comparison

``` r

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
#>   parameter true   egpd gamlss
#> 1     sigma  2.0 2.0957 2.0900
#> 2        xi  0.2 0.1933 0.1943
#> 3     kappa  1.5 1.3904 1.3928
```

Both approaches recover very similar estimates for the continuous model
as well.
