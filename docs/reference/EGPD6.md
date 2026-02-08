# gamlss Family for Continuous EGPD Model 6

The `EGPD6()` function defines a `gamlss.family` distribution for
fitting a continuous Extended Generalized Pareto Distribution with model
6 (truncated beta G-transformation).

The functions `dEGPD6`, `pEGPD6`, `qEGPD6`, and `rEGPD6` define the
density, distribution function, quantile function, and random
generation.

## Usage

``` r
EGPD6(mu.link = "log", sigma.link = "log", nu.link = "log")

dEGPD6(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE)
pEGPD6(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
qEGPD6(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
rEGPD6(n, mu = 1, sigma = 0.5, nu = 1)
```

## Arguments

- mu.link:

  link function for mu, default `"log"`

- sigma.link:

  link function for sigma, default `"log"`

- nu.link:

  link function for nu, default `"log"`

- x, q:

  vector of quantiles

- p:

  vector of probabilities

- n:

  number of observations

- mu:

  GPD scale parameter (sigma), positive

- sigma:

  GPD shape parameter (xi), positive

- nu:

  G-transformation parameter (kappa), positive

- log, log.p:

  logical; if TRUE, probabilities/densities given as log

- lower.tail:

  logical; if TRUE (default), probabilities are P(X \<= x)

## Value

`EGPD6()` returns an object of class `"gamlss.family"`.

`dEGPD6` gives the density, `pEGPD6` gives the distribution function,
`qEGPD6` gives the quantile function, and `rEGPD6` generates random
deviates.

## References

Gamet, P. and Jalbert, J. (2022). A flexible extended generalized Pareto
distribution for tail estimation. *Environmetrics*, 33(6), e2744.

## See also

[`EGPD1`](https://sdwfrost.github.io/egpd/reference/EGPD1.md),
[`EGPD5`](https://sdwfrost.github.io/egpd/reference/EGPD5.md),
[`egpd_bamlss`](https://sdwfrost.github.io/egpd/reference/egpd_bamlss.md),
[`pegpd`](https://sdwfrost.github.io/egpd/reference/pegpd.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  library(gamlss)
  set.seed(1)
  y <- rEGPD6(500, mu = 2, sigma = 0.2, nu = 1.5)
  fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                family = EGPD6())
} # }
```
