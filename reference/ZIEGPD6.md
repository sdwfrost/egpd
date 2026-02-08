# gamlss Family for Zero-Inflated Continuous EGPD Model 6

The `ZIEGPD6()` function defines a `gamlss.family` distribution for
fitting a zero-inflated continuous Extended Generalized Pareto
Distribution with model 6 (truncated beta G-transformation).

The functions `dZIEGPD6`, `pZIEGPD6`, `qZIEGPD6`, and `rZIEGPD6` define
the density, distribution function, quantile function, and random
generation.

## Usage

``` r
ZIEGPD6(mu.link = "log", sigma.link = "log", nu.link = "log",
        tau.link = "logit")

dZIEGPD6(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE)
pZIEGPD6(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
         log.p = FALSE)
qZIEGPD6(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
         log.p = FALSE)
rZIEGPD6(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1)
```

## Arguments

- mu.link:

  link function for mu, default `"log"`

- sigma.link:

  link function for sigma, default `"log"`

- nu.link:

  link function for nu, default `"log"`

- tau.link:

  link function for tau, default `"logit"`

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

- tau:

  zero-inflation probability (pi), in (0,1)

- log, log.p:

  logical; if TRUE, probabilities/densities given as log

- lower.tail:

  logical; if TRUE (default), probabilities are P(X \<= x)

## Value

`ZIEGPD6()` returns an object of class `"gamlss.family"`.

`dZIEGPD6` gives the density, `pZIEGPD6` gives the distribution
function, `qZIEGPD6` gives the quantile function, and `rZIEGPD6`
generates random deviates.

## References

Gamet, P. and Jalbert, J. (2022). A flexible extended generalized Pareto
distribution for tail estimation. *Environmetrics*, 33(6), e2744.

Abbas, A., Ahmad, T. and Ahmad, I. (2025). Modeling zero-inflated
precipitation extremes. *arXiv preprint* arXiv:2504.11058.
<https://arxiv.org/abs/2504.11058>

## See also

[`ZIEGPD1`](https://sdwfrost.github.io/egpd/reference/ZIEGPD1.md),
[`ZIEGPD5`](https://sdwfrost.github.io/egpd/reference/ZIEGPD5.md),
[`ziegpd_bamlss`](https://sdwfrost.github.io/egpd/reference/ziegpd_bamlss.md),
[`dziegpd`](https://sdwfrost.github.io/egpd/reference/dziegpd.md)
