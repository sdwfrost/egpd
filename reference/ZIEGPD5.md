# gamlss Family for Zero-Inflated Continuous EGPD Model 5

The `ZIEGPD5()` function defines a `gamlss.family` distribution for
fitting a zero-inflated continuous Extended Generalized Pareto
Distribution with model 5 (truncated normal G-transformation).

The functions `dZIEGPD5`, `pZIEGPD5`, `qZIEGPD5`, and `rZIEGPD5` define
the density, distribution function, quantile function, and random
generation.

## Usage

``` r
ZIEGPD5(mu.link = "log", sigma.link = "log", nu.link = "log",
        tau.link = "logit")

dZIEGPD5(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE)
pZIEGPD5(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
         log.p = FALSE)
qZIEGPD5(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
         log.p = FALSE)
rZIEGPD5(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1)
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

`ZIEGPD5()` returns an object of class `"gamlss.family"`.

`dZIEGPD5` gives the density, `pZIEGPD5` gives the distribution
function, `qZIEGPD5` gives the quantile function, and `rZIEGPD5`
generates random deviates.

## References

Gamet, P. and Jalbert, J. (2022). A flexible extended generalized Pareto
distribution for tail estimation. *Environmetrics*, 33(6), e2744.

Abbas, A., Ahmad, T. and Ahmad, I. (2025). Modeling zero-inflated
precipitation extremes. *arXiv preprint* arXiv:2504.11058.
<https://arxiv.org/abs/2504.11058>

## See also

[`ZIEGPD1`](https://sdwfrost.github.io/egpd/reference/ZIEGPD1.md),
[`ZIEGPD6`](https://sdwfrost.github.io/egpd/reference/ZIEGPD6.md),
[`ziegpd_bamlss`](https://sdwfrost.github.io/egpd/reference/ziegpd_bamlss.md),
[`dziegpd`](https://sdwfrost.github.io/egpd/reference/dziegpd.md)
