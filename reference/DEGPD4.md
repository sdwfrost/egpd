# gamlss Family for Discrete EGPD Model 4

The `DEGPD4()` function defines a `gamlss.family` distribution for
fitting a discrete Extended Generalized Pareto Distribution with model 4
(power-beta G-transformation). This is a 4-parameter family.

The functions `dDEGPD4`, `pDEGPD4`, `qDEGPD4`, and `rDEGPD4` define the
PMF, distribution function, quantile function, and random generation.

## Usage

``` r
DEGPD4(mu.link = "log", sigma.link = "log", nu.link = "log",
       tau.link = "log")

dDEGPD4(x, mu = 1, sigma = 0.5, nu = 1, tau = 1, log = FALSE)
pDEGPD4(q, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE,
        log.p = FALSE)
qDEGPD4(p, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE,
        log.p = FALSE)
rDEGPD4(n, mu = 1, sigma = 0.5, nu = 1, tau = 1)
```

## Arguments

- mu.link:

  link function for mu, default `"log"`

- sigma.link:

  link function for sigma, default `"log"`

- nu.link:

  link function for nu, default `"log"`

- tau.link:

  link function for tau, default `"log"`

- x, q:

  vector of quantiles (non-negative integers)

- p:

  vector of probabilities

- n:

  number of observations

- mu:

  GPD scale parameter (sigma), positive

- sigma:

  GPD shape parameter (xi), positive

- nu:

  G-transformation parameter (delta), positive

- tau:

  G-transformation parameter (kappa), positive

- log, log.p:

  logical; if TRUE, probabilities/densities given as log

- lower.tail:

  logical; if TRUE (default), probabilities are P(X \<= x)

## Value

`DEGPD4()` returns an object of class `"gamlss.family"`.

`dDEGPD4` gives the PMF, `pDEGPD4` gives the distribution function,
`qDEGPD4` gives the quantile function, and `rDEGPD4` generates random
deviates.

## References

Ahmad, T. and Arshad, I. A. (2024). New flexible versions of extended
generalized Pareto model for count data. *arXiv preprint*
arXiv:2409.18719. <https://arxiv.org/abs/2409.18719>

## See also

[`DEGPD1`](https://sdwfrost.github.io/egpd/reference/DEGPD1.md),
[`DEGPD3`](https://sdwfrost.github.io/egpd/reference/DEGPD3.md),
[`degpd_bamlss`](https://sdwfrost.github.io/egpd/reference/degpd_bamlss.md)
