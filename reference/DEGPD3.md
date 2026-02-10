# gamlss Family for Discrete EGPD Model 3

The `DEGPD3()` function defines a `gamlss.family` distribution for
fitting a discrete Extended Generalized Pareto Distribution with model 3
(incomplete beta G-transformation).

The functions `dDEGPD3`, `pDEGPD3`, `qDEGPD3`, and `rDEGPD3` define the
PMF, distribution function, quantile function, and random generation.

## Usage

``` r
DEGPD3(mu.link = "log", sigma.link = "log", nu.link = "log")

dDEGPD3(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE)
pDEGPD3(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
qDEGPD3(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
rDEGPD3(n, mu = 1, sigma = 0.5, nu = 1)
```

## Arguments

- mu.link:

  link function for mu, default `"log"`

- sigma.link:

  link function for sigma, default `"log"`

- nu.link:

  link function for nu, default `"log"`

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

- log, log.p:

  logical; if TRUE, probabilities/densities given as log

- lower.tail:

  logical; if TRUE (default), probabilities are P(X \<= x)

## Value

`DEGPD3()` returns an object of class `"gamlss.family"`.

`dDEGPD3` gives the PMF, `pDEGPD3` gives the distribution function,
`qDEGPD3` gives the quantile function, and `rDEGPD3` generates random
deviates.

## References

Ahmad, T. and Arshad, I. A. (2024). New flexible versions of extended
generalized Pareto model for count data. *arXiv preprint*
arXiv:2409.18719. <https://arxiv.org/abs/2409.18719>

## See also

[`DEGPD1`](https://sdwfrost.github.io/egpd/reference/DEGPD1.md),
[`DEGPD4`](https://sdwfrost.github.io/egpd/reference/DEGPD4.md),
[`degpd_bamlss`](https://sdwfrost.github.io/egpd/reference/degpd_bamlss.md)
