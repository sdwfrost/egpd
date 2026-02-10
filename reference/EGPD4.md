# gamlss Family for Continuous EGPD Model 4

The `EGPD4()` function defines a `gamlss.family` distribution for
fitting a continuous Extended Generalized Pareto Distribution with model
4 (power-beta G-transformation). This is a 4-parameter family.

The functions `dEGPD4`, `pEGPD4`, `qEGPD4`, and `rEGPD4` define the
density, distribution function, quantile function, and random
generation.

## Usage

``` r
EGPD4(mu.link = "log", sigma.link = "log", nu.link = "log",
      tau.link = "log")

dEGPD4(x, mu = 1, sigma = 0.5, nu = 1, tau = 1, log = FALSE)
pEGPD4(q, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE,
       log.p = FALSE)
qEGPD4(p, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE,
       log.p = FALSE)
rEGPD4(n, mu = 1, sigma = 0.5, nu = 1, tau = 1)
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

  G-transformation parameter (delta), positive

- tau:

  G-transformation parameter (kappa), positive

- log, log.p:

  logical; if TRUE, probabilities/densities given as log

- lower.tail:

  logical; if TRUE (default), probabilities are P(X \<= x)

## Value

`EGPD4()` returns an object of class `"gamlss.family"`.

`dEGPD4` gives the density, `pEGPD4` gives the distribution function,
`qEGPD4` gives the quantile function, and `rEGPD4` generates random
deviates.

## References

Abbas, A., Ahmad, T. and Ahmad, I. (2025). Modeling zero-inflated
precipitation extremes. *arXiv preprint* arXiv:2504.11058.
<https://arxiv.org/abs/2504.11058>

## See also

[`EGPD1`](https://sdwfrost.github.io/egpd/reference/EGPD1.md),
[`EGPD3`](https://sdwfrost.github.io/egpd/reference/EGPD3.md),
[`egpd_bamlss`](https://sdwfrost.github.io/egpd/reference/egpd_bamlss.md)
