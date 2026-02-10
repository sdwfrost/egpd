# gamlss Family for Continuous EGPD Model 3

The `EGPD3()` function defines a `gamlss.family` distribution for
fitting a continuous Extended Generalized Pareto Distribution with model
3 (incomplete beta G-transformation).

The functions `dEGPD3`, `pEGPD3`, `qEGPD3`, and `rEGPD3` define the
density, distribution function, quantile function, and random
generation.

## Usage

``` r
EGPD3(mu.link = "log", sigma.link = "log", nu.link = "log")

dEGPD3(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE)
pEGPD3(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
qEGPD3(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
rEGPD3(n, mu = 1, sigma = 0.5, nu = 1)
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

  G-transformation parameter (delta), positive

- log, log.p:

  logical; if TRUE, probabilities/densities given as log

- lower.tail:

  logical; if TRUE (default), probabilities are P(X \<= x)

## Value

`EGPD3()` returns an object of class `"gamlss.family"`.

`dEGPD3` gives the density, `pEGPD3` gives the distribution function,
`qEGPD3` gives the quantile function, and `rEGPD3` generates random
deviates.

## References

Abbas, A., Ahmad, T. and Ahmad, I. (2025). Modeling zero-inflated
precipitation extremes. *arXiv preprint* arXiv:2504.11058.
<https://arxiv.org/abs/2504.11058>

## See also

[`EGPD1`](https://sdwfrost.github.io/egpd/reference/EGPD1.md),
[`EGPD4`](https://sdwfrost.github.io/egpd/reference/EGPD4.md),
[`egpd_bamlss`](https://sdwfrost.github.io/egpd/reference/egpd_bamlss.md)
