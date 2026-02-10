# gamlss Family for Zero-Inflated Discrete EGPD Model 3

The `ZIDEGPD3()` function defines a `gamlss.family` distribution for
fitting a zero-inflated discrete Extended Generalized Pareto
Distribution with model 3 (incomplete beta G-transformation).

The functions `dZIDEGPD3`, `pZIDEGPD3`, `qZIDEGPD3`, and `rZIDEGPD3`
define the PMF, distribution function, quantile function, and random
generation.

## Usage

``` r
ZIDEGPD3(mu.link = "log", sigma.link = "log", nu.link = "log",
         tau.link = "logit")

dZIDEGPD3(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE)
pZIDEGPD3(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
          log.p = FALSE)
qZIDEGPD3(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
          log.p = FALSE)
rZIDEGPD3(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1)
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

  zero-inflation probability (pi), in (0,1)

- log, log.p:

  logical; if TRUE, probabilities/densities given as log

- lower.tail:

  logical; if TRUE (default), probabilities are P(X \<= x)

## Value

`ZIDEGPD3()` returns an object of class `"gamlss.family"`.

`dZIDEGPD3` gives the PMF, `pZIDEGPD3` gives the distribution function,
`qZIDEGPD3` gives the quantile function, and `rZIDEGPD3` generates
random deviates.

## References

Ahmad, T. and Arshad, I. A. (2024). New flexible versions of extended
generalized Pareto model for count data. *arXiv preprint*
arXiv:2409.18719. <https://arxiv.org/abs/2409.18719>

Ahmad, T. and Hussain, A. (2025). Flexible model for varying levels of
zeros and outliers in count data. *arXiv preprint* arXiv:2510.27365.
<https://arxiv.org/abs/2510.27365>

## See also

[`ZIDEGPD1`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD1.md),
[`zidegpd_bamlss`](https://sdwfrost.github.io/egpd/reference/zidegpd_bamlss.md),
[`dzidiscegpd`](https://sdwfrost.github.io/egpd/reference/dzidiscegpd.md)
