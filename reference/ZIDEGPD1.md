# gamlss Family for Zero-Inflated Discrete EGPD Model 1

The `ZIDEGPD1()` function defines a `gamlss.family` distribution for
fitting a zero-inflated discrete Extended Generalized Pareto
Distribution with model 1 (power G-transformation).

The functions `dZIDEGPD1`, `pZIDEGPD1`, `qZIDEGPD1`, and `rZIDEGPD1`
define the PMF, distribution function, quantile function, and random
generation.

## Usage

``` r
ZIDEGPD1(mu.link = "log", sigma.link = "log", nu.link = "log",
         tau.link = "logit")

dZIDEGPD1(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE)
pZIDEGPD1(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
          log.p = FALSE)
qZIDEGPD1(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
          log.p = FALSE)
rZIDEGPD1(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1)
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

  G-transformation parameter (kappa), positive

- tau:

  zero-inflation probability (pi), in (0,1)

- log, log.p:

  logical; if TRUE, probabilities/densities given as log

- lower.tail:

  logical; if TRUE (default), probabilities are P(X \<= x)

## Value

`ZIDEGPD1()` returns an object of class `"gamlss.family"`.

`dZIDEGPD1` gives the PMF, `pZIDEGPD1` gives the distribution function,
`qZIDEGPD1` gives the quantile function, and `rZIDEGPD1` generates
random deviates.

## References

Ahmad, T. and Arshad, I. A. (2024). New flexible versions of extended
generalized Pareto model for count data. *arXiv preprint*
arXiv:2409.18719. <https://arxiv.org/abs/2409.18719>

Ahmad, T. and Hussain, A. (2025). Flexible model for varying levels of
zeros and outliers in count data. *arXiv preprint* arXiv:2510.27365.
<https://arxiv.org/abs/2510.27365>

## See also

[`ZIDEGPD3`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD3.md),
[`zidegpd_bamlss`](https://sdwfrost.github.io/egpd/reference/zidegpd_bamlss.md),
[`dzidiscegpd`](https://sdwfrost.github.io/egpd/reference/dzidiscegpd.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  library(gamlss)
  set.seed(1)
  y <- rZIDEGPD1(500, mu = 3, sigma = 0.15, nu = 2, tau = 0.25)
  fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                tau.formula = ~ 1, family = ZIDEGPD1())
} # }
```
