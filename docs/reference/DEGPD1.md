# gamlss Family for Discrete EGPD Model 1

The `DEGPD1()` function defines a `gamlss.family` distribution for
fitting a discrete Extended Generalized Pareto Distribution with model 1
(power G-transformation).

The functions `dDEGPD1`, `pDEGPD1`, `qDEGPD1`, and `rDEGPD1` define the
PMF, distribution function, quantile function, and random generation.

## Usage

``` r
DEGPD1(mu.link = "log", sigma.link = "log", nu.link = "log")

dDEGPD1(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE)
pDEGPD1(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
qDEGPD1(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
rDEGPD1(n, mu = 1, sigma = 0.5, nu = 1)
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

  G-transformation parameter (kappa), positive

- log, log.p:

  logical; if TRUE, probabilities/densities given as log

- lower.tail:

  logical; if TRUE (default), probabilities are P(X \<= x)

## Value

`DEGPD1()` returns an object of class `"gamlss.family"`.

`dDEGPD1` gives the PMF, `pDEGPD1` gives the distribution function,
`qDEGPD1` gives the quantile function, and `rDEGPD1` generates random
deviates.

## References

Ahmad, T. and Arshad, I. A. (2024). New flexible versions of extended
generalized Pareto model for count data. *arXiv preprint*
arXiv:2409.18719. <https://arxiv.org/abs/2409.18719>

## See also

[`DEGPD3`](https://sdwfrost.github.io/egpd/reference/DEGPD3.md),
[`DEGPD4`](https://sdwfrost.github.io/egpd/reference/DEGPD4.md),
[`degpd_bamlss`](https://sdwfrost.github.io/egpd/reference/degpd_bamlss.md),
[`ddiscegpd`](https://sdwfrost.github.io/egpd/reference/ddiscegpd.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  library(gamlss)
  set.seed(1)
  y <- rDEGPD1(500, mu = 3, sigma = 0.15, nu = 2)
  fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                family = DEGPD1())
} # }
```
