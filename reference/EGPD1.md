# gamlss Family for Continuous EGPD Model 1

The `EGPD1()` function defines a `gamlss.family` distribution for
fitting a continuous Extended Generalized Pareto Distribution with model
1 (power G-transformation: \\G(u) = u^\kappa\\).

The functions `dEGPD1`, `pEGPD1`, `qEGPD1`, and `rEGPD1` define the
density, distribution function, quantile function, and random generation
for this distribution using gamlss parameter names.

## Usage

``` r
EGPD1(mu.link = "log", sigma.link = "log", nu.link = "log")

dEGPD1(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE)
pEGPD1(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
qEGPD1(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
rEGPD1(n, mu = 1, sigma = 0.5, nu = 1)
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

  GPD scale parameter (sigma in EGPD notation), positive

- sigma:

  GPD shape parameter (xi in EGPD notation), positive

- nu:

  G-transformation parameter (kappa in EGPD notation), positive

- log, log.p:

  logical; if TRUE, probabilities/densities given as log

- lower.tail:

  logical; if TRUE (default), probabilities are P(X \<= x)

## Value

`EGPD1()` returns an object of class `"gamlss.family"`.

`dEGPD1` gives the density, `pEGPD1` gives the distribution function,
`qEGPD1` gives the quantile function, and `rEGPD1` generates random
deviates.

## Details

The gamlss parameters map to EGPD parameters as: `mu` = sigma (GPD
scale), `sigma` = xi (GPD shape), `nu` = kappa (G-transformation power).

First derivatives are computed numerically via central finite
differences; second derivatives use the outer product of gradients (OPG)
approximation.

## References

Abbas, A., Ahmad, T. and Ahmad, I. (2025). Modeling zero-inflated
precipitation extremes. *arXiv preprint* arXiv:2504.11058.
<https://arxiv.org/abs/2504.11058>

## See also

[`egpd_bamlss`](https://sdwfrost.github.io/egpd/reference/egpd_bamlss.md),
[`degpd_density`](https://sdwfrost.github.io/egpd/reference/degpd_density.md),
[`pegpd`](https://sdwfrost.github.io/egpd/reference/pegpd.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  library(gamlss)
  set.seed(1)
  y <- rEGPD1(500, mu = 2, sigma = 0.2, nu = 1.5)
  fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                family = EGPD1())
} # }
```
