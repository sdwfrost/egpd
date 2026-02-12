# \\Experimental\\ Random generation from the zero-inflated MDGPD

Generates random samples from the Zero-Inflated MDGPD (ZIMDGPD). With
probability `pi0`, a row is set to all zeros; otherwise, it is drawn
from [`rmdgpd`](https://sdwfrost.github.io/egpd/reference/rmdgpd.md).

## Usage

``` r
rzimdgpd(n, sigma, xi, lambda, rho, pi0, d = 2L)
```

## Arguments

- n:

  integer; number of observations to generate.

- sigma:

  positive numeric; GPD scale parameter (common across dimensions).

- xi:

  non-negative numeric; GPD shape parameter. When `xi = 0` the marginals
  reduce to scaled geometric.

- lambda:

  positive numeric; Poisson rate for the generator components. Controls
  the spread of the dependence structure.

- rho:

  numeric in \\0, 1); equicorrelation of the `d`-dimensional Poisson
  generator. Higher values give stronger positive dependence: `rho -> 1`
  gives near-perfect dependence, `rho = 0` gives weakest dependence for
  the given `lambda`.

- pi0:

  numeric in (0, 1); joint zero-inflation probability.

- d:

  integer \>= 2; dimension of the multivariate distribution (default
  `2L`).

## Value

An `n` by `d` integer matrix with columns `Y1`, `Y2`, ..., `Yd`.

## Details

This function is **experimental** and its interface may change.

## See also

[`rmdgpd`](https://sdwfrost.github.io/egpd/reference/rmdgpd.md)

## Examples

``` r
Y <- rzimdgpd(1000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = 0.3)
mean(rowSums(Y) == 0)  # approximately 0.3 + natural zeros
#> [1] 0.72

# Trivariate
Y3 <- rzimdgpd(1000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = 0.3, d = 3)
```
