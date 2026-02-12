# \\Experimental\\ Random generation from the zero-inflated bivariate discrete EGPD

Generates random samples from the Zero-Inflated Bivariate Discrete
Extended Generalized Pareto Distribution (BZIDEGPD). With probability
`pi0`, a row is set to `(0, 0)`; otherwise, it is drawn from
[`rbdegpd`](https://sdwfrost.github.io/egpd/reference/rbdegpd.md).

## Usage

``` r
rbzidegpd(n, kappa, sigma, xi, thL, thU, thw, pi0)
```

## Arguments

- n:

  integer; number of observations to generate.

- kappa:

  positive numeric; EGPD shape parameter (power transform).

- sigma:

  positive numeric; GPD scale parameter.

- xi:

  positive numeric; GPD shape parameter.

- thL:

  positive numeric; lower tail dependence parameter (beta shape).

- thU:

  positive numeric; upper tail dependence parameter (beta shape).

- thw:

  numeric in (0, 0.5); weight mixing parameter.

- pi0:

  numeric in (0, 1); joint zero-inflation probability.

## Value

An `n` by 2 integer matrix with columns `Y1` and `Y2`.

## Details

This function is **experimental** and its interface may change.

## See also

[`rbdegpd`](https://sdwfrost.github.io/egpd/reference/rbdegpd.md),
[`rbegpd`](https://sdwfrost.github.io/egpd/reference/rbegpd.md)

## Examples

``` r
Y <- rbzidegpd(1000, kappa = 2, sigma = 1, xi = 0.1,
               thL = 5, thU = 5, thw = 0.2, pi0 = 0.3)
mean(Y[,1] == 0 & Y[,2] == 0)  # approximately 0.3
#> [1] 0.713
```
