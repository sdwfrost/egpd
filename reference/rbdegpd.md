# \\Experimental\\ Random generation from the bivariate discrete EGPD

Generates random samples from the Bivariate Discrete Extended
Generalized Pareto Distribution (BDEGPD) by applying
[`floor()`](https://rdrr.io/r/base/Round.html) to continuous bivariate
BEGPD samples. This function is pure R and does not require Julia.

## Usage

``` r
rbdegpd(n, kappa, sigma, xi, thL, thU, thw)
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

## Value

An `n` by 2 integer matrix with columns `Y1` and `Y2`.

## Details

The BDEGPD constructs bivariate discrete observations by generating
continuous BEGPD samples via
[`rbegpd`](https://sdwfrost.github.io/egpd/reference/rbegpd.md) and
applying [`floor()`](https://rdrr.io/r/base/Round.html) to obtain
non-negative integers.

This function is **experimental** and its interface may change.

## See also

[`rbegpd`](https://sdwfrost.github.io/egpd/reference/rbegpd.md),
[`rbzidegpd`](https://sdwfrost.github.io/egpd/reference/rbzidegpd.md)

## Examples

``` r
Y <- rbdegpd(1000, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
table(Y[,1])
#> 
#>   0   1   2   3   4   5   6 
#> 705 216  51  18   7   2   1 
```
