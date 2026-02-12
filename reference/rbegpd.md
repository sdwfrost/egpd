# Random generation from the bivariate BEGPD

Generates random samples from the bivariate Multivariate Extended
Generalized Pareto Distribution (BEGPD). This function is pure R and
does not require Julia.

## Usage

``` r
rbegpd(n, kappa, sigma, xi, thL, thU, thw)
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

An `n` by 2 numeric matrix with columns `Y1` and `Y2`.

## Details

The BEGPD constructs bivariate observations by:

1.  Generating a radial component \\R\\ from a power-transformed GPD.

2.  Generating lower and upper dependence components from symmetric Beta
    distributions with parameters `thL` and `thU`.

3.  Mixing the components using a weight function based on `thw`.

## References

Alotaibi, N., Sainsbury-Dale, M., Naveau, P., Gaetan, C., and Huser, R.
(2025). Joint modeling of low and high extremes using a multivariate
extended generalized Pareto distribution. *arXiv preprint*
arXiv:2509.05982.

## Examples

``` r
Y <- rbegpd(1000, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
plot(Y[,1], Y[,2], pch = ".", main = "Bivariate BEGPD sample")

```
