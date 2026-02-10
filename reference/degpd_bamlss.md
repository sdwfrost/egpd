# bamlss Family for Discrete EGPD (DEGPD)

Creates a `family.bamlss` object for fitting discrete Extended
Generalized Pareto Distribution (DEGPD) models with `bamlss()`. Uses
[`ddiscegpd`](https://sdwfrost.github.io/egpd/reference/ddiscegpd.md),
[`pdiscegpd`](https://sdwfrost.github.io/egpd/reference/pdiscegpd.md),
[`qdiscegpd`](https://sdwfrost.github.io/egpd/reference/qdiscegpd.md),
and
[`rdiscegpd`](https://sdwfrost.github.io/egpd/reference/rdiscegpd.md)
internally.

## Usage

``` r
degpd_bamlss(m = 1, ...)
```

## Arguments

- m:

  integer 1–4 selecting the G transformation type:

  1

  :   Power: \\G(u) = u^\kappa\\

  2

  :   Mixture: \\G(u) = p u^{\kappa} + (1-p) u^{\delta}\\

  3

  :   Incomplete beta

  4

  :   Power-beta

- ...:

  arguments passed to link specification

## Value

An object of class `family.bamlss`

## References

Ahmad, T. and Arshad, I. A. (2024). New flexible versions of extended
generalized Pareto model for count data. *arXiv preprint*
arXiv:2409.18719. <https://arxiv.org/abs/2409.18719>

## See also

[`DEGPD1`](https://sdwfrost.github.io/egpd/reference/DEGPD1.md),
[`egpd_bamlss`](https://sdwfrost.github.io/egpd/reference/egpd_bamlss.md),
[`zidegpd_bamlss`](https://sdwfrost.github.io/egpd/reference/zidegpd_bamlss.md),
[`ddiscegpd`](https://sdwfrost.github.io/egpd/reference/ddiscegpd.md)
