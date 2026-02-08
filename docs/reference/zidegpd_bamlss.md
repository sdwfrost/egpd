# bamlss Family for Zero-Inflated Discrete EGPD (ZIDEGPD)

Creates a `family.bamlss` object for fitting zero-inflated discrete
Extended Generalized Pareto Distribution (ZIDEGPD) models with
`bamlss()`. Uses
[`dzidiscegpd`](https://sdwfrost.github.io/egpd/reference/dzidiscegpd.md),
[`pzidiscegpd`](https://sdwfrost.github.io/egpd/reference/pzidiscegpd.md),
[`qzidiscegpd`](https://sdwfrost.github.io/egpd/reference/qzidiscegpd.md),
and
[`rzidiscegpd`](https://sdwfrost.github.io/egpd/reference/rzidiscegpd.md)
internally.

## Usage

``` r
zidegpd_bamlss(m = 1, ...)
```

## Arguments

- m:

  integer 1–4 selecting the G transformation type (see
  [`egpd_bamlss`](https://sdwfrost.github.io/egpd/reference/egpd_bamlss.md)
  for details)

- ...:

  arguments passed to link specification

## Value

An object of class `family.bamlss`

## References

Ahmad, T. and Arshad, I. A. (2024). New flexible versions of extended
generalized Pareto model for count data. *arXiv preprint*
arXiv:2409.18719. <https://arxiv.org/abs/2409.18719>

Ahmad, T. and Hussain, A. (2025). Flexible model for varying levels of
zeros and outliers in count data. *arXiv preprint* arXiv:2510.27365.
<https://arxiv.org/abs/2510.27365>

## See also

[`ZIDEGPD1`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD1.md),
[`egpd_bamlss`](https://sdwfrost.github.io/egpd/reference/egpd_bamlss.md),
[`degpd_bamlss`](https://sdwfrost.github.io/egpd/reference/degpd_bamlss.md),
[`ziegpd_bamlss`](https://sdwfrost.github.io/egpd/reference/ziegpd_bamlss.md),
[`dzidiscegpd`](https://sdwfrost.github.io/egpd/reference/dzidiscegpd.md)
