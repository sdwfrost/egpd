# CDF of the Zero-Inflated Extended GPD

CDF of the Zero-Inflated Extended GPD

## Usage

``` r
pziegpd(
  q,
  pi = NA,
  prob = NA,
  kappa = NA,
  delta = NA,
  sigma = NA,
  xi = NA,
  type = 1
)
```

## Arguments

- q:

  non-negative quantiles

- pi:

  zero-inflation probability

- prob:

  mixing probability (type 6)

- kappa:

  shape parameter for G transformation

- delta:

  shape parameter for G transformation (types 4-6)

- sigma:

  GPD scale parameter

- xi:

  GPD shape parameter

- type:

  integer 1-6 specifying G type

## Value

CDF values

## References

Abbas, A., Ahmad, T. and Ahmad, I. (2025). Modeling zero-inflated
precipitation extremes. *arXiv preprint* arXiv:2504.11058.
<https://arxiv.org/abs/2504.11058>
