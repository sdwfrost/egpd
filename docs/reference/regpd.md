# Random generation from the Extended GPD

Random generation from the Extended GPD

## Usage

``` r
regpd(
  n,
  prob = NA,
  kappa = NA,
  delta = NA,
  sigma = NA,
  xi = NA,
  type = 1,
  unifsamp = NULL
)
```

## Arguments

- n:

  number of samples

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

- unifsamp:

  optional uniform samples

## Value

Random samples

## References

Abbas, A., Ahmad, T. and Ahmad, I. (2025). Modeling zero-inflated
precipitation extremes. *arXiv preprint* arXiv:2504.11058.
<https://arxiv.org/abs/2504.11058>
