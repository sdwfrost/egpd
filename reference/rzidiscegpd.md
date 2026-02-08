# Random generation from the Zero-Inflated Discrete Extended GPD

Random generation from the Zero-Inflated Discrete Extended GPD

## Usage

``` r
rzidiscegpd(
  n,
  pi = NA,
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

- unifsamp:

  optional uniform samples

## Value

Random non-negative integer samples

## References

Ahmad, T. and Arshad, I. A. (2024). New flexible versions of extended
generalized Pareto model for count data. *arXiv preprint*
arXiv:2409.18719. <https://arxiv.org/abs/2409.18719>

Ahmad, T. and Hussain, A. (2025). Flexible model for varying levels of
zeros and outliers in count data. *arXiv preprint* arXiv:2510.27365.
<https://arxiv.org/abs/2510.27365>
