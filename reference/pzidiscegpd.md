# CDF of the Zero-Inflated Discrete Extended GPD

CDF of the Zero-Inflated Discrete Extended GPD

## Usage

``` r
pzidiscegpd(
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

  non-negative integer quantiles

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
