# Density of the Zero-Inflated Extended GPD

Continuous zero-inflated EGPD density: \\h(x) = \pi I(x=0) + (1-\pi)
f\_{EGPD}(x)\\ for \\x \ge 0\\.

## Usage

``` r
dziegpd(
  x,
  pi = NA,
  prob = NA,
  kappa = NA,
  delta = NA,
  sigma = NA,
  xi = NA,
  type = 1,
  log = FALSE
)
```

## Arguments

- x:

  non-negative values

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

- log:

  logical: return log-density?

## Value

Density values
