# Distribution function of Compound Poisson-Discrete EGPD

Distribution function of Compound Poisson-Discrete EGPD

## Usage

``` r
pcpdegpd(q, lambda, prob = NA, kappa = NA, delta = NA, sigma, xi, type = 1,
  K = NULL)
```

## Arguments

- q:

  numeric vector of quantiles

- lambda:

  Poisson rate parameter (\> 0)

- prob:

  G-transformation parameter (type 6 only)

- kappa:

  G-transformation parameter

- delta:

  G-transformation parameter

- sigma:

  GPD scale parameter (\> 0)

- xi:

  GPD shape parameter

- type:

  integer 1-6 specifying G-transformation type

- K:

  grid size for Panjer recursion (default: automatically determined from
  data range)

## Value

numeric vector of cumulative probabilities

## See also

[`dcpdegpd`](https://sdwfrost.github.io/egpd/reference/dcpdegpd.md),
[`qcpdegpd`](https://sdwfrost.github.io/egpd/reference/qcpdegpd.md),
[`rcpdegpd`](https://sdwfrost.github.io/egpd/reference/rcpdegpd.md)
