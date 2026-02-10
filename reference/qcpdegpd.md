# Quantile function of Compound Poisson-Discrete EGPD

Quantile function of Compound Poisson-Discrete EGPD

## Usage

``` r
qcpdegpd(p, lambda, prob = NA, kappa = NA, delta = NA, sigma, xi, type = 1,
  K = NULL)
```

## Arguments

- p:

  numeric vector of probabilities

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

  grid size for Panjer recursion (default: 5000)

## Value

numeric vector of quantiles (non-negative integers)

## See also

[`dcpdegpd`](https://sdwfrost.github.io/egpd/reference/dcpdegpd.md),
[`pcpdegpd`](https://sdwfrost.github.io/egpd/reference/pcpdegpd.md),
[`rcpdegpd`](https://sdwfrost.github.io/egpd/reference/rcpdegpd.md)
