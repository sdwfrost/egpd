# Plot diagnostics for a fitegpd object

Produces a 4-panel diagnostic plot. For univariate fits: histogram with
fitted density, empirical vs fitted CDF, Q-Q plot, and P-P plot. For
bivariate BEGPD fits: observed scatter, simulated scatter, radial Q-Q
plot, and posterior marginals or parameter bar chart.

## Usage

``` r
# S3 method for class 'fitegpd'
plot(x, ...)
```

## Arguments

- x:

  a `fitegpd` object

- ...:

  additional graphical parameters
