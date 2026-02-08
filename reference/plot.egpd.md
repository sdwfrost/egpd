# Plot a fitted `egpd` object

Plot a fitted `egpd` object

## Usage

``` r
# S3 method for class 'egpd'
plot(x, onepage = TRUE, which = NULL, main, ask = !onepage, ...)
```

## Arguments

- x:

  a fitted `egpd` object

- onepage:

  logical: should all plots be on one page? Defaults to TRUE

- which:

  a vector of integers identifying which smooths to plot

- main:

  a character string or vector of plot titles

- ask:

  logical: ask to show next plots?

- ...:

  extra arguments to pass to
  [`mgcv::plot.gam`](https://rdrr.io/pkg/mgcv/man/plot.gam.html)

## Value

Plots representing smooth terms
