# bamlss family for discrete EGPD (DEGPD)

Creates a `family.bamlss` object for fitting discrete Extended
Generalized Pareto Distribution models with `bamlss()`.

## Usage

``` r
degpd_bamlss(m = 1, ...)
```

## Arguments

- m:

  integer 1–4 selecting the G transformation (see
  [`egpd_bamlss`](https://sdwfrost.github.io/egpd/reference/egpd_bamlss.md)
  for details)

- ...:

  arguments passed to link specification

## Value

An object of class `family.bamlss`
