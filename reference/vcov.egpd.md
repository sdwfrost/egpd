# Variance-covariance matrix for egpd fits

Variance-covariance matrix for egpd fits

## Usage

``` r
# S3 method for class 'egpd'
vcov(object, ...)
```

## Arguments

- object:

  A fitted `egpd` object.

- ...:

  Not used.

## Value

The variance-covariance matrix `object$Vp` with row and column names
matching the coefficient names.
