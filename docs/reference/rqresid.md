# Randomized Quantile Residuals

Computes randomized quantile residuals (Dunn & Smyth, 1996) for discrete
or continuous EGPD models. For discrete models, a uniform random variate
is drawn between the lower and upper CDF bounds at each observation and
transformed to the normal scale via
[`qnorm`](https://rdrr.io/r/stats/Normal.html). For continuous models,
the CDF value is transformed directly.

## Usage

``` r
rqresid(object, ...)

# S3 method for class 'egpd'
rqresid(object, seed = NULL, ...)
```

## Arguments

- object:

  a fitted [`egpd`](https://sdwfrost.github.io/egpd/reference/egpd.md)
  object

- seed:

  optional random seed for reproducibility

- ...:

  unused

## Value

A numeric vector of randomized quantile residuals on the standard normal
scale. Values that are infinite are replaced with `NA`.

## Details

For a well-specified model the residuals should be approximately
standard normal. The generic `rqresid` dispatches to the `rqresid.egpd`
method for objects of class `"egpd"`.

The method supports all EGPD families (`"egpd"`, `"degpd"`, `"zidegpd"`)
and all model types (1–6).

## References

Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
*Journal of Computational and Graphical Statistics*, 5(3), 236–244.

## See also

[`egpd`](https://sdwfrost.github.io/egpd/reference/egpd.md),
[`predict.egpd`](https://sdwfrost.github.io/egpd/reference/predict.egpd.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  set.seed(42)
  dat <- data.frame(y = regpd(200, sigma = 1, xi = 0.2, kappa = 2, type = 1))
  fit <- egpd(y ~ 1, data = dat, family = "egpd", egpd = list(m = 1))
  r <- rqresid(fit, seed = 1)
  qqnorm(r); qqline(r)
} # }
```
