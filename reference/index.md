# Package index

## Package Overview

- [`egpd-package`](https://sdwfrost.github.io/egpd/reference/egpd-package.md)
  : egpd: Extended Generalized Pareto Distribution GAMs

## Model Fitting

Main function for fitting EGPD, DEGPD, and ZIDEGPD GAMs.

- [`egpd()`](https://sdwfrost.github.io/egpd/reference/egpd.md) : Fit
  Extended Generalized Pareto Distribution GAMs

## Prediction and Methods

S3 methods for fitted egpd objects.

- [`predict(`*`<egpd>`*`)`](https://sdwfrost.github.io/egpd/reference/predict.egpd.md)
  :

  Predictions from a fitted `egpd` object

- [`fitted(`*`<egpd>`*`)`](https://sdwfrost.github.io/egpd/reference/fitted.egpd.md)
  : Extract Model Fitted Values

- [`simulate(`*`<egpd>`*`)`](https://sdwfrost.github.io/egpd/reference/simulate.egpd.md)
  :

  Simulations from a fitted `egpd` object

- [`logLik(`*`<egpd>`*`)`](https://sdwfrost.github.io/egpd/reference/logLik.egpd.md)
  :

  Log-likelihood from a fitted `egpd` object

- [`plot(`*`<egpd>`*`)`](https://sdwfrost.github.io/egpd/reference/plot.egpd.md)
  :

  Plot a fitted `egpd` object

- [`summary(`*`<egpd>`*`)`](https://sdwfrost.github.io/egpd/reference/summary.egpd.md)
  [`print(`*`<summary.egpd>`*`)`](https://sdwfrost.github.io/egpd/reference/summary.egpd.md)
  :

  Summary method for a fitted `egpd` object

- [`print(`*`<egpd>`*`)`](https://sdwfrost.github.io/egpd/reference/print.egpd.md)
  :

  Print a fitted `egpd` object

- [`nobs(`*`<egpd>`*`)`](https://sdwfrost.github.io/egpd/reference/nobs.egpd.md)
  : Number of observations

- [`predictive_coverage()`](https://sdwfrost.github.io/egpd/reference/predictive_coverage.md)
  : Predictive coverage for EGPD models

## EGPD Distributions

Continuous Extended Generalized Pareto Distribution functions.

- [`pegpd()`](https://sdwfrost.github.io/egpd/reference/pegpd.md) : CDF
  of the Extended GPD
- [`degpd_density()`](https://sdwfrost.github.io/egpd/reference/degpd_density.md)
  : Density of the Extended GPD
- [`qegpd()`](https://sdwfrost.github.io/egpd/reference/qegpd.md) :
  Quantile function of the Extended GPD
- [`regpd()`](https://sdwfrost.github.io/egpd/reference/regpd.md) :
  Random generation from the Extended GPD

## Discrete EGPD Distributions

Discrete Extended Generalized Pareto Distribution functions.

- [`ddiscegpd()`](https://sdwfrost.github.io/egpd/reference/ddiscegpd.md)
  : Density of the Discrete Extended GPD
- [`pdiscegpd()`](https://sdwfrost.github.io/egpd/reference/pdiscegpd.md)
  : CDF of the Discrete Extended GPD
- [`qdiscegpd()`](https://sdwfrost.github.io/egpd/reference/qdiscegpd.md)
  : Quantile function of the Discrete Extended GPD
- [`rdiscegpd()`](https://sdwfrost.github.io/egpd/reference/rdiscegpd.md)
  : Random generation from the Discrete Extended GPD

## Zero-Inflated Discrete EGPD Distributions

Zero-inflated discrete EGPD distribution functions.

- [`dzidiscegpd()`](https://sdwfrost.github.io/egpd/reference/dzidiscegpd.md)
  : Density of the Zero-Inflated Discrete Extended GPD
- [`pzidiscegpd()`](https://sdwfrost.github.io/egpd/reference/pzidiscegpd.md)
  : CDF of the Zero-Inflated Discrete Extended GPD
- [`qzidiscegpd()`](https://sdwfrost.github.io/egpd/reference/qzidiscegpd.md)
  : Quantile function of the Zero-Inflated Discrete Extended GPD
- [`rzidiscegpd()`](https://sdwfrost.github.io/egpd/reference/rzidiscegpd.md)
  : Random generation from the Zero-Inflated Discrete Extended GPD

## Zero-Inflated Continuous EGPD Distributions

Zero-inflated continuous EGPD distribution functions.

- [`dziegpd()`](https://sdwfrost.github.io/egpd/reference/dziegpd.md) :
  Density of the Zero-Inflated Extended GPD
- [`pziegpd()`](https://sdwfrost.github.io/egpd/reference/pziegpd.md) :
  CDF of the Zero-Inflated Extended GPD
- [`qziegpd()`](https://sdwfrost.github.io/egpd/reference/qziegpd.md) :
  Quantile function of the Zero-Inflated Extended GPD
- [`rziegpd()`](https://sdwfrost.github.io/egpd/reference/rziegpd.md) :
  Random generation from the Zero-Inflated Extended GPD

## Diagnostics

Model diagnostic tools.

- [`rqresid()`](https://sdwfrost.github.io/egpd/reference/rqresid.md) :
  Randomized Quantile Residuals

## gamlss Families

Family constructors and d/p/q/r wrappers for fitting EGPD models with
gamlss.

- [`EGPD1()`](https://sdwfrost.github.io/egpd/reference/EGPD1.md)
  [`dEGPD1()`](https://sdwfrost.github.io/egpd/reference/EGPD1.md)
  [`pEGPD1()`](https://sdwfrost.github.io/egpd/reference/EGPD1.md)
  [`qEGPD1()`](https://sdwfrost.github.io/egpd/reference/EGPD1.md)
  [`rEGPD1()`](https://sdwfrost.github.io/egpd/reference/EGPD1.md) :
  gamlss Family for Continuous EGPD Model 1
- [`EGPD3()`](https://sdwfrost.github.io/egpd/reference/EGPD3.md)
  [`dEGPD3()`](https://sdwfrost.github.io/egpd/reference/EGPD3.md)
  [`pEGPD3()`](https://sdwfrost.github.io/egpd/reference/EGPD3.md)
  [`qEGPD3()`](https://sdwfrost.github.io/egpd/reference/EGPD3.md)
  [`rEGPD3()`](https://sdwfrost.github.io/egpd/reference/EGPD3.md) :
  gamlss Family for Continuous EGPD Model 3
- [`EGPD4()`](https://sdwfrost.github.io/egpd/reference/EGPD4.md)
  [`dEGPD4()`](https://sdwfrost.github.io/egpd/reference/EGPD4.md)
  [`pEGPD4()`](https://sdwfrost.github.io/egpd/reference/EGPD4.md)
  [`qEGPD4()`](https://sdwfrost.github.io/egpd/reference/EGPD4.md)
  [`rEGPD4()`](https://sdwfrost.github.io/egpd/reference/EGPD4.md) :
  gamlss Family for Continuous EGPD Model 4
- [`EGPD5()`](https://sdwfrost.github.io/egpd/reference/EGPD5.md)
  [`dEGPD5()`](https://sdwfrost.github.io/egpd/reference/EGPD5.md)
  [`pEGPD5()`](https://sdwfrost.github.io/egpd/reference/EGPD5.md)
  [`qEGPD5()`](https://sdwfrost.github.io/egpd/reference/EGPD5.md)
  [`rEGPD5()`](https://sdwfrost.github.io/egpd/reference/EGPD5.md) :
  gamlss Family for Continuous EGPD Model 5
- [`EGPD6()`](https://sdwfrost.github.io/egpd/reference/EGPD6.md)
  [`dEGPD6()`](https://sdwfrost.github.io/egpd/reference/EGPD6.md)
  [`pEGPD6()`](https://sdwfrost.github.io/egpd/reference/EGPD6.md)
  [`qEGPD6()`](https://sdwfrost.github.io/egpd/reference/EGPD6.md)
  [`rEGPD6()`](https://sdwfrost.github.io/egpd/reference/EGPD6.md) :
  gamlss Family for Continuous EGPD Model 6
- [`DEGPD1()`](https://sdwfrost.github.io/egpd/reference/DEGPD1.md)
  [`dDEGPD1()`](https://sdwfrost.github.io/egpd/reference/DEGPD1.md)
  [`pDEGPD1()`](https://sdwfrost.github.io/egpd/reference/DEGPD1.md)
  [`qDEGPD1()`](https://sdwfrost.github.io/egpd/reference/DEGPD1.md)
  [`rDEGPD1()`](https://sdwfrost.github.io/egpd/reference/DEGPD1.md) :
  gamlss Family for Discrete EGPD Model 1
- [`DEGPD3()`](https://sdwfrost.github.io/egpd/reference/DEGPD3.md)
  [`dDEGPD3()`](https://sdwfrost.github.io/egpd/reference/DEGPD3.md)
  [`pDEGPD3()`](https://sdwfrost.github.io/egpd/reference/DEGPD3.md)
  [`qDEGPD3()`](https://sdwfrost.github.io/egpd/reference/DEGPD3.md)
  [`rDEGPD3()`](https://sdwfrost.github.io/egpd/reference/DEGPD3.md) :
  gamlss Family for Discrete EGPD Model 3
- [`DEGPD4()`](https://sdwfrost.github.io/egpd/reference/DEGPD4.md)
  [`dDEGPD4()`](https://sdwfrost.github.io/egpd/reference/DEGPD4.md)
  [`pDEGPD4()`](https://sdwfrost.github.io/egpd/reference/DEGPD4.md)
  [`qDEGPD4()`](https://sdwfrost.github.io/egpd/reference/DEGPD4.md)
  [`rDEGPD4()`](https://sdwfrost.github.io/egpd/reference/DEGPD4.md) :
  gamlss Family for Discrete EGPD Model 4
- [`DEGPD5()`](https://sdwfrost.github.io/egpd/reference/DEGPD5.md)
  [`dDEGPD5()`](https://sdwfrost.github.io/egpd/reference/DEGPD5.md)
  [`pDEGPD5()`](https://sdwfrost.github.io/egpd/reference/DEGPD5.md)
  [`qDEGPD5()`](https://sdwfrost.github.io/egpd/reference/DEGPD5.md)
  [`rDEGPD5()`](https://sdwfrost.github.io/egpd/reference/DEGPD5.md) :
  gamlss Family for Discrete EGPD Model 5
- [`DEGPD6()`](https://sdwfrost.github.io/egpd/reference/DEGPD6.md)
  [`dDEGPD6()`](https://sdwfrost.github.io/egpd/reference/DEGPD6.md)
  [`pDEGPD6()`](https://sdwfrost.github.io/egpd/reference/DEGPD6.md)
  [`qDEGPD6()`](https://sdwfrost.github.io/egpd/reference/DEGPD6.md)
  [`rDEGPD6()`](https://sdwfrost.github.io/egpd/reference/DEGPD6.md) :
  gamlss Family for Discrete EGPD Model 6
- [`ZIEGPD1()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD1.md)
  [`dZIEGPD1()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD1.md)
  [`pZIEGPD1()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD1.md)
  [`qZIEGPD1()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD1.md)
  [`rZIEGPD1()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD1.md) :
  gamlss Family for Zero-Inflated Continuous EGPD Model 1
- [`ZIEGPD3()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD3.md)
  [`dZIEGPD3()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD3.md)
  [`pZIEGPD3()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD3.md)
  [`qZIEGPD3()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD3.md)
  [`rZIEGPD3()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD3.md) :
  gamlss Family for Zero-Inflated Continuous EGPD Model 3
- [`ZIEGPD5()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD5.md)
  [`dZIEGPD5()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD5.md)
  [`pZIEGPD5()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD5.md)
  [`qZIEGPD5()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD5.md)
  [`rZIEGPD5()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD5.md) :
  gamlss Family for Zero-Inflated Continuous EGPD Model 5
- [`ZIEGPD6()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD6.md)
  [`dZIEGPD6()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD6.md)
  [`pZIEGPD6()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD6.md)
  [`qZIEGPD6()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD6.md)
  [`rZIEGPD6()`](https://sdwfrost.github.io/egpd/reference/ZIEGPD6.md) :
  gamlss Family for Zero-Inflated Continuous EGPD Model 6
- [`ZIDEGPD1()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD1.md)
  [`dZIDEGPD1()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD1.md)
  [`pZIDEGPD1()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD1.md)
  [`qZIDEGPD1()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD1.md)
  [`rZIDEGPD1()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD1.md)
  : gamlss Family for Zero-Inflated Discrete EGPD Model 1
- [`ZIDEGPD3()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD3.md)
  [`dZIDEGPD3()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD3.md)
  [`pZIDEGPD3()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD3.md)
  [`qZIDEGPD3()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD3.md)
  [`rZIDEGPD3()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD3.md)
  : gamlss Family for Zero-Inflated Discrete EGPD Model 3
- [`ZIDEGPD5()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD5.md)
  [`dZIDEGPD5()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD5.md)
  [`pZIDEGPD5()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD5.md)
  [`qZIDEGPD5()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD5.md)
  [`rZIDEGPD5()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD5.md)
  : gamlss Family for Zero-Inflated Discrete EGPD Model 5
- [`ZIDEGPD6()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD6.md)
  [`dZIDEGPD6()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD6.md)
  [`pZIDEGPD6()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD6.md)
  [`qZIDEGPD6()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD6.md)
  [`rZIDEGPD6()`](https://sdwfrost.github.io/egpd/reference/ZIDEGPD6.md)
  : gamlss Family for Zero-Inflated Discrete EGPD Model 6

## bamlss Families

Family constructors for fitting EGPD models with bamlss.

- [`egpd_bamlss()`](https://sdwfrost.github.io/egpd/reference/egpd_bamlss.md)
  : bamlss Family for Continuous EGPD
- [`ziegpd_bamlss()`](https://sdwfrost.github.io/egpd/reference/ziegpd_bamlss.md)
  : bamlss Family for Zero-Inflated Continuous EGPD
- [`degpd_bamlss()`](https://sdwfrost.github.io/egpd/reference/degpd_bamlss.md)
  : bamlss Family for Discrete EGPD (DEGPD)
- [`zidegpd_bamlss()`](https://sdwfrost.github.io/egpd/reference/zidegpd_bamlss.md)
  : bamlss Family for Zero-Inflated Discrete EGPD (ZIDEGPD)

## G-function Utilities

Transformation functions for the EGPD G-function (types 1–6).

- [`p.G()`](https://sdwfrost.github.io/egpd/reference/p.G.md) :
  Transformation CDF for EGPD
- [`d.G()`](https://sdwfrost.github.io/egpd/reference/d.G.md) :
  Transformation density for EGPD
- [`q.G()`](https://sdwfrost.github.io/egpd/reference/q.G.md) : Inverse
  transformation (quantile) for EGPD
- [`r.G()`](https://sdwfrost.github.io/egpd/reference/r.G.md) : Random
  generation from transformation G

## Helpers

Utility functions.

- [`dfbind()`](https://sdwfrost.github.io/egpd/reference/dfbind.md) :
  Bind a list of data frames
- [`pinv()`](https://sdwfrost.github.io/egpd/reference/pinv.md)
  [`ginv.egpd()`](https://sdwfrost.github.io/egpd/reference/pinv.md) :
  Moore-Penrose pseudo-inverse of a matrix
- [`seq_between()`](https://sdwfrost.github.io/egpd/reference/seq_between.md)
  : Generate a sequence between a range

## Datasets

Example datasets included in the package.

- [`ny_complaints`](https://sdwfrost.github.io/egpd/reference/ny_complaints.md)
  : New York Auto Insurance Complaint Rankings
- [`docvisits`](https://sdwfrost.github.io/egpd/reference/docvisits.md)
  : Doctor Visits
- [`nsw_offenses`](https://sdwfrost.github.io/egpd/reference/nsw_offenses.md)
  : Gaming and Betting Offenses in New South Wales
