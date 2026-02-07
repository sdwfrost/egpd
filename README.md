# egpd: Extended Generalized Pareto Distribution GAMs

The **egpd** package fits Extended Generalized Pareto Distribution (EGPD),
Discrete EGPD (DEGPD), and Zero-Inflated Discrete EGPD (ZIDEGPD) models
within a GAM (Generalized Additive Model) framework, with additional
support for fitting via
[gamlss](https://CRAN.R-project.org/package=gamlss) and
[bamlss](https://CRAN.R-project.org/package=bamlss). The fitting
infrastructure is adapted from the
[evgam](https://CRAN.R-project.org/package=evgam) package by Ben Youngman.

Models 1--4 are supported for each distribution family, following the
parameterizations of Naveau et al. (2016).

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("sdwfrost/egpd")
```

## Quick start

```r
library(egpd)

# Fit a DEGPD model to insurance complaint counts
data(ny_complaints)
fit <- egpd(
  list(upheld ~ s(year), ~ 1, ~ 1),
  data = ny_complaints,
  family = "degpd",
  degpd.args = list(m = 1)
)
summary(fit)
```

## Vignettes

- [Discrete EGPD Models for Insurance Complaint Counts](articles/insurance-complaints.html)
- [Zero-Inflated Discrete EGPD Models for Doctor Visit Counts](articles/doctor-visits.html)
- [Threshold Exceedance Modeling with DEGPD](articles/gaming-offenses.html)
- [Continuous EGPD Models for Temperature Extremes](articles/temperature-extremes.html)
- [Simulation Examples](articles/simulation.html)
- [Comparing egpd and bamlss Fits](articles/bamlss-comparison.html)
- [Comparing egpd and gamlss Fits](articles/gamlss-comparison.html)
- [Predictive Coverage Assessment](articles/predictive-coverage.html)
- [Parameter Estimate Coverage Assessment](articles/parameter-coverage.html)

## References

Abbas, A., Ahmad, T. and Ahmad, I. (2025). Modeling zero-inflated
precipitation extremes. *arXiv preprint* arXiv:2504.11058.
<https://arxiv.org/abs/2504.11058>

Ahmad, T. and Arshad, I. A. (2024). New flexible versions of extended
generalized Pareto model for count data. *arXiv preprint*
arXiv:2409.18719. <https://arxiv.org/abs/2409.18719>

Ahmad, T. and Hussain, A. (2025). Flexible model for varying levels of
zeros and outliers in count data. *arXiv preprint* arXiv:2510.27365.
<https://arxiv.org/abs/2510.27365>

Naveau, P., Huser, R., Ribereau, P., and Hannart, A. (2016). Modeling
jointly low, moderate, and heavy rainfall intensities without a
threshold selection. *Water Resources Research*, 52(4), 2897--2911.
